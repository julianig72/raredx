"""
raredx web — FastAPI backend for clinical VCF analysis.

A clinician uploads a VCF, optionally pastes a clinical note and/or HPO terms, picks the
genome build and which AI layers to enable, and gets back a prioritized candidate report.

Jobs run in a background thread with live progress; the browser polls /api/status.
Results (HTML report + CSV) are served from a per-job temp directory.

Run:
    uvicorn web.server:app --host 0.0.0.0 --port 8000
    # or: python -m web.server

SECURITY / DEPLOYMENT NOTES (read before exposing to real clinicians):
  * No PHI is persisted beyond the job's temp dir; set RAREDX_DATA_DIR to a volume with
    an appropriate retention/wipe policy, or leave it in the OS temp dir (cleared on reboot).
  * This dev server has NO authentication. Put it behind an authenticating reverse proxy
    (OAuth2/OIDC) and TLS before any clinical use. Do not expose port 8000 directly.
  * Uploads are size-capped (RAREDX_MAX_MB, default 50). VCF is parsed, never executed.
  * The tool is DECISION SUPPORT, not a diagnosis (surfaced in the UI and every report).
"""
import asyncio, datetime, json, os, re, sys, uuid, threading, shutil, tempfile, time, traceback
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager, suppress
from pathlib import Path

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

# import the pipeline (parent dir on path)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raredx_pipeline as rx

DATA_DIR = Path(os.environ.get("RAREDX_DATA_DIR", tempfile.gettempdir())) / "raredx_jobs"
DATA_DIR.mkdir(parents=True, exist_ok=True)
# Persistent patient-session store (repo-local by default). Unlike DATA_DIR this is NEVER
# TTL-cleaned, so past results stay recoverable. Holds PHI (VCF + report) → kept out of git
# via .gitignore; point RAREDX_SESSIONS_DIR at an encrypted volume for real clinical use.
SESSIONS_DIR = Path(os.environ.get(
    "RAREDX_SESSIONS_DIR", Path(__file__).resolve().parent.parent / "patient_sessions"))
SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
MAX_MB = int(os.environ.get("RAREDX_MAX_MB", "50"))
CONTACT_EMAIL = os.environ.get("RAREDX_EMAIL")  # optional NCBI contact
MAX_ACTIVE_JOBS = max(1, int(os.environ.get("RAREDX_MAX_ACTIVE_JOBS", "4")))
JOB_TTL_SECONDS = max(3600, int(os.environ.get("RAREDX_JOB_TTL_HOURS", "24")) * 3600)
MAX_NOTE_CHARS = max(1000, int(os.environ.get("RAREDX_MAX_NOTE_CHARS", "20000")))
MAX_HPO_CHARS = max(1000, int(os.environ.get("RAREDX_MAX_HPO_CHARS", "10000")))
MAX_SAMPLE_CHARS = max(32, int(os.environ.get("RAREDX_MAX_SAMPLE_CHARS", "200")))
MAX_LLM_REQUESTS = max(1, int(os.environ.get("RAREDX_MAX_LLM_REQUESTS", "2")))
MAX_VARIANTS = max(1, int(os.environ.get("RAREDX_MAX_VARIANTS", "50000")))
JOB_TIMEOUT_SECONDS = max(60, int(os.environ.get("RAREDX_JOB_TIMEOUT_SECONDS", "7200")))


@asynccontextmanager
async def lifespan(app):
    _cleanup_jobs()
    async def cleanup_loop():
        interval=max(60,min(JOB_TTL_SECONDS//4,3600))
        while True:
            await asyncio.sleep(interval)
            await asyncio.to_thread(_cleanup_jobs)
    cleanup_task=asyncio.create_task(cleanup_loop())
    try:
        yield
    finally:
        cleanup_task.cancel()
        with suppress(asyncio.CancelledError):
            await cleanup_task
        await asyncio.to_thread(_cleanup_jobs)


app = FastAPI(
    title="raredx",
    description="AI-assisted rare-disease VCF analysis (decision support)",
    lifespan=lifespan,
)

# in-memory job registry: job_id -> {status, done, total, message, error, prefix, sample, started}
JOBS = {}
JOBS_LOCK = threading.Lock()
JOB_SLOTS = threading.BoundedSemaphore(MAX_ACTIVE_JOBS)
EXECUTOR = ThreadPoolExecutor(max_workers=MAX_ACTIVE_JOBS, thread_name_prefix="raredx")
LLM_SLOTS = threading.BoundedSemaphore(MAX_LLM_REQUESTS)

def _set(job_id, **kw):
    with JOBS_LOCK:
        JOBS.setdefault(job_id, {}).update(kw)


def _cleanup_jobs():
    cutoff=time.time()-JOB_TTL_SECONDS
    with JOBS_LOCK:
        expired=[
            job_id for job_id,j in JOBS.items()
            if j.get("status") in ("done","error","cancelled")
            and j.get("finished",j.get("started",0))<cutoff
        ]
        for job_id in expired:
            JOBS.pop(job_id,None)
        active=set(JOBS)
    for path in DATA_DIR.iterdir():
        try:
            if path.is_dir() and path.name not in active and path.stat().st_mtime<cutoff:
               shutil.rmtree(path,ignore_errors=True)
        except OSError:
            continue


JOB_ID_RE = re.compile(r"^[0-9a-fA-F]{6,32}$")


def _job_dir(job_id):
    """Resolve a job's on-disk directory, guarding against path traversal."""
    if not JOB_ID_RE.match(job_id or ""):
        return None
    root=DATA_DIR.resolve()
    candidate=(root / job_id).resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        return None
    return candidate if candidate.is_dir() else None


def _session_dir(job_id):
    """Resolve a persisted patient-session directory, guarding against path traversal."""
    if not JOB_ID_RE.match(job_id or ""):
        return None
    root=SESSIONS_DIR.resolve()
    candidate=(root / job_id).resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        return None
    return candidate if candidate.is_dir() else None


def _job_file(job_id, in_memory_path, disk_name):
    """Return an existing result file: prefer the in-memory path, then the live job dir,
    then the persistent session store (which survives DATA_DIR TTL cleanup)."""
    if in_memory_path and os.path.exists(in_memory_path):
        return in_memory_path
    for directory in (_job_dir(job_id), _session_dir(job_id)):
        if directory is not None:
            candidate=directory / disk_name
            if candidate.exists():
                return str(candidate)
    return None


def _persist_session(job_id, meta):
    """Copy a finished analysis into the persistent patient-session store so results stay
    recoverable after DATA_DIR is TTL-cleaned. Best-effort: never breaks the running job."""
    if not JOB_ID_RE.match(job_id or ""):
        return
    try:
        dest=SESSIONS_DIR / job_id
        dest.mkdir(parents=True, exist_ok=True)
        src=DATA_DIR / job_id
        for name in ("result_report.html","result_annotated.csv","input.vcf"):
            source=src / name
            if source.exists():
                shutil.copy2(source, dest / name)
        (dest / "metadata.json").write_text(
            json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError:
        pass


def _run_job(job_id, vcf_path, sample, hpo, note, assembly, use_esm, use_am, use_agentic, reflect_k,
             father_vcf=None, mother_vcf=None, expanded_hpo=None, cancel_event=None):
    outdir = DATA_DIR / job_id
    prefix = str(outdir / "result")
    rx.set_request_deadline(time.monotonic()+JOB_TIMEOUT_SECONDS)
    rx.set_request_cancel_event(cancel_event)
    try:
        def progress(done, total, message):
            _set(job_id,done=done,total=total,message=message,updated=time.time())
        _set(job_id,status="running",done=0,total=0,message="Iniciando…",
             run_started=time.time(),updated=time.time())
        result = rx.run_pipeline(
            vcf_path, sample=sample, hpo=hpo, clinical_note_text=note or None,
            assembly=assembly, use_esm=use_esm, use_am=use_am, agentic=use_agentic,
            reflect_k=reflect_k, father_vcf=father_vcf, mother_vcf=mother_vcf,
            email=CONTACT_EMAIL, progress=progress, max_variants=MAX_VARIANTS,
            reviewed_hpo_expansion=expanded_hpo,
        )
        if cancel_event is not None and cancel_event.is_set():
            raise rx.AnalysisCancelled("analysis cancelled by user")
        rx.write_outputs(result, prefix)
        variants = result["variants"]
        top = [
            {k: v.get(k) for k in ("rank","gene","consequence","protein","af","clinvar",
                                   "af_status","call","variant_score","pheno_score","pheno_disease",
                                   "pheno_shared","combined","am_pathogenicity","esm2_llr","acmg_tags")}
            for v in variants[:15]
        ]
        # surface the agentic differential (disease-level hypotheses) to the UI, if computed
        ag = result.get("agentic") or {}
        differential = [
            {k: d.get(k) for k in ("disease","genes","inheritance","likelihood",
                                   "supporting_variants","rationale","evidence","next_steps")}
            for d in ag.get("differential", [])
        ]
        _set(job_id, status="done", message="Análisis completo",
             n_variants=len(variants), top=top, differential=differential,
             warnings=result.get("warnings",[]),
             llm_providers=result.get("llm_providers",[]),
             report=f"{prefix}_report.html", csv=f"{prefix}_annotated.csv",
             finished=time.time(),updated=time.time())
        ai_layers=[]
        if use_esm: ai_layers.append("ESM-2")
        if use_am: ai_layers.append("AlphaMissense")
        if use_agentic: ai_layers.append("Agentic")
        if father_vcf or mother_vcf: ai_layers.append("Trio")
        top_variant = top[0] if top else {}
        _persist_session(job_id, {
            "job_id": job_id,
            "sample": sample or "",
            "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "created_ts": time.time(),
            "assembly": result.get("assembly") or assembly,
            "hpo": hpo or "",
            "expanded_hpo": bool(expanded_hpo),
            "n_variants": len(variants),
            "top_gene": top_variant.get("gene"),
            "top_call": top_variant.get("call"),
            "ai_layers": ai_layers,
            "n_warnings": len(result.get("warnings",[])),
        })
    except rx.AnalysisCancelled:
        shutil.rmtree(outdir,ignore_errors=True)
        _set(job_id,status="cancelled",message="Análisis detenido por el usuario",
             finished=time.time(),updated=time.time())
    except Exception as e:
        _set(job_id,status="error",message=str(e),error=traceback.format_exc(),
             finished=time.time(),updated=time.time())
    finally:
        rx.set_request_deadline(None)
        rx.set_request_cancel_event(None)
        JOB_SLOTS.release()


@app.post("/api/extract-hpo")
def extract_hpo(clinical_note: str = Form(...)):
    note=clinical_note.strip()
    if not note:
        raise HTTPException(400,"clinical note is empty")
    if len(note)>MAX_NOTE_CHARS:
        raise HTTPException(413,f"clinical note exceeds {MAX_NOTE_CHARS} characters")
    if not LLM_SLOTS.acquire(blocking=False):
        raise HTTPException(503,"LLM extraction capacity reached; retry later")
    try:
        rx._reset_llm_state()
        terms=rx.extract_hpo_from_note(note)
        diagnostics=rx.llm_diagnostics()
        if not terms and diagnostics["errors"]:
            raise HTTPException(503,"HPO extraction failed; verify GitHub Copilot authentication")
        return {"terms":terms,"llm_providers":diagnostics["providers"],
                "warnings":diagnostics["errors"]}
    finally:
        LLM_SLOTS.release()


@app.post("/api/expand-hpo")
def expand_hpo(hpo: str = Form(...)):
    if len(hpo)>MAX_HPO_CHARS:
        raise HTTPException(413,f"HPO input exceeds {MAX_HPO_CHARS} characters")
    tokens=[token for token in re.split(r"[,\n]",hpo) if token.strip()]
    if not tokens:
        return {"terms":[],"warnings":[]}
    terms,available=rx.expand_hpo_profile(tokens,return_status=True)
    if not terms:
        raise HTTPException(400,"no HPO terms could be resolved")
    warnings=[] if available else ["HPO expansion was partially unavailable"]
    return {"terms":terms,"warnings":warnings}


@app.post("/api/analyze")
async def analyze(
    vcf: UploadFile = File(...),
    sample: str = Form("SAMPLE"),
    hpo: str = Form(""),
    expanded_hpo: str | None = Form(None),
    clinical_note: str = Form(""),
    assembly: str = Form("auto"),
    esm: str = Form("false"),
    alphamissense: str = Form("false"),
    agentic: str = Form("false"),
    reflect_k: str = Form("8"),
    father: UploadFile = File(None),
    mother: UploadFile = File(None),
):
    _cleanup_jobs()
    if assembly not in ("auto","GRCh38","GRCh37"):
        raise HTTPException(400,"assembly must be auto, GRCh38, or GRCh37")
    sample=sample.strip() or "SAMPLE"
    if len(sample)>MAX_SAMPLE_CHARS:
        raise HTTPException(413,f"sample ID exceeds {MAX_SAMPLE_CHARS} characters")
    if len(hpo)>MAX_HPO_CHARS:
        raise HTTPException(413,f"HPO input exceeds {MAX_HPO_CHARS} characters")
    if expanded_hpo is not None and len(expanded_hpo)>MAX_HPO_CHARS:
        raise HTTPException(413,f"expanded HPO input exceeds {MAX_HPO_CHARS} characters")
    if len(clinical_note)>MAX_NOTE_CHARS:
        raise HTTPException(413,f"clinical note exceeds {MAX_NOTE_CHARS} characters")
    bool_values={"esm":esm,"alphamissense":alphamissense,"agentic":agentic}
    invalid=[name for name,value in bool_values.items() if value.lower() not in ("true","false")]
    if invalid:
        raise HTTPException(400,f"{', '.join(invalid)} must be true or false")
    try:
        reflect_k_i = max(1, min(int(reflect_k), 50))
    except (TypeError, ValueError):
        raise HTTPException(400, "reflect_k must be an integer")
    job_id = uuid.uuid4().hex[:12]
    outdir = DATA_DIR / job_id
    outdir.mkdir(parents=True, exist_ok=True)
    vcf_path = outdir / "input.vcf"

    async def _stream(upload, dest):
        """Stream an uploaded file to dest with the size cap; return bytes written."""
        n = 0
        exceeded=False
        with open(dest, "wb") as fh:
            while chunk := await upload.read(1 << 20):
                n += len(chunk)
                if n > MAX_MB * (1 << 20):
                    exceeded=True
                    break
                fh.write(chunk)
        if exceeded:
            shutil.rmtree(outdir,ignore_errors=True)
            raise HTTPException(413,f"VCF exceeds {MAX_MB} MB limit")
        return n

    if await _stream(vcf, vcf_path) == 0:
        shutil.rmtree(outdir, ignore_errors=True)
        raise HTTPException(400, "empty upload")

    # optional trio parents (de novo detection → PS2)
    father_path = mother_path = None
    if father is not None and getattr(father, "filename", None):
        father_path = str(outdir / "father.vcf"); await _stream(father, father_path)
    if mother is not None and getattr(mother, "filename", None):
        mother_path = str(outdir / "mother.vcf"); await _stream(mother, mother_path)

    def _validate_and_detect():
        rx.parse_vcf(
            vcf_path,sample=sample,allow_empty=False,max_variants=MAX_VARIANTS,called_only=True
        )
        resolved = rx.detect_vcf_assembly(vcf_path) if assembly=="auto" else assembly
        if father_path:
            rx.parse_vcf(father_path)
        if mother_path:
            rx.parse_vcf(mother_path)
        return resolved

    try:
        # Parsing a large VCF is blocking + CPU-bound; run it off the event loop so
        # status polling, /api/jobs and static assets stay responsive meanwhile.
        assembly = await asyncio.to_thread(_validate_and_detect)
    except (OSError,rx.VCFParseError) as e:
        shutil.rmtree(outdir,ignore_errors=True)
        raise HTTPException(400,str(e))
    if not JOB_SLOTS.acquire(blocking=False):
        shutil.rmtree(outdir,ignore_errors=True)
        raise HTTPException(503,"server is at analysis capacity; retry later")

    layers=["Ensembl VEP","ClinVar","gnomAD"]
    if hpo.strip():
        layers.append("Open Targets / HPO")
    if alphamissense.lower()=="true":
        layers.append("AlphaMissense")
    if esm.lower()=="true":
        layers.append("ESM-2")
    if agentic.lower()=="true":
        layers.append("Razonamiento agéntico")
    if father_path or mother_path:
        layers.append("Herencia del trío")
    now=time.time()
    cancel_event=threading.Event()
    _set(job_id,status="queued",sample=sample,assembly=assembly,started=now,updated=now,
         layers=layers,cancel_event=cancel_event)
    try:
        EXECUTOR.submit(
            _run_job,job_id,str(vcf_path),sample,hpo,clinical_note,assembly,
            esm.lower()=="true",alphamissense.lower()=="true",
            agentic.lower()=="true",reflect_k_i,father_path,mother_path,
            expanded_hpo,cancel_event,
        )
    except Exception:
        JOB_SLOTS.release()
        with JOBS_LOCK:
            JOBS.pop(job_id,None)
        shutil.rmtree(outdir,ignore_errors=True)
        raise
    return {"job_id":job_id,"assembly":assembly}


@app.post("/api/cancel/{job_id}")
def cancel(job_id: str):
    with JOBS_LOCK:
        job=JOBS.get(job_id)
        if not job:
            raise HTTPException(404,"unknown job")
        if job.get("status") not in ("queued","running","cancelling"):
            raise HTTPException(409,"job is not running")
        event=job.get("cancel_event")
        if event is None:
            raise HTTPException(409,"job cannot be cancelled")
        event.set()
        job.update(status="cancelling",message="Deteniendo análisis…",updated=time.time())
    return {"status":"cancelling"}


@app.get("/api/status/{job_id}")
def status(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    if not j:
        raise HTTPException(404, "unknown job")
    now=time.time()
    run_started=j.get("run_started") or j.get("started") or now
    ended=j.get("finished") if j.get("status") in ("done","error","cancelled") else now
    elapsed=max(0.0,ended-run_started)
    done=j.get("done",0)
    total=j.get("total",0)
    rate_per_minute=(done/elapsed*60) if done and elapsed>0 else None
    eta_seconds=(
        max(0.0,(total-done)/(rate_per_minute/60))
        if rate_per_minute and total and done<total else
        0.0 if j.get("status")=="done" else None
    )
    pct = int(100 * j.get("done", 0) / j["total"]) if j.get("total") else 0
    return JSONResponse({
        "status": j.get("status"), "message": j.get("message"), "percent": pct,
        "assembly": j.get("assembly"),
        "done": j.get("done", 0), "total": j.get("total", 0),
        "n_variants": j.get("n_variants"), "top": j.get("top"),
        "differential": j.get("differential"),
        "warnings": j.get("warnings",[]),
        "llm_providers": j.get("llm_providers",[]),
        "layers":j.get("layers",[]),
        "elapsed_seconds":round(elapsed,1),
        "rate_per_minute":round(rate_per_minute,1) if rate_per_minute else None,
        "eta_seconds":round(eta_seconds,1) if eta_seconds is not None else None,
        "last_update_seconds":round(max(0.0,now-j.get("updated",now)),1),
        "has_report": bool(j.get("report")), "has_csv": bool(j.get("csv")),
    })

@app.get("/api/jobs")
def jobs():
    """List persisted patient sessions (newest first) so past results stay recoverable
    across restarts and after the ephemeral job dir is TTL-cleaned."""
    with JOBS_LOCK:
        mem={jid:{"n_variants":j.get("n_variants"),"sample":j.get("sample"),
                  "assembly":j.get("assembly"),"status":j.get("status")}
             for jid,j in JOBS.items()}
    out=[]
    for directory in SESSIONS_DIR.iterdir():
        if not directory.is_dir() or not JOB_ID_RE.match(directory.name):
            continue
        report_path=directory / "result_report.html"
        if not report_path.exists():
            continue
        csv_path=directory / "result_annotated.csv"
        meta={}
        meta_path=directory / "metadata.json"
        if meta_path.exists():
            try:
                meta=json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError,ValueError):
                meta={}
        info=mem.get(directory.name,{})
        n_variants=meta.get("n_variants")
        if n_variants is None:
            n_variants=info.get("n_variants")
        if n_variants is None and csv_path.exists():
            try:
                import csv as _csv
                with open(csv_path,encoding="utf-8",newline="") as fh:
                    n_variants=max(0,sum(1 for _ in _csv.reader(fh))-1)
            except OSError:
                n_variants=None
        out.append({
            "job_id":directory.name,
            "modified":meta.get("created_ts") or report_path.stat().st_mtime,
            "created":meta.get("created"),
            "n_variants":n_variants,
            "sample":meta.get("sample") or info.get("sample"),
            "assembly":meta.get("assembly") or info.get("assembly"),
            "hpo":meta.get("hpo"),
            "expanded_hpo":meta.get("expanded_hpo"),
            "top_gene":meta.get("top_gene"),
            "top_call":meta.get("top_call"),
            "ai_layers":meta.get("ai_layers",[]),
            "n_warnings":meta.get("n_warnings"),
            "status":info.get("status") or "done",
            "has_csv":csv_path.exists(),
        })
    out.sort(key=lambda item:item["modified"],reverse=True)
    return {"jobs":out[:100]}


@app.get("/api/report/{job_id}")
def report(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    path=_job_file(job_id,(j or {}).get("report"),"result_report.html")
    if not path:
        raise HTTPException(404, "report not ready")
    return FileResponse(path, media_type="text/html")

@app.get("/api/csv/{job_id}")
def csv(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    path=_job_file(job_id,(j or {}).get("csv"),"result_annotated.csv")
    if not path:
        raise HTTPException(404, "csv not ready")
    return FileResponse(path, media_type="text/csv",
                        filename=f"raredx_{job_id}_candidates.csv")

# serve the static SPA (index.html + assets) at /
_static = Path(__file__).resolve().parent / "static"
app.mount("/", StaticFiles(directory=str(_static), html=True), name="static")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("PORT", "8000")))
