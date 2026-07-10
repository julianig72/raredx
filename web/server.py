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
import os, sys, uuid, threading, shutil, tempfile, time, traceback
from pathlib import Path

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

# import the pipeline (parent dir on path)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raredx_pipeline as rx

DATA_DIR = Path(os.environ.get("RAREDX_DATA_DIR", tempfile.gettempdir())) / "raredx_jobs"
DATA_DIR.mkdir(parents=True, exist_ok=True)
MAX_MB = int(os.environ.get("RAREDX_MAX_MB", "50"))
CONTACT_EMAIL = os.environ.get("RAREDX_EMAIL")  # optional NCBI contact

app = FastAPI(title="raredx", description="AI-assisted rare-disease VCF analysis (decision support)")

# in-memory job registry: job_id -> {status, done, total, message, error, prefix, sample, started}
JOBS = {}
JOBS_LOCK = threading.Lock()

def _set(job_id, **kw):
    with JOBS_LOCK:
        JOBS.setdefault(job_id, {}).update(kw)

def _run_job(job_id, vcf_path, sample, hpo, note, assembly, use_esm, use_am, use_agentic):
    outdir = DATA_DIR / job_id
    prefix = str(outdir / "result")
    try:
        def progress(done, total, message):
            _set(job_id, done=done, total=total, message=message)
        _set(job_id, status="running", done=0, total=0, message="Iniciando…")
        result = rx.run_pipeline(
            vcf_path, sample=sample, hpo=hpo, clinical_note_text=note or None,
            assembly=assembly, use_esm=use_esm, use_am=use_am, agentic=use_agentic,
            email=CONTACT_EMAIL, progress=progress,
        )
        rx.write_outputs(result, prefix)
        variants = result["variants"]
        top = [
            {k: v.get(k) for k in ("rank","gene","consequence","protein","af","clinvar",
                                   "call","combined","am_pathogenicity","esm2_llr","acmg_tags")}
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
             report=f"{prefix}_report.html", csv=f"{prefix}_annotated.csv")
    except Exception as e:
        _set(job_id, status="error", message=str(e), error=traceback.format_exc())

@app.post("/api/analyze")
async def analyze(
    vcf: UploadFile = File(...),
    sample: str = Form("SAMPLE"),
    hpo: str = Form(""),
    clinical_note: str = Form(""),
    assembly: str = Form("GRCh38"),
    esm: str = Form("false"),
    alphamissense: str = Form("false"),
    agentic: str = Form("false"),
):
    if assembly not in ("GRCh38", "GRCh37"):
        raise HTTPException(400, "assembly must be GRCh38 or GRCh37")
    job_id = uuid.uuid4().hex[:12]
    outdir = DATA_DIR / job_id
    outdir.mkdir(parents=True, exist_ok=True)
    vcf_path = outdir / "input.vcf"

    # stream upload to disk with a size cap
    size = 0
    with open(vcf_path, "wb") as fh:
        while chunk := await vcf.read(1 << 20):
            size += len(chunk)
            if size > MAX_MB * (1 << 20):
                shutil.rmtree(outdir, ignore_errors=True)
                raise HTTPException(413, f"VCF exceeds {MAX_MB} MB limit")
            fh.write(chunk)
    if size == 0:
        shutil.rmtree(outdir, ignore_errors=True)
        raise HTTPException(400, "empty upload")

    _set(job_id, status="queued", sample=sample, started=time.time())
    t = threading.Thread(
        target=_run_job,
        args=(job_id, str(vcf_path), sample, hpo, clinical_note, assembly,
              esm.lower() == "true", alphamissense.lower() == "true",
              agentic.lower() == "true"),
        daemon=True,
    )
    t.start()
    return {"job_id": job_id}

@app.get("/api/status/{job_id}")
def status(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    if not j:
        raise HTTPException(404, "unknown job")
    pct = int(100 * j.get("done", 0) / j["total"]) if j.get("total") else 0
    return JSONResponse({
        "status": j.get("status"), "message": j.get("message"), "percent": pct,
        "done": j.get("done", 0), "total": j.get("total", 0),
        "n_variants": j.get("n_variants"), "top": j.get("top"),
        "differential": j.get("differential"),
        "has_report": bool(j.get("report")), "has_csv": bool(j.get("csv")),
    })

@app.get("/api/report/{job_id}")
def report(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    if not j or not j.get("report") or not os.path.exists(j["report"]):
        raise HTTPException(404, "report not ready")
    return FileResponse(j["report"], media_type="text/html")

@app.get("/api/csv/{job_id}")
def csv(job_id: str):
    with JOBS_LOCK:
        j = JOBS.get(job_id)
    if not j or not j.get("csv") or not os.path.exists(j["csv"]):
        raise HTTPException(404, "csv not ready")
    return FileResponse(j["csv"], media_type="text/csv",
                        filename=f"raredx_{job_id}_candidates.csv")

# serve the static SPA (index.html + assets) at /
_static = Path(__file__).resolve().parent / "static"
app.mount("/", StaticFiles(directory=str(_static), html=True), name="static")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("PORT", "8000")))
