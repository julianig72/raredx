"""Smoke tests for the raredx web backend (FastAPI TestClient — no socket bind needed).

These hit live REST APIs (Ensembl/gnomAD/ClinVar), so they need network access and take
~10-30 s. Run from the repo root:  python -m pytest web/test_server.py -v
"""
import os, sys, threading, time
from pathlib import Path

os.environ.setdefault("RAREDX_DATA_DIR", "/tmp")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest
from fastapi.testclient import TestClient
import web.server as server
from web.server import app

client = TestClient(app)

MINI_VCF = (
    "##fileformat=VCFv4.2\n##reference=GRCh37\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t16-DR636\n"
    "chr2\t166859140\t.\tA\tG\t500\tPASS\t.\tGT\t0/1\n"       # SCN1A C1376R (pathogenic missense)
    "chr19\t13318785\t.\tT\tG\t500\tPASS\t.\tGT\t0/1\n"       # CACNA1A Q2288P (benign missense)
)


def test_static_index_served():
    r = client.get("/")
    assert r.status_code == 200
    assert "Analizar VCF" in r.text  # the upload UI
    assert "Extraer HPO con Copilot" in r.text
    assert "Ensamblaje del genoma" in r.text
    assert "Expandir y revisar HPO" in r.text
    assert "Detener análisis" in r.text


def test_empty_upload_rejected():
    r = client.post("/api/analyze", files={"vcf": ("empty.vcf", b"", "text/plain")},
                    data={"assembly": "GRCh38"})
    assert r.status_code == 400


def test_bad_assembly_rejected():
    r = client.post("/api/analyze", files={"vcf": ("x.vcf", MINI_VCF.encode(), "text/plain")},
                    data={"assembly": "hg99"})
    assert r.status_code == 400


def test_invalid_boolean_rejected():
    r = client.post(
        "/api/analyze",
        files={"vcf":("x.vcf",MINI_VCF.encode(),"text/plain")},
        data={"assembly":"GRCh37","alphamissense":"yes"},
    )
    assert r.status_code == 400


def test_oversized_clinical_note_rejected(monkeypatch):
    monkeypatch.setattr(server,"MAX_NOTE_CHARS",10)
    r = client.post(
        "/api/analyze",
        files={"vcf":("x.vcf",MINI_VCF.encode(),"text/plain")},
        data={"assembly":"GRCh37","clinical_note":"x"*11},
    )
    assert r.status_code == 413


def test_malformed_vcf_rejected_before_job_creation():
    r = client.post(
        "/api/analyze",
        files={"vcf": ("bad.vcf", b"not a VCF\n", "text/plain")},
        data={"assembly": "GRCh38"},
    )
    assert r.status_code == 400


def test_extract_hpo_endpoint_returns_reviewable_terms(monkeypatch):
    monkeypatch.setattr(
        "web.server.rx.extract_hpo_from_note",
        lambda note: [
            {
                "hpo_id": "HP:0001250",
                "label": "Seizure",
                "note_evidence": "convulsiones",
            }
        ],
    )
    monkeypatch.setattr(
        "web.server.rx.llm_diagnostics",
        lambda: {"providers": ["GitHub Copilot"], "errors": []},
    )

    r = client.post("/api/extract-hpo", data={"clinical_note": "Presenta convulsiones"})

    assert r.status_code == 200
    assert r.json()["terms"][0]["hpo_id"] == "HP:0001250"
    assert r.json()["llm_providers"] == ["GitHub Copilot"]


def test_expand_hpo_endpoint_returns_reviewable_ancestors(monkeypatch):
    monkeypatch.setattr(
        server.rx,
        "expand_hpo_profile",
        lambda tokens,return_status=False: ([
            {"hpo_id":"HP:0001250","label":"Seizure","kind":"direct",
             "source_hpo_ids":["HP:0001250"]},
            {"hpo_id":"HP:0000707","label":"Abnormality of the nervous system",
             "kind":"ancestor","source_hpo_ids":["HP:0001250"]},
        ],True),
    )

    response=client.post("/api/expand-hpo",data={"hpo":"HP:0001250"})

    assert response.status_code == 200
    assert [term["kind"] for term in response.json()["terms"]] == ["direct","ancestor"]


def test_status_returns_progress_telemetry():
    job_id="telemetry-test"
    now=time.time()
    with server.JOBS_LOCK:
        server.JOBS[job_id]={
            "status":"running","message":"Anotando variantes: 60/120",
            "done":60,"total":120,"started":now-120,"run_started":now-120,
            "updated":now-2,"layers":["Ensembl VEP","ClinVar"],
        }
    try:
        payload=client.get(f"/api/status/{job_id}").json()
    finally:
        with server.JOBS_LOCK:
            server.JOBS.pop(job_id,None)

    assert payload["percent"] == 50
    assert 29 <= payload["rate_per_minute"] <= 31
    assert 119 <= payload["eta_seconds"] <= 121
    assert payload["layers"] == ["Ensembl VEP","ClinVar"]
    assert payload["last_update_seconds"] >= 2


def test_cancel_endpoint_signals_running_job():
    job_id="cancel-test"
    event=threading.Event()
    with server.JOBS_LOCK:
        server.JOBS[job_id]={
            "status":"running","message":"running","started":time.time(),
            "updated":time.time(),"cancel_event":event,
        }
    try:
        response=client.post(f"/api/cancel/{job_id}")
        with server.JOBS_LOCK:
            status=server.JOBS[job_id]["status"]
    finally:
        with server.JOBS_LOCK:
            server.JOBS.pop(job_id,None)

    assert response.status_code == 200
    assert event.is_set()
    assert status == "cancelling"


def test_auto_assembly_is_returned_before_job_starts(monkeypatch):
    class NoopExecutor:
        def submit(self, *args, **kwargs):
            return None

    monkeypatch.setattr(server, "EXECUTOR", NoopExecutor())
    r = client.post(
        "/api/analyze",
        files={"vcf": ("mini37.vcf", MINI_VCF.encode(), "text/plain")},
        data={"sample": "16-DR636", "assembly": "auto"},
    )

    assert r.status_code == 200
    assert r.json()["assembly"] == "GRCh37"
    job = r.json()["job_id"]
    with server.JOBS_LOCK:
        server.JOBS.pop(job, None)
    server.JOB_SLOTS.release()
    import shutil
    shutil.rmtree(server.DATA_DIR / job, ignore_errors=True)


def test_oversized_upload_is_removed_after_file_close(monkeypatch):
    before = {p.name for p in server.DATA_DIR.iterdir()}
    monkeypatch.setattr(server, "MAX_MB", 0)

    r = client.post(
        "/api/analyze",
        files={"vcf": ("too-large.vcf", MINI_VCF.encode(), "text/plain")},
        data={"assembly": "auto"},
    )

    assert r.status_code == 413
    assert {p.name for p in server.DATA_DIR.iterdir()} == before


def _run(**data):
    r = client.post("/api/analyze", files={"vcf": ("mini37.vcf", MINI_VCF.encode(), "text/plain")},
                    data={"sample": "16-DR636", "assembly": "GRCh37", **data})
    assert r.status_code == 200
    job = r.json()["job_id"]
    for _ in range(180):
        s = client.get(f"/api/status/{job}").json()
        if s["status"] in ("done", "error", "cancelled"):
            break
        time.sleep(1)
    return job, s


@pytest.mark.network
def test_full_analysis_alphamissense():
    job, s = _run(alphamissense="true", hpo="HP:0001250")
    assert s["status"] == "done", s.get("message")
    assert s["n_variants"] == 2
    top = {v["gene"]: v for v in s["top"]}
    # SCN1A missense should be ranked #1 with its precomputed AlphaMissense score.
    # AlphaMissense is supporting computational evidence, not a clinical verdict by itself.
    assert top["SCN1A"]["rank"] == 1
    assert float(top["SCN1A"]["am_pathogenicity"]) >= 0.9
    assert "benign" not in (top["SCN1A"]["call"] or "").lower()
    # report + csv retrievable
    assert client.get(f"/api/report/{job}").status_code == 200
    assert client.get(f"/api/csv/{job}").status_code == 200
    assert "AlphaMissense" in client.get(f"/api/report/{job}").text


if __name__ == "__main__":
    # runnable without pytest: quick smoke
    test_static_index_served(); print("static index: OK")
    test_empty_upload_rejected(); print("empty guard: OK")
    test_bad_assembly_rejected(); print("assembly guard: OK")
    job, s = _run(alphamissense="true", hpo="HP:0001250")
    print("full analysis:", s["status"], "| SCN1A rank", {v["gene"]: v["rank"] for v in s["top"]}.get("SCN1A"))
