"""Smoke tests for the raredx web backend (FastAPI TestClient — no socket bind needed).

These hit live REST APIs (Ensembl/gnomAD/ClinVar), so they need network access and take
~10-30 s. Run from the repo root:  python -m pytest web/test_server.py -v
"""
import os, sys, time
from pathlib import Path

os.environ.setdefault("RAREDX_DATA_DIR", "/tmp")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest
from fastapi.testclient import TestClient
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


def test_empty_upload_rejected():
    r = client.post("/api/analyze", files={"vcf": ("empty.vcf", b"", "text/plain")},
                    data={"assembly": "GRCh38"})
    assert r.status_code == 400


def test_bad_assembly_rejected():
    r = client.post("/api/analyze", files={"vcf": ("x.vcf", MINI_VCF.encode(), "text/plain")},
                    data={"assembly": "hg99"})
    assert r.status_code == 400


def _run(**data):
    r = client.post("/api/analyze", files={"vcf": ("mini37.vcf", MINI_VCF.encode(), "text/plain")},
                    data={"sample": "16-DR636", "assembly": "GRCh37", **data})
    assert r.status_code == 200
    job = r.json()["job_id"]
    for _ in range(180):
        s = client.get(f"/api/status/{job}").json()
        if s["status"] in ("done", "error"):
            break
        time.sleep(1)
    return job, s


@pytest.mark.network
def test_full_analysis_alphamissense():
    job, s = _run(alphamissense="true", hpo="HP:0001250")
    assert s["status"] == "done", s.get("message")
    assert s["n_variants"] == 2
    top = {v["gene"]: v for v in s["top"]}
    # SCN1A missense should be ranked #1 and flagged pathogenic by AlphaMissense
    assert top["SCN1A"]["rank"] == 1
    assert float(top["SCN1A"]["am_pathogenicity"]) >= 0.9
    assert "pathogenic" in (top["SCN1A"]["call"] or "").lower()
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
