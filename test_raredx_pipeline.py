import json
import threading
from pathlib import Path

import pytest

import raredx_pipeline as rx


HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFIRST\tREQUESTED\n"
)


def write_vcf(path, records, header=HEADER):
    path.write_text(header + "".join(records), encoding="utf-8")
    return path


def test_parse_vcf_splits_alts_and_selects_requested_sample(tmp_path):
    path = write_vcf(
        tmp_path / "multi.vcf",
        ["chr1\t10\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/0\t1/2\n"],
    )

    variants = rx.parse_vcf(path, sample="REQUESTED", allow_empty=False)

    assert [v["alt"] for v in variants] == ["G", "T"]
    assert [v["sample"]["GT"] for v in variants] == ["1/0", "0/1"]


def test_parse_vcf_rejects_unknown_sample_in_multi_sample_input(tmp_path):
    path = write_vcf(
        tmp_path / "multi.vcf",
        ["chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"],
    )
    with pytest.raises(rx.VCFParseError, match="not found"):
        rx.parse_vcf(path, sample="MISSING", allow_empty=False)


def test_parse_vcf_requires_sample_in_multi_sample_input(tmp_path):
    path = write_vcf(
        tmp_path / "multi.vcf",
        ["chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"],
    )
    with pytest.raises(rx.VCFParseError, match="requires an explicit sample ID"):
        rx.parse_vcf(path, allow_empty=False)


def test_http_get_retries_transient_gateway_errors(monkeypatch):
    statuses=iter([502,504,200])

    class Response:
        def __init__(self,status):
            self.status_code=status

        def json(self):
            return {"ok":True}

    monkeypatch.setattr(rx.requests,"get",lambda *args,**kwargs: Response(next(statuses)))
    monkeypatch.setattr(rx.time,"sleep",lambda *args: None)

    assert rx._get("https://example.test") == {"ok":True}


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        ("##reference=GRCh38\n", "GRCh38"),
        ("##reference=file:///refs/human_g1k_v37.fasta\n", "GRCh37"),
        ("##reference=file:///refs/hs37d5.fa\n", "GRCh37"),
        ("##contig=<ID=chr1,length=248956422>\n", "GRCh38"),
        ("##contig=<ID=1,length=249250621>\n", "GRCh37"),
        ("##contig=<ID=1,assembly=GRCh37>\n", "GRCh37"),
    ],
)
def test_detect_vcf_assembly_from_metadata(tmp_path, metadata, expected):
    path = tmp_path / "assembly.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + metadata
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "1\t10\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    assert rx.detect_vcf_assembly(path) == expected


def test_detect_vcf_assembly_rejects_ambiguous_metadata(tmp_path):
    path = tmp_path / "assembly.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t10\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    with pytest.raises(rx.VCFParseError, match="could not detect"):
        rx.detect_vcf_assembly(path)


def test_request_cancellation_interrupts_analysis():
    event=threading.Event()
    event.set()
    rx.set_request_cancel_event(event)
    try:
        with pytest.raises(rx.AnalysisCancelled):
            rx._request_timeout()
    finally:
        rx.set_request_cancel_event(None)


def test_expand_hpo_profile_labels_direct_terms_and_ancestors(monkeypatch):
    def fake_get(url,**kwargs):
        if url.endswith("/search"):
            return {"response":{"docs":[{"obo_id":"HP:0001250","label":"Seizure"}]}}
        return {"_embedded":{"terms":[
            {"obo_id":"HP:0000707","label":"Abnormality of the nervous system"},
            {"obo_id":"HP:0000118","label":"Phenotypic abnormality"},
        ]}}

    monkeypatch.setattr(rx,"_get",fake_get)

    terms,available=rx.expand_hpo_profile(["HP:0001250"],return_status=True)

    assert available
    assert terms == [
        {"hpo_id":"HP:0001250","label":"Seizure","kind":"direct",
         "source_hpo_ids":["HP:0001250"]},
        {"hpo_id":"HP:0000707","label":"Abnormality of the nervous system",
         "kind":"ancestor","source_hpo_ids":["HP:0001250"]},
    ]


def test_parse_vcf_rejects_non_vcf_and_empty_variant_set(tmp_path):
    invalid = tmp_path / "invalid.vcf"
    invalid.write_text("not a VCF\n", encoding="utf-8")
    with pytest.raises(rx.VCFParseError):
        rx.parse_vcf(invalid, allow_empty=False)

    empty = tmp_path / "empty.vcf"
    empty.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(rx.VCFParseError, match="no supported variant"):
        rx.parse_vcf(empty, allow_empty=False)


def test_parse_vcf_enforces_variant_limit(tmp_path):
    path = write_vcf(
        tmp_path / "many.vcf",
        [
            "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n",
            "1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\n",
        ],
    )
    with pytest.raises(rx.VCFParseError, match="variant analysis limit"):
        rx.parse_vcf(path, sample="REQUESTED", max_variants=1)


def test_called_only_excludes_reference_and_uncalled_genotypes(tmp_path):
    path = write_vcf(
        tmp_path / "calls.vcf",
        [
            "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\n",
            "1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\n",
            "1\t12\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\n",
        ],
    )

    variants = rx.parse_vcf(path, sample="REQUESTED", called_only=True)

    assert [(v["pos"], v["sample"]["GT"]) for v in variants] == [(12, "0/1")]


def test_trio_requires_explicit_parent_calls_for_de_novo(tmp_path):
    trio_header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    child = write_vcf(
        tmp_path / "child.vcf",
        ["1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"],
        trio_header,
    )
    father = write_vcf(
        tmp_path / "father.vcf",
        ["1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"],
        trio_header,
    )
    mother = write_vcf(
        tmp_path / "mother.vcf",
        ["1\t30\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"],
        trio_header,
    )

    variants = rx.parse_vcf(child)
    rx.trio_inheritance(variants, father, mother)
    assert variants[0]["inheritance"] == "absent_in_parents"

    reference_record = ["1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n"]
    write_vcf(father, reference_record, trio_header)
    write_vcf(mother, reference_record, trio_header)
    variants = rx.parse_vcf(child)
    rx.trio_inheritance(variants, father, mother)
    assert variants[0]["inheritance"] == "de_novo"

    write_vcf(father, ["1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/.\n"], trio_header)
    variants = rx.parse_vcf(child)
    rx.trio_inheritance(variants, father, mother)
    assert variants[0]["inheritance"] == "absent_in_parents"


def test_ba1_dominates_pathogenic_evidence_and_conflicts_stay_vus():
    call, tags = rx.classify(0.10, "stop_gained", 0.99, 0.2, None, 0, None, None)
    assert call == "Benign"
    assert {t[0] for t in tags} >= {"BA1", "PVS1"}

    call, tags = rx.classify(
        None,
        "missense_variant",
        None,
        None,
        "Conflicting classifications of pathogenicity",
        0,
        None,
        None,
    )
    assert call == "Uncertain significance (VUS)"
    assert "ClinVar_conflicting" in {t[0] for t in tags}


def test_unavailable_gnomad_does_not_assign_pm2():
    _, tags = rx.classify(
        None, "missense_variant", None, None, None, 0, None, None, af_available=False
    )
    codes = {t[0] for t in tags}
    assert "gnomAD_unavailable" in codes
    assert "PM2" not in codes


def test_gnomad_graphql_errors_with_null_variant_are_unavailable(monkeypatch):
    monkeypatch.setattr(
        rx,
        "_gql",
        lambda *args, **kwargs: {
            "data": {"variant": None},
            "errors": [{"message": "backend unavailable"}],
        },
    )
    assert rx.gnomad_af(
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
        return_status=True,
    ) == (None, False)


def test_grch37_frequency_is_bound_to_alt_allele(monkeypatch):
    monkeypatch.setattr(
        rx,
        "_get",
        lambda *args, **kwargs: [
            {
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [
                    {
                        "canonical": 1,
                        "gene_symbol": "GENE",
                        "gene_id": "ENSG",
                        "amino_acids": "A/G",
                    }
                ],
                "colocated_variants": [
                    {
                        "frequencies": {
                            "A": {"gnomade": 0.42},
                            "G": {"gnomade": 0.00001},
                        }
                    }
                ],
            }
        ],
    )

    annotation = rx.vep({"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}, "GRCh37")

    assert annotation["gnomad_af"] == 0.00001


def test_clinvar_selects_only_exact_allele(monkeypatch):
    wrong = {
        "variation_set": [
            {
                "canonical_spdi": "NC_000001.11:99:A:T",
                "variation_loc": [
                    {"assembly_name": "GRCh38", "chr": "1", "start": "100"}
                ],
            }
        ],
        "germline_classification": {"description": "Pathogenic", "review_status": "expert panel"},
    }
    exact = {
        "variation_set": [
            {
                "canonical_spdi": "NC_000001.11:99:A:G",
                "variation_loc": [
                    {"assembly_name": "GRCh38", "chr": "1", "start": "100"}
                ],
            }
        ],
        "germline_classification": {
            "description": "Likely benign",
            "review_status": "criteria provided, single submitter",
        },
    }

    def fake_get(url, **kwargs):
        if "esearch" in url:
            return {"esearchresult": {"idlist": ["wrong", "exact"]}}
        return {"result": {"wrong": wrong, "exact": exact}}

    monkeypatch.setattr(rx, "_get", fake_get)
    result = rx.clinvar(
        {"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "rsid": None},
        assembly="GRCh38",
    )

    assert result["significance"] == "Likely benign"


def test_clinvar_matches_normalized_vcf_indel():
    doc = {
        "variation_set": [
            {
                "canonical_spdi": "NC_000001.11:100:T:",
                "variation_loc": [
                    {"assembly_name": "GRCh38", "chr": "1", "start": "101", "stop": "101"}
                ],
            }
        ]
    }
    rec = {"chrom": "1", "pos": 100, "ref": "AT", "alt": "A"}

    assert rx._clinvar_matches_allele(doc, rec, "GRCh38")


def test_computational_evidence_alone_stays_vus_and_ps2_can_reclassify():
    ai_call, _ = rx.classify(
        None,
        "missense_variant",
        None,
        None,
        None,
        0,
        None,
        None,
        extra_tags=[("PP3_AM", "deleterious")],
    )
    assert ai_call == "Uncertain significance (VUS)"

    call, tags = rx.classify(
        None, "missense_variant", None, None, None, 0, None, None
    )
    assert call == "Uncertain significance (VUS)"
    tags.append(("PS2", "de novo"))
    assert rx._call_from_tags(tags) == "Likely pathogenic"


def test_liftover_transforms_reverse_strand_alleles_for_alphamissense(monkeypatch):
    requested=[]

    def fake_get(url,**kwargs):
        requested.append(url)
        if "/map/" in url:
            return {"mappings":[{"mapped":{
                "seq_region_name":"2","start":200,"end":200,"strand":-1
            }}]}
        return [{"allele_string":"T/C","transcript_consequences":[{
            "gene_symbol":"GENE",
            "alphamissense":{"am_pathogenicity":0.9,"am_class":"likely_pathogenic"},
        }]}]

    monkeypatch.setattr(rx,"_get",fake_get)

    result=rx.alphamissense_score("2",100,"A","G",assembly="GRCh37",gene="GENE")

    assert result["am_pathogenicity"] == 0.9
    assert any("/2:200-200/C" in url for url in requested)


def test_alphamissense_rejects_reference_mismatch(monkeypatch):
    monkeypatch.setattr(
        rx,"_get",lambda *args,**kwargs: [{
            "allele_string":"T/G",
            "transcript_consequences":[{
                "alphamissense":{"am_pathogenicity":0.99,"am_class":"likely_pathogenic"}
            }],
        }]
    )

    assert rx.alphamissense_score("1",100,"A","G") is None


def test_absolute_score_does_not_make_single_common_variant_ambiguous():
    call, tags = rx.classify(0.10, "missense_variant", None, None, None, 0, None, None)
    assert rx._variant_evidence_score(tags, call) == 0.0


def test_agentic_prompt_contains_trio_evidence(monkeypatch):
    prompts = []

    def fake_llm(system, user, *args, **kwargs):
        prompts.append(user)
        if len(prompts) == 1:
            return {
                "candidates": [
                    {"rank": 1, "gene": "GENE", "verdict": "support", "reasoning": "x"}
                ],
                "all_refuted": False,
            }
        return {"differential": []}

    monkeypatch.setattr(rx, "_llm_json", fake_llm)
    monkeypatch.setattr(rx, "verify_links", lambda urls: {u: True for u in urls})
    variant = {
        "rank": 1,
        "gene": "GENE",
        "chrom": "1",
        "pos": 1,
        "ref": "A",
        "alt": "G",
        "consequence": "missense_variant",
        "call": "Likely pathogenic",
        "combined": 0.9,
        "inheritance": "de_novo",
        "father_gt": "0/0",
        "mother_gt": "0/0",
        "sample": {"GT": "0/1"},
    }

    rx.agentic_diagnosis([variant], [], sample="S")
    candidate_json = prompts[0].split("Candidates (already ranked):\n", 1)[1].split(
        "\n\nReturn JSON:", 1
    )[0]
    candidate = json.loads(candidate_json)[0]

    assert candidate["inheritance"] == "de_novo"
    assert candidate["father_gt"] == candidate["mother_gt"] == "0/0"


def test_llm_auto_prefers_copilot(monkeypatch):
    monkeypatch.delenv("RAREDX_LLM_PROVIDER", raising=False)
    monkeypatch.setattr(rx, "_copilot_raw", lambda system, user: '{"ok":true}')
    monkeypatch.setattr(
        rx,
        "_anthropic_raw",
        lambda *args, **kwargs: pytest.fail("Anthropic fallback should not be called"),
    )
    rx._reset_llm_state()

    assert rx._llm_raw("system", "user") == '{"ok":true}'
    assert rx._LLM_STATE.providers == {"GitHub Copilot"}


def test_llm_auto_falls_back_to_anthropic(monkeypatch):
    monkeypatch.delenv("RAREDX_LLM_PROVIDER", raising=False)
    monkeypatch.setattr(
        rx, "_copilot_raw", lambda *args: (_ for _ in ()).throw(RuntimeError("offline"))
    )
    monkeypatch.setattr(rx, "_anthropic_raw", lambda *args: '{"fallback":true}')
    rx._reset_llm_state()

    assert rx._llm_raw("system", "user") == '{"fallback":true}'
    assert rx._LLM_STATE.providers == {"Anthropic"}
    assert any("GitHub Copilot" in error for error in rx._LLM_STATE.errors)


def test_hpo_extraction_uses_common_llm_provider(monkeypatch):
    monkeypatch.setattr(
        rx,
        "_llm_json",
        lambda *args, **kwargs: {
            "phenotypes": [
                {"note_evidence": "convulsiones", "hpo_phrase": "Seizure"}
            ]
        },
    )
    monkeypatch.setattr(
        rx,
        "_get",
        lambda *args, **kwargs: {
            "response": {"docs": [{"obo_id": "HP:0001250", "label": "Seizure"}]}
        },
    )
    monkeypatch.setattr(rx,"_verify_hpo_candidates",lambda note,candidates,api_key=None: candidates)

    assert rx.extract_hpo_from_note("Paciente con convulsiones") == [
        {
            "hpo_id": "HP:0001250",
            "label": "Seizure",
            "note_evidence": "convulsiones",
        }
    ]


def test_hpo_extraction_ignores_invalid_items_and_unverified_quotes(monkeypatch):
    monkeypatch.setattr(
        rx,
        "_llm_json",
        lambda *args,**kwargs: {
            "phenotypes":[
                "invalid",
                {"note_evidence":"hallucinated quote","hpo_phrase":"Seizure"},
            ]
        },
    )
    monkeypatch.setattr(
        rx,
        "_get",
        lambda *args,**kwargs: {
            "response":{"docs":[{"obo_id":"HP:0001250","label":"Seizure"}]}
        },
    )
    monkeypatch.setattr(rx,"_verify_hpo_candidates",lambda note,candidates,api_key=None: candidates)

    assert rx.extract_hpo_from_note("Paciente con convulsiones") == [{
        "hpo_id":"HP:0001250","label":"Seizure","note_evidence":""
    }]


def test_hpo_extraction_broadens_unsupported_seizure_subtype(monkeypatch):
    responses=iter([
        {
            "phenotypes":[{
                "note_evidence":"crisis epilepticas",
                "hpo_phrase":"Focal-onset seizure evolving into bilateral convulsive status epilepticus",
            }]
        },
        {
            "decisions":[{
                "candidate_hpo_id":"HP:0032662",
                "selected_hpo_id":"HP:0001250",
                "reason":"The note does not state the subtype qualifiers",
            }]
        },
    ])
    monkeypatch.setattr(
        rx,
        "_llm_json",
        lambda *args,**kwargs: next(responses),
    )

    monkeypatch.setattr(
        rx,
        "_get",
        lambda *args,**kwargs: {"response":{"docs":[{
            "obo_id":"HP:0032662",
            "label":"Focal-onset seizure evolving into bilateral convulsive status epilepticus",
        }]}},
    )
    monkeypatch.setattr(
        rx,
        "_hpo_ancestor_options",
        lambda candidate: ([
            {"hpo_id":"HP:0032662","label":candidate["label"]},
            {"hpo_id":"HP:0001250","label":"Seizure"},
        ],True),
    )

    terms=rx.extract_hpo_from_note(
        "varon de 10 meses de edad con crisis epilepticas y farmacoresistente"
    )

    assert terms == [{
        "hpo_id":"HP:0001250",
        "label":"Seizure",
        "note_evidence":"crisis epilepticas",
    }]


def test_hpo_extraction_keeps_explicit_focal_specificity(monkeypatch):
    responses=iter([
        {
            "phenotypes":[{
                "note_evidence":"crisis epilepticas focales",
                "hpo_phrase":"Focal-onset seizure",
            }]
        },
        {
            "decisions":[{
                "candidate_hpo_id":"HP:0007359",
                "selected_hpo_id":"HP:0007359",
                "reason":"Focal is explicit",
            }]
        },
    ])
    monkeypatch.setattr(
        rx,
        "_llm_json",
        lambda *args,**kwargs: next(responses),
    )
    monkeypatch.setattr(
        rx,
        "_get",
        lambda *args,**kwargs: {
            "response":{"docs":[{"obo_id":"HP:0007359","label":"Focal-onset seizure"}]}
        },
    )
    monkeypatch.setattr(
        rx,
        "_hpo_ancestor_options",
        lambda candidate: ([{"hpo_id":"HP:0007359","label":"Focal-onset seizure"},
                            {"hpo_id":"HP:0001250","label":"Seizure"}],True),
    )

    assert rx.extract_hpo_from_note("Presenta crisis epilepticas focales")[0]["hpo_id"] == "HP:0007359"


def test_hpo_verifier_generalizes_beyond_seizures(monkeypatch):
    monkeypatch.setattr(
        rx,
        "_llm_json",
        lambda *args,**kwargs: {
            "decisions":[{
                "candidate_hpo_id":"HP:0008619",
                "selected_hpo_id":"HP:0000365",
                "reason":"The note does not state bilateral involvement",
            }]
        },
    )
    monkeypatch.setattr(
        rx,
        "_hpo_ancestor_options",
        lambda candidate: ([
            {"hpo_id":"HP:0008619","label":"Bilateral sensorineural hearing impairment"},
            {"hpo_id":"HP:0000365","label":"Hearing impairment"},
        ],True),
    )
    candidate={
        "hpo_id":"HP:0008619",
        "label":"Bilateral sensorineural hearing impairment",
        "note_evidence":"hipoacusia",
    }

    assert rx._verify_hpo_candidates("Paciente con hipoacusia",[candidate]) == [{
        "hpo_id":"HP:0000365",
        "label":"Hearing impairment",
        "note_evidence":"hipoacusia",
    }]


def test_hpo_verifier_fails_closed_on_invalid_response(monkeypatch):
    monkeypatch.setattr(rx,"_llm_json",lambda *args,**kwargs: {"unexpected":[]})
    monkeypatch.setattr(
        rx,"_hpo_ancestor_options",
        lambda candidate: ([{"hpo_id":candidate["hpo_id"],"label":candidate["label"]}],True),
    )
    rx._reset_llm_state()

    assert rx._verify_hpo_candidates("note",[{
        "hpo_id":"HP:1","label":"Unsupported","note_evidence":""
    }]) == []
    assert any("HPO verification" in error for error in rx._LLM_STATE.errors)


def test_agentic_payload_sanitizers_reject_malformed_shapes():
    assert rx._sanitize_reflection({"candidates":["bad"]}) is None
    assert rx._sanitize_differential({"differential":[{"disease":"x","genes":"GENE"}]},set()) == [{
        "disease":"x",
        "genes":[],
        "inheritance":"",
        "supporting_variants":[],
        "likelihood":"low",
        "rationale":"",
        "evidence":[],
        "next_steps":"",
    }]


def test_agentic_verdicts_reorder_variants_and_update_reflection_ranks():
    variants = [
        {"rank": 1, "gene": "REFUTED", "combined": 0.9, "reflect_verdict": "refute"},
        {"rank": 2, "gene": "SUPPORTED", "combined": 0.7, "reflect_verdict": "support"},
        {"rank": 3, "gene": "UNEVALUATED", "combined": 0.8},
    ]
    result = {
        "reflection": [
            {"rank": 1, "gene": "REFUTED"},
            {"rank": 2, "gene": "SUPPORTED"},
        ]
    }

    rx._rerank_agentic(variants, result)

    assert [v["gene"] for v in variants] == ["SUPPORTED", "UNEVALUATED", "REFUTED"]
    assert [v["rank"] for v in variants] == [1, 2, 3]
    assert result["reflection"][0]["rank"] == 3
    assert result["reflection"][1]["rank"] == 1


def test_html_report_escapes_vcf_controlled_fields(tmp_path):
    variant = {
        "rank": 1,
        "gene": "GENE",
        "chrom": "1<script>alert(1)</script>",
        "pos": 1,
        "ref": "A",
        "alt": "G",
        "af": None,
        "clinvar": None,
        "stars": 0,
        "call": "Uncertain significance (VUS)",
        "variant_score": 0.5,
        "pheno_score": 0.0,
        "combined": 0.5,
        "filter": "PASS<img src=x onerror=alert(2)>",
    }
    prefix = str(tmp_path / "result")

    rx.write_html([variant], [], prefix)
    report = Path(prefix + "_report.html").read_text(encoding="utf-8")

    assert "<script>alert(1)</script>" not in report
    assert "<img src=x onerror=alert(2)>" not in report
    assert "&lt;script&gt;" in report


def test_html_report_distinguishes_unavailable_frequency(tmp_path):
    variant = {
        "rank": 1,
        "gene": "GENE",
        "chrom": "1",
        "pos": 1,
        "ref": "A",
        "alt": "G",
        "af": None,
        "af_status": "unavailable",
        "clinvar": None,
        "stars": 0,
        "call": "Uncertain significance (VUS)",
        "variant_score": 0.5,
        "pheno_score": 0.0,
        "combined": 0.5,
        "filter": "PASS",
    }
    prefix = str(tmp_path / "result")

    rx.write_html([variant], [], prefix)
    report = Path(prefix + "_report.html").read_text(encoding="utf-8")

    assert "no disponible" in report
