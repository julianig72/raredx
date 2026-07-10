#!/usr/bin/env python3
"""
raredx - Modern VCF annotation & phenotype-driven prioritization for rare-disease diagnosis.

End-to-end, self-contained: annotates each variant against LIVE public APIs and
ranks candidates Exomiser-style by combining an ACMG/AMP-inspired variant score
with an HPO phenotype-match score.

Data sources (all public REST, no API key needed):
  - Ensembl VEP     https://rest.ensembl.org/vep/human/...   (consequence, SIFT/PolyPhen)
  - gnomAD          https://gnomad.broadinstitute.org/api    (GraphQL: AF, pLI/LOEUF)
  - ClinVar         NCBI E-utilities esearch/esummary        (clinical significance, stars)
  - Open Targets    https://api.platform.opentargets.org     (gene->disease->HPO phenotypes)
  - HPO/OLS4        https://www.ebi.ac.uk/ols4/api           (resolve & expand patient terms)

Usage:
  # variant annotation + ACMG ranking only:
  python raredx_pipeline.py input.vcf --sample PATIENT_001 --out-prefix out/patient

  # add phenotype-driven prioritization (HPO IDs or free-text, comma-separated or @file):
  python raredx_pipeline.py input.vcf --hpo "HP:0002205,HP:0001738,Bronchiectasis" \
         --out-prefix out/patient

Outputs:  <prefix>_annotated.csv   <prefix>_report.html
Requires: requests  (pip install requests)

NOTE: run respectfully - the script paces requests. NCBI asks for a contact email;
pass --email you@inst.org to be a good citizen (optional).
"""
import argparse, json, sys, time, html, datetime, re
import requests

ENSEMBL="https://rest.ensembl.org"
GNOMAD="https://gnomad.broadinstitute.org/api"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OT="https://api.platform.opentargets.org/api/v4/graphql"
OLS="https://www.ebi.ac.uk/ols4/api"
UA={"User-Agent":"raredx-pipeline","Accept":"application/json"}
LOF_TERMS={"stop_gained","frameshift_variant","splice_acceptor_variant",
           "splice_donor_variant","start_lost","stop_lost","transcript_ablation"}

def _get(url, **kw):
    for i in range(4):
        try:
            r=requests.get(url, headers=UA, timeout=30, **kw)
            if r.status_code==200: return r.json()
            if r.status_code in (429,500,503): time.sleep(1.5*(i+1)); continue
            return None
        except requests.RequestException:
            time.sleep(1.0*(i+1))
    return None

def _gql(url, query, variables=None):
    for i in range(4):
        try:
            r=requests.post(url, json={"query":query,"variables":variables or {}}, headers=UA, timeout=30)
            if r.status_code==200: return r.json()
            time.sleep(1.5*(i+1))
        except requests.RequestException:
            time.sleep(1.0*(i+1))
    return None

# ---------- VCF ----------
def parse_vcf(path):
    out=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<8: continue
            rec=dict(chrom=f[0].replace("chr",""),pos=int(f[1]),rsid=f[2] if f[2]!="." else None,
                     ref=f[3],alt=f[4].split(",")[0],qual=f[5],filter=f[6] or "PASS")
            if len(f)>9:
                rec["sample"]=dict(zip(f[8].split(":"),f[9].split(":")))
            out.append(rec)
    return out

# ---------- annotation ----------
def vep(rec, assembly="GRCh38"):
    """Ensembl VEP by region+allele (works without an rsID). Assembly-aware: for GRCh37 it
    uses the GRCh37 REST endpoint and also harvests gnomAD AF from colocated_variants, since
    the GraphQL gnomAD API (gnomad_af) is GRCh38-only and would silently miss GRCh37 coords."""
    base = ENSEMBL if assembly == "GRCh38" else ENSEMBL_GRCH37
    reg=f'{rec["chrom"]}:{rec["pos"]}-{rec["pos"]+len(rec["ref"])-1}'
    url=f'{base}/vep/human/region/{reg}/{rec["alt"]}'
    params={"content-type":"application/json"}
    if assembly=="GRCh37": params.update({"AF_gnomade":1,"AF_gnomadg":1,"sift":1,"polyphen":1})
    j=_get(url, params=params)
    if not j: return {}
    r=j[0]
    tcs=r.get("transcript_consequences") or []
    sift=poly=None; aa=None; gene=None; gid=None
    # canonical first, else most severe
    ranked=sorted(tcs, key=lambda t:(-int(t.get("canonical",0)),))
    for t in ranked:
        gene=gene or t.get("gene_symbol"); gid=gid or t.get("gene_id")
        if t.get("sift_prediction") and not sift: sift=t["sift_prediction"]
        if t.get("polyphen_prediction") and not poly: poly=t["polyphen_prediction"]
        if t.get("amino_acids") and not aa: aa=t["amino_acids"]
    # gnomAD AF from colocated frequencies (populated on the GRCh37 endpoint)
    gaf=None
    for cv in (r.get("colocated_variants") or []):
        fr=cv.get("frequencies")
        if fr:
            for _al,d in fr.items():
                for k in ("gnomade","gnomadg","gnomad"):
                    if d.get(k) is not None and gaf is None: gaf=d[k]
    return dict(most_severe=r.get("most_severe_consequence"),gene=gene,gene_id=gid,
                amino_acids=aa,sift=sift,polyphen=poly,gnomad_af=gaf,
                colocated=[c.get("id") for c in (r.get("colocated_variants") or []) if str(c.get("id","")).startswith("rs")])

GNOMAD_VAR_Q="""query($vid:String!,$ds:DatasetId!){
 variant(variantId:$vid, dataset:$ds){ exome{af homozygote_count} genome{af homozygote_count} } }"""
def gnomad_af(rec):
    vid=f'{rec["chrom"]}-{rec["pos"]}-{rec["ref"]}-{rec["alt"]}'
    j=_gql(GNOMAD, GNOMAD_VAR_Q, {"vid":vid,"ds":"gnomad_r4"})
    v=((j or {}).get("data") or {}).get("variant")
    if not v: return None
    afs=[x.get("af") for x in (v.get("exome"),v.get("genome")) if x and x.get("af") is not None]
    return max(afs) if afs else None

GNOMAD_CONS_Q="""query($sym:String!){ gene(gene_symbol:$sym, reference_genome:GRCh38){
  gnomad_constraint{ pli oe_lof_upper } } }"""
def gnomad_constraint(sym, cache):
    if sym in cache: return cache[sym]
    j=_gql(GNOMAD, GNOMAD_CONS_Q, {"sym":sym})
    c=(((j or {}).get("data") or {}).get("gene") or {}).get("gnomad_constraint") or {}
    cache[sym]={"pli":c.get("pli"),"loeuf":c.get("oe_lof_upper")}
    return cache[sym]

def clinvar(rec, email=None):
    """ClinVar via E-utilities: esearch by rsid or chrom/pos, esummary for significance+stars."""
    params={"db":"clinvar","retmode":"json"}
    if email: params["email"]=email
    term=f'{rec["rsid"]}' if rec.get("rsid") else f'{rec["chrom"]}[chr] AND {rec["pos"]}[chrpos37] OR {rec["pos"]}[chrpos38]'
    es=_get(f"{EUTILS}/esearch.fcgi", params={**params,"term":term})
    ids=(((es or {}).get("esearchresult") or {}).get("idlist")) or []
    if not ids: return {"significance":None,"stars":0}
    su=_get(f"{EUTILS}/esummary.fcgi", params={**params,"id":ids[0]})
    doc=(((su or {}).get("result") or {}).get(ids[0])) or {}
    germ=doc.get("germline_classification") or {}
    desc=germ.get("description")
    rs=(germ.get("review_status") or "").lower()
    stars=(4 if "practice guideline" in rs else 3 if "expert panel" in rs
           else 2 if "multiple" in rs and "conflict" not in rs else 1 if "single" in rs else 0)
    conds=[t.get("trait_name") for t in (germ.get("trait_set") or []) if t.get("trait_name")]
    return {"significance":desc,"stars":stars,"conditions":conds[:3],
            "protein_change":doc.get("protein_change")}

# ---------- classification ----------
def classify(af, cons, pli, loeuf, sig, stars, sift, poly):
    tags=[]
    if af is not None and af>=0.05: tags.append(("BA1",f"gnomAD AF {af:.3f} >=5%"))
    elif af is not None and af>=0.01: tags.append(("BS1",f"gnomAD AF {af:.3f} >=1%"))
    if af is None: tags.append(("PM2","absent from gnomAD r4"))
    elif af<1e-4: tags.append(("PM2",f"gnomAD AF {af:.2e} <0.01%"))
    lof_intol=(pli is not None and pli>=0.9) or (loeuf is not None and loeuf<0.6)
    if cons in LOF_TERMS:
        tags.append(("PVS1" if lof_intol else "PVS1_mod",
                     f"{cons}"+(" in LoF-intolerant gene" if lof_intol else " (constraint modest)")))
    if sift and "deleterious" in sift and poly and "damaging" in poly:
        tags.append(("PP3",f"SIFT {sift} & PolyPhen {poly}"))
    elif sift=="tolerated" and poly and "benign" in poly:
        tags.append(("BP4",f"SIFT {sift} & PolyPhen {poly}"))
    sig=sig or ""
    if "athogenic" in sig: tags.append((("PS_ClinVar" if stars>=2 else "PP5_ClinVar"),f"ClinVar {sig} ({stars} stars)"))
    elif "enign" in sig: tags.append((("BS_ClinVar" if stars>=2 else "BP6_ClinVar"),f"ClinVar {sig} ({stars} stars)"))
    codes={t[0] for t in tags}
    ps={"PVS1","PS_ClinVar"}&codes; psup={"PM2","PP3","PP5_ClinVar","PVS1_mod"}&codes
    bs={"BA1","BS1","BS_ClinVar"}&codes; bsup={"BP4","BP6_ClinVar"}&codes
    if bs and not ps: call="Benign / Likely benign"
    elif ps and "PVS1" in codes and ({"PM2","PS_ClinVar"}&codes): call="Pathogenic"
    elif ps: call="Likely pathogenic"
    elif psup and not bs: call="Likely pathogenic" if len(psup)>=2 else "Uncertain significance (VUS)"
    elif bsup: call="Likely benign"
    else: call="Uncertain significance (VUS)"
    if "drug response" in sig.lower(): call+=" | Pharmacogenomic (drug response)"
    return call, tags

# ---------- phenotype (HPO) ----------
def resolve_hpo(tokens):
    """Accept HPO IDs (HP:xxxxxxx) or free text; resolve text via OLS4."""
    ids=[]
    for tok in tokens:
        tok=tok.strip()
        if re.fullmatch(r"HP:\d+", tok, re.I): ids.append({"hpo_id":tok.upper(),"label":tok}); continue
        j=_get(f"{OLS}/search", params={"q":tok,"ontology":"hp","rows":3})
        docs=(((j or {}).get("response") or {}).get("docs")) or []
        pick=next((d for d in docs if d.get("label","").lower()==tok.lower()), docs[0] if docs else None)
        if pick: ids.append({"hpo_id":pick["obo_id"],"label":pick["label"]})
    return ids

def expand_hpo(hpo_ids):
    full=set(hpo_ids)
    for hid in list(hpo_ids):
        iri=f"http://purl.obolibrary.org/obo/{hid.replace(':','_')}"
        j=_get(f"{OLS}/ontologies/hp/ancestors", params={"id":iri,"size":200})
        for t in (((j or {}).get("_embedded") or {}).get("terms")) or []:
            oid=t.get("obo_id")
            if oid and oid.startswith("HP:"): full.add(oid)
    full-={"HP:0000001","HP:0000118"}
    return full

OT_Q="""query($id:String!){ target(ensemblId:$id){ associatedDiseases(page:{size:6,index:0}){
  rows{ score disease{ id name phenotypes(page:{size:80,index:0}){ rows{ phenotypeHPO{ id name } } } } } } } }"""
def gene_pheno_score(ensg, patient_ids, patient_full):
    if not ensg: return 0.0,0,0,"",""
    j=_gql(OT, OT_Q, {"id":ensg})
    t=(((j or {}).get("data") or {}).get("target")) or {}
    hpo={}; best=None
    for row in ((t.get("associatedDiseases") or {}).get("rows") or []):
        d=row["disease"]
        dh={p["phenotypeHPO"]["id"].replace("_",":"):p["phenotypeHPO"]["name"]
            for p in ((d.get("phenotypes") or {}).get("rows") or []) if p.get("phenotypeHPO")}
        shared=[h for h in dh if h in patient_full]
        directsh=[h for h in dh if h in patient_ids]
        if shared and (best is None or len(directsh)>best[1]):
            best=(d["name"], len(directsh), [dh[h] for h in shared][:5])
        hpo.update(dh)
    direct=set(hpo)&patient_ids; matched=set(hpo)&patient_full
    score=(len(direct)*1.0+(len(matched)-len(direct))*0.4)/max(len(patient_ids),1)
    return round(min(score,1.0),3), len(direct), len(matched), (best[0] if best else ""), ("; ".join(best[2]) if best else "")


# ---------- AI module A: ESM-2 missense pathogenicity (optional) ----------
_ESM_CACHE = {}
def esm_score_missense(seq, pos1, wt_aa, mut_aa, model_name="esm2_t6_8M_UR50D", window=511):
    """Masked-marginal log-likelihood ratio log P(mut)/P(wt) from ESM-2.
    Negative => mutation less likely than WT under the protein LM (deleterious).
    Requires: pip install torch fair-esm  (weights auto-download on first use).
    Returns None if esm/torch unavailable or WT residue mismatches."""
    try:
        import torch, esm
    except ImportError:
        return None
    if not seq or not (1 <= pos1 <= len(seq)) or seq[pos1-1] != wt_aa:
        return None
    if model_name not in _ESM_CACHE:
        model, alphabet = getattr(esm.pretrained, model_name)()
        _ESM_CACHE[model_name] = (model.eval(), alphabet, alphabet.get_batch_converter())
    model, alphabet, bc = _ESM_CACHE[model_name]
    half = window // 2
    start = max(0, pos1-1-half); end = min(len(seq), start+window); start = max(0, end-window)
    sub = seq[start:end]; rel = pos1-1-start
    toks = bc([("v", sub)])[2]; mask_i = rel+1
    toks_m = toks.clone(); toks_m[0, mask_i] = alphabet.mask_idx
    import torch as _t
    with _t.no_grad():
        logits = model(toks_m)["logits"][0, mask_i]
    lp = _t.log_softmax(logits, dim=-1)
    return round(float(lp[alphabet.get_idx(mut_aa)] - lp[alphabet.get_idx(wt_aa)]), 3)

def esm_call(llr, del_thr=-3.0, tol_thr=-0.5):
    if llr is None: return ""
    return "deleterious" if llr <= del_thr else ("tolerated" if llr >= tol_thr else "ambiguous")


# ---------- AI module C: AlphaMissense pathogenicity (precomputed, via Ensembl VEP) ----------
# AlphaMissense (Cheng et al. 2023, Science) is NOT executable — DeepMind released only
# precomputed pathogenicity scores for ~71M human missense variants under a non-commercial
# licence. Ensembl VEP serves them (flag AlphaMissense=1). Scores are GRCh38-only, so a
# GRCh37 variant is lifted over to GRCh38 first.
ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org"

def liftover_37_to_38(chrom, pos):
    """Map a GRCh37 position to GRCh38 via the Ensembl assembly-mapping REST endpoint."""
    c = str(chrom).replace("chr", "")
    j = _get(f"{ENSEMBL_GRCH37}/map/human/GRCh37/{c}:{pos}..{pos}/GRCh38")
    if not j:
        return None
    mp = j.get("mappings") or []
    return mp[0]["mapped"]["start"] if mp else None

def alphamissense_score(chrom, pos, ref, alt, assembly="GRCh38", gene=None):
    """Return {'am_pathogenicity': float, 'am_class': str} for a missense variant, or None.
    Queries Ensembl VEP GRCh38 with AlphaMissense=1. If the input is GRCh37, lifts over first."""
    c = str(chrom).replace("chr", "")
    pos38 = pos if assembly == "GRCh38" else liftover_37_to_38(c, pos)
    if not pos38:
        return None
    j = _get(f"{ENSEMBL}/vep/human/region/{c}:{pos38}-{pos38}/{alt}", params={"AlphaMissense": 1})
    if not j:
        return None
    tcs = j[0].get("transcript_consequences") or []
    # prefer the transcript for the annotated gene, else any transcript carrying an AM score
    withscore = [t for t in tcs if t.get("alphamissense")]
    if not withscore:
        return None
    same = [t for t in withscore if gene and t.get("gene_symbol") == gene]
    t = (same or withscore)[0]
    am = t["alphamissense"]
    return {"am_pathogenicity": am.get("am_pathogenicity"), "am_class": am.get("am_class")}

def am_call(am):
    """Normalize AlphaMissense class to the pipeline's deleterious/tolerated/ambiguous vocabulary."""
    if not am or am.get("am_class") is None:
        return ""
    cls = am["am_class"]
    return {"likely_pathogenic": "deleterious", "likely_benign": "tolerated",
            "ambiguous": "ambiguous"}.get(cls, "ambiguous")

def protein_context(rsid, chrom, pos, ref, alt, gene=None):
    """Get protein sequence + AA position for a missense variant via Ensembl VEP + sequence.
    Always queries by REGION+ALLELE (not rsID): an rsID can carry several alt alleles, and the
    VEP-by-id route may return the amino-acid change for a DIFFERENT allele than the one in the
    VCF. Region+allele guarantees the scored mutant residue matches the patient's actual alt."""
    url = f"{ENSEMBL}/vep/human/region/{chrom}:{pos}-{pos+len(ref)-1}/{alt}"
    j = _get(url, params={"content-type":"application/json"})
    if not j: return None
    tcs = j[0].get("transcript_consequences") or []
    cand = [t for t in tcs if t.get("amino_acids") and "/" in t.get("amino_acids","") and t.get("protein_start")]
    if not cand: return None
    # prefer the transcript for the gene we annotated, then MANE/canonical if flagged, else first
    same_gene = [t for t in cand if gene and t.get("gene_symbol")==gene]
    pool = same_gene or cand
    t = sorted(pool, key=lambda x:(-int(x.get("canonical",0) or 0),))[0]
    wt, mut = t["amino_acids"].split("/")
    pid = t.get("protein_id") or t["transcript_id"]
    s = _get(f"{ENSEMBL}/sequence/id/{pid}", params={"type":"protein","content-type":"application/json"})
    seq = (s or {}).get("seq")
    if not seq: return None
    return {"seq":seq, "pos":t["protein_start"], "wt":wt, "mut":mut, "protein_change":f"{wt}{t['protein_start']}{mut}"}


# ---------- AI module B: extract HPO terms from a free-text clinical note ----------
def extract_hpo_from_note(note_text, api_key=None):
    """Use an LLM to extract present phenotypes from a clinical note, then ground each
    to an official HPO ID via OLS4. Requires ANTHROPIC_API_KEY (env or arg) and the
    `anthropic` package (pip install anthropic). Returns [{hpo_id,label,note_evidence}]."""
    import os
    key = api_key or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        print("[raredx] --clinical-note needs ANTHROPIC_API_KEY (or use --hpo directly)", file=sys.stderr)
        return []
    try:
        import anthropic
    except ImportError:
        print("[raredx] pip install anthropic to use --clinical-note", file=sys.stderr)
        return []
    client = anthropic.Anthropic(api_key=key)
    sys_p = ("You extract phenotypic abnormalities from a clinical note and normalize each to a "
             "concise English HPO-style phrase. Extract ONLY explicitly present abnormal phenotypes "
             "(not negated/absent findings, not family history, not treatments). Return STRICT JSON.")
    user_p = (f"Clinical note:\n{note_text}\n\nReturn JSON: "
              '{"phenotypes":[{"note_evidence":"quote","hpo_phrase":"English phenotype term"}]}')
    msg = client.messages.create(model="claude-sonnet-4-5", max_tokens=1500,
                                 system=sys_p, messages=[{"role":"user","content":user_p}])
    txt = re.sub(r"^```[a-z]*|```$", "", msg.content[0].text.strip(), flags=re.M).strip()
    try:
        phenos = json.loads(txt).get("phenotypes", [])
    except json.JSONDecodeError:
        return []
    # ground each phrase to an HPO ID via OLS4
    seen=set(); out=[]
    for p in phenos:
        phrase = p.get("hpo_phrase","")
        for q in [phrase] + ([] if phrase else []):
            j = _get(f"{OLS}/search", params={"q":q,"ontology":"hp","rows":3})
            docs = (((j or {}).get("response") or {}).get("docs")) or []
            pick = next((d for d in docs if d.get("label","").lower()==q.lower()), docs[0] if docs else None)
            if pick and pick["obo_id"] not in seen:
                seen.add(pick["obo_id"])
                out.append({"hpo_id":pick["obo_id"],"label":pick["label"],"note_evidence":p.get("note_evidence","")})
                break
    return out

# ---------- report ----------
def write_html(variants, patient_hpo, prefix, assembly="GRCh38"):
    TC={"Pathogenic":"#c0392b","Likely pathogenic":"#e67e22","Uncertain significance (VUS)":"#7f8c8d",
        "Likely benign":"#27ae60","Benign":"#2ecc71"}
    tier=lambda c:"Benign" if c.split(" |")[0].startswith("Benign") else c.split(" |")[0]
    col=lambda c:TC.get(tier(c),"#7f8c8d")
    faf=lambda a:"absent" if a is None else (f"{a:.2e}" if a<1e-3 else f"{a:.3f}")
    def ai_cell(v):
        # AlphaMissense pathogenicity (0-1) and/or ESM-2 LLR, whichever were computed
        parts=[]
        am=v.get("am_pathogenicity")
        if am is not None:
            amc="#c0392b" if float(am)>=0.564 else ("#27ae60" if float(am)<=0.34 else "#7f8c8d")
            parts.append(f'<span style="color:{amc};font-weight:600" title="AlphaMissense">AM {float(am):.2f}</span>')
        llr=v.get("esm2_llr")
        if llr is not None:
            ec="#c0392b" if float(llr)<=-3 else ("#27ae60" if float(llr)>=-0.5 else "#7f8c8d")
            parts.append(f'<span style="color:{ec}" title="ESM-2 LLR">ESM {float(llr):.1f}</span>')
        return " ".join(parts) or "—"
    chips="".join(f'<span class="hpo">{html.escape(t["label"])} <code>{t["hpo_id"]}</code></span>' for t in patient_hpo)
    rows="".join(f'<tr><td>{v["rank"]}</td><td><b>{html.escape(v.get("gene") or "")}</b></td>'
                 f'<td>{v["chrom"]}:{v["pos"]} {html.escape(v["ref"])}>{html.escape(v["alt"])}</td>'
                 f'<td>{faf(v.get("af"))}</td><td>{html.escape(str(v.get("clinvar") or ""))} {"*"*int(v.get("stars") or 0)}</td>'
                 f'<td>{ai_cell(v)}</td>'
                 f'<td style="color:{col(v["call"])};font-weight:600">{html.escape(tier(v["call"]))}</td>'
                 f'<td>{v.get("variant_score",0):.2f}</td><td style="color:#8e44ad">{v.get("pheno_score",0):.2f}</td>'
                 f'<td><b>{v.get("combined",0):.2f}</b></td><td>{v.get("filter")}</td></tr>' for v in variants)
    now=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    doc=f"""<!doctype html><html lang="es"><head><meta charset="utf-8"><title>raredx report</title>
<style>body{{font-family:-apple-system,Segoe UI,sans-serif;max-width:1040px;margin:auto;padding:24px;background:#f5f6fa;color:#2c3e50}}
header{{background:#1e2a38;color:#fff;padding:20px 28px;border-radius:10px}} h1{{font-size:20px;margin:0}}
.hpo{{display:inline-block;background:#f3eafc;color:#6c3483;border-radius:14px;padding:3px 11px;font-size:12px;margin:3px}}
.hpo code{{color:#8e44ad;font-size:10px}} .hpobar{{background:#fff;border-radius:10px;padding:12px 16px;margin:14px 0}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:10px;overflow:hidden;font-size:12.5px}}
th{{background:#eef1f5;text-align:left;padding:8px;font-size:10.5px;text-transform:uppercase;color:#7f8c8d}}
td{{padding:7px 8px;border-top:1px solid #eef1f5}} .note{{background:#fff8e1;border:1px solid #ffe082;border-radius:8px;padding:12px 16px;font-size:12px;color:#795548;margin-top:20px}}</style></head><body>
<header><h1>raredx - Informe genomico priorizado por fenotipo</h1>
<div style="color:#9fb3c8;font-size:12px">{assembly} | {len(variants)} variantes | {now}</div></header>
<div class="hpobar"><b>Perfil HPO del paciente ({len(patient_hpo)}):</b><br>{chips or "(ninguno - solo ranking por variante)"}</div>
<table><thead><tr><th>#</th><th>Gen</th><th>Variante</th><th>gnomAD</th><th>ClinVar</th><th>IA</th><th>Clase</th><th>Var</th><th>Feno</th><th>Comb</th><th>QC</th></tr></thead>
<tbody>{rows}</tbody></table>
<div class="note"><b>Sistema de apoyo a la decision, no diagnostico.</b> VEP + gnomAD + ClinVar (score de variante ACMG-lite) x fenotipos HPO de Open Targets (score fenotipico). Confirmar por metodo ortogonal e interpretar en contexto clinico.</div>
</body></html>"""
    open(f"{prefix}_report.html","w").write(doc)

# ---------- main ----------
CSV_COLS=["rank","gene","rsid","chrom","pos","ref","alt","consequence","protein","af",
          "clinvar","stars","pli","loeuf","sift","polyphen","call","variant_score",
          "pheno_score","pheno_direct","pheno_disease","pheno_shared","combined","filter",
          "acmg_tags","esm2_llr","esm2_call","am_pathogenicity","am_call"]

def run_pipeline(vcf_path, sample="SAMPLE", hpo="", clinical_note_text=None, assembly="GRCh38",
                 use_esm=False, use_am=False, email=None, anthropic_key=None, progress=None):
    """Annotate a VCF and prioritize variants. Reusable entry point for CLI and web server.

    Args:
        vcf_path: path to the input VCF.
        sample: sample id (label only).
        hpo: comma/newline-separated HPO IDs or free-text terms.
        clinical_note_text: raw clinical-note text (LLM -> HPO), or None.
        assembly: 'GRCh38' or 'GRCh37' (GRCh37 lifts over to GRCh38 for AlphaMissense).
        use_esm / use_am: enable the ESM-2 / AlphaMissense missense layers.
        email: contact email for NCBI E-utilities.
        anthropic_key: key for the clinical-note LLM layer (else ANTHROPIC_API_KEY env).
        progress: optional callable(done:int, total:int, message:str) for live status.

    Returns dict: {variants, patient_hpo, csv_cols, n_input}. Pass the result to
    write_outputs(result, prefix) to emit <prefix>_annotated.csv and <prefix>_report.html.
    """
    def _prog(done, total, msg):
        if progress:
            try: progress(done, total, msg)
            except Exception: pass

    variants=parse_vcf(vcf_path)
    _prog(0, len(variants), f"VCF leído: {len(variants)} variantes")

    # patient HPO profile — from clinical note (LLM) and/or explicit hpo string
    patient_hpo=[]
    if clinical_note_text:
        patient_hpo=extract_hpo_from_note(clinical_note_text, anthropic_key)
        _prog(0, len(variants), f"{len(patient_hpo)} fenotipos extraídos de la nota clínica")
    tokens=[t for t in re.split(r"[,\n]", hpo or "") if t.strip()]
    if tokens:
        existing={t["hpo_id"] for t in patient_hpo}
        for t in resolve_hpo(tokens):
            if t["hpo_id"] not in existing: patient_hpo.append(t); existing.add(t["hpo_id"])
    patient_ids={t["hpo_id"] for t in patient_hpo}
    patient_full=expand_hpo(patient_ids) if patient_ids else set()
    if patient_hpo: _prog(0, len(variants), f"{len(patient_ids)} términos HPO (expandidos a {len(patient_full)})")

    cons_cache={}
    for i,v in enumerate(variants,1):
        ve=vep(v, assembly); v.update(consequence=ve.get("most_severe"), gene=ve.get("gene"),
                            gene_id=ve.get("gene_id"), protein=ve.get("amino_acids"),
                            sift=ve.get("sift"), polyphen=ve.get("polyphen"))
        v["af"]=gnomad_af(v) if assembly=="GRCh38" else ve.get("gnomad_af")
        con=gnomad_constraint(v["gene"], cons_cache) if v.get("gene") else {"pli":None,"loeuf":None}
        v["pli"],v["loeuf"]=con.get("pli"),con.get("loeuf")
        cv=clinvar(v, email); v.update(clinvar=cv.get("significance"), stars=cv.get("stars"),
                                       conditions=cv.get("conditions"), protein=v.get("protein") or cv.get("protein_change"))
        v["esm2_llr"]=None; v["esm2_call"]=""
        if use_esm and v["consequence"]=="missense_variant":
            ctx=protein_context(v.get("rsid"), v["chrom"], v["pos"], v["ref"], v["alt"], v.get("gene"))
            if ctx:
                v["esm2_llr"]=esm_score_missense(ctx["seq"], ctx["pos"], ctx["wt"], ctx["mut"])
                v["esm2_call"]=esm_call(v["esm2_llr"])
                v["protein"]=v.get("protein") or ctx["protein_change"]
        v["am_pathogenicity"]=None; v["am_call"]=""
        if use_am and v["consequence"]=="missense_variant":
            am=alphamissense_score(v["chrom"], v["pos"], v["ref"], v["alt"], assembly=assembly, gene=v.get("gene"))
            if am:
                v["am_pathogenicity"]=am["am_pathogenicity"]; v["am_call"]=am_call(am)
        call,tags=classify(v["af"],v["consequence"],v["pli"],v["loeuf"],v["clinvar"],v["stars"],v["sift"],v["polyphen"])
        if v["esm2_call"]=="deleterious": tags.append(("PP3_ESM",f"ESM-2 LLR {v['esm2_llr']} (deleterious)"))
        elif v["esm2_call"]=="tolerated": tags.append(("BP4_ESM",f"ESM-2 LLR {v['esm2_llr']} (tolerated)"))
        if v["am_call"]=="deleterious": tags.append(("PP3_AM",f"AlphaMissense {v['am_pathogenicity']} (likely_pathogenic)"))
        elif v["am_call"]=="tolerated": tags.append(("BP4_AM",f"AlphaMissense {v['am_pathogenicity']} (likely_benign)"))
        v["call"]=call; v["acmg_tags"]=",".join(t[0] for t in tags); v["evidence"]="; ".join(f"{t[0]}: {t[1]}" for t in tags)
        codes={t[0] for t in tags}; s=0
        s+=50*len({"PVS1","PS_ClinVar"}&codes)+20*len({"PM2","PP3","PP5_ClinVar","PVS1_mod"}&codes)
        s-=40*len({"BA1","BS_ClinVar"}&codes)+15*len({"BS1","BP4","BP6_ClinVar"}&codes)
        if call.split(" |")[0]=="Pathogenic": s+=30
        if v["filter"]!="PASS": s-=25
        v["_raw"]=s
        if patient_ids:
            ps,pd_,pm_,pdis,psh=gene_pheno_score(v.get("gene_id"),patient_ids,patient_full)
            v.update(pheno_score=ps,pheno_direct=pd_,pheno_matched=pm_,pheno_disease=pdis,pheno_shared=psh)
        else:
            v.update(pheno_score=0.0,pheno_direct=0,pheno_matched=0,pheno_disease="",pheno_shared="")
        if i % 10 == 0 or i == len(variants):
            _prog(i, len(variants), f"Anotando variantes: {i}/{len(variants)}")

    raws=[v["_raw"] for v in variants] or [0]; lo,hi=min(raws),max(raws)
    for v in variants:
        vs=round((v["_raw"]-lo)/(hi-lo),3) if hi>lo else 0.5
        if v.get("esm2_call")=="deleterious": vs=min(vs+0.08,1.0)
        elif v.get("esm2_call")=="tolerated": vs=max(vs-0.08,0.0)
        if v.get("am_call")=="deleterious": vs=min(vs+0.10,1.0)
        elif v.get("am_call")=="tolerated": vs=max(vs-0.10,0.0)
        v["variant_score"]=round(vs,3)
        c=0.55*v["variant_score"]+0.45*v["pheno_score"] if patient_ids else v["variant_score"]
        if v["filter"]!="PASS": c*=0.5
        v["combined"]=round(c,3)
    variants.sort(key=lambda v:-v["combined"])
    for i,v in enumerate(variants,1): v["rank"]=i
    _prog(len(variants), len(variants), "Priorización completa")

    return {"variants":variants, "patient_hpo":patient_hpo, "csv_cols":CSV_COLS,
            "n_input":len(variants), "assembly":assembly}

def write_outputs(result, out_prefix):
    """Write <prefix>_annotated.csv and <prefix>_report.html from a run_pipeline() result."""
    import csv
    with open(f"{out_prefix}_annotated.csv","w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=result["csv_cols"], extrasaction="ignore"); w.writeheader()
        for v in result["variants"]: w.writerow(v)
    write_html(result["variants"], result["patient_hpo"], out_prefix, result.get("assembly","GRCh38"))

def main():
    ap=argparse.ArgumentParser(description="raredx VCF annotation & phenotype prioritization")
    ap.add_argument("vcf")
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--hpo", default="", help="Comma-separated HPO IDs or free-text terms, or @file.txt")
    ap.add_argument("--out-prefix", default="raredx_out")
    ap.add_argument("--email", default=None, help="Contact email for NCBI E-utilities (optional)")
    ap.add_argument("--esm", action="store_true", help="Score missense with ESM-2 (needs torch+fair-esm)")
    ap.add_argument("--alphamissense", action="store_true", help="Score missense with AlphaMissense (precomputed, via Ensembl VEP; no GPU)")
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38","GRCh37"], help="Genome build of the input VCF (AlphaMissense lifts GRCh37->GRCh38)")
    ap.add_argument("--clinical-note", default=None, help="Path to free-text clinical note; extracts HPO via LLM (needs ANTHROPIC_API_KEY)")
    ap.add_argument("--anthropic-key", default=None, help="Anthropic API key (else uses ANTHROPIC_API_KEY env)")
    a=ap.parse_args()

    hpo_arg=a.hpo
    if hpo_arg.startswith("@"): hpo_arg=open(hpo_arg[1:]).read()
    note=open(a.clinical_note).read() if a.clinical_note else None

    def cli_progress(done, total, msg): print(f"[raredx] {msg}", file=sys.stderr)
    result=run_pipeline(a.vcf, sample=a.sample, hpo=hpo_arg, clinical_note_text=note,
                        assembly=a.assembly, use_esm=a.esm, use_am=a.alphamissense,
                        email=a.email, anthropic_key=a.anthropic_key, progress=cli_progress)
    write_outputs(result, a.out_prefix)
    print(f"[raredx] wrote {a.out_prefix}_annotated.csv and {a.out_prefix}_report.html", file=sys.stderr)

if __name__=="__main__":
    main()
