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

# ---------- agentic reasoning layer (inspired by DeepRare, Nature 2025) ----------
LLM_MODEL = "claude-sonnet-4-5"

def _llm_raw(system_p, user_p, api_key=None, max_tokens=2500):
    """Return raw LLM text via the anthropic SDK (ANTHROPIC_API_KEY) or, if unavailable,
    a `host.llm` accessor when running inside a Claude Science kernel. None if neither."""
    import os
    key = api_key or os.environ.get("ANTHROPIC_API_KEY")
    if key:
        try:
            import anthropic
            msg = anthropic.Anthropic(api_key=key).messages.create(
                model=LLM_MODEL, max_tokens=max_tokens, system=system_p,
                messages=[{"role":"user","content":user_p}])
            return msg.content[0].text
        except Exception as e:
            print(f"[raredx] anthropic error: {e}", file=sys.stderr)
    _host = globals().get("host")  # injected in Claude Science kernels
    if _host is not None and hasattr(_host, "llm"):
        try:
            return _host.llm(user_p, system=system_p, max_tokens=max_tokens)["text"]
        except Exception as e:
            print(f"[raredx] host.llm error: {e}", file=sys.stderr)
    print("[raredx] agentic layer needs ANTHROPIC_API_KEY (or a host.llm accessor)", file=sys.stderr)
    return None

def _extract_json(raw):
    """Extract a JSON object from an LLM reply that may wrap it in a ```json fence and/or
    append prose commentary after it. Tries: fenced block -> whole string -> first {...last }."""
    for cand in (
        (re.search(r"```(?:json)?\s*(\{.*?\})\s*```", raw, re.S) or [None, None])[1],
        raw.strip(),
        raw[raw.find("{"): raw.rfind("}") + 1] if "{" in raw and "}" in raw else None,
    ):
        if not cand: continue
        try:
            return json.loads(cand)
        except json.JSONDecodeError:
            continue
    return None

def _llm_json(system_p, user_p, api_key=None, max_tokens=2500):
    """Call the LLM and parse a STRICT-JSON reply. Returns dict/list or None on failure."""
    raw = _llm_raw(system_p, user_p, api_key, max_tokens)
    if raw is None: return None
    obj = _extract_json(raw)
    if obj is None:
        print("[raredx] could not parse JSON from LLM reply", file=sys.stderr)
    return obj

def verify_links(urls, timeout=6):
    """Reference-link verification (DeepRare anti-hallucination step): keep only URLs that
    resolve (HTTP < 400). Returns {url: bool}. Best-effort; network failures count as invalid."""
    out = {}
    for u in dict.fromkeys(urls):
        ok = False
        try:
            r = requests.head(u, timeout=timeout, allow_redirects=True, headers=UA)
            if r.status_code >= 400 or r.status_code == 405:  # some servers reject HEAD
                r = requests.get(u, timeout=timeout, allow_redirects=True, headers=UA, stream=True)
            ok = r.status_code < 400
        except Exception:
            ok = False
        out[u] = ok
    return out

def build_evidence(v):
    """Build STRUCTURED evidence URLs for a variant (deterministic — no LLM-invented links).
    Feeding the LLM a fixed menu of real URLs prevents hallucinated citations entirely."""
    ev = {}
    gene = v.get("gene"); rsid = v.get("rsid"); gid = v.get("gene_id")
    if rsid and str(rsid).startswith("rs"):
        ev["ClinVar (rsID)"] = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}"
        ev["dbSNP"] = f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
    if gene:
        ev["ClinVar (gene)"] = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}%5Bgene%5D"
        ev["OMIM (gene)"] = f"https://www.omim.org/search?search={gene}"
        ev["GeneReviews/PubMed"] = f"https://pubmed.ncbi.nlm.nih.gov/?term={gene}+{v.get('consequence','').replace('_variant','')}"
    if gid:
        ev["Open Targets"] = f"https://platform.opentargets.org/target/{gid}"
        ev["Ensembl gene"] = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gid}"
    return ev

def _zygosity(v):
    gt = ((v.get("sample") or {}).get("GT") or "").replace("|", "/")
    if gt in ("1/1","1"): return "homozygous"
    if gt in ("0/1","1/0","0/2","1/2"): return "heterozygous"
    return gt or "unknown"

def _gt_alleles(sample):
    """Return the set of allele indices in a genotype call, or None if missing/uncalled."""
    gt = (sample or {}).get("GT")
    if not gt or gt in (".", "./.", ".|."): return None
    parts = gt.replace("|", "/").split("/")
    try:
        return {int(p) for p in parts if p != "."}
    except ValueError:
        return None

def _has_alt(sample):
    """True if the genotype carries the alt allele (any index >= 1)."""
    al = _gt_alleles(sample)
    return bool(al and any(a >= 1 for a in al))

def trio_inheritance(variants, father_vcf=None, mother_vcf=None, progress=None):
    """Classify the inheritance of each candidate variant against parental genotypes.

    For each proband variant, look up the same (chrom,pos,ref,alt) in each parent's VCF and
    assign an inheritance mode. A *de novo* call (alt present in the child, absent in both
    parents, with both parents genotyped) adds the ACMG PS2 supporting-strong criterion.

    Returns the number of variants annotated. Sets per-variant keys:
      inheritance      one of: de_novo, paternal, maternal, biparental, homozygous_recessive,
                       absent_in_parents (incomplete), or "" if no parental data
      father_gt/mother_gt   the parental genotype strings (or 'absent' / 'NA')
    """
    if not father_vcf and not mother_vcf:
        return 0
    def _index(path):
        idx = {}
        if not path: return None
        for r in parse_vcf(path):
            idx[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r.get("sample") or {}
        return idx
    fidx = _index(father_vcf); midx = _index(mother_vcf)
    if progress:
        try: progress(0, 0, "Analizando herencia (trío)…")
        except Exception: pass
    n = 0
    for v in variants:
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        fs = fidx.get(key) if fidx is not None else None
        ms = midx.get(key) if midx is not None else None
        f_alt = _has_alt(fs) if fs is not None else None
        m_alt = _has_alt(ms) if ms is not None else None
        v["father_gt"] = (fs or {}).get("GT", "absent" if fidx is not None else "NA")
        v["mother_gt"] = (ms or {}).get("GT", "absent" if midx is not None else "NA")
        child_hom = _zygosity(v) == "homozygous"
        both_typed = fidx is not None and midx is not None
        inh = ""
        if both_typed:
            if not f_alt and not m_alt:
                inh = "de_novo"
            elif f_alt and m_alt:
                inh = "homozygous_recessive" if child_hom else "biparental"
            elif f_alt:
                inh = "paternal"
            elif m_alt:
                inh = "maternal"
        else:  # only one parent available — partial call
            present = (f_alt if fidx is not None else m_alt)
            inh = ("paternal" if fidx is not None else "maternal") if present else "absent_in_parents"
        v["inheritance"] = inh
        n += 1
    return n

def agentic_diagnosis(variants, patient_hpo, api_key=None, sample="", progress=None,
                      start_k=8, relax_step=8, max_rounds=3):
    """DeepRare-style agentic layer over the ranked variant list:
      1) SELF-REFLECTION: an LLM reviews the top-K candidates against the phenotype profile
         and inheritance/zygosity, returning support|uncertain|refute + reasoning per candidate.
         If ALL top-K are refuted, it deepens the candidate window (K += relax_step) and
         re-reflects — the faithful analog of DeepRare increasing search depth N and re-iterating.
      2) DIFFERENTIAL SYNTHESIS: groups the surviving candidates into a disease-level
         differential with traceable reasoning, citing only STRUCTURED (pre-built) evidence URLs.
      3) LINK VERIFICATION: every cited URL is checked to resolve; dead links are dropped.
    Returns {differential, reflection, rounds, k_considered} or None if the LLM is unavailable."""
    def _p(msg):
        if progress:
            try: progress(0, 0, msg)
            except Exception: pass
    if not variants: return None
    hpo_str = ", ".join(f"{t['label']} ({t['hpo_id']})" for t in patient_hpo) or "(no HPO provided)"

    def describe(v):
        return {
            "rank": v.get("rank"), "gene": v.get("gene"), "variant": f"{v.get('chrom')}:{v.get('pos')} {v.get('ref')}>{v.get('alt')}",
            "consequence": v.get("consequence"), "protein": v.get("protein"),
            "gnomad_af": v.get("af"), "clinvar": v.get("clinvar"), "clinvar_stars": v.get("stars"),
            "acmg": v.get("acmg_tags"), "classification": (v.get("call") or "").split(" |")[0],
            "alphamissense": v.get("am_pathogenicity"), "esm2_llr": v.get("esm2_llr"),
            "zygosity": _zygosity(v), "chrom": v.get("chrom"),
            "phenotype_match_disease": v.get("pheno_disease"), "phenotype_shared_terms": v.get("pheno_shared"),
            "combined_score": v.get("combined"),
        }

    k = min(start_k, len(variants)); rounds = 0; reflection = None; supported = []
    reflect_sys = (
        "You are a clinical geneticist reviewing an automatically-prioritized variant list for a "
        "rare-disease case. For EACH candidate, judge whether it plausibly explains the patient's "
        "phenotype, weighing: gene-disease/phenotype fit, variant classification and ACMG evidence, "
        "population frequency (common variants rarely cause rare disease), in-silico predictors "
        "(AlphaMissense/ESM-2), and CRUCIALLY the inheritance pattern vs the observed zygosity "
        "(e.g. a single heterozygous variant in a recessive gene with no second hit is a carrier, "
        "not a cause; an X-linked heterozygous call in a male is suspicious). "
        "Verdict must be one of: support, uncertain, refute. Return STRICT JSON.")
    while rounds < max_rounds:
        rounds += 1
        _p(f"Autorreflexión (ronda {rounds}, {k} candidatos)…")
        cands = [describe(v) for v in variants[:k]]
        user = (f"Patient sample: {sample or 'n/a'}\nPhenotype (HPO): {hpo_str}\n\n"
                f"Candidates (already ranked):\n{json.dumps(cands, ensure_ascii=False, indent=1)}\n\n"
                'Return JSON: {"candidates":[{"rank":int,"gene":str,"verdict":"support|uncertain|refute",'
                '"reasoning":str}],"all_refuted":bool}')
        reflection = _llm_json(reflect_sys, user, api_key)
        if not reflection: return None
        verdicts = reflection.get("candidates", [])
        supported = [c for c in verdicts if c.get("verdict") in ("support","uncertain")]
        if supported or k >= len(variants):
            break
        k = min(k + relax_step, len(variants))  # deepen window and re-reflect

    # attach verdicts back to variants by rank
    vmap = {v.get("rank"): v for v in variants}
    for c in reflection.get("candidates", []):
        if c.get("rank") in vmap:
            vmap[c["rank"]]["reflect_verdict"] = c.get("verdict")
            vmap[c["rank"]]["reflect_reasoning"] = c.get("reasoning")

    # build verified evidence menu for the supported candidates
    _p("Verificando enlaces de evidencia…")
    focus = supported or reflection.get("candidates", [])[:3]
    focus_ranks = [c["rank"] for c in focus if c.get("rank") in vmap]
    evidence = {}   # rank -> {label:url} verified
    all_urls = []
    per_rank_ev = {}
    for rk in focus_ranks:
        ev = build_evidence(vmap[rk]); per_rank_ev[rk] = ev; all_urls += list(ev.values())
    valid = verify_links(all_urls)
    for rk in focus_ranks:
        evidence[rk] = {lab: u for lab, u in per_rank_ev[rk].items() if valid.get(u)}

    # differential synthesis grounded on verified links only
    _p("Sintetizando diagnóstico diferencial…")
    diff_sys = (
        "You are a clinical geneticist writing a differential diagnosis from the SUPPORTED variant "
        "candidates. Group them into disease-level hypotheses (a disease, its gene(s), the supporting "
        "variant(s), the inheritance pattern, and a concise evidence-grounded rationale). Rank the "
        "differential by likelihood. Cite ONLY URLs from the provided verified-evidence menu — never "
        "invent links. This is decision support, not a diagnosis. Return STRICT JSON.")
    focus_desc = [dict(describe(vmap[rk]),
                       verified_evidence=evidence.get(rk, {}),
                       verdict=next((c["verdict"] for c in focus if c.get("rank")==rk), None),
                       reflection=vmap[rk].get("reflect_reasoning"))
                  for rk in focus_ranks]
    duser = (f"Phenotype (HPO): {hpo_str}\n\nSupported candidates with verified evidence links:\n"
             f"{json.dumps(focus_desc, ensure_ascii=False, indent=1)}\n\n"
             'Return JSON: {"differential":[{"disease":str,"genes":[str],"inheritance":str,'
             '"supporting_variants":[str],"likelihood":"high|moderate|low","rationale":str,'
             '"evidence":[{"label":str,"url":str}],"next_steps":str}]}')
    diff = _llm_json(diff_sys, duser, api_key)
    differential = (diff or {}).get("differential", []) if diff else []
    # final guard: strip any URL the LLM emitted that isn't in the verified set
    verified_set = {u for u,ok in valid.items() if ok}
    for d in differential:
        d["evidence"] = [e for e in d.get("evidence", []) if e.get("url") in verified_set]

    return {"differential": differential,
            "reflection": reflection.get("candidates", []),
            "rounds": rounds, "k_considered": k,
            "n_verified_links": len(verified_set)}

# ---------- report ----------
def _differential_html(agentic):
    """Render the DeepRare-style differential + self-reflection block, or '' if absent."""
    if not agentic or not agentic.get("differential"):
        return ""
    LK={"high":"#c0392b","moderate":"#e67e22","low":"#7f8c8d"}
    cards=[]
    for i,d in enumerate(agentic["differential"],1):
        lk=(d.get("likelihood") or "low").lower()
        links="".join(f'<a href="{html.escape(e["url"])}" target="_blank" style="font-size:11px;margin-right:10px;color:#2980b9">{html.escape(e["label"])} ↗</a>'
                      for e in d.get("evidence",[]))
        genes=", ".join(d.get("genes",[])); svar="; ".join(d.get("supporting_variants",[]))
        cards.append(
            f'<div style="background:#fff;border-left:4px solid {LK.get(lk,"#7f8c8d")};border-radius:8px;padding:14px 18px;margin:10px 0">'
            f'<div style="font-size:15px;font-weight:700">{i}. {html.escape(d.get("disease","?"))} '
            f'<span style="font-size:11px;font-weight:600;color:{LK.get(lk,"#7f8c8d")};text-transform:uppercase">· {html.escape(lk)}</span></div>'
            f'<div style="font-size:12px;color:#7f8c8d;margin:3px 0"><b>Gen(es):</b> {html.escape(genes)} &nbsp;·&nbsp; '
            f'<b>Herencia:</b> {html.escape(d.get("inheritance","?"))} &nbsp;·&nbsp; <b>Variante(s):</b> {html.escape(svar)}</div>'
            f'<div style="font-size:12.5px;margin:6px 0">{html.escape(d.get("rationale",""))}</div>'
            f'<div style="font-size:12px;color:#34495e;margin-top:4px"><b>Siguiente paso:</b> {html.escape(d.get("next_steps",""))}</div>'
            f'<div style="margin-top:6px">{links or "<span style=\'font-size:11px;color:#95a5a6\'>(sin enlaces verificados)</span>"}</div></div>')
    ref=agentic.get("reflection",[])
    refuted=[c for c in ref if c.get("verdict")=="refute"]
    ref_html=""
    if refuted:
        items="".join(f'<li style="margin:4px 0"><b>{html.escape(str(c.get("gene") or "?"))}</b> (#{c.get("rank")}): {html.escape(c.get("reasoning",""))}</li>' for c in refuted[:8])
        ref_html=(f'<details style="margin:10px 0"><summary style="cursor:pointer;font-size:12.5px;color:#7f8c8d">'
                  f'Autorreflexión: {len(refuted)} candidato(s) descartado(s) por incoherencia clínica</summary>'
                  f'<ul style="font-size:12px;color:#555">{items}</ul></details>')
    meta=(f'<div style="font-size:11px;color:#95a5a6;margin-top:4px">Capa agéntica: {agentic.get("rounds",1)} ronda(s) de '
          f'autorreflexión · {agentic.get("k_considered","?")} candidatos evaluados · '
          f'{agentic.get("n_verified_links",0)} enlaces de evidencia verificados</div>')
    return (f'<div style="margin:16px 0"><h2 style="font-size:16px;color:#1e2a38;margin:0 0 4px">'
            f'🧬 Diagnóstico diferencial (razonamiento agéntico)</h2>{meta}{"".join(cards)}{ref_html}</div>')

def write_html(variants, patient_hpo, prefix, assembly="GRCh38", agentic=None):
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
    has_trio=any(v.get("inheritance") for v in variants)
    INH_LABEL={"de_novo":"de novo","paternal":"paterna","maternal":"materna","biparental":"biparental",
               "homozygous_recessive":"hom. recesiva","absent_in_parents":"ausente en padres"}
    def inh_cell(v):
        inh=v.get("inheritance") or ""
        if not inh: return "<td>—</td>"
        color="#c0392b" if inh=="de_novo" else ("#16a085" if inh=="homozygous_recessive" else "#34495e")
        wt=";font-weight:700" if inh=="de_novo" else ""
        return f'<td style="color:{color}{wt}">{INH_LABEL.get(inh,inh)}</td>'
    rows="".join(f'<tr><td>{v["rank"]}</td><td><b>{html.escape(v.get("gene") or "")}</b></td>'
                 f'<td>{v["chrom"]}:{v["pos"]} {html.escape(v["ref"])}>{html.escape(v["alt"])}</td>'
                 f'<td>{faf(v.get("af"))}</td><td>{html.escape(str(v.get("clinvar") or ""))} {"*"*int(v.get("stars") or 0)}</td>'
                 f'<td>{ai_cell(v)}</td>'
                 f'<td style="color:{col(v["call"])};font-weight:600">{html.escape(tier(v["call"]))}</td>'
                 + (inh_cell(v) if has_trio else "")
                 + f'<td>{v.get("variant_score",0):.2f}</td><td style="color:#8e44ad">{v.get("pheno_score",0):.2f}</td>'
                 f'<td><b>{v.get("combined",0):.2f}</b></td><td>{v.get("filter")}</td></tr>' for v in variants)
    inh_th="<th>Herencia</th>" if has_trio else ""
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
{_differential_html(agentic)}
<h2 style="font-size:14px;color:#1e2a38;margin:18px 0 4px">Variantes candidatas priorizadas</h2>
<table><thead><tr><th>#</th><th>Gen</th><th>Variante</th><th>gnomAD</th><th>ClinVar</th><th>IA</th><th>Clase</th>{inh_th}<th>Var</th><th>Feno</th><th>Comb</th><th>QC</th></tr></thead>
<tbody>{rows}</tbody></table>
<div class="note"><b>Sistema de apoyo a la decision, no diagnostico.</b> VEP + gnomAD + ClinVar (score de variante ACMG-lite) x fenotipos HPO de Open Targets (score fenotipico). Confirmar por metodo ortogonal e interpretar en contexto clinico.</div>
</body></html>"""
    open(f"{prefix}_report.html","w").write(doc)

# ---------- main ----------
CSV_COLS=["rank","gene","rsid","chrom","pos","ref","alt","consequence","protein","af",
          "clinvar","stars","pli","loeuf","sift","polyphen","call","variant_score",
          "pheno_score","pheno_direct","pheno_disease","pheno_shared","combined","filter",
          "acmg_tags","esm2_llr","esm2_call","am_pathogenicity","am_call",
          "agentic_evaluated","reflect_verdict","inheritance","father_gt","mother_gt"]

def run_pipeline(vcf_path, sample="SAMPLE", hpo="", clinical_note_text=None, assembly="GRCh38",
                 use_esm=False, use_am=False, agentic=False, reflect_k=8,
                 father_vcf=None, mother_vcf=None, email=None,
                 anthropic_key=None, progress=None):
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

    Returns dict: {variants, patient_hpo, csv_cols, n_input, assembly}. Pass the result to
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

    # trio inheritance (optional) — assigns v["inheritance"]; de novo adds PS2 (supporting strong)
    if father_vcf or mother_vcf:
        n_trio=trio_inheritance(variants, father_vcf, mother_vcf, progress)
        for v in variants:
            if v.get("inheritance")=="de_novo":
                v["_raw"]=v.get("_raw",0)+30  # PS2 bonus
                v["acmg_tags"]=(v.get("acmg_tags","")+(",PS2" if v.get("acmg_tags") else "PS2"))
                v["evidence"]=(v.get("evidence","")+"; PS2: de novo (absent in both genotyped parents)").lstrip("; ")
        _prog(len(variants), len(variants), f"Herencia analizada en {n_trio} variantes")

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

    # agentic reasoning layer (self-reflection + traceable differential) — optional, needs LLM
    agentic_result=None
    if agentic:
        agentic_result=agentic_diagnosis(variants, patient_hpo, api_key=anthropic_key,
                                          sample=sample, progress=progress, start_k=reflect_k)
        # flag which rows the LLM self-reflection actually examined (only the top-K window)
        for v in variants:
            v["agentic_evaluated"]="yes" if v.get("reflect_verdict") else "no"

    return {"variants":variants, "patient_hpo":patient_hpo, "csv_cols":CSV_COLS,
            "n_input":len(variants), "assembly":assembly, "agentic":agentic_result}

def write_outputs(result, out_prefix):
    """Write <prefix>_annotated.csv and <prefix>_report.html from a run_pipeline() result."""
    import csv
    with open(f"{out_prefix}_annotated.csv","w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=result["csv_cols"], extrasaction="ignore"); w.writeheader()
        for v in result["variants"]: w.writerow(v)
    write_html(result["variants"], result["patient_hpo"], out_prefix, result.get("assembly","GRCh38"),
               result.get("agentic"))

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
    ap.add_argument("--agentic", action="store_true", help="Agentic layer: LLM self-reflection over candidates + traceable differential diagnosis with link verification (needs ANTHROPIC_API_KEY)")
    ap.add_argument("--reflect-k", type=int, default=8, help="Agentic layer: number of top candidates the LLM self-reflection examines (default 8; window widens only if all are refuted)")
    ap.add_argument("--father", default=None, help="Father VCF for trio analysis (de novo detection → ACMG PS2)")
    ap.add_argument("--mother", default=None, help="Mother VCF for trio analysis (de novo detection → ACMG PS2)")
    ap.add_argument("--anthropic-key", default=None, help="Anthropic API key (else uses ANTHROPIC_API_KEY env)")
    a=ap.parse_args()

    hpo_arg=a.hpo
    if hpo_arg.startswith("@"): hpo_arg=open(hpo_arg[1:]).read()
    note=open(a.clinical_note).read() if a.clinical_note else None

    def cli_progress(done, total, msg): print(f"[raredx] {msg}", file=sys.stderr)
    result=run_pipeline(a.vcf, sample=a.sample, hpo=hpo_arg, clinical_note_text=note,
                        assembly=a.assembly, use_esm=a.esm, use_am=a.alphamissense,
                        agentic=a.agentic, reflect_k=a.reflect_k,
                        father_vcf=a.father, mother_vcf=a.mother, email=a.email,
                        anthropic_key=a.anthropic_key, progress=cli_progress)
    write_outputs(result, a.out_prefix)
    print(f"[raredx] wrote {a.out_prefix}_annotated.csv and {a.out_prefix}_report.html", file=sys.stderr)

if __name__=="__main__":
    main()
