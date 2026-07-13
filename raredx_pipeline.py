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
import argparse, asyncio, json, os, sys, time, html, datetime, re, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

ENSEMBL="https://rest.ensembl.org"
GNOMAD="https://gnomad.broadinstitute.org/api"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OT="https://api.platform.opentargets.org/api/v4/graphql"
OLS="https://www.ebi.ac.uk/ols4/api"
UA={"User-Agent":"raredx-pipeline","Accept":"application/json"}
LOF_TERMS={"stop_gained","frameshift_variant","splice_acceptor_variant",
           "splice_donor_variant","start_lost","stop_lost","transcript_ablation"}
_REQUEST_STATE=threading.local()


class AnalysisCancelled(RuntimeError):
    pass


def set_request_deadline(deadline=None):
    _REQUEST_STATE.deadline=deadline


def set_request_cancel_event(event=None):
    _REQUEST_STATE.cancel_event=event


def _check_cancelled():
    event=getattr(_REQUEST_STATE,"cancel_event",None)
    if event is not None and event.is_set():
        raise AnalysisCancelled("analysis cancelled by user")


def _request_timeout():
    _check_cancelled()
    deadline=getattr(_REQUEST_STATE,"deadline",None)
    if deadline is None:
        return 30
    remaining=deadline-time.monotonic()
    if remaining<=0:
        raise TimeoutError("analysis exceeded its configured wall-clock deadline")
    return max(0.1,min(30,remaining))


def _retry_wait(attempt):
    _check_cancelled()
    if attempt >= 3:
        return
    delay=1.5*(attempt+1)
    deadline=getattr(_REQUEST_STATE,"deadline",None)
    if deadline is not None:
        remaining=deadline-time.monotonic()
        if remaining<=0:
            raise TimeoutError("analysis exceeded its configured wall-clock deadline")
        delay=min(delay,remaining)
    event=getattr(_REQUEST_STATE,"cancel_event",None)
    if event is not None:
        if event.wait(delay):
            raise AnalysisCancelled("analysis cancelled by user")
    else:
        time.sleep(delay)


def _annotation_workers():
    """Concurrency for live-API annotation (env RAREDX_ANNOTATION_WORKERS, default 8)."""
    try:
        return max(1,int(os.environ.get("RAREDX_ANNOTATION_WORKERS","8")))
    except (TypeError,ValueError):
        return 8


def _prefilter_af_cutoff():
    """Population AF at/above which a variant is Benign by BA1 (env RAREDX_PREFILTER_AF)."""
    try:
        return float(os.environ.get("RAREDX_PREFILTER_AF","0.05"))
    except (TypeError,ValueError):
        return 0.05


def _map_concurrent(items, fn, workers, on_progress=None):
    """Apply fn to each item across a thread pool, preserving input order.

    fn must re-inject the per-request deadline/cancel into its own worker thread
    (thread-local state does not propagate automatically). The first cancellation,
    timeout, or error is re-raised after the remaining futures are cancelled. Falls back
    to a serial loop for a single worker or a single item.
    """
    n=len(items)
    if n==0:
        return []
    if workers<=1 or n==1:
        out=[]
        for k,it in enumerate(items,1):
            out.append(fn(it))
            if on_progress: on_progress(k)
        return out
    out=[None]*n
    done=0
    with ThreadPoolExecutor(max_workers=min(workers,n)) as ex:
        futures={ex.submit(fn,it):idx for idx,it in enumerate(items)}
        try:
            for fut in as_completed(futures):
                out[futures[fut]]=fut.result()
                done+=1
                if on_progress: on_progress(done)
        except BaseException:
            ex.shutdown(wait=False,cancel_futures=True)
            raise
    return out


def _get(url, **kw):
    for i in range(4):
        try:
            r=requests.get(url, headers=UA, timeout=_request_timeout(), **kw)
            if r.status_code==200:
                return r.json()
            if r.status_code in (408,425,429) or 500<=r.status_code<600:
                _retry_wait(i)
                continue
            return None
        except requests.RequestException:
            _retry_wait(i)
    return None

def _gql(url, query, variables=None):
    for i in range(4):
        try:
            r=requests.post(url, json={"query":query,"variables":variables or {}}, headers=UA,
                            timeout=_request_timeout())
            if r.status_code==200:
                return r.json()
            if r.status_code not in (408,425,429) and not 500<=r.status_code<600:
                return None
            _retry_wait(i)
        except requests.RequestException:
            _retry_wait(i)
    return None

# ---------- VCF ----------
class VCFParseError(ValueError):
    pass


ASSEMBLY_CONTIG_LENGTHS={
    "GRCh38":{"1":248956422,"2":242193529,"X":156040895},
    "GRCh37":{"1":249250621,"2":243199373,"X":155270560},
}


def detect_vcf_assembly(path):
    """Detect GRCh37/GRCh38 from ##reference or assembly-specific contig lengths."""
    references=[]
    contigs={}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##reference="):
                references.append(line.split("=",1)[1].strip().lower())
            elif line.startswith("##contig=<"):
                id_match=re.search(r"(?:^|,)ID=([^,>]+)",line[10:])
                length_match=re.search(r"(?:^|,)length=(\d+)",line[10:],re.I)
                assembly_match=re.search(r"(?:^|,)assembly=([^,>]+)",line[10:],re.I)
                if id_match and length_match:
                    chrom=id_match.group(1)
                    chrom=chrom[3:] if chrom.lower().startswith("chr") else chrom
                    contigs[chrom.upper()]=int(length_match.group(1))
                if assembly_match:
                    references.append(assembly_match.group(1).strip().lower())
            elif line.startswith("#CHROM"):
                break
    reference=" ".join(references)
    aliases={
        "GRCh38":("grch38","hg38","b38","hs38","hs38dh","gcf_000001405.38"),
        "GRCh37":("grch37","hg19","b37","hs37d5","human_g1k_v37","gcf_000001405.25"),
    }
    matches={assembly for assembly,terms in aliases.items() if any(term in reference for term in terms)}
    for assembly,lengths in ASSEMBLY_CONTIG_LENGTHS.items():
        if any(contigs.get(chrom)==length for chrom,length in lengths.items()):
            matches.add(assembly)
    if len(matches)==1:
        return matches.pop()
    if len(matches)>1:
        raise VCFParseError("VCF metadata contains conflicting genome assembly markers")
    raise VCFParseError(
        "could not detect genome assembly; add ##reference=GRCh38/GRCh37 or standard contig lengths"
    )


def _sample_for_alt(format_keys, sample_value, alt_index):
    sample = dict(zip(format_keys, sample_value.split(":")))
    gt = sample.get("GT")
    if not gt:
        return sample
    sep = "|" if "|" in gt else "/"
    normalized = []
    for allele in gt.replace("|", "/").split("/"):
        if allele == ".":
            normalized.append(".")
            continue
        try:
            normalized.append("1" if int(allele) == alt_index else "0")
        except ValueError:
            normalized.append(".")
    sample["GT"] = sep.join(normalized)
    return sample


def _called_alt(sample):
    gt=(sample or {}).get("GT")
    if not gt:
        return None
    parts=gt.replace("|","/").split("/")
    if any(part=="." for part in parts):
        return None
    try:
        return any(int(part)>=1 for part in parts)
    except ValueError:
        return None


def _numeric(value, default=-1.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _genotype_signature(sample):
    gt=(sample or {}).get("GT")
    if not gt:
        return None
    parts=gt.replace("|","/").split("/")
    if any(part=="." for part in parts):
        return None
    try:
        return tuple(sorted(int(part) for part in parts))
    except ValueError:
        return None


def _vcf_record_rank(record):
    """Prefer the strongest call when the same exact allele occurs more than once."""
    sample=record.get("sample") or {}
    filter_value=str(record.get("filter") or "").upper()
    filter_rank=2 if filter_value=="PASS" else 1 if filter_value in {"",".","UNFILTERED"} else 0
    return (
        filter_rank,
        int(_genotype_signature(sample) is not None),
        _numeric(sample.get("GQ")),
        _numeric(record.get("qual")),
        _numeric(sample.get("DP")),
        int(record.get("rsid") is not None),
        sum(str(value) not in {"","."} for key,value in sample.items() if key!="GT"),
    )


def _merge_vcf_records(existing, candidate):
    winner=candidate if _vcf_record_rank(candidate)>_vcf_record_rank(existing) else existing
    existing_gt=_genotype_signature(existing.get("sample"))
    candidate_gt=_genotype_signature(candidate.get("sample"))
    conflict=bool(existing.get("genotype_conflict")) or (
        existing_gt is not None and candidate_gt is not None and existing_gt!=candidate_gt
    )
    if not conflict:
        return winner
    winner=dict(winner)
    winner["sample"]=dict(winner.get("sample") or {})
    winner["sample"]["GT"]="./."
    winner["genotype_conflict"]=True
    return winner


def parse_vcf(path, sample=None, allow_empty=True, max_variants=None, called_only=False):
    """Parse one record per ALT allele and select the requested sample when present."""
    out=[]
    allele_indexes={}
    excluded_records={}
    header_seen=False
    sample_names=[]
    sample_col=None
    with open(path) as fh:
        for line_no,line in enumerate(fh,1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header=line.rstrip("\n").split("\t")
                if len(header)<8:
                    raise VCFParseError("invalid VCF header")
                header_seen=True
                sample_names=header[9:]
                if sample_names:
                    if len(sample_names)>1:
                        if sample in (None,"","SAMPLE"):
                            names=", ".join(sample_names[:10])
                            suffix="..." if len(sample_names)>10 else ""
                            raise VCFParseError(
                                f"multi-sample VCF requires an explicit sample ID; available: {names}{suffix}"
                            )
                        if sample not in sample_names:
                            raise VCFParseError(
                                f"sample {sample!r} not found in multi-sample VCF"
                            )
                    sample_col=(sample_names.index(sample) if sample in sample_names else 0)+9
                continue
            if line.startswith("#") or not line.strip():
                continue
            if not header_seen:
                raise VCFParseError("invalid VCF: missing #CHROM header")
            f=line.rstrip("\n").split("\t")
            if len(f)<8:
                raise VCFParseError(f"invalid VCF record at line {line_no}: expected at least 8 columns")
            try:
                pos=int(f[1])
            except ValueError as exc:
                raise VCFParseError(f"invalid VCF position at line {line_no}") from exc
            if pos<1 or not f[3] or f[3]=="." or not f[4] or f[4]==".":
                continue
            chrom=f[0][3:] if f[0].lower().startswith("chr") else f[0]
            for alt_index,alt in enumerate(f[4].split(","),1):
                if not alt or alt=="." or alt=="*" or alt.startswith("<") or "[" in alt or "]" in alt:
                    continue
                rec=dict(chrom=chrom,pos=pos,rsid=f[2] if f[2]!="." else None,
                         ref=f[3],alt=alt,qual=f[5],filter=f[6] or "PASS")
                if sample_col is not None and len(f)>sample_col:
                    rec["sample"]=_sample_for_alt(f[8].split(":"),f[sample_col],alt_index)
                allele_key=(chrom,pos,f[3].upper(),alt.upper())
                if called_only and _called_alt(rec.get("sample")) is not True:
                    existing_index=allele_indexes.get(allele_key)
                    if existing_index is not None:
                        out[existing_index]=_merge_vcf_records(out[existing_index],rec)
                    else:
                        previous=excluded_records.get(allele_key)
                        excluded_records[allele_key]=(
                            _merge_vcf_records(previous,rec) if previous is not None else rec
                        )
                    continue
                if called_only and allele_key in excluded_records:
                    rec=_merge_vcf_records(excluded_records.pop(allele_key),rec)
                    if _called_alt(rec.get("sample")) is not True:
                        continue
                existing_index=allele_indexes.get(allele_key)
                if existing_index is not None:
                    out[existing_index]=_merge_vcf_records(out[existing_index],rec)
                    continue
                allele_indexes[allele_key]=len(out)
                out.append(rec)
                if max_variants is not None and len(out)>max_variants:
                    raise VCFParseError(f"VCF exceeds the {max_variants} variant analysis limit")
    if not header_seen:
        raise VCFParseError("invalid VCF: missing #CHROM header")
    if called_only:
        out=[record for record in out if _called_alt(record.get("sample")) is True]
    if not out and not allow_empty:
        raise VCFParseError("VCF contains no supported variant records")
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
    if not j: return {"annotation_available":False}
    r=j[0]
    tcs=r.get("transcript_consequences") or []
    most_severe=r.get("most_severe_consequence")
    severe=[t for t in tcs if most_severe in (t.get("consequence_terms") or [])]
    ranked=sorted(severe or tcs,key=lambda t:-int(t.get("canonical",0) or 0))
    representative=ranked[0] if ranked else {}
    gene=representative.get("gene_symbol")
    gid=representative.get("gene_id")
    sift=representative.get("sift_prediction")
    poly=representative.get("polyphen_prediction")
    aa=representative.get("amino_acids")
    # gnomAD AF from colocated frequencies (populated on the GRCh37 endpoint)
    gaf=None
    for cv in (r.get("colocated_variants") or []):
        fr=cv.get("frequencies") or {}
        d=fr.get(rec["alt"]) or fr.get(rec["alt"].upper())
        if d:
            vals=[d.get(k) for k in ("gnomade","gnomadg","gnomad") if d.get(k) is not None]
            if vals:
                gaf=max(vals) if gaf is None else max(gaf,*vals)
    return dict(most_severe=most_severe,gene=gene,gene_id=gid,
                amino_acids=aa,sift=sift,polyphen=poly,gnomad_af=gaf,
                annotation_available=True,
                colocated=[c.get("id") for c in (r.get("colocated_variants") or []) if str(c.get("id","")).startswith("rs")])

GNOMAD_VAR_Q="""query($vid:String!,$ds:DatasetId!){
 variant(variantId:$vid, dataset:$ds){ exome{af homozygote_count} genome{af homozygote_count} } }"""
def gnomad_af(rec, return_status=False):
    vid=f'{rec["chrom"]}-{rec["pos"]}-{rec["ref"]}-{rec["alt"]}'
    j=_gql(GNOMAD, GNOMAD_VAR_Q, {"vid":vid,"ds":"gnomad_r4"})
    if not isinstance(j,dict):
        return (None,False) if return_status else None
    v=(j.get("data") or {}).get("variant")
    if j.get("errors") and not v:
        return (None,False) if return_status else None
    if not v:
        return (None,True) if return_status else None
    afs=[x.get("af") for x in (v.get("exome"),v.get("genome")) if x and x.get("af") is not None]
    af=max(afs) if afs else None
    return (af,True) if return_status else af

GNOMAD_CONS_Q="""query($sym:String!){ gene(gene_symbol:$sym, reference_genome:GRCh38){
  gnomad_constraint{ pli oe_lof_upper } } }"""
def gnomad_constraint(sym, cache):
    if sym in cache: return cache[sym]
    j=_gql(GNOMAD, GNOMAD_CONS_Q, {"sym":sym})
    gene=((j or {}).get("data") or {}).get("gene") if isinstance(j,dict) else None
    available=isinstance(j,dict) and not (j.get("errors") and not gene)
    c=(((j or {}).get("data") or {}).get("gene") or {}).get("gnomad_constraint") or {}
    cache[sym]={"pli":c.get("pli"),"loeuf":c.get("oe_lof_upper"),"available":available}
    return cache[sym]

def _clinvar_stars(doc):
    rs=((doc.get("germline_classification") or {}).get("review_status") or "").lower()
    return (4 if "practice guideline" in rs else 3 if "expert panel" in rs
            else 2 if "multiple" in rs and "conflict" not in rs else 1 if "single" in rs else 0)


def _normalized_vcf_allele(rec):
    pos0=int(rec["pos"])-1
    ref=str(rec["ref"]).upper()
    alt=str(rec["alt"]).upper()
    while ref and alt and ref[0]==alt[0]:
        ref=ref[1:]
        alt=alt[1:]
        pos0+=1
    while ref and alt and ref[-1]==alt[-1]:
        ref=ref[:-1]
        alt=alt[:-1]
    return pos0,ref,alt


def _clinvar_matches_allele(doc, rec, assembly):
    vcf_pos0,vcf_ref,vcf_alt=_normalized_vcf_allele(rec)
    for variation in doc.get("variation_set") or []:
        locations=variation.get("variation_loc") or []
        location_match=any(
            loc.get("assembly_name")==assembly
            and str(loc.get("chr","")).removeprefix("chr")==str(rec["chrom"]).removeprefix("chr")
            and min(int(loc.get("start") or 0),int(loc.get("stop") or loc.get("start") or 0))
                <= vcf_pos0+1
                <= max(int(loc.get("start") or 0),int(loc.get("stop") or loc.get("start") or 0))
            for loc in locations
        )
        spdi=variation.get("canonical_spdi") or ""
        parts=spdi.rsplit(":",3)
        if len(parts)!=4:
            continue
        try:
            spdi_pos0=int(parts[-3])
        except ValueError:
            continue
        allele_match=parts[-2].upper()==vcf_ref and parts[-1].upper()==vcf_alt
        position_match=(spdi_pos0==vcf_pos0 if assembly=="GRCh38" else location_match)
        if position_match and allele_match:
            return True
    return False


def clinvar(rec, email=None, assembly="GRCh38"):
    """Return ClinVar evidence only when its genomic placement and allele match the VCF."""
    params={"db":"clinvar","retmode":"json"}
    if email: params["email"]=email
    term=(f'{rec["rsid"]}[All Fields]' if rec.get("rsid")
          else f'{rec["chrom"]}[chr] AND ({rec["pos"]}[chrpos37] OR {rec["pos"]}[chrpos38])')
    es=_get(f"{EUTILS}/esearch.fcgi", params={**params,"term":term,"retmax":100})
    if es is None:
        return {"significance":None,"stars":0,"available":False}
    ids=(((es or {}).get("esearchresult") or {}).get("idlist")) or []
    if not ids:
        return {"significance":None,"stars":0,"available":True}
    su=_get(f"{EUTILS}/esummary.fcgi", params={**params,"id":",".join(ids)})
    if su is None:
        return {"significance":None,"stars":0,"available":False}
    result=(su or {}).get("result") or {}
    docs=[result.get(i) or {} for i in ids]
    matches=[d for d in docs if _clinvar_matches_allele(d,rec,assembly)]
    if not matches:
        return {"significance":None,"stars":0,"available":True}
    doc=max(matches,key=_clinvar_stars)
    germ=doc.get("germline_classification") or {}
    desc=germ.get("description")
    stars=_clinvar_stars(doc)
    conds=[t.get("trait_name") for t in (germ.get("trait_set") or []) if t.get("trait_name")]
    return {"significance":desc,"stars":stars,"conditions":conds[:3],
            "protein_change":doc.get("protein_change"),"available":True}

# ---------- classification ----------
def _call_from_tags(tags, sig=""):
    codes={t[0] for t in tags}
    if "BA1" in codes:
        call="Benign"
    elif "ClinVar_conflicting" in codes:
        call="Uncertain significance (VUS)"
    else:
        very_strong=int("PVS1" in codes)
        strong=len({"PS_ClinVar","PS2"}&codes)
        moderate=len({"PM2","PVS1_mod"}&codes)
        supporting=int("PP5_ClinVar" in codes)+int(bool({"PP3","PP3_ESM","PP3_AM"}&codes))
        benign_strong=len({"BS1","BS_ClinVar"}&codes)
        benign_supporting=int("BP6_ClinVar" in codes)+int(
            bool({"BP4","BP4_ESM","BP4_AM"}&codes)
        )
        pathogenic_evidence=very_strong+strong+moderate+supporting
        benign_evidence=benign_strong+benign_supporting
        if pathogenic_evidence and benign_evidence:
            call="Uncertain significance (VUS)"
        elif benign_strong>=2:
            call="Benign"
        elif (benign_strong>=1 and benign_supporting>=1) or benign_supporting>=2:
            call="Likely benign"
        elif (
            (very_strong and strong>=1)
            or (very_strong and moderate>=2)
            or (very_strong and moderate>=1 and supporting>=1)
            or (very_strong and supporting>=2)
            or strong>=2
            or (strong>=1 and moderate>=3)
            or (strong>=1 and moderate>=2 and supporting>=2)
            or (strong>=1 and moderate>=1 and supporting>=4)
        ):
            call="Pathogenic"
        elif (
            (very_strong and moderate>=1)
            or (strong>=1 and 1<=moderate<=2)
            or (strong>=1 and supporting>=2)
            or moderate>=3
            or (moderate>=2 and supporting>=2)
            or (moderate>=1 and supporting>=4)
        ):
            call="Likely pathogenic"
        else:
            call="Uncertain significance (VUS)"
    if "drug response" in (sig or "").lower():
        call+=" | Pharmacogenomic (drug response)"
    return call


def classify(af, cons, pli, loeuf, sig, stars, sift, poly, af_available=True, extra_tags=None):
    tags=[]
    if af is not None and af>=0.05: tags.append(("BA1",f"gnomAD AF {af:.3f} >=5%"))
    elif af is not None and af>=0.01: tags.append(("BS1",f"gnomAD AF {af:.3f} >=1%"))
    if not af_available: tags.append(("gnomAD_unavailable","gnomAD lookup unavailable; absence not assessed"))
    elif af is None: tags.append(("PM2","absent from gnomAD r4"))
    elif af<1e-4: tags.append(("PM2",f"gnomAD AF {af:.2e} <0.01%"))
    if cons in LOF_TERMS:
        constraint=(
            "gene is LoF-constrained"
            if (pli is not None and pli>=0.9) or (loeuf is not None and loeuf<0.6)
            else "gene constraint is modest or unavailable"
        )
        tags.append((
            "LoF_predicted",
            f"{cons}; {constraint}; PVS1 requires a confirmed disease mechanism, "
            "transcript relevance, and NMD assessment",
        ))
    if sift and "deleterious" in sift and poly and "damaging" in poly:
        tags.append(("PP3",f"SIFT {sift} & PolyPhen {poly}"))
    elif sift=="tolerated" and poly and "benign" in poly:
        tags.append(("BP4",f"SIFT {sift} & PolyPhen {poly}"))
    sig=sig or ""; sig_l=sig.lower()
    if "conflict" in sig_l:
        tags.append(("ClinVar_conflicting",f"ClinVar {sig} ({stars} stars)"))
    elif "pathogenic" in sig_l:
        tags.append((("PS_ClinVar" if stars>=2 else "PP5_ClinVar"),f"ClinVar {sig} ({stars} stars)"))
    elif "benign" in sig_l:
        tags.append((("BS_ClinVar" if stars>=2 else "BP6_ClinVar"),f"ClinVar {sig} ({stars} stars)"))
    tags.extend(extra_tags or [])
    return _call_from_tags(tags,sig),tags

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

def expand_hpo(hpo_ids, return_status=False):
    full=set(hpo_ids)
    available=True
    for hid in list(hpo_ids):
        iri=f"http://purl.obolibrary.org/obo/{hid.replace(':','_')}"
        j=_get(f"{OLS}/ontologies/hp/ancestors", params={"id":iri,"size":200})
        if not isinstance(j,dict):
            available=False
        for t in (((j or {}).get("_embedded") or {}).get("terms")) or []:
            oid=t.get("obo_id")
            if oid and oid.startswith("HP:"): full.add(oid)
    full-={"HP:0000001","HP:0000118"}
    return (full,available) if return_status else full


def expand_hpo_profile(tokens, return_status=False):
    """Resolve direct terms and return their labeled HPO ancestor graph for human review."""
    direct=resolve_hpo(tokens)
    terms={}
    available=True
    for item in direct:
        hpo_id=item["hpo_id"]
        label=item.get("label") or hpo_id
        if label==hpo_id:
            j=_get(f"{OLS}/search",params={"q":hpo_id,"ontology":"hp","rows":10})
            if j is None:
                available=False
            docs=(((j or {}).get("response") or {}).get("docs")) or []
            exact=next((d for d in docs if d.get("obo_id")==hpo_id),None)
            if exact and exact.get("label"):
                label=exact["label"]
        terms[hpo_id]={
            "hpo_id":hpo_id,"label":label,"kind":"direct","source_hpo_ids":[hpo_id]
        }
        iri=f"http://purl.obolibrary.org/obo/{hpo_id.replace(':','_')}"
        j=_get(f"{OLS}/ontologies/hp/ancestors",params={"id":iri,"size":200})
        if j is None:
            available=False
        for ancestor in (((j or {}).get("_embedded") or {}).get("terms")) or []:
            ancestor_id=ancestor.get("obo_id")
            if not ancestor_id or ancestor_id in {"HP:0000001","HP:0000118"}:
                continue
            existing=terms.get(ancestor_id)
            if existing:
                if hpo_id not in existing["source_hpo_ids"]:
                    existing["source_hpo_ids"].append(hpo_id)
                continue
            terms[ancestor_id]={
                "hpo_id":ancestor_id,
                "label":ancestor.get("label") or ancestor_id,
                "kind":"ancestor",
                "source_hpo_ids":[hpo_id],
            }
    result=sorted(
        terms.values(),
        key=lambda term:(term["kind"]!="direct",str(term["label"]).casefold(),term["hpo_id"]),
    )
    return (result,available) if return_status else result

OT_Q="""query($id:String!){ target(ensemblId:$id){ associatedDiseases(page:{size:25,index:0}){
  rows{ score disease{ id name phenotypes(page:{size:80,index:0}){ rows{ phenotypeHPO{ id name } } } } } } } }"""
def gene_pheno_score(ensg, patient_ids, patient_full, return_status=False):
    empty=(0.0,0,0,"","")
    if not ensg:
        return (*empty,True) if return_status else empty
    j=_gql(OT, OT_Q, {"id":ensg})
    available=isinstance(j,dict) and not (j.get("errors") and not (j.get("data") or {}).get("target"))
    if not available:
        return (*empty,False) if return_status else empty
    t=(((j or {}).get("data") or {}).get("target")) or {}
    rows=(t.get("associatedDiseases") or {}).get("rows") or []
    best=None
    valid_associations=0
    missing_associations=False
    for row in rows:
        association=_numeric(row.get("score"),None)
        if association is None:
            missing_associations=True
            continue
        valid_associations+=1
        d=row["disease"]
        dh={p["phenotypeHPO"]["id"].replace("_",":"):p["phenotypeHPO"]["name"]
            for p in ((d.get("phenotypes") or {}).get("rows") or []) if p.get("phenotypeHPO")}
        direct=set(dh)&patient_ids
        matched=set(dh)&patient_full
        ancestor=matched-direct
        coverage=(len(direct)+0.4*len(ancestor))/max(len(patient_ids),1)
        association=min(1.0,max(0.0,association))
        score=min(coverage,1.0)*association
        candidate=(
            score,
            len(direct),
            len(matched),
            association,
            d.get("name") or "",
            [dh[h] for h in sorted(matched)][:5],
        )
        if matched and (best is None or candidate[:4]>best[:4]):
            best=candidate
    result=(
        round(best[0],3) if best else 0.0,
        best[1] if best else 0,
        best[2] if best else 0,
        best[4] if best else "",
        "; ".join(best[5]) if best else "",
    )
    if rows and missing_associations and valid_associations==0:
        available=False
    return (*result,available) if return_status else result


# ---------- AI module A: ESM-2 missense pathogenicity (optional) ----------
_ESM_CACHE = {}
_ESM_LOCK = threading.Lock()
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
        with _ESM_LOCK:
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


def _reverse_complement(allele):
    try:
        return allele.upper().translate(str.maketrans("ACGTN","TGCAN"))[::-1]
    except AttributeError:
        return None


def _liftover_allele_37_to_38(chrom, pos, ref, alt):
    """Map a complete GRCh37 allele to GRCh38, including strand and chromosome changes."""
    c=str(chrom).removeprefix("chr")
    end=int(pos)+len(ref)-1
    j=_get(f"{ENSEMBL_GRCH37}/map/human/GRCh37/{c}:{pos}..{end}/GRCh38")
    mappings=(j or {}).get("mappings") or []
    if len(mappings)!=1:
        return None
    mapped=mappings[0].get("mapped") or {}
    try:
        start=int(mapped["start"])
        mapped_end=int(mapped["end"])
        strand=int(mapped.get("strand",1))
    except (KeyError,TypeError,ValueError):
        return None
    if abs(mapped_end-start)+1!=len(ref):
        return None
    mapped_ref=ref.upper()
    mapped_alt=alt.upper()
    if strand==-1:
        mapped_ref=_reverse_complement(mapped_ref)
        mapped_alt=_reverse_complement(mapped_alt)
    if not mapped_ref or not mapped_alt:
        return None
    return {
        "chrom":str(mapped.get("seq_region_name") or c).removeprefix("chr"),
        "pos":min(start,mapped_end),
        "ref":mapped_ref,
        "alt":mapped_alt,
        "strand":strand,
    }


def liftover_37_to_38(chrom, pos):
    """Map a GRCh37 position to GRCh38 via the Ensembl assembly-mapping REST endpoint."""
    mapped=_liftover_allele_37_to_38(chrom,pos,"N","N")
    return mapped["pos"] if mapped else None


def _vep_snv_matches(j, ref, alt):
    if not j or len(ref)!=1 or len(alt)!=1:
        return False
    alleles=str(j[0].get("allele_string") or "").upper().split("/")
    return len(alleles)>=2 and alleles[0]==ref.upper() and alt.upper() in alleles[1:]

def alphamissense_score(chrom, pos, ref, alt, assembly="GRCh38", gene=None):
    """Return {'am_pathogenicity': float, 'am_class': str} for a missense variant, or None.
    Queries Ensembl VEP GRCh38 with AlphaMissense=1. If the input is GRCh37, lifts over first."""
    if len(ref)!=1 or len(alt)!=1:
        return None
    mapped={"chrom":str(chrom).removeprefix("chr"),"pos":int(pos),"ref":ref.upper(),"alt":alt.upper()}
    if assembly=="GRCh37":
        mapped=_liftover_allele_37_to_38(chrom,pos,ref,alt)
    if not mapped:
        return None
    j=_get(
        f"{ENSEMBL}/vep/human/region/{mapped['chrom']}:{mapped['pos']}-{mapped['pos']}/{mapped['alt']}",
        params={"AlphaMissense":1},
    )
    if not _vep_snv_matches(j,mapped["ref"],mapped["alt"]):
        hgvs=f"{mapped['chrom']}:g.{mapped['pos']}{mapped['ref']}>{mapped['alt']}"
        j=_get(f"{ENSEMBL}/vep/human/hgvs/{hgvs}",params={"AlphaMissense":1})
    if not _vep_snv_matches(j,mapped["ref"],mapped["alt"]):
        return None
    tcs = j[0].get("transcript_consequences") or []
    # prefer the transcript for the annotated gene, else any transcript carrying an AM score
    withscore = [t for t in tcs if t.get("alphamissense")]
    if not withscore:
        return None
    same = [t for t in withscore if gene and t.get("gene_symbol") == gene]
    t = (same or withscore)[0]
    am = t["alphamissense"]
    try:
        score=float(am.get("am_pathogenicity"))
    except (TypeError,ValueError):
        return None
    am_class=am.get("am_class")
    if not 0<=score<=1 or am_class not in {"likely_pathogenic","likely_benign","ambiguous"}:
        return None
    return {"am_pathogenicity":score,"am_class":am_class}

def am_call(am):
    """Normalize AlphaMissense class to the pipeline's deleterious/tolerated/ambiguous vocabulary."""
    if not am or am.get("am_class") is None:
        return ""
    cls = am["am_class"]
    return {"likely_pathogenic": "deleterious", "likely_benign": "tolerated",
            "ambiguous": "ambiguous"}.get(cls, "ambiguous")

def protein_context(rsid, chrom, pos, ref, alt, gene=None, assembly="GRCh38"):
    """Get protein sequence + AA position for a missense variant via Ensembl VEP + sequence.
    Always queries by REGION+ALLELE (not rsID): an rsID can carry several alt alleles, and the
    VEP-by-id route may return the amino-acid change for a DIFFERENT allele than the one in the
    VCF. Region+allele guarantees the scored mutant residue matches the patient's actual alt."""
    mapped={"chrom":str(chrom).removeprefix("chr"),"pos":int(pos),"ref":ref.upper(),"alt":alt.upper()}
    if assembly=="GRCh37":
        mapped=_liftover_allele_37_to_38(chrom,pos,ref,alt)
    if not mapped:
        return None
    url = (
        f"{ENSEMBL}/vep/human/region/{mapped['chrom']}:{mapped['pos']}-"
        f"{mapped['pos']+len(mapped['ref'])-1}/{mapped['alt']}"
    )
    j = _get(url, params={"content-type":"application/json"})
    if not j:
        return None
    if len(mapped["ref"])==len(mapped["alt"])==1 and not _vep_snv_matches(
            j,mapped["ref"],mapped["alt"]):
        return None
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


# ---------- LLM providers ----------
LLM_MODEL = "claude-sonnet-4-5"
_LLM_STATE = threading.local()


def _reset_llm_state():
    _LLM_STATE.providers=set()
    _LLM_STATE.errors=[]


def _record_llm_provider(provider):
    if not hasattr(_LLM_STATE,"providers"):
        _reset_llm_state()
    _LLM_STATE.providers.add(provider)


def _record_llm_error(provider,error):
    if not hasattr(_LLM_STATE,"errors"):
        _reset_llm_state()
    _LLM_STATE.errors.append(f"{provider}: {error}")


async def _copilot_invoke(system_p,user_p,model,timeout,isolated_home):
    from copilot import CopilotClient
    from copilot.generated.rpc import PermissionDecisionReject
    from copilot.session_events import AssistantMessageData

    def deny_tools(request,invocation):
        return PermissionDecisionReject(feedback="Tools are disabled for this inference-only session.")

    client=CopilotClient(
        log_level="error",
        base_directory=isolated_home,
        working_directory=isolated_home,
    )
    session=None
    try:
        await asyncio.wait_for(client.start(),timeout=20)
        session=await asyncio.wait_for(client.create_session(
            model=model,
            on_permission_request=deny_tools,
            available_tools=[],
            system_message={"mode":"replace","content":system_p+"\nDo not use tools. Return only the requested output."},
            skip_custom_instructions=True,
            enable_config_discovery=False,
            enable_on_demand_instruction_discovery=False,
            enable_file_hooks=False,
            enable_host_git_operations=False,
            enable_session_store=False,
            enable_skills=False,
            hooks={},
            config_directory=isolated_home,
            memory={"enabled":False},
            infinite_sessions={"enabled":False},
        ),timeout=20)
        reply=await asyncio.wait_for(session.send_and_wait(user_p),timeout=timeout)
        if not reply or not isinstance(reply.data,AssistantMessageData):
            raise RuntimeError("Copilot returned no assistant message")
        return reply.data.content
    finally:
        if session is not None:
            try:
               await asyncio.wait_for(client.delete_session(session.session_id),timeout=5)
            except Exception:
               pass
        try:
            await asyncio.wait_for(client.stop(),timeout=5)
        except Exception:
            try:
               await asyncio.wait_for(client.force_stop(),timeout=5)
            except Exception:
               pass


def _copilot_process_worker(result_queue,system_p,user_p,model,timeout,isolated_home):
    try:
        raw=asyncio.run(_copilot_invoke(system_p,user_p,model,timeout,isolated_home))
        result_queue.put(("ok",raw))
    except BaseException as exc:
        result_queue.put(("error",f"{type(exc).__name__}: {exc}"))


def _copilot_raw(system_p,user_p):
    import multiprocessing,os,queue,shutil,tempfile
    try:
        import copilot
    except ImportError as exc:
        raise RuntimeError("install github-copilot-sdk to use the Copilot provider") from exc
    model=os.environ.get("RAREDX_COPILOT_MODEL","gpt-5-mini")
    timeout=max(10,int(os.environ.get("RAREDX_LLM_TIMEOUT_SECONDS","120")))
    isolated_home=tempfile.mkdtemp(prefix="raredx-copilot-")
    context=multiprocessing.get_context("spawn")
    result_queue=context.Queue(maxsize=1)
    process=context.Process(
        target=_copilot_process_worker,
        args=(result_queue,system_p,user_p,model,timeout,isolated_home),
    )
    try:
        process.start()
        end=time.monotonic()+timeout+55
        while process.is_alive() and time.monotonic()<end:
            process.join(0.5)
            _check_cancelled()
        if process.is_alive():
            process.terminate()
            process.join(5)
            if process.is_alive():
               process.kill()
               process.join(5)
            raise TimeoutError(f"Copilot inference exceeded {timeout} seconds")
        try:
            status,payload=result_queue.get(timeout=2)
        except queue.Empty as exc:
            raise RuntimeError(f"Copilot worker exited with code {process.exitcode} without a result") from exc
        if status!="ok":
            raise RuntimeError(payload)
        return payload
    finally:
        if process.is_alive():
            process.kill()
            process.join(5)
        result_queue.close()
        result_queue.join_thread()
        shutil.rmtree(isolated_home,ignore_errors=True)


def llm_diagnostics():
    return {
        "providers":sorted(getattr(_LLM_STATE,"providers",set())),
        "errors":list(getattr(_LLM_STATE,"errors",[])),
    }


def _anthropic_raw(system_p,user_p,api_key=None,max_tokens=2500):
    import os
    key=api_key or os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        raise RuntimeError("ANTHROPIC_API_KEY is not configured")
    try:
        import anthropic
    except ImportError as exc:
        raise RuntimeError("install anthropic to use the Anthropic provider") from exc
    msg=anthropic.Anthropic(api_key=key).messages.create(
        model=LLM_MODEL,max_tokens=max_tokens,system=system_p,
        messages=[{"role":"user","content":user_p}])
    return msg.content[0].text


# ---------- AI module B: extract HPO terms from a free-text clinical note ----------
def _ground_hpo_phrase(phrase):
    j=_get(f"{OLS}/search",params={"q":phrase,"ontology":"hp","rows":5})
    if j is None:
        return None,False
    docs=(((j or {}).get("response") or {}).get("docs")) or []
    exact=next(
        (d for d in docs if str(d.get("label") or "").casefold()==phrase.casefold()),
        None,
    )
    pick=exact or (docs[0] if docs else None)
    if not pick or not str(pick.get("obo_id") or "").startswith("HP:"):
        return None,True
    return {"obo_id":pick["obo_id"],"label":str(pick.get("label") or phrase)},True


def _hpo_ancestor_options(candidate):
    iri=f"http://purl.obolibrary.org/obo/{candidate['hpo_id'].replace(':','_')}"
    j=_get(f"{OLS}/ontologies/hp/ancestors",params={"id":iri,"size":200})
    options={candidate["hpo_id"]:{"hpo_id":candidate["hpo_id"],"label":candidate["label"]}}
    if j is None:
        return list(options.values()),False
    for term in (((j or {}).get("_embedded") or {}).get("terms")) or []:
        hpo_id=str(term.get("obo_id") or "")
        label=str(term.get("label") or "")
        if hpo_id.startswith("HP:") and hpo_id not in {"HP:0000001","HP:0000118"} and label:
            options.setdefault(hpo_id,{"hpo_id":hpo_id,"label":label})
    return list(options.values())[:100],True


def _verify_hpo_candidates(note_text,candidates,api_key=None):
    """Select the most specific supported term from each candidate's HPO ancestor chain."""
    if not candidates:
        return []
    verification_items=[]; allowed_by_candidate={}; ancestor_failures=0
    for candidate in candidates:
        allowed,available=_hpo_ancestor_options(candidate)
        if not available:
            ancestor_failures+=1
        allowed_by_candidate[candidate["hpo_id"]]={term["hpo_id"]:term for term in allowed}
        verification_items.append({"candidate":candidate,"allowed_terms":allowed})
    if ancestor_failures:
        _record_llm_error("OLS HPO verification",f"{ancestor_failures} ancestor lookup(s) failed")
    system_p=(
        "You are an independent clinical-note entailment verifier for HPO coding. For every candidate, "
        "compare its official label and allowed HPO ancestor terms with the complete note. Select the "
        "MOST SPECIFIC allowed term whose phenotype and every qualifier are explicitly stated. Never "
        "infer subtype, anatomy, "
        "laterality, severity, frequency, duration, age of onset, chronicity, cause, inheritance, "
        "treatment response, or associated findings. The selected_hpo_id MUST be one of allowed_terms; "
        "use an empty string when none is supported. Evaluate all medical domains identically. Return "
        "STRICT JSON."
    )
    user_p=(
        f"Clinical note:\n{note_text}\n\nVerification items:\n"
        f"{json.dumps(verification_items,ensure_ascii=False)}\n\n"
        'Return JSON: {"decisions":[{"candidate_hpo_id":"HP:...",'
        '"selected_hpo_id":"HP:... or empty","reason":"brief"}]}'
    )
    parsed=_llm_json(system_p,user_p,api_key,max_tokens=2500)
    decisions=(parsed or {}).get("decisions") if isinstance(parsed,dict) else None
    if not isinstance(decisions,list):
        _record_llm_error("HPO verification","invalid verifier response")
        return []
    by_id={candidate["hpo_id"]:candidate for candidate in candidates}
    verified=[]; seen=set()
    for decision in decisions[:len(candidates)*2]:
        if not isinstance(decision,dict):
            continue
        candidate_id=str(decision.get("candidate_hpo_id") or "")
        candidate=by_id.get(candidate_id)
        if not candidate:
            continue
        selected_id=str(decision.get("selected_hpo_id") or "")
        selected_term=allowed_by_candidate[candidate_id].get(selected_id)
        selected=(
            {"hpo_id":selected_term["hpo_id"],"label":selected_term["label"],
             "note_evidence":candidate["note_evidence"]}
            if selected_term else None
        )
        if selected and selected["hpo_id"] not in seen:
            seen.add(selected["hpo_id"])
            verified.append(selected)
    return verified


def extract_hpo_from_note(note_text, api_key=None):
    """Use an LLM to extract present phenotypes from a clinical note, then ground each
    to an official HPO ID via OLS4. Uses GitHub Copilot by default, with Anthropic fallback."""
    sys_p = (
        "You extract phenotypic abnormalities from a clinical note and normalize each to a concise "
        "English HPO-style phrase. Extract ONLY explicitly present abnormalities: not negated findings, "
        "family history, treatments, or inferred subtypes. Every modifier in hpo_phrase must be stated "
        "in the quoted evidence. Never add specificity about subtype, anatomy, laterality, severity, "
        "frequency, duration, onset, chronicity, cause, inheritance, or treatment response. Return "
        "STRICT JSON."
    )
    user_p = (f"Clinical note:\n{note_text}\n\nReturn JSON: "
              '{"phenotypes":[{"note_evidence":"quote","hpo_phrase":"English phenotype term"}]}')
    parsed=_llm_json(sys_p,user_p,api_key,max_tokens=1500)
    if not isinstance(parsed,dict):
        return []
    phenos=parsed.get("phenotypes",[])
    if not isinstance(phenos,list):
        return []
    # Ground candidates first, then verify every official label against the complete note.
    seen=set(); candidates=[]; grounding_failures=0
    note_fold=note_text.casefold()
    for p in phenos[:50]:
        if not isinstance(p,dict):
            continue
        phrase=str(p.get("hpo_phrase") or "").strip()[:200]
        if not phrase:
            continue
        pick,available=_ground_hpo_phrase(phrase)
        if not available:
            grounding_failures+=1
        if not pick:
            continue
        evidence=str(p.get("note_evidence") or "").strip()[:500]
        if evidence and evidence.casefold() not in note_fold:
            evidence=""
        if pick["obo_id"] not in seen:
            seen.add(pick["obo_id"])
            candidates.append({"hpo_id":pick["obo_id"],"label":pick["label"],
                               "note_evidence":evidence})
    if grounding_failures:
        _record_llm_error("OLS HPO grounding",f"{grounding_failures} lookup(s) failed")
    return _verify_hpo_candidates(note_text,candidates,api_key)

def _llm_raw(system_p, user_p, api_key=None, max_tokens=2500):
    """Return LLM text through Copilot, Anthropic, or a host accessor."""
    import os
    configured=os.environ.get("RAREDX_LLM_PROVIDER")
    provider=(configured or ("anthropic" if api_key else "auto")).lower()
    if provider not in {"auto","copilot","anthropic","host"}:
        raise ValueError("RAREDX_LLM_PROVIDER must be auto, copilot, anthropic, or host")
    if provider in {"auto","copilot"}:
        try:
            raw=_copilot_raw(system_p,user_p)
            _record_llm_provider("GitHub Copilot")
            return raw
        except Exception as e:
            _record_llm_error("GitHub Copilot",e)
            print(f"[raredx] GitHub Copilot error: {e}",file=sys.stderr)
            if provider=="copilot":
                return None
    if provider in {"auto","anthropic"}:
        try:
            raw=_anthropic_raw(system_p,user_p,api_key,max_tokens)
            _record_llm_provider("Anthropic")
            return raw
        except Exception as e:
            _record_llm_error("Anthropic",e)
            print(f"[raredx] Anthropic error: {e}",file=sys.stderr)
            if provider=="anthropic":
                return None
    _host=globals().get("host")
    if provider in {"auto","host"} and _host is not None and hasattr(_host,"llm"):
        try:
            raw=_host.llm(user_p,system=system_p,max_tokens=max_tokens)["text"]
            _record_llm_provider("host.llm")
            return raw
        except Exception as e:
            _record_llm_error("host.llm",e)
            print(f"[raredx] host.llm error: {e}",file=sys.stderr)
    print("[raredx] no configured LLM provider is available",file=sys.stderr)
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
    if any(p=="." for p in parts):
        return None
    try:
        return {int(p) for p in parts}
    except ValueError:
        return None

def _has_alt(sample):
    """Return alt presence, or None when the genotype is not fully callable."""
    al = _gt_alleles(sample)
    return None if al is None else any(a >= 1 for a in al)

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
        v["father_gt"] = ((fs or {}).get("GT") or ("not_present" if fidx is not None else "NA"))
        v["mother_gt"] = ((ms or {}).get("GT") or ("not_present" if midx is not None else "NA"))
        child_alt = _has_alt(v.get("sample"))
        child_hom = _zygosity(v) == "homozygous"
        both_typed = f_alt is not None and m_alt is not None
        inh = ""
        if child_alt is not True:
            inh = ""
        elif both_typed:
            if not f_alt and not m_alt:
                inh = "de_novo"
            elif f_alt and m_alt:
                inh = "homozygous_recessive" if child_hom else "biparental"
            elif f_alt:
                inh = "paternal"
            elif m_alt:
                inh = "maternal"
        elif f_alt is True:
            inh = "paternal"
        elif m_alt is True:
            inh = "maternal"
        elif fidx is not None or midx is not None:
            inh = "absent_in_parents"
        v["inheritance"] = inh
        n += 1
    return n


def _bounded_text(value,limit):
    return str(value or "").strip()[:limit]


def _sanitize_reflection(value):
    if not isinstance(value,dict) or not isinstance(value.get("candidates"),list):
        return None
    candidates=[]
    for item in value["candidates"][:100]:
        if not isinstance(item,dict):
            continue
        try:
            rank=int(item.get("rank"))
        except (TypeError,ValueError):
            continue
        verdict=str(item.get("verdict") or "").lower()
        if rank<1 or verdict not in {"support","uncertain","refute"}:
            continue
        candidates.append({
            "rank":rank,
            "gene":_bounded_text(item.get("gene"),100),
            "verdict":verdict,
            "reasoning":_bounded_text(item.get("reasoning"),2000),
        })
    if not candidates:
        return None
    return {"candidates":candidates,
            "all_refuted":all(c["verdict"]=="refute" for c in candidates)}


def _sanitize_differential(value,verified_urls):
    if not isinstance(value,dict) or not isinstance(value.get("differential"),list):
        return []
    out=[]
    for item in value["differential"][:20]:
        if not isinstance(item,dict):
            continue
        likelihood=str(item.get("likelihood") or "low").lower()
        if likelihood not in {"high","moderate","low"}:
            likelihood="low"
        raw_genes=item.get("genes")
        raw_variants=item.get("supporting_variants")
        raw_evidence=item.get("evidence")
        genes=[_bounded_text(g,100) for g in (raw_genes if isinstance(raw_genes,list) else [])[:20]
               if isinstance(g,(str,int,float)) and _bounded_text(g,100)]
        variants=[_bounded_text(v,200) for v in (raw_variants if isinstance(raw_variants,list) else [])[:50]
                  if isinstance(v,(str,int,float)) and _bounded_text(v,200)]
        evidence=[]
        for link in (raw_evidence if isinstance(raw_evidence,list) else [])[:30]:
            if not isinstance(link,dict):
                continue
            url=str(link.get("url") or "")
            if url in verified_urls:
                evidence.append({"label":_bounded_text(link.get("label"),200),"url":url})
        disease=_bounded_text(item.get("disease"),300)
        if not disease:
            continue
        out.append({
            "disease":disease,
            "genes":genes,
            "inheritance":_bounded_text(item.get("inheritance"),200),
            "supporting_variants":variants,
            "likelihood":likelihood,
            "rationale":_bounded_text(item.get("rationale"),4000),
            "evidence":evidence,
            "next_steps":_bounded_text(item.get("next_steps"),2000),
        })
    return out


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
            "inheritance": v.get("inheritance"), "father_gt": v.get("father_gt"),
            "mother_gt": v.get("mother_gt"),
            "phenotype_match_disease": v.get("pheno_disease"), "phenotype_shared_terms": v.get("pheno_shared"),
            "combined_score": v.get("combined"),
        }

    k = max(1,min(start_k,len(variants))); rounds = 0; reflection = None; supported = []
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
        reflection = _sanitize_reflection(_llm_json(reflect_sys,user,api_key))
        if not reflection:
            return None
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
    verified_set = {u for u,ok in valid.items() if ok}
    differential=_sanitize_differential(_llm_json(diff_sys,duser,api_key),verified_set)

    return {"differential": differential,
            "reflection": reflection.get("candidates", []),
            "rounds": rounds, "k_considered": k,
            "n_verified_links": len(verified_set)}

# ---------- report ----------
def _differential_html(agentic):
    """Render the DeepRare-style differential + self-reflection block, or '' if absent."""
    if not agentic:
        return ""
    LK={"high":"#c0392b","moderate":"#e67e22","low":"#7f8c8d"}
    cards=[]
    for i,d in enumerate(agentic["differential"],1):
        lk=(d.get("likelihood") or "low").lower()
        links="".join(f'<a href="{html.escape(e["url"])}" target="_blank" rel="noopener noreferrer" style="font-size:11px;margin-right:10px;color:#2980b9">{html.escape(e["label"])} ↗</a>'
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
    heading=(
        '<h2 style="font-size:16px;color:#1e2a38;margin:0 0 4px">'
        '🧬 Diagnóstico diferencial (razonamiento agéntico)</h2>'
        if cards else
        '<h2 style="font-size:16px;color:#1e2a38;margin:0 0 4px">'
        '🧬 Autorreflexión agéntica</h2>'
    )
    return f'<div style="margin:16px 0">{heading}{meta}{"".join(cards)}{ref_html}</div>'

def write_html(variants, patient_hpo, prefix, assembly="GRCh38", agentic=None, warnings=None):
    TC={"Pathogenic":"#c0392b","Likely pathogenic":"#e67e22","Uncertain significance (VUS)":"#7f8c8d",
        "Likely benign":"#27ae60","Benign":"#2ecc71"}
    tier=lambda c:"Benign" if c.split(" |")[0].startswith("Benign") else c.split(" |")[0]
    col=lambda c:TC.get(tier(c),"#7f8c8d")
    def faf(v):
        if v.get("af_status")=="unavailable":
            return "no disponible"
        a=v.get("af")
        return "ausente" if a is None else (f"{a:.2e}" if a<1e-3 else f"{a:.3f}")
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
    def pheno_cell(v):
        score=float(v.get("pheno_score") or 0)
        disease=html.escape(str(v.get("pheno_disease") or "sin enfermedad coincidente"))
        shared=html.escape(str(v.get("pheno_shared") or "sin HPO compartidos"))
        return (
            f'<td style="color:#8e44ad"><b>{score:.2f}</b><br>'
            f'<span style="font-size:10px;color:#34495e">{disease}</span><br>'
            f'<span style="font-size:9px;color:#7f8c8d">{shared}</span></td>'
        )
    chips="".join(f'<span class="hpo">{html.escape(t["label"])} <code>{html.escape(t["hpo_id"])}</code></span>' for t in patient_hpo)
    warning_html=""
    if warnings:
        warning_html=(
            '<div class="note"><b>Advertencias de datos:</b><ul>'
            + "".join(f"<li>{html.escape(str(w))}</li>" for w in warnings)
            + "</ul></div>"
        )
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
                 f'<td>{html.escape(str(v["chrom"]))}:{v["pos"]} {html.escape(v["ref"])}>{html.escape(v["alt"])}</td>'
                 f'<td>{faf(v)}</td><td>{html.escape(str(v.get("clinvar") or ""))} {"*"*int(v.get("stars") or 0)}</td>'
                 f'<td>{ai_cell(v)}</td>'
                 f'<td style="color:{col(v["call"])};font-weight:600">{html.escape(tier(v["call"]))}</td>'
                 + (inh_cell(v) if has_trio else "")
                 + f'<td>{v.get("variant_score",0):.2f}</td>{pheno_cell(v)}'
                 f'<td><b>{v.get("combined",0):.2f}</b></td><td>{html.escape(str(v.get("filter") or ""))}</td></tr>' for v in variants)
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
{warning_html}
{_differential_html(agentic)}
<h2 style="font-size:14px;color:#1e2a38;margin:18px 0 4px">Variantes candidatas priorizadas</h2>
<table><thead><tr><th>#</th><th>Gen</th><th>Variante</th><th>gnomAD</th><th>ClinVar</th><th>IA</th><th>Clase</th>{inh_th}<th>Var</th><th>Encaje HPO</th><th>Comb</th><th>QC</th></tr></thead>
<tbody>{rows}</tbody></table>
<div class="note"><b>Sistema de apoyo a la decision, no diagnostico.</b> VEP + gnomAD + ClinVar (score de variante ACMG-lite) x fenotipos HPO de Open Targets (score fenotipico). Confirmar por metodo ortogonal e interpretar en contexto clinico.</div>
</body></html>"""
    with open(f"{prefix}_report.html","w",encoding="utf-8") as fh:
        fh.write(doc)

# ---------- main ----------
CSV_COLS=["rank","gene","rsid","chrom","pos","ref","alt","consequence","protein","af",
          "af_status","clinvar","stars","pli","loeuf","sift","polyphen","call","variant_score",
          "pheno_score","pheno_direct","pheno_disease","pheno_shared","combined","filter",
          "acmg_tags","esm2_llr","esm2_call","am_pathogenicity","am_call",
          "agentic_evaluated","reflect_verdict","inheritance","father_gt","mother_gt"]


def _variant_evidence_score(tags, call):
    """Absolute 0-1 evidence score; unlike min-max scaling, it is stable across VCFs."""
    codes={t[0] for t in tags}
    if "BA1" in codes:
        return 0.0
    score=0.5
    score+=0.30 if "PVS1" in codes else 0
    score+=0.25 if "PS_ClinVar" in codes else 0
    score+=0.25 if "PS2" in codes else 0
    score+=0.15 if "PM2" in codes else 0
    score+=0.12 if "PVS1_mod" in codes else 0
    score+=0.12 if "LoF_predicted" in codes else 0
    score+=(0.10 if "PP3_AM" in codes else 0.08 if {"PP3","PP3_ESM"}&codes else 0)
    score+=0.06 if "PP5_ClinVar" in codes else 0
    score-=0.30 if "BS_ClinVar" in codes else 0
    score-=0.25 if "BS1" in codes else 0
    score-=(0.10 if "BP4_AM" in codes else 0.08 if {"BP4","BP4_ESM"}&codes else 0)
    score-=0.06 if "BP6_ClinVar" in codes else 0
    if call.split(" |")[0].startswith("Benign"):
        score=min(score,0.2)
    if "ClinVar_conflicting" in codes:
        score=0.5
    elif call.split(" |")[0]=="Uncertain significance (VUS)" and (
            {"BS1","BS_ClinVar"}&codes and {"PVS1","PS_ClinVar","PS2"}&codes):
        score=0.5
    return round(max(0.0,min(score,1.0)),3)


def _rerank_agentic(variants, agentic_result):
    original_to_variant={v["rank"]:v for v in variants}
    verdict_order={"support":0,"uncertain":1,"refute":3}
    variants.sort(key=lambda v:(verdict_order.get(v.get("reflect_verdict"),2),-v["combined"]))
    old_to_new={}
    for new_rank,v in enumerate(variants,1):
        old_to_new[v["rank"]]=new_rank
        v["rank"]=new_rank
    for reflected in agentic_result.get("reflection",[]):
        old_rank=reflected.get("rank")
        if old_rank in original_to_variant:
            reflected["original_rank"]=old_rank
            reflected["rank"]=old_to_new[old_rank]


def run_pipeline(vcf_path, sample="SAMPLE", hpo="", clinical_note_text=None, assembly="auto",
                 use_esm=False, use_am=False, agentic=False, reflect_k=8,
                 father_vcf=None, mother_vcf=None, email=None,
                 anthropic_key=None, progress=None, max_variants=None,
                 reviewed_hpo_expansion=None, expand_hpo_terms=False):
    """Annotate a VCF and prioritize variants. Reusable entry point for CLI and web server.

    Args:
        vcf_path: path to the input VCF.
        sample: sample ID to select in a multi-sample VCF; single-sample VCFs use their only sample.
        hpo: comma/newline-separated HPO IDs or free-text terms.
        clinical_note_text: raw clinical-note text (LLM -> HPO), or None.
        assembly: 'auto', 'GRCh38', or 'GRCh37' (auto reads VCF metadata).
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

    _reset_llm_state()
    _check_cancelled()
    _prog(0, 0, "Leyendo y filtrando VCF…")
    if assembly in (None,"auto","AUTO"):
        assembly=detect_vcf_assembly(vcf_path)
    variants=parse_vcf(
        vcf_path,sample=sample,allow_empty=False,max_variants=max_variants,called_only=True
    )
    _prog(0, len(variants), f"VCF leído: {len(variants)} variantes")
    warnings=[]

    # patient HPO profile — from clinical note (LLM) and/or explicit hpo string
    patient_hpo=[]
    if clinical_note_text:
        patient_hpo=extract_hpo_from_note(clinical_note_text, anthropic_key)
        llm_errors=list(getattr(_LLM_STATE,"errors",[]))
        if not patient_hpo:
            warnings.append("Clinical note did not yield HPO terms; check LLM configuration or provide HPO IDs")
        if llm_errors:
            warnings.extend(f"Clinical-note extraction warning: {error}" for error in llm_errors)
        providers=sorted(getattr(_LLM_STATE,"providers",set()))
        provider_note=f" con {', '.join(providers)}" if providers else ""
        _prog(0, len(variants), f"{len(patient_hpo)} fenotipos extraídos de la nota clínica{provider_note}")
    tokens=[t for t in re.split(r"[,\n]", hpo or "") if t.strip()]
    if tokens:
        existing={t["hpo_id"] for t in patient_hpo}
        resolved=resolve_hpo(tokens)
        for t in resolved:
            if t["hpo_id"] not in existing: patient_hpo.append(t); existing.add(t["hpo_id"])
        if len(resolved)<len(tokens):
            warnings.append(f"{len(tokens)-len(resolved)} HPO term(s) could not be resolved")
    patient_ids={t["hpo_id"] for t in patient_hpo}
    if patient_ids and reviewed_hpo_expansion is not None:
        reviewed_tokens=[
            token.strip() for token in re.split(r"[,\n]",reviewed_hpo_expansion)
            if token.strip()
        ]
        reviewed_terms=resolve_hpo(reviewed_tokens)
        patient_full={term["hpo_id"] for term in reviewed_terms}|patient_ids
        hpo_expansion_available=len(reviewed_terms)==len(reviewed_tokens)
        if not hpo_expansion_available:
            warnings.append("Some reviewed expanded HPO terms could not be resolved")
    elif patient_ids and expand_hpo_terms:
        patient_full,hpo_expansion_available=expand_hpo(patient_ids,return_status=True)
        if not hpo_expansion_available:
            warnings.append("HPO ancestor expansion was partially unavailable; phenotype matching may be incomplete")
    elif patient_ids:
        # Default: match on the directly observed HPO terms only. Ancestor expansion is an
        # opt-in step (reviewed_hpo_expansion from the UI, or expand_hpo_terms from the CLI),
        # not something done automatically at the start of every analysis.
        patient_full=set(patient_ids)
    else:
        patient_full=set()
    if len(patient_ids)==1:
        only_hpo=next(iter(patient_ids))
        warnings.append(
            f"Phenotype profile contains only one direct HPO term ({only_hpo}); "
            "phenotype ranking has low specificity. Add other observed clinical "
            "features before interpreting gene order."
        )
    if patient_hpo:
        if len(patient_full)>len(patient_ids):
            _prog(0, len(variants), f"{len(patient_ids)} términos HPO (expandidos a {len(patient_full)})")
        else:
            _prog(0, len(variants), f"{len(patient_ids)} términos HPO directos (sin expansión)")

    cons_cache={}
    pheno_cache={}
    AF_COMMON=_prefilter_af_cutoff()
    n_workers=_annotation_workers()
    _deadline=getattr(_REQUEST_STATE,"deadline",None)
    _cancel=getattr(_REQUEST_STATE,"cancel_event",None)

    def _annotate_frequency(v):
        """Phase 1 (cheap): consequence + population allele frequency only."""
        set_request_deadline(_deadline); set_request_cancel_event(_cancel)
        _request_timeout()
        ve=vep(v, assembly); v.update(consequence=ve.get("most_severe"), gene=ve.get("gene"),
                            gene_id=ve.get("gene_id"), protein=ve.get("amino_acids"),
                            sift=ve.get("sift"), polyphen=ve.get("polyphen"))
        vep_ok=bool(ve.get("annotation_available"))
        if assembly=="GRCh38":
            v["af"],af_ok=gnomad_af(v,return_status=True)
        else:
            v["af"],af_ok=ve.get("gnomad_af"),vep_ok
        v["af_status"]=("observed" if v["af"] is not None else "absent" if af_ok else "unavailable")
        v["_af_ok"]=af_ok
        return {"vep_ok":vep_ok,"af_ok":af_ok}

    def _default_deep_fields(v):
        v["pli"]=None; v["loeuf"]=None
        v["clinvar"]=None; v["stars"]=None; v["conditions"]=None
        v["esm2_llr"]=None; v["esm2_call"]=""
        v["am_pathogenicity"]=None; v["am_call"]=""

    def _finalize(v, do_pheno):
        """Build ACMG tags from the evidence on v, classify, and optionally score phenotype."""
        extra_tags=[]
        if v["esm2_call"]=="deleterious": extra_tags.append(("PP3_ESM",f"ESM-2 LLR {v['esm2_llr']} (deleterious)"))
        elif v["esm2_call"]=="tolerated": extra_tags.append(("BP4_ESM",f"ESM-2 LLR {v['esm2_llr']} (tolerated)"))
        if v["am_call"]=="deleterious": extra_tags.append(("PP3_AM",f"AlphaMissense {v['am_pathogenicity']} (likely_pathogenic)"))
        elif v["am_call"]=="tolerated": extra_tags.append(("BP4_AM",f"AlphaMissense {v['am_pathogenicity']} (likely_benign)"))
        call,tags=classify(v["af"],v["consequence"],v["pli"],v["loeuf"],v["clinvar"],v["stars"],
                           v["sift"],v["polyphen"],af_available=v.get("_af_ok",True),extra_tags=extra_tags)
        v["call"]=call
        v["_tags"]=tags
        if do_pheno and patient_ids:
            gene_id=v.get("gene_id")
            if gene_id not in pheno_cache:
                pheno_cache[gene_id]=gene_pheno_score(
                    gene_id,patient_ids,patient_full,return_status=True
                )
            ps,pd_,pm_,pdis,psh,pheno_available=pheno_cache[gene_id]
            v.update(pheno_score=ps,pheno_direct=pd_,pheno_matched=pm_,pheno_disease=pdis,pheno_shared=psh)
            return bool(gene_id) and not pheno_available
        v.update(pheno_score=0.0,pheno_direct=0,pheno_matched=0,pheno_disease="",pheno_shared="")
        return False

    def _annotate_deep(v):
        """Phase 2 (expensive): constraint, ClinVar, ESM/AM, classification, phenotype."""
        set_request_deadline(_deadline); set_request_cancel_event(_cancel)
        flags={"clinvar_ok":True,"constraint_ok":True,"am_attempt":False,"am_ok":True,
               "esm_attempt":False,"esm_ok":True,"pheno_failed":False}
        _request_timeout()
        con=gnomad_constraint(v["gene"], cons_cache) if v.get("gene") else {"pli":None,"loeuf":None}
        v["pli"],v["loeuf"]=con.get("pli"),con.get("loeuf")
        if v.get("gene") and not con.get("available"):
            flags["constraint_ok"]=False
        cv=clinvar(v,email,assembly); v.update(clinvar=cv.get("significance"), stars=cv.get("stars"),
                                      conditions=cv.get("conditions"), protein=v.get("protein") or cv.get("protein_change"))
        if not cv.get("available"):
            flags["clinvar_ok"]=False
        v["esm2_llr"]=None; v["esm2_call"]=""
        if use_esm and v["consequence"]=="missense_variant":
            flags["esm_attempt"]=True
            ctx=protein_context(v.get("rsid"),v["chrom"],v["pos"],v["ref"],v["alt"],
                                v.get("gene"),assembly=assembly)
            if ctx:
                v["esm2_llr"]=esm_score_missense(ctx["seq"], ctx["pos"], ctx["wt"], ctx["mut"])
                v["esm2_call"]=esm_call(v["esm2_llr"])
                v["protein"]=v.get("protein") or ctx["protein_change"]
            if v["esm2_llr"] is None:
                flags["esm_ok"]=False
        v["am_pathogenicity"]=None; v["am_call"]=""
        if use_am and v["consequence"]=="missense_variant":
            flags["am_attempt"]=True
            am=alphamissense_score(v["chrom"], v["pos"], v["ref"], v["alt"], assembly=assembly, gene=v.get("gene"))
            if am:
                v["am_pathogenicity"]=am["am_pathogenicity"]; v["am_call"]=am_call(am)
            else:
                flags["am_ok"]=False
        flags["pheno_failed"]=_finalize(v, do_pheno=True)
        return flags

    # Phase 1: consequence + population allele frequency for every called variant (concurrent).
    total=len(variants)
    freq_flags=_map_concurrent(
        variants, _annotate_frequency, n_workers,
        lambda done: (done % 20 == 0 or done == total)
                     and _prog(done, total, f"Anotando frecuencias: {done}/{total}"))
    vep_failures=sum(1 for f in freq_flags if not f["vep_ok"])
    gnomad_failures=sum(1 for f in freq_flags if not f["af_ok"])

    if vep_failures==total:
        raise RuntimeError("Ensembl VEP unavailable for every variant; analysis aborted")

    # Prefilter: variants at population AF >= 5% are Benign by BA1 (see classify/_call_from_tags,
    # where BA1 short-circuits every other criterion). Skip the expensive evidence layers
    # (constraint/ClinVar/ESM/AlphaMissense/phenotype) for them and classify directly — the
    # verdict is identical, and the deep budget is spent only on plausible candidates.
    common=[v for v in variants if v.get("af") is not None and v["af"]>=AF_COMMON]
    candidates=[v for v in variants if not (v.get("af") is not None and v["af"]>=AF_COMMON)]
    for v in common:
        _default_deep_fields(v)
        _finalize(v, do_pheno=False)
    _prog(0, total,
          f"Prefiltrado: {len(candidates)} candidatas para análisis profundo "
          f"({len(common)} comunes AF≥{AF_COMMON:.0%} clasificadas como benignas)")

    # Phase 2: deep evidence only for the plausible candidates (concurrent).
    n_cand=len(candidates)
    deep_flags=_map_concurrent(
        candidates, _annotate_deep, n_workers,
        lambda done: (done % 10 == 0 or done == n_cand)
                     and _prog(done, n_cand, f"Análisis profundo de candidatas: {done}/{n_cand}"))
    clinvar_failures=sum(1 for f in deep_flags if not f["clinvar_ok"])
    constraint_failures=sum(1 for f in deep_flags if not f["constraint_ok"])
    am_attempts=sum(1 for f in deep_flags if f["am_attempt"])
    am_failures=sum(1 for f in deep_flags if f["am_attempt"] and not f["am_ok"])
    esm_attempts=sum(1 for f in deep_flags if f["esm_attempt"])
    esm_failures=sum(1 for f in deep_flags if f["esm_attempt"] and not f["esm_ok"])
    pheno_failures=sum(1 for f in deep_flags if f["pheno_failed"])
    for v in variants:
        v.pop("_af_ok",None)
    if vep_failures:
        warnings.append(f"Ensembl VEP unavailable for {vep_failures} variant(s)")
    if gnomad_failures:
        warnings.append(f"gnomAD frequency unavailable for {gnomad_failures} variant(s); PM2 was not assigned")
    if clinvar_failures:
        warnings.append(f"ClinVar unavailable for {clinvar_failures} variant(s)")
    if constraint_failures:
        warnings.append(f"gnomAD constraint unavailable for {constraint_failures} variant(s)")
    if am_attempts and am_failures:
        warnings.append(f"AlphaMissense unavailable for {am_failures} of {am_attempts} missense variant(s)")
    if esm_attempts and esm_failures:
        warnings.append(f"ESM-2 unavailable for {esm_failures} of {esm_attempts} missense variant(s)")
    if pheno_failures:
        warnings.append(f"Open Targets phenotype data unavailable for {pheno_failures} variant(s)")

    # trio inheritance (optional) — assigns v["inheritance"]; de novo adds PS2 (supporting strong)
    if father_vcf or mother_vcf:
        n_trio=trio_inheritance(variants, father_vcf, mother_vcf, progress)
        for v in variants:
            if v.get("inheritance")=="de_novo":
                v["_tags"].append(("PS2","de novo (absent in both explicitly genotyped parents)"))
                v["call"]=_call_from_tags(v["_tags"],v.get("clinvar"))
        _prog(len(variants), len(variants), f"Herencia analizada en {n_trio} variantes")

    for v in variants:
        tags=v.pop("_tags")
        v["acmg_tags"]=",".join(t[0] for t in tags)
        v["evidence"]="; ".join(f"{t[0]}: {t[1]}" for t in tags)
        v["variant_score"]=_variant_evidence_score(tags,v["call"])
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
        if agentic_result is None:
            warnings.append("Agentic reasoning was requested but unavailable; numeric ranking was retained")
        # flag which rows the LLM self-reflection actually examined (only the top-K window)
        for v in variants:
            v["agentic_evaluated"]="yes" if v.get("reflect_verdict") else "no"
        if agentic_result:
            _rerank_agentic(variants,agentic_result)

    return {"variants":variants, "patient_hpo":patient_hpo, "csv_cols":CSV_COLS,
            "n_input":len(variants), "assembly":assembly, "agentic":agentic_result,
            "warnings":warnings,
            "llm_providers":sorted(getattr(_LLM_STATE,"providers",set()))}

def write_outputs(result, out_prefix):
    """Write <prefix>_annotated.csv and <prefix>_report.html from a run_pipeline() result."""
    import csv
    from pathlib import Path
    Path(out_prefix).parent.mkdir(parents=True,exist_ok=True)
    with open(f"{out_prefix}_annotated.csv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh, fieldnames=result["csv_cols"], extrasaction="ignore"); w.writeheader()
        for v in result["variants"]: w.writerow(v)
    write_html(result["variants"], result["patient_hpo"], out_prefix, result.get("assembly","GRCh38"),
               result.get("agentic"),result.get("warnings"))

def main():
    ap=argparse.ArgumentParser(description="raredx VCF annotation & phenotype prioritization")
    ap.add_argument("vcf")
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--hpo", default="", help="Comma-separated HPO IDs or free-text terms, or @file.txt")
    ap.add_argument("--out-prefix", default="raredx_out")
    ap.add_argument("--email", default=None, help="Contact email for NCBI E-utilities (optional)")
    ap.add_argument("--esm", action="store_true", help="Score missense with ESM-2 (needs torch+fair-esm)")
    ap.add_argument("--alphamissense", action="store_true", help="Score missense with AlphaMissense (precomputed, via Ensembl VEP; no GPU)")
    ap.add_argument("--assembly", default="auto", choices=["auto","GRCh38","GRCh37"], help="Genome build (default: detect from VCF metadata)")
    ap.add_argument("--clinical-note", default=None, help="Path to free-text clinical note; extracts HPO via configured LLM")
    ap.add_argument("--agentic", action="store_true", help="Agentic layer: LLM self-reflection over candidates + traceable differential diagnosis with link verification (needs ANTHROPIC_API_KEY)")
    ap.add_argument("--reflect-k", type=int, default=8, help="Agentic layer: number of top candidates the LLM self-reflection examines (default 8; window widens only if all are refuted)")
    ap.add_argument("--father", default=None, help="Father VCF for trio analysis (de novo detection → ACMG PS2)")
    ap.add_argument("--mother", default=None, help="Mother VCF for trio analysis (de novo detection → ACMG PS2)")
    ap.add_argument("--anthropic-key", default=None, help="Anthropic API key (else uses ANTHROPIC_API_KEY env)")
    ap.add_argument("--expand-hpo", action="store_true", help="Also match on HPO ancestor terms (opt-in; by default only the directly observed terms are used)")
    a=ap.parse_args()

    hpo_arg=a.hpo
    if hpo_arg.startswith("@"): hpo_arg=open(hpo_arg[1:]).read()
    note=open(a.clinical_note).read() if a.clinical_note else None

    def cli_progress(done, total, msg): print(f"[raredx] {msg}", file=sys.stderr)
    result=run_pipeline(a.vcf, sample=a.sample, hpo=hpo_arg, clinical_note_text=note,
                        assembly=a.assembly, use_esm=a.esm, use_am=a.alphamissense,
                        agentic=a.agentic, reflect_k=a.reflect_k,
                        father_vcf=a.father, mother_vcf=a.mother, email=a.email,
                        anthropic_key=a.anthropic_key, progress=cli_progress,
                        expand_hpo_terms=a.expand_hpo)
    write_outputs(result, a.out_prefix)
    print(f"[raredx] wrote {a.out_prefix}_annotated.csv and {a.out_prefix}_report.html", file=sys.stderr)

if __name__=="__main__":
    main()
