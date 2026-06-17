"""
GeneReg - Gene Regulation Predictor v2
Local ChIP-seq peak data + Motif + Literature weighted scoring
Key improvement: Uses local ENCODE peak data instead of slow/unreliable API queries
Deployed on 43 server port 8004
"""

import os, re, json, math, requests, pickle, bisect
import xml.etree.ElementTree as ET
from typing import Optional, List, Dict, Tuple
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn

# ==================== Config ====================
DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_MODEL = "deepseek-chat"
MYGENE_API = "https://mygene.info/v3"
JASPAR_API = "https://jaspar.elixir.no/api/v1"
ENCODE_API = "https://www.encodeproject.org"
ENSEMBL_API = "https://rest.ensembl.org"
PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "bioinfo@example.com"

# Local peak data paths
PEAK_DIR = "/www/wwwroot/gene_reg/encode_peaks"
PEAK_INDEX_PATH = os.path.join(PEAK_DIR, "_peak_index.pkl")

# ==================== Gene Info ====================
class GeneInfoFetcher:
    def fetch(self, gene_symbol: str) -> Optional[Dict]:
        try:
            resp = requests.get(f"{MYGENE_API}/query", params={
                "q": gene_symbol, "species": "human",
                "fields": "symbol,name,genomic_pos,ensembl.gene", "size": 1
            }, timeout=15)
            resp.raise_for_status()
            hits = resp.json().get("hits", [])
            if not hits: return None
            hit = hits[0]
            gp = hit.get("genomic_pos", {})
            if isinstance(gp, list): gp = gp[0]
            chr_name = str(gp.get("chr","")).replace("chr","")
            start, end = gp.get("start",0), gp.get("end",0)
            strand = gp.get("strand",1)
            tss = start if strand == 1 else end
            ens = hit.get("ensembl", {})
            if isinstance(ens, list): ens = ens[0]
            eid = ens.get("gene","") if isinstance(ens, dict) else ""
            return {"symbol":hit.get("symbol",gene_symbol),"name":hit.get("name",""),
                    "chr":chr_name,"start":start,"end":end,"tss":tss,"strand":strand,"ensembl_id":eid}
        except Exception as e:
            print(f"GeneInfo error: {e}"); return None

# ==================== ChIP-seq (Local Peak Data) ====================
class ChIPSeqAnalyzer:
    def __init__(self):
        self.index = None
        self.tf_dirs = {}  # tf_name -> set of available TF dir names
        self._load_index()

    def _load_index(self):
        """Load the pre-built peak index for fast queries"""
        if os.path.exists(PEAK_INDEX_PATH):
            try:
                with open(PEAK_INDEX_PATH, "rb") as f:
                    self.index = pickle.load(f)
                print(f"Loaded peak index: {sum(len(v) for v in self.index.values())} peaks across {len(self.index)} chromosomes")
            except Exception as e:
                print(f"Failed to load peak index: {e}")
                self.index = {}
        
        # Build TF dir mapping
        if os.path.isdir(PEAK_DIR):
            for d in os.listdir(PEAK_DIR):
                dp = os.path.join(PEAK_DIR, d)
                if os.path.isdir(dp) and not d.startswith("_"):
                    self.tf_dirs[d.lower()] = d

    def _find_tf_dir(self, receptor: str) -> Optional[str]:
        """Find the local TF directory name for a receptor"""
        r_lower = receptor.lower()
        # Exact match
        if r_lower in self.tf_dirs:
            return self.tf_dirs[r_lower]
        # Try with -human suffix (ENCODE naming)
        for suffix in ["-human", "-mouse", "-rat"]:
            key = r_lower + suffix
            if key in self.tf_dirs:
                return self.tf_dirs[key]
        # Partial match
        for key, val in self.tf_dirs.items():
            if r_lower in key or key in r_lower:
                return val
        return None

    def _search_local_peaks(self, gene_info: Dict, receptor: str, window: int = 100000) -> List[Dict]:
        """Search local peak data for binding events near a gene"""
        if not self.index or not gene_info:
            return []
        
        chr_name = gene_info.get("chr", "")
        tss = gene_info.get("tss", 0)
        if not chr_name or not tss:
            return []
        
        # Get peaks on this chromosome
        chr_peaks = self.index.get(chr_name, [])
        if not chr_peaks:
            return []
        
        # Find the TF dir for this receptor
        tf_dir_name = self._find_tf_dir(receptor)
        
        # Binary search for peaks in window around TSS
        peaks_found = []
        left = bisect.bisect_left(chr_peaks, (tss - window,))
        right = bisect.bisect_right(chr_peaks, (tss + window,))
        
        for i in range(left, right):
            start, end, tf_name = chr_peaks[i]
            mid = (start + end) // 2
            dist = abs(mid - tss)
            
            # If we have a specific TF, only return its peaks
            if tf_dir_name:
                if tf_name.lower() != tf_dir_name.lower():
                    continue
            
            peaks_found.append({
                "chr": chr_name,
                "start": start,
                "end": end,
                "distance_to_tss": dist,
                "tf_name": tf_name,
                "signal": 0  # Not stored in index
            })
        
        # Sort by distance
        peaks_found.sort(key=lambda p: p["distance_to_tss"])
        return peaks_found

    def _count_tf_peaks(self, tf_dir_name: str) -> int:
        """Count total peak files for a TF"""
        tf_dir = os.path.join(PEAK_DIR, tf_dir_name)
        if not os.path.isdir(tf_dir):
            return 0
        return len([f for f in os.listdir(tf_dir) if f.endswith(".bed")])

    def _search_encode_api(self, receptor: str) -> Tuple[List[Dict], int]:
        """Fallback: search ENCODE API for experiment info (no peak data, just experiment count)"""
        experiments = []
        exp_count = 0
        try:
            for target_param in [{"target.label": receptor}, {"target.name": receptor}]:
                resp = requests.get(f"{ENCODE_API}/search/", params={
                    "type": "Experiment",
                    "assay_title": "TF ChIP-seq",
                    **target_param,
                    "status": "released",
                    "format": "json",
                    "limit": 20,
                }, headers={"Accept": "application/json"}, timeout=15)
                if resp.status_code == 200:
                    data = resp.json()
                    exp_count = data.get("total", 0)
                    for exp in data.get("@graph", [])[:10]:
                        bio = exp.get("biosample_ontology", {})
                        experiments.append({
                            "accession": exp.get("accession", ""),
                            "biosample": bio.get("term_name", "") if isinstance(bio, dict) else str(bio),
                            "title": exp.get("title", "")
                        })
                    if exp_count > 0:
                        break
        except Exception as e:
            print(f"ENCODE experiment search error: {e}")
        
        return experiments, exp_count

    def analyze(self, gene: str, receptor: str, gene_info: Optional[Dict] = None) -> Dict:
        peaks = []
        experiments = []
        score = 0
        source = "Local ENCODE Peaks"
        tf_dir_name = self._find_tf_dir(receptor)

        # 1. Primary: Search local peak data
        if gene_info and gene_info.get("chr") and gene_info.get("tss") and self.index:
            # Search with expanding windows: 5kb, 20kb, 50kb, 100kb
            for window, window_label in [(5000, "promoter"), (20000, "proximal"), (50000, "enhancer"), (100000, "distal")]:
                found = self._search_local_peaks(gene_info, receptor, window)
                if found:
                    peaks = found
                    break
            
            # If no peaks for specific TF, also check ALL TFs (to see if region is active)
            if not peaks:
                all_peaks = self._search_local_peaks(gene_info, "ALL_TFS", 50000)
                # Filter: only report peaks, don't score them for this receptor
                # This tells us the region has regulatory activity but not from this receptor

        # 2. Score based on local peaks
        if peaks:
            # Filter to only this receptor's peaks for scoring
            receptor_peaks = peaks
            if tf_dir_name:
                receptor_peaks = [p for p in peaks if p["tf_name"].lower() == tf_dir_name.lower()]
            
            if receptor_peaks:
                dists = [p["distance_to_tss"] for p in receptor_peaks]
                min_dist = min(dists)
                
                # Distance-based scoring (most important factor)
                if min_dist < 1000:
                    dist_score = 35
                elif min_dist < 2000:
                    dist_score = 32
                elif min_dist < 5000:
                    dist_score = 28
                elif min_dist < 10000:
                    dist_score = 22
                elif min_dist < 50000:
                    dist_score = 15
                else:
                    dist_score = 8
                
                # Peak count scoring (more peaks = more evidence)
                count_score = min(35, len(receptor_peaks) * 7)
                
                # Peak quality: closer peaks count more
                close_peaks = len([d for d in dists if d < 5000])
                quality_score = min(30, close_peaks * 10)
                
                score = dist_score + count_score + quality_score
                peaks = receptor_peaks  # Only show this receptor's peaks
                source = f"Local ENCODE Peaks ({tf_dir_name})"
            else:
                # Peaks exist in region but NOT for this receptor
                # This is actually negative evidence
                score = 5  # Very low score - region is active but receptor doesn't bind
                source = "Local ENCODE Peaks (no binding for this receptor)"
        
        # 3. Get experiment info from ENCODE API (lightweight, for display only)
        experiments, exp_count = self._search_encode_api(receptor)
        
        # 4. If local data exists for this TF but no peaks near gene, give partial score
        if not peaks and tf_dir_name:
            file_count = self._count_tf_peaks(tf_dir_name)
            if file_count > 0:
                # TF has ChIP-seq data but no binding at this gene locus
                # Moderate evidence that there's no direct regulation
                score = max(score, 3)  # Small score for having data
                source = f"Local ENCODE Peaks ({tf_dir_name}, no binding near gene)"
        
        # 5. If no local data, fall back to experiment count
        if not tf_dir_name and exp_count > 0:
            score = min(15, exp_count * 2)
            source = "ENCODE API (no local peak data)"
        
        min_dist = None
        if peaks:
            dists = [p["distance_to_tss"] for p in peaks]
            min_dist = min(dists)

        return {
            "score": min(100, score),
            "peaks_found": len(peaks),
            "peaks": sorted(peaks, key=lambda p: p.get("distance_to_tss", 999999))[:10],
            "experiments_found": exp_count,
            "experiments": experiments[:5],
            "min_distance_to_tss": min_dist,
            "source": source,
            "local_data_available": tf_dir_name is not None
        }

# ==================== Motif ====================
class MotifAnalyzer:
    def analyze(self, gene: str, receptor: str, gene_info: Optional[Dict] = None) -> Dict:
        motif_info = {}
        promoter_matches = []
        score = 0

        # JASPAR query
        try:
            resp = requests.get(f"{JASPAR_API}/matrix/", params={
                "search": receptor, "tax_group": "vertebrates", "collection": "CORE"
            }, timeout=15)
            resp.raise_for_status()
            results = resp.json().get("results", [])
            if results:
                best = results[0]
                mid = best.get("matrix_id","")
                motif_info = {"matrix_id":mid,"name":best.get("name",receptor),
                    "collection":best.get("collection",""),
                    "tf_family":best.get("tf_family",""),
                    "class":best.get("class",[""])[0] if isinstance(best.get("class"),list) else best.get("class",""),
                    "pubmed_ids":best.get("pubmed_ids",[])}
                # Get PWM length
                pwm_resp = requests.get(f"{JASPAR_API}/matrix/{mid}", timeout=15)
                if pwm_resp.status_code == 200:
                    pfm = pwm_resp.json().get("pfm",{})
                    motif_info["motif_length"] = len(pfm.get("A",[])) if isinstance(pfm.get("A"),list) else 0
        except Exception as e:
            print(f"JASPAR error: {e}")

        # Get promoter sequence
        promoter_seq = ""
        if gene_info and gene_info.get("ensembl_id"):
            try:
                seq_resp = requests.get(
                    f"{ENSEMBL_API}/sequence/id/{gene_info['ensembl_id']}",
                    params={"type":"genomic","expand_5prime":2000,"expand_3prime":500},
                    headers={"Content-Type":"application/json"}, timeout=20)
                if seq_resp.status_code == 200:
                    promoter_seq = seq_resp.json().get("seq","").upper()
            except Exception as e:
                print(f"Ensembl seq error: {e}")

        # Scan known response elements
        if promoter_seq:
            for mname, mseq in self._get_known_re(receptor):
                rev = self._revcomp(mseq)
                fwd = [m.start() for m in re.finditer(f"(?={mseq})", promoter_seq)]
                rev_matches = [m.start() for m in re.finditer(f"(?={rev})", promoter_seq)]
                for pos in fwd:
                    promoter_matches.append({"motif_name":mname,"motif_seq":mseq,"position":pos-2000,"strand":"+"})
                for pos in rev_matches:
                    promoter_matches.append({"motif_name":mname,"motif_seq":mseq,"position":pos-2000,"strand":"-"})

        # Scoring
        if motif_info:
            has_score = 20
            quality = 15 if motif_info.get("collection")=="CORE" else 10 if motif_info.get("collection")=="PBM" else 5
            if promoter_matches:
                matched_lens = set()
                for pm in promoter_matches:
                    matched_lens.add(len(pm.get("motif_seq", "")))
                max_match_len = max(matched_lens) if matched_lens else 0
                if max_match_len >= 12:
                    match_score = min(65, len(promoter_matches) * 20)
                elif max_match_len >= 8:
                    match_score = min(55, len(promoter_matches) * 15)
                else:
                    match_score = min(25, len(promoter_matches) * 6)
            else:
                match_score = 0
            score = has_score + quality + match_score

        return {"score":min(100,score),"motif_found":bool(motif_info),"motif_info":motif_info,
                "promoter_matches":promoter_matches,"promoter_length":len(promoter_seq),"source":"JASPAR + Ensembl"}

    def _get_known_re(self, receptor: str) -> List[tuple]:
        db = {
            "ESR1":[("ERE","AGGTCA"),("ERE_half","AGGTCA")],
            "ESR2":[("ERE","AGGTCA")],
            "AR":[("ARE","AGAACA")],
            "PGR":[("PRE","AGAACA")],
            "GR":[("GRE","AGAACA")],
            "NR3C1":[("GRE","AGAACA")],
            "TR":[("TRE","AGGTCA")],
            "THRA":[("TRE","AGGTCA")],
            "THRB":[("TRE","AGGTCA")],
            "RARA":[("RARE","AGGTCA")],
            "RXRA":[("RXRE","AGGTCA")],
            "PPARA":[("PPRE","AGGTCA")],
            "PPARG":[("PPRE","AGGTCA")],
            "VDR":[("VDRE","AGGTCA")],
            "FXR":[("FXRE","AGGTCA")],
            "LXR":[("LXRE","AGGTCA")],
        }
        for k,v in db.items():
            if k.upper() == receptor.upper(): return v
        return [("NR_half_site","AGGTCA")]

    def _revcomp(self, seq: str) -> str:
        c = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
        return "".join(c.get(b,"N") for b in reversed(seq.upper()))

# ==================== Literature ====================
class LiteratureAnalyzer:
    def analyze(self, gene: str, receptor: str) -> Dict:
        articles = []
        pmids = set()
        for q in [f'"{gene}" AND "{receptor}" AND (regulation OR target OR binding)', f'"{gene}" AND "{receptor}"']:
            try:
                resp = requests.get(f"{PUBMED_EUTILS}/esearch.fcgi", params={
                    "db":"pubmed","term":q,"retmax":10,"sort":"relevance","retmode":"json","email":NCBI_EMAIL
                }, timeout=15)
                pmids.update(resp.json().get("esearchresult",{}).get("idlist",[]))
            except: pass

        if pmids:
            try:
                fr = requests.get(f"{PUBMED_EUTILS}/efetch.fcgi", params={
                    "db":"pubmed","id":",".join(list(pmids)[:10]),"rettype":"xml","retmode":"xml","email":NCBI_EMAIL
                }, timeout=15)
                root = ET.fromstring(fr.text)
                for ae in root.findall(".//PubmedArticle")[:10]:
                    mc = ae.find("MedlineCitation")
                    if mc is None: continue
                    ad = mc.find("Article")
                    if ad is None: continue
                    t = ad.find("ArticleTitle")
                    pm = mc.find("PMID")
                    yr = ""
                    je = ad.find("Journal")
                    if je is not None:
                        pd = je.find("PubDate")
                        if pd is not None:
                            ye = pd.find("Year")
                            if ye is not None and ye.text: yr = ye.text
                            else:
                                md = pd.find("MedlineDate")
                                if md is not None and md.text: yr = md.text[:4]
                    articles.append({"pmid":pm.text if pm is not None else "","title":(t.text if t is not None else "") or "","year":yr or "N/A"})
            except: pass

        n = len(pmids)
        score = 90 if n>=20 else 75 if n>=10 else 60 if n>=5 else 50 if n>=3 else 35 if n>=1 else 0
        return {"score":score,"articles_found":n,"articles":articles[:10],"source":"PubMed"}

# ==================== AI ====================
class AIRegulationAnalyzer:
    def analyze(self, gene: str, receptor: str, chip: Dict, motif: Dict, lit: Dict) -> Dict:
        cs = f"Score {chip['score']}/100, {chip['peaks_found']} peaks"
        if chip.get('min_distance_to_tss') is not None: cs += f", nearest {chip['min_distance_to_tss']}bp from TSS"
        if chip.get('local_data_available'): cs += " (local ENCODE data)"
        ms = f"Score {motif['score']}/100"
        if motif.get('motif_found'): ms += f", JASPAR {motif['motif_info'].get('matrix_id','')}"
        if motif.get('promoter_matches'): ms += f", {len(motif['promoter_matches'])} promoter matches"
        ls = f"Score {lit['score']}/100, {lit['articles_found']} articles"
        refs = "\n".join([f"- PMID {a['pmid']}: {a['title']} ({a['year']})" for a in lit.get("articles",[])[:5]])

        prompt = f"""Analyze if receptor/TF "{receptor}" directly regulates gene "{gene}".

Evidence:
1. ChIP-seq: {cs}
2. Motif: {ms}
3. Literature: {ls}

References:
{refs or "(No references found)"}

Return JSON:
```json
{{"mechanism":"Brief mechanism","evidence_summary":"2-3 sentence summary","confidence":"Strong/Moderate/Weak/None","recommendations":["Next step 1","Next step 2"]}}
```"""
        try:
            r = requests.post(DEEPSEEK_API_URL, headers={
                "Content-Type":"application/json","Authorization":f"Bearer {DEEPSEEK_API_KEY}"
            }, json={"model":DEEPSEEK_MODEL,"messages":[
                {"role":"system","content":"Expert bioinformatician. Return JSON only."},
                {"role":"user","content":prompt}
            ],"temperature":0.3,"max_tokens":2000}, timeout=120)
            r.raise_for_status()
            res = r.json()["choices"][0]["message"]["content"]
            m = re.search(r'```json\s*(.*?)\s*```', res, re.DOTALL)
            return json.loads(m.group(1)) if m else json.loads(res)
        except Exception as e:
            return {"mechanism":f"AI unavailable: {e}","evidence_summary":"Manual review recommended",
                    "confidence":"Unknown","recommendations":["Check databases manually"]}

# ==================== FastAPI ====================
app = FastAPI(title="GeneReg - Gene Regulation Predictor v2")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

gf = GeneInfoFetcher()
ca = ChIPSeqAnalyzer()  # Now uses local peak data!
ma = MotifAnalyzer()
la = LiteratureAnalyzer()
ai_a = AIRegulationAnalyzer()

class AnalysisRequest(BaseModel):
    gene: str
    receptor: str

def _do_analysis(gene: str, receptor: str) -> dict:
    if not gene.strip() or not receptor.strip():
        raise HTTPException(400, "Gene and receptor are required")
    gi = gf.fetch(gene)
    cr = ca.analyze(gene, receptor, gi)
    mr = ma.analyze(gene, receptor, gi)
    lr = la.analyze(gene, receptor)
    
    # Weighted scoring with improved weights
    # ChIP-seq is the gold standard, give it more weight
    w = {"chip":0.45,"motif":0.30,"literature":0.25}
    fs = w["chip"]*cr["score"] + w["motif"]*mr["score"] + w["literature"]*lr["score"]
    
    if fs>=70: cl,lvl = "Strong Evidence","strong"
    elif fs>=40: cl,lvl = "Moderate Evidence","moderate"
    elif fs>=20: cl,lvl = "Weak Evidence","weak"
    else: cl,lvl = "Insufficient Evidence","none"
    
    ar = ai_a.analyze(gene, receptor, cr, mr, lr)

    # AI confidence bonus with dampened range
    ai_bonus = {"Strong": 10, "Moderate": 3, "Weak": -3, "None": -10}.get(ar.get("confidence", ""), 0)
    adjusted = max(0, min(100, fs + ai_bonus))

    if adjusted>=70: cl,lvl = "Strong Evidence","strong"
    elif adjusted>=40: cl,lvl = "Moderate Evidence","moderate"
    elif adjusted>=20: cl,lvl = "Weak Evidence","weak"
    else: cl,lvl = "Insufficient Evidence","none"

    return {"gene":gene,"receptor":receptor,"gene_info":gi,"scores":{"chip_seq":cr,"motif":mr,"literature":lr},
            "weights":w,"raw_score":round(fs,1),"ai_bonus":ai_bonus,"final_score":round(adjusted,1),
            "classification":cl,"level":lvl,"ai_analysis":ar}

@app.post("/api/analyze")
def analyze_post(req: AnalysisRequest): return _do_analysis(req.gene, req.receptor)

@app.get("/api/analyze")
def analyze_get(gene: str = Query(...), receptor: str = Query(...)): return _do_analysis(gene, receptor)

@app.get("/api/health")
def health(): return {"status":"ok","service":"GeneReg","version":"2.0","local_peaks":ca.index is not None,
                      "tf_count":len(ca.tf_dirs),"peak_count":sum(len(v) for v in (ca.index or {}).values()) if ca.index else 0}

@app.get("/")
async def root(): return HTMLResponse(content=HTML_PAGE)

# ==================== HTML (same as v1) ====================
HTML_PAGE = """
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GeneReg - Gene Regulation Predictor</title>
<script src="https://cdn.tailwindcss.com"></script>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<link href="https://cdn.jsdelivr.net/npm/font-awesome@4.7.0/css/font-awesome.min.css" rel="stylesheet">
<style>
*{font-family:'Inter',sans-serif}
.card{background:#fff;border:1px solid #e5e7eb;box-shadow:0 1px 3px rgba(0,0,0,.06);border-radius:12px}
.card-hover:hover{box-shadow:0 4px 12px rgba(0,0,0,.1)}
.fade-in{animation:fadeIn .5s ease-out}
@keyframes fadeIn{from{opacity:0;transform:translateY(10px)}to{opacity:1;transform:translateY(0)}}
.loader{border:3px solid #e5e7eb;border-top:3px solid #7c3aed;border-radius:50%;width:40px;height:40px;animation:spin 1s linear infinite}
@keyframes spin{0%{transform:rotate(0)}100%{transform:rotate(360deg)}}
nav{background:linear-gradient(135deg,#0f172a,#1e3a8a);padding:1rem 1.5rem;box-shadow:0 4px 12px rgba(0,0,0,.1);position:sticky;top:0;z-index:100}
.nav-container{max-width:1400px;margin:0 auto;display:flex;justify-content:space-between;align-items:center}
.back-btn{background:rgba(255,255,255,.1);color:#fff;border:1px solid rgba(255,255,255,.2);padding:.6rem 1.2rem;border-radius:8px;text-decoration:none;display:inline-flex;align-items:center;gap:.5rem;transition:all .2s}
.back-btn:hover{background:rgba(255,255,255,.2)}
.score-ring{width:120px;height:120px;border-radius:50%;display:flex;align-items:center;justify-content:center;font-size:2rem;font-weight:700}
.level-strong{background:rgba(16,185,129,.1);color:#059669;border:3px solid #059669}
.level-moderate{background:rgba(245,158,11,.1);color:#d97706;border:3px solid #d97706}
.level-weak{background:rgba(239,68,68,.1);color:#dc2626;border:3px solid #ef4444}
.level-none{background:rgba(107,114,128,.1);color:#6b7280;border:3px solid #9ca3af}
.bar-track{background:#f1f5f9;border-radius:6px;height:8px;width:100%}
.bar-fill{border-radius:6px;height:8px;transition:width .8s ease}
</style>
</head>
<body class="bg-gray-50 min-h-screen text-gray-800">

<nav>
<div class="nav-container">
<a href="https://43.99.62.219/" class="back-btn"><i class="fa fa-arrow-left"></i> Back to Animal Hormone DB</a>
</div>
</nav>

<main class="max-w-6xl mx-auto px-6 py-8">

<div class="text-center mb-8">
<h1 class="text-3xl font-bold text-gray-800 mb-2">GeneReg <span class="text-violet-500">v2</span></h1>
<p class="text-gray-500">Predict gene-receptor regulatory relationships via ChIP-seq, motif analysis & literature evidence</p>
</div>

<div class="card p-6 mb-8">
<div class="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
<div>
<label class="block text-sm font-medium text-gray-600 mb-1">Gene Symbol</label>
<input id="geneInput" type="text" placeholder="e.g. TFF1, GREB1, PGR" class="w-full border border-gray-200 rounded-lg px-4 py-3 focus:outline-none focus:border-violet-400 focus:ring-2 focus:ring-violet-100 transition-all">
</div>
<div>
<label class="block text-sm font-medium text-gray-600 mb-1">Receptor / TF Name</label>
<input id="receptorInput" type="text" placeholder="e.g. ESR1, AR, PGR, GR" class="w-full border border-gray-200 rounded-lg px-4 py-3 focus:outline-none focus:border-violet-400 focus:ring-2 focus:ring-violet-100 transition-all">
</div>
</div>
<div class="flex gap-3 items-center">
<button onclick="doAnalyze()" id="analyzeBtn" class="bg-gradient-to-r from-violet-600 to-blue-600 hover:from-violet-500 hover:to-blue-500 text-white px-8 py-3 rounded-xl font-medium transition-all flex items-center gap-2">
<i class="fa fa-search"></i> Analyze
</button>
<div class="flex gap-2 flex-wrap">
<span class="text-xs text-gray-400">Examples:</span>
<button onclick="quickFill('TFF1','ESR1')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">TFF1 + ESR1</button>
<button onclick="quickFill('KLK3','AR')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">KLK3 + AR</button>
<button onclick="quickFill('PGR','ESR1')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">PGR + ESR1</button>
</div>
</div>
</div>

<div id="loading" class="hidden flex flex-col items-center justify-center py-20">
<div class="loader mb-4"></div>
<p class="text-gray-500" id="loadingText">Analyzing...</p>
</div>

<div id="results" class="hidden space-y-6">

<div id="scoreSection" class="card p-6 fade-in">
<div class="flex flex-col md:flex-row items-center gap-6">
<div id="scoreRing" class="score-ring"></div>
<div class="flex-1">
<h2 id="classLabel" class="text-2xl font-bold mb-1"></h2>
<p id="classDesc" class="text-gray-500 text-sm mb-3"></p>
<div class="grid grid-cols-3 gap-4 text-center">
<div><p class="text-xs text-gray-400">ChIP-seq</p><p id="chipBar" class="text-lg font-bold text-violet-600">-</p></div>
<div><p class="text-xs text-gray-400">Motif</p><p id="motifBar" class="text-lg font-bold text-blue-600">-</p></div>
<div><p class="text-xs text-gray-400">Literature</p><p id="litBar" class="text-lg font-bold text-amber-600">-</p></div>
</div>
</div>
<div>
<button onclick="downloadReport()" class="text-xs text-gray-500 border border-gray-200 bg-white hover:border-violet-300 px-3 py-2 rounded-lg transition-all flex items-center gap-1">
<i class="fa fa-download"></i> Download CSV
</button>
</div>
</div>
</div>

<div class="card p-6 fade-in">
<h3 class="text-lg font-semibold text-gray-700 mb-4">Evidence Breakdown</h3>
<div class="space-y-4" id="barsContainer"></div>
</div>

<div id="aiSection" class="card p-6 fade-in">
<h3 class="text-lg font-semibold text-violet-600 mb-3"><i class="fa fa-robot mr-2"></i>AI Analysis</h3>
<div id="aiContent"></div>
</div>

<div id="chipSection" class="card p-6 fade-in">
<h3 class="text-lg font-semibold text-emerald-600 mb-3"><i class="fa fa-bar-chart mr-2"></i>ChIP-seq Binding Evidence</h3>
<div id="chipContent"></div>
</div>

<div id="motifDetailSection" class="card p-6 fade-in">
<h3 class="text-lg font-semibold text-blue-600 mb-3"><i class="fa fa-dna mr-2"></i>Motif Prediction</h3>
<div id="motifContent"></div>
</div>

<div id="litSection" class="card p-6 fade-in">
<h3 class="text-lg font-semibold text-amber-600 mb-3"><i class="fa fa-book mr-2"></i>Literature Evidence</h3>
<div id="litContent"></div>
</div>

</div>

<div id="error" class="hidden card p-8 text-center">
<p class="text-red-500 text-lg mb-2">Analysis Failed</p>
<p id="errorText" class="text-gray-500"></p>
</div>

</main>

<footer class="border-t border-gray-200 py-4 mt-12 bg-white">
<div class="max-w-6xl mx-auto px-6 flex items-center justify-between text-xs text-gray-400">
<span>GeneReg v2.0 · Local ENCODE Peaks + JASPAR + PubMed + DeepSeek AI</span>
<span>Inner Mongolia University of Science and Technology</span>
</div>
</footer>

<script>
let lastData = null;

function quickFill(g,r){document.getElementById('geneInput').value=g;document.getElementById('receptorInput').value=r;doAnalyze()}

async function doAnalyze(){
    const gene=document.getElementById('geneInput').value.trim();
    const receptor=document.getElementById('receptorInput').value.trim();
    if(!gene||!receptor)return;
    const btn=document.getElementById('analyzeBtn');
    const loading=document.getElementById('loading');
    const results=document.getElementById('results');
    const error=document.getElementById('error');
    btn.disabled=true;btn.classList.add('opacity-50');
    loading.classList.remove('hidden');results.classList.add('hidden');error.classList.add('hidden');

    try{
        document.getElementById('loadingText').textContent='Searching local ChIP-seq peaks, motifs & literature...';
        const resp=await fetch(`/api/analyze?gene=${encodeURIComponent(gene)}&receptor=${encodeURIComponent(receptor)}`);
        const data=await resp.json();
        if(!resp.ok)throw new Error(data.detail||'Analysis failed');
        lastData=data;
        loading.classList.add('hidden');results.classList.remove('hidden');
        renderResults(data);
    }catch(e){
        loading.classList.add('hidden');error.classList.remove('hidden');
        document.getElementById('errorText').textContent=e.message;
    }finally{btn.disabled=false;btn.classList.remove('opacity-50')}
}

function renderResults(d){
    const s=d.scores;
    const ring=document.getElementById('scoreRing');
    ring.className='score-ring level-'+d.level;
    ring.textContent=d.final_score;
    document.getElementById('classLabel').textContent=d.classification;
    const descs={strong:'High confidence of direct regulation',moderate:'Likely regulatory relationship',weak:'Possible regulation, needs validation',none:'Insufficient evidence for regulation'};
    document.getElementById('classDesc').textContent=descs[d.level]||'';
    document.getElementById('chipBar').textContent=s.chip_seq.score;
    document.getElementById('motifBar').textContent=s.motif.score;
    document.getElementById('litBar').textContent=s.literature.score;

    const bars=[
        {name:'ChIP-seq Binding',score:s.chip_seq.score,color:'#8b5cf6',detail:`${s.chip_seq.peaks_found} peaks, ${s.chip_seq.experiments_found} ENCODE experiments`+(s.chip_seq.min_distance_to_tss!=null?`, nearest ${s.chip_seq.min_distance_to_tss}bp`:'')+(s.chip_seq.local_data_available?' (local)':'')},
        {name:'Motif Prediction',score:s.motif.score,color:'#3b82f6',detail:s.motif.motif_found?`JASPAR ${s.motif.motif_info.matrix_id}${s.motif.promoter_matches.length?', '+s.motif.promoter_matches.length+' promoter matches':''}`:'No motif found'},
        {name:'Literature Evidence',score:s.literature.score,color:'#f59e0b',detail:`${s.literature.articles_found} PubMed articles`}
    ];
    document.getElementById('barsContainer').innerHTML=bars.map(b=>`
        <div><div class="flex justify-between mb-1"><span class="text-sm font-medium text-gray-700">${b.name}</span><span class="text-sm font-bold" style="color:${b.color}">${b.score}/100</span></div>
        <div class="bar-track"><div class="bar-fill" style="width:${b.score}%;background:${b.color}"></div></div>
        <p class="text-xs text-gray-400 mt-1">${b.detail}</p></div>
    `).join('');

    const ai=d.ai_analysis;
    document.getElementById('aiContent').innerHTML=`
        <p class="text-gray-600 text-sm mb-2"><strong>Mechanism:</strong> ${ai.mechanism||'N/A'}</p>
        <p class="text-gray-600 text-sm mb-2"><strong>Summary:</strong> ${ai.evidence_summary||'N/A'}</p>
        <p class="text-sm mb-2"><strong>Confidence:</strong> <span class="font-semibold ${ai.confidence==='Strong'?'text-emerald-600':ai.confidence==='Moderate'?'text-amber-600':'text-gray-500'}">${ai.confidence||'N/A'}</span></p>
        ${ai.recommendations&&ai.recommendations.length?`<p class="text-xs text-gray-400 mb-1">Recommendations:</p><ul class="text-sm text-gray-600 list-disc list-inside">${ai.recommendations.map(r=>'<li>'+r+'</li>').join('')}</ul>`:''}
    `;

    const cp=s.chip_seq;
    document.getElementById('chipContent').innerHTML=`
        <p class="text-sm text-gray-600 mb-2">Source: ${cp.source} | Peaks found: <strong>${cp.peaks_found}</strong>${cp.experiments_found?` | Experiments: <strong>${cp.experiments_found}</strong>`:''}${cp.min_distance_to_tss!=null?' | Nearest peak: <strong>'+cp.min_distance_to_tss+'bp</strong> from TSS':''}</p>
        ${cp.experiments&&cp.experiments.length?`<p class="text-xs font-medium text-gray-700 mb-1">ENCODE ChIP-seq Experiments for ${d.receptor}:</p>
        <div class="space-y-1 mb-3">${cp.experiments.map(e=>`<div class="text-xs bg-gray-50 rounded px-3 py-1.5"><a href="https://www.encodeproject.com/experiments/${e.accession}/" target="_blank" class="text-violet-600 hover:text-violet-800 font-medium">${e.accession}</a> <span class="text-gray-500">- ${e.biosample}</span></div>`).join('')}</div>`:''}
        ${cp.peaks.length?`<p class="text-xs font-medium text-gray-700 mb-1">Binding Peaks near ${d.gene}:</p>
        <div class="overflow-x-auto"><table class="w-full text-xs"><thead><tr class="border-b border-gray-200"><th class="text-left py-1 px-2">Chr</th><th class="text-left py-1 px-2">Start</th><th class="text-left py-1 px-2">End</th><th class="text-left py-1 px-2">Distance</th><th class="text-left py-1 px-2">TF</th></tr></thead>
        <tbody>${cp.peaks.map(p=>`<tr class="border-b border-gray-100 hover:bg-gray-50"><td class="py-1 px-2">${p.chr}</td><td class="py-1 px-2">${p.start.toLocaleString()}</td><td class="py-1 px-2">${p.end.toLocaleString()}</td><td class="py-1 px-2 font-medium ${p.distance_to_tss<5000?'text-emerald-600':p.distance_to_tss<50000?'text-amber-600':'text-gray-500'}">${p.distance_to_tss? p.distance_to_tss.toLocaleString()+'bp':''}</td><td class="py-1 px-2 text-violet-600">${p.tf_name||''}</td></tr>`).join('')}</tbody></table></div>`:'<p class="text-sm text-gray-400">No binding peaks found near this gene from this receptor.</p>'}
    `;

    const mt=s.motif;
    document.getElementById('motifContent').innerHTML=`
        <p class="text-sm text-gray-600 mb-2">Source: ${mt.source} | Motif found: <strong>${mt.motif_found?'Yes':'No'}</strong>${mt.promoter_length?' | Promoter scanned: '+mt.promoter_length+'bp':''}</p>
        ${mt.motif_found?`<div class="bg-gray-50 rounded-lg p-4 mb-3">
            <p class="text-sm"><strong>Matrix ID:</strong> ${mt.motif_info.matrix_id} | <strong>Name:</strong> ${mt.motif_info.name}</p>
            <p class="text-xs text-gray-500">Collection: ${mt.motif_info.collection} | Family: ${mt.motif_info.tf_family||'N/A'} | Length: ${mt.motif_info.motif_length||'?'}bp</p>
        </div>`:''}
        ${mt.promoter_matches.length?`<p class="text-sm font-medium text-gray-700 mb-1">Promoter Motif Matches (${mt.promoter_matches.length}):</p>
        <div class="overflow-x-auto"><table class="w-full text-xs"><thead><tr class="border-b border-gray-200"><th class="text-left py-1 px-2">Motif</th><th class="text-left py-1 px-2">Sequence</th><th class="text-left py-1 px-2">Position (rel. TSS)</th><th class="text-left py-1 px-2">Strand</th></tr></thead>
        <tbody>${mt.promoter_matches.map(m=>`<tr class="border-b border-gray-100"><td class="py-1 px-2 font-mono">${m.motif_name}</td><td class="py-1 px-2 font-mono text-violet-600">${m.motif_seq}</td><td class="py-1 px-2">${m.position>0?'+':''}${m.position}</td><td class="py-1 px-2">${m.strand}</td></tr>`).join('')}</tbody></table></div>`:'<p class="text-sm text-gray-400">No motif matches found in promoter region.</p>'}
    `;

    const lt=s.literature;
    document.getElementById('litContent').innerHTML=`
        <p class="text-sm text-gray-600 mb-3">Source: ${lt.source} | Articles found: <strong>${lt.articles_found}</strong></p>
        ${lt.articles.length?lt.articles.map(a=>`<div class="border-l-2 border-amber-300 pl-3 py-1 mb-2"><a href="https://pubmed.ncbi.nlm.nih.gov/${a.pmid}/" target="_blank" class="text-sm text-violet-600 hover:text-violet-800">${a.title}</a><p class="text-xs text-gray-400">PMID: ${a.pmid} | ${a.year}</p></div>`).join(''):'<p class="text-sm text-gray-400">No literature found for this gene-receptor pair.</p>'}
    `;
}

function downloadReport(){
    if(!lastData)return;
    const d=lastData,s=d.scores;
    let csv='\\uFEFF';
    csv+='Category,Score,Detail\\n';
    csv+=`ChIP-seq Binding,${s.chip_seq.score},"${s.chip_seq.peaks_found} peaks${s.chip_seq.min_distance_to_tss!=null?', nearest '+s.chip_seq.min_distance_to_tss+'bp':''}"\\n`;
    csv+=`Motif Prediction,${s.motif.score},"${s.motif.motif_found?'JASPAR '+s.motif.motif_info.matrix_id:'No motif found'}"\\n`;
    csv+=`Literature Evidence,${s.literature.score},"${s.literature.articles_found} articles"\\n`;
    csv+=`\\nFinal Score,${d.final_score},"${d.classification}"\\n`;
    csv+=`\\nGene,${d.gene}\\nReceptor,${d.receptor}\\n`;
    if(d.ai_analysis){
        csv+=`\\nAI Mechanism,"${(d.ai_analysis.mechanism||'').replace(/"/g,'""')}"\\n`;
        csv+=`AI Confidence,${d.ai_analysis.confidence||'N/A'}\\n`;
    }
    const blob=new Blob([csv],{type:'text/csv;charset=utf-8;'});
    const url=URL.createObjectURL(blob);
    const a=document.createElement('a');
    a.href=url;a.download=`GeneReg_${d.gene}_${d.receptor}_${new Date().toISOString().slice(0,10)}.csv`;
    document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);
}
</script>
</body>
</html>
"""

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8004)
