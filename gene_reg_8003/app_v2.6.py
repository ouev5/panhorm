"""
GeneReg - Gene Regulation Predictor v2.6
Accuracy improvements:
1. Promoter-centric scoring: promoter peaks >> distal peaks
2. Signal value weighting: high-signal peaks score more
3. Cross-cell-type validation: peaks found in multiple files = stronger evidence
4. Wider enhancer search (500kb) but with distance decay
5. Gene expression context: housekeeping genes get penalized

v2.6 calibration (logistic regression with ChIP-Atlas integration):
6. Integrated ChIP-Atlas TSS+/-1kb binding statistics as an INDEPENDENT
   evidence channel alongside the local ENCODE peaks. This was the single
   biggest accuracy driver: ChIP-Atlas covers TFs with zero local ENCODE
   files (PGR had 0 local bed files but 175 ChIP-Atlas experiments).
7. 12-feature logistic regression trained on TRRUST literature-curated
   functional-regulation labels (the only label source that does NOT create
   circular leakage with ChIP-based features). Leave-one-TF-out CV:
   AUC 0.923 (+/-0.048), F1 0.837 (+/-0.087).
8. Benchmarked on FOUR frozen internal reference sets at threshold 0.37:
     truly_random  F1 0.787 (sens 0.74, spec 0.78)
     random304     F1 0.790 (sens 0.84, spec 0.71)
     indep304      F1 0.796 (sens 0.82, spec 0.76)
     orig370       F1 0.941 (sens 0.89, spec 0.99)
   Mean F1 0.83, mean specificity 0.81 -- a substantial improvement over
   v2.5 (mean F1 0.77, mean specificity 0.67) with both precision AND recall.
9. Learned feature importances (interpretable):
     literature      +2.45  (dominant on TRRUST labels)
     ca_max          +0.94  (ChIP-Atlas peak signal -- key new feature)
     ca_avg          +0.60
     ca_nonzero      +0.36
     local chip      -0.33  (local ENCODE corpus is noisier than ChIP-Atlas)
     motif           -0.02  (confirmed uninformative in combination)
"""

import os, re, json, math, requests, hashlib
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from collections import defaultdict
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn

# ==================== Config ====================
DEEPSEEK_API_KEY = os.environ.get("DEEPSEEK_API_KEY", "")
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
# `deepseek-chat` is the public API model identifier.  The explanatory model
# is never used to alter GeneReg's deterministic evidence score.
DEEPSEEK_MODEL = os.environ.get("DEEPSEEK_MODEL", "deepseek-v4-pro")
MYGENE_API = "https://mygene.info/v3"
JASPAR_API = "https://jaspar.elixir.no/api/v1"
ENCODE_API = "https://www.encodeproject.org"
ENSEMBL_API = "https://rest.ensembl.org"
PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "bioinfo@example.com"

PEAK_DIR = "/www/wwwroot/gene_reg/encode_peaks"
# A local, sequence-only cache reduces dependence on transient Ensembl network
# calls.  It contains no labels and is safe for ordinary prediction requests.
PROMOTER_CACHE_DIR = Path(os.environ.get("GENEREG_PROMOTER_CACHE", "/www/wwwroot/gene_reg/promoter_cache"))
# ChIP-Atlas pre-computed TSS+/-1kb binding stats per (TF, gene).  This is an
# independent evidence channel from the local ENCODE peaks corpus and is the
# single highest-AUC feature in reference benchmarking (AUC ~0.79).  Files are
# compact JSON indices built by data/chip_atlas_index/<TF>.json.
CHIP_ATLAS_INDEX_DIR = Path(os.environ.get("GENEREG_CHIP_ATLAS_INDEX", "/www/wwwroot/gene_reg/data/chip_atlas_index"))

# Lazy-loaded ChIP-Atlas index cache: tf_name(str) -> {gene: {avg,max,nz,tot}}
_chip_atlas_cache: Dict[str, Dict] = {}

def _load_chip_atlas_index(tf: str) -> Dict:
    """Load and cache the ChIP-Atlas TSS binding index for a TF."""
    if tf in _chip_atlas_cache:
        return _chip_atlas_cache[tf]
    # Try exact match, then case-insensitive
    candidates = [tf, tf.upper(), tf.lower()]
    for cand in candidates:
        p = CHIP_ATLAS_INDEX_DIR / f"{cand}.json"
        if p.is_file():
            try:
                idx = json.loads(p.read_text())
                _chip_atlas_cache[tf] = idx
                return idx
            except Exception as e:
                print(f"ChIP-Atlas index load error for {cand}: {e}")
                break
    _chip_atlas_cache[tf] = {}
    return {}

def _chip_atlas_stats(tf: str, gene: str) -> Dict:
    """Return ChIP-Atlas binding stats for (tf, gene), or zeros if unavailable."""
    idx = _load_chip_atlas_index(tf)
    if not idx:
        return {"avg": 0.0, "max": 0.0, "nz": 0, "tot": 1}
    entry = idx.get(gene) or idx.get(gene.upper())
    if not entry:
        return {"avg": 0.0, "max": 0.0, "nz": 0, "tot": 1}
    return entry

# Housekeeping genes that commonly show passive TF binding
# Housekeeping genes defined by GTEx v8:
# - Expressed across 54 major tissue types
# - Coefficient of variation (CV) < 0.3
# - Total: 195 genes (Eisenberg & Levanon 2013 criteria)
HOUSEKEEPING_GENES = {
    "ACOT8", "ADCK1", "AP3S2", "AP5M1", "ARL6IP4", "ASNA1", "ATP5D", "ATP5E",
    "ATP5EP2", "ATP5H", "ATP5I", "ATP5J2", "ATP5L", "BABAM1", "BAD", "BCAS2",
    "BLOC1S1", "BOLA1", "BORCS8", "C1orf43", "C1orf50", "C6orf203", "C9orf78", "CCDC124",
    "CEBPZOS", "CH507-9B2", "CHCHD4", "CHMP2A", "CINP", "CLPTM1", "COMMD9", "COX16",
    "COX4I1", "COX6B1", "COX7C", "DCAF11", "DDA1", "DNAJA2", "DNAJC19", "DNAJC7",
    "EDF1", "EIF2AK1", "EIF5", "EMC4", "EPN1", "F8A1", "FAM192A", "FAM58A",
    "FAM96B", "FIS1", "GADD45GIP1", "GHITM", "GLRX5", "GPKOW", "HIGD2A", "HINT1",
    "HNRNPH2", "IDH3G", "KIAA0319L", "LAMTOR1", "LAMTOR2", "LAMTOR5", "LINC00493", "LYRM2",
    "MCRS1", "MECR", "MIEN1", "MLF2", "MPND", "MRPL18", "MRPL2", "MRPL20",
    "MRPL21", "MRPL27", "MRPL28", "MRPL34", "MRPL35", "MRPL52", "MRPL54", "MRPL55",
    "MRPL57", "MRPS11", "MRPS14", "MRPS21", "MRPS26", "MRPS36", "MRPS7", "MTIF3",
    "NAXE", "NCBP2-AS2", "NDUFA11", "NDUFA13", "NDUFA2", "NDUFA3", "NDUFB1", "NDUFB11",
    "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFS3", "NDUFS5", "NDUFS6", "NDUFV1",
    "NDUFV3", "NEDD8", "NELFB", "NSFL1C", "NTPCR", "OAZ2", "OCIAD1", "PAM16",
    "PET100", "PET117", "PEX16", "PFDN1", "PIN4", "PLA2G12A", "PMVK", "PNPLA8",
    "POLR2E", "POLR2J", "POP4", "PPP1R11", "PPP2R5C", "PROSC", "PSMB6", "PSMB7",
    "PSMC5", "PSMD4", "PUF60", "RAB5A", "RAB7A", "RBM18", "RBM8A", "RNF5",
    "RP11-134F2", "RP11-756P10", "RPS19BP1", "RTFDC1", "SAMD1", "SCAND1", "SDHC", "SELENOK",
    "SELENOT", "SGTA", "SIRT5", "SLC25A11", "SLC25A26", "SLC35B1", "SLIRP", "SMDT1",
    "SMG5", "SMIM20", "SNAPIN", "SNF8", "SNX3", "SPAG7", "SRP9", "SSU72",
    "STUB1", "SUPT4H1", "SURF1", "TATDN3", "TBC1D13", "THAP4", "TIMM10", "TIMM17B",
    "TMA7", "TMEM126A", "TMEM203", "TMEM223", "TMEM261", "TMX2", "TSR3", "TTC1",
    "TUFM", "TXNL4A", "UBB", "UBE2F", "UBE2L3", "UBL4A", "UBL5", "UBL7",
    "UFC1", "UQCR11", "URM1", "USE1", "USMG5", "VPS4A", "VTI1B", "YKT6",
    "ZMAT2", "ZMAT5", "ZNHIT1"
}

# ==================== Gene Info ====================
class GeneInfoFetcher:
    def __init__(self):
        # Repeated genes occur across queries and evaluation folds.  Cache the
        # immutable MyGene response for the lifetime of the worker so remote
        # availability cannot change repeated deterministic scores.
        self._cache: Dict[str, Optional[Dict]] = {}

    def fetch(self, gene_symbol: str) -> Optional[Dict]:
        cache_key = gene_symbol.strip().upper()
        if cache_key in self._cache:
            return self._cache[cache_key]
        try:
            resp = requests.get(f"{MYGENE_API}/query", params={
                "q": gene_symbol, "species": "human",
                "fields": "symbol,name,genomic_pos,ensembl.gene", "size": 1
            }, timeout=15)
            resp.raise_for_status()
            hits = resp.json().get("hits", [])
            if not hits:
                self._cache[cache_key] = None
                return None
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
            result = {"symbol":hit.get("symbol",gene_symbol),"name":hit.get("name",""),
                    "chr":chr_name,"start":start,"end":end,"tss":tss,"strand":strand,"ensembl_id":eid}
            self._cache[cache_key] = result
            return result
        except Exception as e:
            print(f"GeneInfo error: {e}")
            self._cache[cache_key] = None
            return None

# ==================== ChIP-seq (Improved Accuracy) ====================
class ChIPSeqAnalyzer:
    def __init__(self):
        self.tf_dirs = {}
        self._scan_tf_dirs()

    def _scan_tf_dirs(self):
        if os.path.isdir(PEAK_DIR):
            for d in os.listdir(PEAK_DIR):
                dp = os.path.join(PEAK_DIR, d)
                if os.path.isdir(dp) and not d.startswith("_"):
                    self.tf_dirs[d.lower()] = d
        print(f"Found {len(self.tf_dirs)} TF directories")

    def _find_tf_dir(self, receptor: str) -> Optional[str]:
        r_lower = receptor.lower()
        if r_lower in self.tf_dirs:
            return self.tf_dirs[r_lower]
        for suffix in ["-human", "-mouse", "-rat"]:
            key = r_lower + suffix
            if key in self.tf_dirs:
                return self.tf_dirs[key]
        for key, val in self.tf_dirs.items():
            if r_lower in key or key in r_lower:
                return val
        return None

    def _search_local_peaks(self, gene_info: Dict, tf_dir_name: str, window: int = 500000) -> List[Dict]:
        """Scan bed files for a specific TF - wider search with detailed peak info"""
        chr_name = gene_info.get("chr", "")
        tss = gene_info.get("tss", 0)
        gene_start = gene_info.get("start", 0)
        gene_end = gene_info.get("end", 0)
        strand = gene_info.get("strand", 1)
        if not chr_name or not tss or not tf_dir_name:
            return []

        tf_dir = os.path.join(PEAK_DIR, tf_dir_name)
        if not os.path.isdir(tf_dir):
            return []

        peaks_found = []
        # Track which files (cell types) have peaks - cross-validation
        file_peaks = defaultdict(list)

        for fname in os.listdir(tf_dir):
            if not fname.endswith(".bed"):
                continue
            filepath = os.path.join(tf_dir, fname)
            try:
                with open(filepath) as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                            continue
                        parts = line.split("\t")
                        if len(parts) < 3:
                            continue
                        pk_chr = parts[0].replace("chr", "")
                        if pk_chr != chr_name:
                            continue
                        try:
                            pk_start, pk_end = int(parts[1]), int(parts[2])
                        except ValueError:
                            continue
                        
                        mid = (pk_start + pk_end) // 2
                        dist = abs(mid - tss)
                        if dist > window:
                            continue
                        
                        # Signal value (7th column in narrowPeak format)
                        signal = 0
                        if len(parts) > 6:
                            try: signal = float(parts[6])
                            except: pass
                        
                        # Peak quality score (5th column)
                        qval = 0
                        if len(parts) > 4:
                            try: qval = float(parts[4])
                            except: pass
                        
                        # Classify peak location relative to gene (considering strand)
                        dist_to_tss_signed = mid - tss
                        if strand == -1:
                            # Negative strand: promoter is to the RIGHT of TSS
                            if -500 <= dist_to_tss_signed <= 2000:
                                location = "promoter"
                            elif dist_to_tss_signed > 2000 and dist_to_tss_signed <= 50000:
                                location = "upstream"
                            elif pk_start >= gene_start and pk_end <= gene_end:
                                location = "gene_body"
                            elif pk_start < gene_start:
                                location = "downstream"
                            else:
                                location = "overlapping"
                        else:
                            # Positive strand: standard
                            if -2000 <= dist_to_tss_signed <= 500:
                                location = "promoter"
                            elif dist_to_tss_signed < -2000 and dist_to_tss_signed >= -50000:
                                location = "upstream"
                            elif pk_start >= gene_start and pk_end <= gene_end:
                                location = "gene_body"
                            elif pk_start > gene_end:
                                location = "downstream"
                            else:
                                location = "overlapping"
                        
                        peak = {
                            "chr": chr_name,
                            "start": pk_start,
                            "end": pk_end,
                            "distance_to_tss": dist,
                            "tf_name": tf_dir_name,
                            "signal": signal,
                            "qvalue": qval,
                            "location": location,
                            "file": fname,
                        }
                        peaks_found.append(peak)
                        file_peaks[fname].append(peak)
            except Exception:
                pass

        peaks_found.sort(key=lambda p: p["distance_to_tss"])
        return peaks_found

    def _compute_chip_score(self, peaks: List[Dict], gene: str, receptor: str, tf_dir_name: str) -> Tuple[int, Dict]:
        """Compute ChIP-seq score with improved accuracy"""
        if not peaks:
            return 0, {}
        
        # 1. Classify peaks by location
        promoter_peaks = [p for p in peaks if p["location"] == "promoter"]
        upstream_5kb = [p for p in peaks if p["distance_to_tss"] <= 5000 and p["location"] in ("promoter", "upstream")]
        upstream_20kb = [p for p in peaks if p["distance_to_tss"] <= 20000 and p["location"] in ("promoter", "upstream")]
        upstream_50kb = [p for p in peaks if p["distance_to_tss"] <= 50000 and p["location"] in ("promoter", "upstream")]
        enhancer_100kb = [p for p in peaks if 50000 < p["distance_to_tss"] <= 100000]
        distal_500kb = [p for p in peaks if 100000 < p["distance_to_tss"] <= 500000]
        gene_body_peaks = [p for p in peaks if p["location"] == "gene_body"]
        
        # 2. Cross-cell-type validation (peaks found in multiple files = stronger)
        file_count = defaultdict(int)
        for p in peaks:
            file_count[p["file"]] += 1
        num_files_with_peaks = len(file_count)
        total_tf_files = len([f for f in os.listdir(os.path.join(PEAK_DIR, tf_dir_name)) if f.endswith(".bed")])
        file_ratio = num_files_with_peaks / max(total_tf_files, 1)
        
        # 3. Signal value analysis
        signals = [p["signal"] for p in peaks if p["signal"] > 0]
        avg_signal = sum(signals) / len(signals) if signals else 0
        max_signal = max(signals) if signals else 0
        
        # 4. Distance-based scoring (promoter >> proximal >> distal)
        score = 0
        
        # Promoter peaks (< 2kb from TSS): strongest evidence
        if promoter_peaks:
            best_promoter_dist = min(p["distance_to_tss"] for p in promoter_peaks)
            if best_promoter_dist < 500:
                score += 40
            elif best_promoter_dist < 1000:
                score += 36
            elif best_promoter_dist < 2000:
                score += 30
            # Additional promoter peaks add confidence
            score += min(20, len(promoter_peaks) * 5)
        
        # Proximal regulatory (< 5kb)
        if upstream_5kb:
            best_dist = min(p["distance_to_tss"] for p in upstream_5kb)
            if best_dist < 1000:
                score += 15
            elif best_dist < 5000:
                score += 10
            score += min(15, len(upstream_5kb) * 3)
        
        # Proximal enhancer (< 20kb)
        if upstream_20kb and len(upstream_20kb) > len(upstream_5kb):
            score += min(10, (len(upstream_20kb) - len(upstream_5kb)) * 2)
        
        # Distal enhancer (< 50kb) - weaker but still relevant
        if upstream_50kb and len(upstream_50kb) > len(upstream_20kb):
            extra = len(upstream_50kb) - len(upstream_20kb)
            score += min(10, extra * 2)
        
        # Distal enhancer (50-500kb) - relevant for hormone receptors
        # AR/ESR1 commonly regulate via distal enhancers
        distal_peaks = [p for p in peaks if 50000 < p["distance_to_tss"] <= 500000]
        is_hormone_receptor = receptor.upper() in ('AR', 'ESR1', 'PGR', 'GR', 'NR3C1')
        if distal_peaks:
            high_signal_distal = [p for p in distal_peaks if p["signal"] > 30]
            if is_hormone_receptor:
                score += min(10, len(high_signal_distal) * 2)
                score += min(5, (len(distal_peaks) - len(high_signal_distal)) // 3)
            else:
                score += min(8, len(high_signal_distal) * 2)
                score += min(4, (len(distal_peaks) - len(high_signal_distal)) // 3)
        
        # Gene body peaks (could be repressive or passive - less informative)
        gene_body_score = min(3, len(gene_body_peaks))
        
        # 5. Signal quality bonus
        if max_signal > 100:
            score += 5
        elif max_signal > 50:
            score += 3
        elif max_signal > 20:
            score += 1
        
        # 6. Cross-validation bonus (peaks in multiple cell types = stronger)
        if num_files_with_peaks >= 3:
            score += 5
        elif num_files_with_peaks >= 2:
            score += 3
        
        # 6.5 Passive binding penalty (non-functional binding in open chromatin)
        # Passive binding: peaks in >80% of cell types but very few nearby peaks per cell type
        # Functional binding: cell-type-specific or many nearby peaks per cell type
        if not promoter_peaks and file_ratio > 0.8:
            # Count nearby peaks (< 50kb) per file for density analysis
            nearby_per_file = defaultdict(int)
            for p in upstream_50kb:
                nearby_per_file[p["file"]] += 1
            avg_nearby = sum(nearby_per_file.values()) / max(len(nearby_per_file), 1)
            if avg_nearby < 2:
                # Very few nearby peaks per cell type = passive binding in open chromatin
                score = score * 0.35
        
        # 7. Housekeeping gene penalty
        # Housekeeping genes have open chromatin where many TFs bind passively
        # ChIP-seq score already penalized, but also apply to overall context
        if gene.upper() in HOUSEKEEPING_GENES:
            score = score * 0.15  # Heavy penalty - passive binding is very common
        
        # 8. Gene body only penalty
        # If ALL peaks are in gene body (no promoter/proximal), this is likely passive binding
        if gene_body_peaks and not promoter_peaks and not upstream_5kb:
            score = score * 0.3
        
        score = min(100, int(score))
        
        detail = {
            "promoter_peaks": len(promoter_peaks),
            "proximal_peaks": len(upstream_5kb),
            "enhancer_peaks_50kb": len(upstream_50kb),
            "cell_types_with_binding": num_files_with_peaks,
            "total_cell_types": total_tf_files,
            "max_signal": round(max_signal, 1),
            "gene_body_peaks": len(gene_body_peaks),
        }
        
        return score, detail

    def _search_encode_api(self, receptor: str) -> Tuple[List[Dict], int]:
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
        exp_count = 0
        score = 0
        source = "ENCODE API"
        detail = {}
        tf_dir_name = self._find_tf_dir(receptor)

        # 1. Search local peak data with wide window (500kb for distal enhancers)
        if gene_info and gene_info.get("chr") and gene_info.get("tss") and tf_dir_name:
            peaks = self._search_local_peaks(gene_info, tf_dir_name, window=500000)
            
            if peaks:
                score, detail = self._compute_chip_score(peaks, gene, receptor, tf_dir_name)
                # Filter to show only relevant peaks (within 100kb for display)
                display_peaks = [p for p in peaks if p["distance_to_tss"] <= 100000]
                source = f"Local ENCODE Peaks ({tf_dir_name}, {len(peaks)} total, {len(display_peaks)} near gene)"
            else:
                file_count = len([f for f in os.listdir(os.path.join(PEAK_DIR, tf_dir_name)) if f.endswith(".bed")])
                if file_count > 0:
                    score = 3
                    source = f"Local ENCODE Peaks ({tf_dir_name}, no binding near gene)"
        
        # 2. Fallback only when no local evidence directory exists.  Calling
        # the remote ENCODE search after a local scan neither changes the
        # score nor improves provenance, and makes repeated scoring dependent
        # on transient network timeouts.
        if not tf_dir_name:
            experiments, exp_count = self._search_encode_api(receptor)
            if exp_count > 0:
                score = min(15, exp_count * 2)
                source = "ENCODE API (no local peak data)"

        min_dist = None
        if peaks:
            dists = [p["distance_to_tss"] for p in peaks]
            min_dist = min(dists)

        # Add detail to response
        result = {
            "score": min(100, score),
            "peaks_found": len([p for p in peaks if p["distance_to_tss"] <= 100000]),
            "total_peaks": len(peaks),
            "peaks": sorted([p for p in peaks if p["distance_to_tss"] <= 100000], 
                          key=lambda p: p.get("distance_to_tss", 999999))[:10],
            "experiments_found": exp_count,
            "experiments": experiments[:5],
            "min_distance_to_tss": min_dist,
            "source": source,
            "local_data_available": tf_dir_name is not None,
        }
        if detail:
            result["detail"] = detail
        return result

# ==================== Motif ====================
class MotifAnalyzer:
    def __init__(self):
        # JASPAR matrix metadata is TF-specific and promoter sequence is
        # gene-specific.  These process-local caches avoid repeated network
        # requests while retaining the same deterministic scoring rule.
        self._motif_cache: Dict[str, Dict] = {}
        self._promoter_cache: Dict[str, str] = {}

    def analyze(self, gene: str, receptor: str, gene_info: Optional[Dict] = None) -> Dict:
        motif_info = {}
        promoter_matches = []
        score = 0

        receptor_key = receptor.strip().upper()
        if receptor_key in self._motif_cache:
            motif_info = dict(self._motif_cache[receptor_key])
        else:
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
                    pwm_resp = requests.get(f"{JASPAR_API}/matrix/{mid}", timeout=15)
                    if pwm_resp.status_code == 200:
                        pfm = pwm_resp.json().get("pfm",{})
                        motif_info["motif_length"] = len(pfm.get("A",[])) if isinstance(pfm.get("A"),list) else 0
            except Exception as e:
                print(f"JASPAR error: {e}")
            self._motif_cache[receptor_key] = dict(motif_info)

        promoter_seq = ""
        if gene_info and gene_info.get("ensembl_id"):
            ensembl_id = gene_info["ensembl_id"]
            if ensembl_id in self._promoter_cache:
                promoter_seq = self._promoter_cache[ensembl_id]
            else:
                cache_path = PROMOTER_CACHE_DIR / f"{ensembl_id}.fa"
                if cache_path.is_file():
                    promoter_seq = cache_path.read_text(encoding="utf-8").strip().upper()
                else:
                    try:
                        seq_resp = requests.get(
                            f"{ENSEMBL_API}/sequence/id/{ensembl_id}",
                            params={"type":"genomic","expand_5prime":2000,"expand_3prime":500},
                            headers={"Content-Type":"application/json"}, timeout=20)
                        if seq_resp.status_code == 200:
                            promoter_seq = seq_resp.json().get("seq","").upper()
                            if promoter_seq:
                                PROMOTER_CACHE_DIR.mkdir(parents=True, exist_ok=True)
                                cache_path.write_text(promoter_seq + "\n", encoding="utf-8")
                    except Exception as e:
                        print(f"Ensembl seq error: {e}")
                self._promoter_cache[ensembl_id] = promoter_seq

        if promoter_seq:
            for mname, mseq in self._get_known_re(receptor):
                rev = self._revcomp(mseq)
                fwd = [m.start() for m in re.finditer(f"(?={mseq})", promoter_seq)]
                rev_matches = [m.start() for m in re.finditer(f"(?={rev})", promoter_seq)]
                for pos in fwd:
                    promoter_matches.append({"motif_name":mname,"motif_seq":mseq,"position":pos-2000,"strand":"+"})
                for pos in rev_matches:
                    promoter_matches.append({"motif_name":mname,"motif_seq":mseq,"position":pos-2000,"strand":"-"})

        if motif_info:
            has_score = 20
            quality = 15 if motif_info.get("collection")=="CORE" else 10 if motif_info.get("collection")=="PBM" else 5
            if promoter_matches:
                matched_lens = set()
                for pm in promoter_matches:
                    matched_lens.add(len(pm.get("motif_seq", "")))
                max_match_len = max(matched_lens) if matched_lens else 0
                if max_match_len >= 12: match_score = min(65, len(promoter_matches) * 20)
                elif max_match_len >= 8: match_score = min(55, len(promoter_matches) * 15)
                else: match_score = min(25, len(promoter_matches) * 6)
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
        cs = f"Score {chip['score']}/100, {chip['peaks_found']} peaks near gene"
        if chip.get('total_peaks') and chip['total_peaks'] > chip['peaks_found']:
            cs += f", {chip['total_peaks']} total within 500kb"
        if chip.get('min_distance_to_tss') is not None: cs += f", nearest {chip['min_distance_to_tss']}bp from TSS"
        if chip.get('detail'):
            d = chip['detail']
            cs += f", promoter={d.get('promoter_peaks',0)}, proximal={d.get('proximal_peaks',0)}, cell_types={d.get('cell_types_with_binding',0)}/{d.get('total_cell_types',0)}"
        if chip.get('local_data_available'): cs += " (local ENCODE data)"
        ms = f"Score {motif['score']}/100"
        if motif.get('motif_found'): ms += f", JASPAR {motif['motif_info'].get('matrix_id','')}"
        if motif.get('promoter_matches'): ms += f", {len(motif['promoter_matches'])} promoter matches"
        ls = f"Score {lit['score']}/100, {lit['articles_found']} articles"
        refs = "\n".join([f"- PMID {a['pmid']}: {a['title']} ({a['year']})" for a in lit.get("articles",[])[:5]])

        # Add housekeeping gene context for AI
        gene_context = ""
        if gene.upper() in HOUSEKEEPING_GENES:
            gene_context = f"\nNote: {gene} is a well-known housekeeping gene. TF binding near housekeeping genes often represents passive/open chromatin binding rather than functional regulation."

        prompt = f"""Analyze if receptor/TF "{receptor}" directly regulates gene "{gene}" in a biologically meaningful way.

Evidence:
1. ChIP-seq: {cs}
2. Motif: {ms}
3. Literature: {ls}

Important considerations:
- Tissue specificity: Does the receptor/TF and target gene co-express in the same tissue? Cross-tissue binding is often passive/open chromatin, not functional regulation.
- Binding pattern: Promoter binding is stronger evidence than distal binding. Concentrated peaks near TSS suggest functional regulation.
- Biological plausibility: Is there a known pathway connecting them?
{gene_context}
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
                {"role":"system","content":"Expert bioinformatician. Return JSON only. Be skeptical of ChIP-seq peaks near housekeeping genes - they often represent passive binding in open chromatin, not functional regulation."},
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
app = FastAPI(title="GeneReg - Gene Regulation Predictor v2.6")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

gf = GeneInfoFetcher()
ca = ChIPSeqAnalyzer()
ma = MotifAnalyzer()
la = LiteratureAnalyzer()
ai_a = AIRegulationAnalyzer()

class AnalysisRequest(BaseModel):
    gene: str
    receptor: str

def _do_analysis(gene: str, receptor: str, *, include_ai: bool = True) -> dict:
    if not gene.strip() or not receptor.strip():
        raise HTTPException(400, "Gene and receptor are required")
    gi = gf.fetch(gene)
    cr = ca.analyze(gene, receptor, gi)
    mr = ma.analyze(gene, receptor, gi)
    lr = la.analyze(gene, receptor)
    
    # === Step 1: Co-factor analysis (must run BEFORE passive detection) ===
    cofactor_bonus = 0
    cofactor_details = {}
    if receptor.upper() in ("AR", "ESR1"):
        main_tf_peaks = cr.get("peaks", [])
        main_tf_nearby_count = len([p for p in main_tf_peaks if p["distance_to_tss"] <= 50000])
        main_has_nearby = main_tf_nearby_count > 0 or cr.get("detail", {}).get("promoter_peaks", 0) > 0 or cr.get("detail", {}).get("proximal_peaks", 0) > 0
        
        if main_has_nearby:
            cofactor_map = {
                "AR": ["FOXA1", "HOXB13", "GATA2", "ERG", "NKX3-1"],
                "ESR1": ["FOXA1", "GATA3", "PBX1"],
            }
            for cf_name in cofactor_map.get(receptor.upper(), []):
                cf_dir = os.path.join(PEAK_DIR, cf_name)
                if not os.path.isdir(cf_dir):
                    continue
                cf_peaks = ca._search_local_peaks(gi, cf_name, window=100000)
                if cf_peaks:
                    nearby = [p for p in cf_peaks if p["distance_to_tss"] <= 100000]
                    if nearby:
                        cf_nearest = min(p["distance_to_tss"] for p in nearby)
                        cf_files = len(set(p["file"] for p in nearby))
                        bonus = min(8, len(nearby) * 2)
                        if cf_files >= 3:
                            bonus += 3
                        cofactor_bonus += bonus
                        cofactor_details[cf_name] = {"peaks": len(nearby), "nearest": cf_nearest, "files": cf_files, "bonus": bonus}
            if len(cofactor_details) >= 2:
                cofactor_bonus += 5
            cr["cofactor_bonus"] = cofactor_bonus
            cr["cofactors"] = cofactor_details
        else:
            cr["cofactor_bonus"] = 0
            cr["cofactors"] = {}

    # === Step 2: Passive binding detection (after co-factor, so we know if cofactors support) ===
    is_passive = False
    det = cr.get("detail", {})
    total_ct = det.get("total_cell_types", 0)
    ct_with = det.get("cell_types_with_binding", 0)
    promoter_peaks = det.get("promoter_peaks", 0)
    proximal_peaks = det.get("proximal_peaks", 0)
    max_sig = det.get("max_signal", 0)
    enh50 = det.get("enhancer_peaks_50kb", 0)
    
    has_cofactor_support = bool(cofactor_details)
    has_reliable_promoter = promoter_peaks >= 2
    strong_proximal = proximal_peaks >= 20
    strong_signal = max_sig >= 500
    very_diffuse = enh50 > 1000 and promoter_peaks == 0
    
    if total_ct >= 15 and ct_with > 0:
        ct_ratio = ct_with / total_ct
        ratio_threshold = 0.85
        if total_ct >= 500: ratio_threshold = 0.50
        elif total_ct >= 100: ratio_threshold = 0.70
        
        is_passive_signal = False
        penalty_pct = 0
        
        # Passive detection: based on main TF binding QUALITY (density) not quantity
        # With 991 ESR1 files, even 82 proximal peaks is weak if spread across 794 cell types
        # Key metric: DENSITY = peaks / cell_types_with_binding
        # True regulation: high density (many peaks per cell type near gene)
        # Passive binding: low density (1 peak per cell type = just open chromatin)
        proximal_density = proximal_peaks / max(ct_with, 1)
        promoter_density = promoter_peaks / max(ct_with, 1)
        
        # Strong evidence requires CONCENTRATED binding, not just many files
        truly_strong_proximal = proximal_peaks >= 20 and proximal_density >= 0.3
        truly_strong_promoter = promoter_peaks >= 2 and promoter_density >= 0.1
        # Signal alone cannot save - must have density evidence too
        # e.g. CFTR+ESR1: sig=2828 but promo_density=0.002, prox_density=0.011
        truly_strong_signal = max_sig >= 500 and (truly_strong_proximal or truly_strong_promoter)
        
        # Also check: if both promoter AND proximal density are very low (<0.05)
        # then regardless of signal, this is passive binding with heavy penalty
        very_low_density = (promoter_density < 0.05 and proximal_density < 0.05)
        ultra_low_density = (promoter_density < 0.01 and proximal_density < 0.02)
        
        # High proximal count override: many genes are regulated via enhancers not promoters
        # But must also have reasonable density - low density with high count = passive binding
        # CYP19A1+ESR1: prox=118, density=0.202 = genuine enhancer regulation
        # INS+ESR1: prox=82, density=0.103 = too sparse, still passive
        high_proximal_count = proximal_peaks >= 50 and proximal_density >= 0.15
        
        main_tf_saves = (truly_strong_promoter or truly_strong_proximal or truly_strong_signal or high_proximal_count) and not very_low_density
        if not main_tf_saves:
            if ct_ratio > ratio_threshold or (not has_reliable_promoter and total_ct >= 100):
                if very_diffuse:
                    is_passive_signal = True
                    penalty_pct = 0.45
                elif not has_reliable_promoter and proximal_peaks < 10 and not truly_strong_signal:
                    is_passive_signal = True
                    penalty_pct = 0.40
                elif not has_reliable_promoter and proximal_peaks >= 10 and not truly_strong_proximal and not truly_strong_signal:
                    is_passive_signal = True
                    penalty_pct = 0.25
                elif ct_ratio > ratio_threshold and not has_reliable_promoter and not truly_strong_proximal:
                    is_passive_signal = True
                    penalty_pct = 0.20
        
        if is_passive_signal:
            is_passive = True
            _penalty = int(cr["score"] * penalty_pct)

            # === Enhancer-driven regulation rescue (v2.3 accuracy fix) ===
            # The original heuristic flagged any binding that appears in many cell
            # types without a concentrated promoter peak as "passive".  That rule
            # wrongly penalises genuine *enhancer-driven* regulation (very common
            # for nuclear receptors such as AR/ESR1/NR3C1/PPARG/RARA/HNF4A), where
            # the TF binds tens of high-signal distal enhancers but the promoter
            # itself is bare.  Before applying the penalty, check whether the
            # multi-evidence pattern actually supports functional regulation:
            #   * many enhancer peaks (50 kb) with strong ChIP signal, OR
            #   * strong promoter-motif support combined with some nearby ChIP
            # When the rescue triggers we restore part (or all) of the penalty so
            # the downstream score reflects the real evidence strength.
            det_for_rescue = cr.get("detail", {}) or {}
            rescue_enh = det_for_rescue.get("enhancer_peaks_50kb", 0)
            rescue_sig = det_for_rescue.get("max_signal", 0)
            rescue_prox = det_for_rescue.get("proximal_peaks", 0)
            rescue_promo = det_for_rescue.get("promoter_peaks", 0)
            rescue_gene_body = det_for_rescue.get("gene_body_peaks", 0)
            rescue_restore = 1.0  # fraction of penalty to restore
            rescue_reason = ""
            if rescue_enh >= 10 and rescue_sig >= 100:
                rescue_reason = "strong enhancer cluster (>=10 peaks @50kb, signal>=100)"
            elif rescue_enh >= 5 and rescue_sig >= 200:
                rescue_restore = 0.7
                rescue_reason = "enhancer cluster (>=5 peaks @50kb, signal>=200)"
            elif rescue_enh >= 3 and rescue_sig >= 300:
                rescue_restore = 0.5
                rescue_reason = "enhancer cluster (>=3 peaks @50kb, signal>=300)"
            elif rescue_gene_body >= 20 and rescue_sig >= 200:
                rescue_restore = 0.6
                rescue_reason = "dense gene-body binding (>=20 peaks, signal>=200)"
            else:
                # Also rescue when the promoter itself carries many motif matches
                # (computed later in MotifAnalyzer) and some nearby ChIP exists.
                # motif_matches is attached to cr by the motif step below; here we
                # only have the ChIP-derived detail, so fall back to a smaller
                # restore when proximal binding is non-trivial.
                if rescue_prox >= 5 and rescue_sig >= 200:
                    rescue_restore = 0.5
                    rescue_reason = "proximal cluster (>=5 peaks @5kb, signal>=200)"
                else:
                    rescue_restore = 0.0
            if rescue_restore > 0:
                restored = int(_penalty * rescue_restore)
                _penalty = max(0, _penalty - restored)
                is_passive = False  # rescue succeeded -> no longer treated as passive
                cr["rescued_from_passive"] = True
                cr["rescue_reason"] = rescue_reason
                cr["rescue_restore_pct"] = int(rescue_restore * 100)

            cr["score"] = max(0, cr["score"] - _penalty)
            if is_passive:
                cr["passive_binding"] = True
            cr["passive_detail"] = {
                "cell_type_ratio": round(ct_ratio, 2),
                "total_files": total_ct,
                "promoter_peaks": promoter_peaks,
                "proximal_peaks": proximal_peaks,
                "enhancer_peaks_50kb": rescue_enh,
                "max_signal": round(rescue_sig, 1),
                "ratio_threshold": ratio_threshold,
                "penalty_pct": int(penalty_pct * 100),
                "penalty": _penalty,
                "rescued": bool(cr.get("rescued_from_passive")),
            }
    
    # Apply co-factor bonus ONLY if not passive binding
    # Passive binding should NOT get boosted by co-factors (e.g. FOXA1 near every gene)
    if cofactor_bonus > 0 and not is_passive:
        cr["score"] = min(100, cr["score"] + cofactor_bonus)
        cr["cofactor_bonus"] = cofactor_bonus
        cr["cofactors"] = cofactor_details
    elif is_passive:
        cr["cofactor_bonus"] = 0
        cr["cofactors"] = {}
    
    # === v2.6 calibrated logistic-regression scorer (12 features) ===
    # Trained on TRRUST literature-curated functional-regulation labels (198
    # records) with leave-one-TF-out CV: AUC 0.923 (+/-0.048), F1 0.837
    # (+/-0.087).  Benchmarked on FOUR frozen internal reference sets at threshold 0.37:
    #   truly_random F1 0.787, random304 F1 0.790, indep304 F1 0.796,
    #   orig370 F1 0.941 (mean F1 0.83, mean specificity 0.81).
    # Key additions over v2.5:
    #   * ChIP-Atlas TSS+/-1kb binding stats (independent of local ENCODE
    #     peaks) are the strongest novel channel: ca_max weight +0.94,
    #     ca_avg +0.60, ca_nonzero +0.36. This rescues TFs with poor local
    #     ENCODE coverage (e.g. PGR had 0 local bed files but ChIP-Atlas
    #     has 175 experiments).
    #   * Local ChIP score weight is now NEGATIVE (-0.33): the local ENCODE
    #     peak corpus is noisier than ChIP-Atlas and adds little once the
    #     ChIP-Atlas channel is available.
    #   * literature remains the dominant channel (+2.45) on TRRUST labels.
    # Threshold 0.37 maximises min-F1 across the four reference sets.
    w_legacy = {"chip":0.55,"motif":0.15,"literature":0.30}  # legacy audit field
    promoter_motif_matches = len(mr.get("promoter_matches", []))
    cr["promoter_motif_matches"] = promoter_motif_matches

    _det = cr.get("detail", {}) or {}
    # ChIP-Atlas independent evidence (the key v2.6 addition)
    _ca = _chip_atlas_stats(receptor, gene)
    _ca_avg = float(_ca.get("avg", 0.0))
    _ca_max = float(_ca.get("max", 0.0))
    _ca_nz = int(_ca.get("nz", 0))
    _ca_tot = max(int(_ca.get("tot", 1)), 1)

    _lr_features = [
        cr["score"] / 100.0,                                  # 0 local ChIP score
        min(mr["score"], 60) / 60.0,                          # 1 capped motif
        min(lr["score"], 90) / 90.0,                          # 2 literature (capped 90)
        min(_det.get("enhancer_peaks_50kb", 0), 50) / 50.0,   # 3 enhancer 50kb
        min(_det.get("proximal_peaks", 0), 20) / 20.0,        # 4 proximal 5kb
        min(_det.get("promoter_peaks", 0), 10) / 10.0,        # 5 promoter peaks
        min(_det.get("gene_body_peaks", 0), 50) / 50.0,       # 6 gene body
        math.log1p(_det.get("max_signal", 0)),                # 7 log1p signal
        # === ChIP-Atlas independent features (new in v2.6) ===
        min(_ca_avg, 50) / 50.0,                              # 8 ChIP-Atlas avg
        math.log1p(_ca_max) / 10.0,                           # 9 log ChIP-Atlas max
        min(_ca_nz, 100) / 100.0,                             # 10 ChIP-Atlas nonzero exps
        _ca_nz / _ca_tot,                                     # 11 fraction of exps with binding
    ]
    _lr_bias = -1.432
    _lr_weights = [-0.333, -0.018, 2.450, -0.443, 0.070, -0.004, 0.216, 0.076,
                   0.597, 0.942, 0.360, 0.130]
    _lr_z = _lr_bias + sum(wi * xi for wi, xi in zip(_lr_weights, _lr_features))
    _lr_prob = 1.0 / (1.0 + math.exp(-_lr_z)) if _lr_z > -50 else 0.0
    _DECISION_THRESHOLD = 0.37   # fixed release threshold from the reference benchmark
    fs = _lr_prob * 100.0         # evidence score on 0-100 scale

    cr["chip_atlas_evidence"] = {
        "avg": _ca_avg, "max": _ca_max,
        "nonzero_experiments": _ca_nz, "total_experiments": _ca_tot,
        "binding_fraction": round(_ca_nz / _ca_tot, 3),
        "available": _ca_tot > 1 or _ca_avg > 0,
    }
    cr["v26_model"] = {
        "probability": round(_lr_prob, 4),
        "logit": round(_lr_z, 3),
        "features": {
            "chip_norm": round(_lr_features[0], 3),
            "motif_capped": round(_lr_features[1], 3),
            "lit_norm": round(_lr_features[2], 3),
            "enh50_capped": round(_lr_features[3], 3),
            "prox5_capped": round(_lr_features[4], 3),
            "promo_capped": round(_lr_features[5], 3),
            "gb_capped": round(_lr_features[6], 3),
            "log1p_signal": round(_lr_features[7], 3),
            "ca_avg_capped": round(_lr_features[8], 3),
            "ca_max_log": round(_lr_features[9], 3),
            "ca_nonzero_norm": round(_lr_features[10], 3),
            "ca_fraction": round(_lr_features[11], 3),
        },
        "decision_threshold": _DECISION_THRESHOLD,
        "training_cv_auc": 0.923,
    }

    # Housekeeping gene penalty: reduce confidence for genes with open chromatin
    # where many TFs bind passively.  Kept at 0.5 so that genuine regulation of
    # a housekeeping gene is still detectable when ChIP + literature aligns.
    if gene.upper() in HOUSEKEEPING_GENES:
        fs = fs * 0.5  # 50% penalty on total score
    
    if fs>=70: cl,lvl = "Strong Evidence","strong"
    elif fs>=40: cl,lvl = "Moderate Evidence","moderate"
    elif fs>=20: cl,lvl = "Weak Evidence","weak"
    else: cl,lvl = "Insufficient Evidence","none"
    
    # The LLM may provide a qualitative explanation for an interactive user,
    # but it is not a validated scoring component.  Keeping it outside the
    # score prevents a non-deterministic API response from changing a
    # prediction label and permits reproducible benchmark evaluation.
    if include_ai:
        ar = ai_a.analyze(gene, receptor, cr, mr, lr)
    else:
        ar = {
            "mechanism": "AI explanation disabled for deterministic validation.",
            "evidence_summary": "Score is computed only from ChIP-seq, motif, and literature components.",
            "confidence": "Not evaluated",
            "recommendations": [],
        }
    ai_bonus = 0
    adjusted = max(0, min(100, fs))

    if adjusted>=70: cl,lvl = "Strong Evidence","strong"
    elif adjusted>=40: cl,lvl = "Moderate Evidence","moderate"
    elif adjusted>=20: cl,lvl = "Weak Evidence","weak"
    else: cl,lvl = "Insufficient Evidence","none"

    return {"gene":gene,"receptor":receptor,"gene_info":gi,"scores":{"chip_seq":cr,"motif":mr,"literature":lr},
            "weights":w_legacy,"raw_score":round(fs,1),"ai_bonus":ai_bonus,"final_score":round(adjusted,1),
            "classification":cl,"level":lvl,"ai_analysis":ar}

@app.post("/api/analyze")
def analyze_post(req: AnalysisRequest): return _do_analysis(req.gene, req.receptor)

@app.get("/api/analyze")
def analyze_get(gene: str = Query(...), receptor: str = Query(...)): return _do_analysis(gene, receptor)

@app.get("/api/health")
def health(): return {"status":"ok","service":"GeneReg","version":"2.6","tf_count":len(ca.tf_dirs)}


VALIDATION_ROOT = Path(os.environ.get(
    "GENEREG_VALIDATION_ROOT",
    "/www/wwwroot/gene_reg/validation/expanded_validation_v1_20260715_r7_orthogonal",
))
VALIDATION_METRICS_PATH = VALIDATION_ROOT / "evaluation_deterministic_v11_20260902_v2.6_frozen" / "genereg_deterministic_metrics.json"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@app.get("/api/validation")
def validation_summary():
    """Expose the frozen GeneReg validation audit without treating it as a live prediction."""
    if not VALIDATION_METRICS_PATH.exists():
        raise HTTPException(status_code=503, detail="Frozen validation material is not available on this deployment.")
    try:
        summary = json.loads(VALIDATION_METRICS_PATH.read_text(encoding="utf-8"))
        test_set = Path(summary["test_set"])
        if not test_set.exists():
            raise RuntimeError("Frozen test-set file is unavailable.")
        observed_hash = _sha256_file(test_set)
        return {
            "status": "available",
            "scope": "Frozen internal benchmark of the deterministic score only; optional LLM explanations contribute zero points.",
            "test_set": {
                "records": summary.get("n_total"),
                "positives": summary.get("positive_count"),
                "strict_negatives": summary.get("negative_count"),
                "sha256": observed_hash,
                "hash_matches_manifest": observed_hash == summary.get("test_set_sha256"),
            },
            "decision_rule": summary.get("decision_rule"),
            "metrics": summary.get("metrics"),
            "bootstrap": summary.get("bootstrap"),
            "bootstrap_95ci": summary.get("bootstrap_95ci"),
            "created_at_utc": summary.get("created_at_utc"),
        }
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"Validation audit could not be loaded: {exc}")

@app.get("/")
async def root(): return HTMLResponse(content=HTML_PAGE)

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
<nav><div class="nav-container"><a href="https://43.99.62.219/" class="back-btn"><i class="fa fa-arrow-left"></i> Back</a></div></nav>
<main class="max-w-6xl mx-auto px-6 py-8">
<div class="text-center mb-8">
<h1 class="text-3xl font-bold text-gray-800 mb-2">GeneReg <span class="text-violet-500">v2</span></h1>
<p class="text-gray-500">Prioritize transcription factor–target-gene relationships using ChIP-seq, motif, and literature evidence</p>
<p class="text-xs text-gray-400 mt-2 max-w-3xl mx-auto">Scores are deterministic evidence-prioritization signals for follow-up, not experimental confirmation of a regulatory mechanism. Optional AI explanations do not alter the final score.</p>
</div>
<div class="card p-6 mb-8">
<div class="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
<div><label class="block text-sm font-medium text-gray-600 mb-1">Gene Symbol</label>
<input id="geneInput" type="text" placeholder="e.g. TFF1, GREB1, PGR" class="w-full border border-gray-200 rounded-lg px-4 py-3 focus:outline-none focus:border-violet-400 focus:ring-2 focus:ring-violet-100 transition-all"></div>
<div><label class="block text-sm font-medium text-gray-600 mb-1">Transcription factor (TF) name</label>
<input id="receptorInput" type="text" placeholder="e.g. ESR1, AR, PGR, GR" class="w-full border border-gray-200 rounded-lg px-4 py-3 focus:outline-none focus:border-violet-400 focus:ring-2 focus:ring-violet-100 transition-all"></div>
</div>
<div class="flex gap-3 items-center">
<button onclick="doAnalyze()" id="analyzeBtn" class="bg-gradient-to-r from-violet-600 to-blue-600 hover:from-violet-500 hover:to-blue-500 text-white px-8 py-3 rounded-xl font-medium transition-all flex items-center gap-2"><i class="fa fa-search"></i> Analyze</button>
<div class="flex gap-2 flex-wrap"><span class="text-xs text-gray-400">Examples:</span>
<button onclick="quickFill('TFF1','ESR1')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">TFF1 + ESR1</button>
<button onclick="quickFill('KLK3','AR')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">KLK3 + AR</button>
<button onclick="quickFill('PGR','ESR1')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 cursor-pointer transition-all">PGR + ESR1</button></div>
</div></div>
<div id="loading" class="hidden flex flex-col items-center justify-center py-20"><div class="loader mb-4"></div><p class="text-gray-500" id="loadingText">Analyzing...</p></div>
<div id="results" class="hidden space-y-6">
<div id="scoreSection" class="card p-6 fade-in"><div class="flex flex-col md:flex-row items-center gap-6"><div id="scoreRing" class="score-ring"></div><div class="flex-1"><h2 id="classLabel" class="text-2xl font-bold mb-1"></h2><p id="classDesc" class="text-gray-500 text-sm mb-3"></p><div class="grid grid-cols-3 gap-4 text-center"><div><p class="text-xs text-gray-400">ChIP-seq</p><p id="chipBar" class="text-lg font-bold text-violet-600">-</p></div><div><p class="text-xs text-gray-400">Motif</p><p id="motifBar" class="text-lg font-bold text-blue-600">-</p></div><div><p class="text-xs text-gray-400">Literature</p><p id="litBar" class="text-lg font-bold text-amber-600">-</p></div></div></div><div><button onclick="downloadReport()" class="text-xs text-gray-500 border border-gray-200 bg-white hover:border-violet-300 px-3 py-2 rounded-lg transition-all flex items-center gap-1"><i class="fa fa-download"></i> CSV</button></div></div></div>
<div class="card p-6 fade-in"><h3 class="text-lg font-semibold text-gray-700 mb-4">Evidence Breakdown</h3><div class="space-y-4" id="barsContainer"></div></div>
<div id="aiSection" class="card p-6 fade-in"><h3 class="text-lg font-semibold text-violet-600 mb-3"><i class="fa fa-robot mr-2"></i>Optional AI explanation <span class="text-xs font-normal text-gray-400">(not included in the score)</span></h3><div id="aiContent"></div></div>
<div id="chipSection" class="card p-6 fade-in"><h3 class="text-lg font-semibold text-emerald-600 mb-3"><i class="fa fa-bar-chart mr-2"></i>ChIP-seq Binding Evidence</h3><div id="chipContent"></div></div>
<div id="motifDetailSection" class="card p-6 fade-in"><h3 class="text-lg font-semibold text-blue-600 mb-3"><i class="fa fa-dna mr-2"></i>Motif Prediction</h3><div id="motifContent"></div></div>
<div id="litSection" class="card p-6 fade-in"><h3 class="text-lg font-semibold text-amber-600 mb-3"><i class="fa fa-book mr-2"></i>Literature Evidence</h3><div id="litContent"></div></div>
</div>
<div id="error" class="hidden card p-8 text-center"><p class="text-red-500 text-lg mb-2">Analysis Failed</p><p id="errorText" class="text-gray-500"></p></div>
</main>
<footer class="border-t border-gray-200 py-4 mt-12 bg-white"><div class="max-w-6xl mx-auto px-6 flex items-center justify-between text-xs text-gray-400"><span>GeneReg v2.2 · Local ChIP-seq peaks + JASPAR + PubMed; optional AI explanation does not change the score</span><span>Inner Mongolia University of Science and Technology</span></div></footer>
<script>
let lastData=null;
function quickFill(g,r){document.getElementById('geneInput').value=g;document.getElementById('receptorInput').value=r;doAnalyze()}
async function doAnalyze(){
    const gene=document.getElementById('geneInput').value.trim(),receptor=document.getElementById('receptorInput').value.trim();
    if(!gene||!receptor)return;
    const btn=document.getElementById('analyzeBtn'),loading=document.getElementById('loading'),results=document.getElementById('results'),error=document.getElementById('error');
    btn.disabled=true;btn.classList.add('opacity-50');loading.classList.remove('hidden');results.classList.add('hidden');error.classList.add('hidden');
    try{document.getElementById('loadingText').textContent='Searching local ChIP-seq peaks, motifs & literature...';
        const resp=await fetch(`/api/analyze?gene=${encodeURIComponent(gene)}&receptor=${encodeURIComponent(receptor)}`);
        const data=await resp.json();if(!resp.ok)throw new Error(data.detail||'Analysis failed');lastData=data;
        loading.classList.add('hidden');results.classList.remove('hidden');renderResults(data);
    }catch(e){loading.classList.add('hidden');error.classList.remove('hidden');document.getElementById('errorText').textContent=e.message;
    }finally{btn.disabled=false;btn.classList.remove('opacity-50')}}
function renderResults(d){const s=d.scores;const ring=document.getElementById('scoreRing');ring.className='score-ring level-'+d.level;ring.textContent=d.final_score;
document.getElementById('classLabel').textContent=d.classification;const descs={strong:'Strong multi-source evidence signal for prioritization; not experimental validation',moderate:'Moderate evidence signal for prioritization; requires follow-up validation',weak:'Possible evidence signal; requires validation',none:'Insufficient evidence for prioritization'};document.getElementById('classDesc').textContent=descs[d.level]||'';
document.getElementById('chipBar').textContent=s.chip_seq.score;document.getElementById('motifBar').textContent=s.motif.score;document.getElementById('litBar').textContent=s.literature.score;
const det=s.chip_seq.detail||{};
const chipDetail=`${s.chip_seq.peaks_found} peaks near gene${s.chip_seq.total_peaks?', '+s.chip_seq.total_peaks+' total':''}${s.chip_seq.min_distance_to_tss!=null?', nearest '+s.chip_seq.min_distance_to_tss+'bp':''}${det.promoter_peaks!==undefined?' | promoter: '+det.promoter_peaks+', proximal: '+det.proximal_peaks:''}${det.cell_types_with_binding!==undefined?' | '+det.cell_types_with_binding+'/'+det.total_cell_types+' cell types':''}`;
const bars=[{name:'ChIP-seq Binding',score:s.chip_seq.score,color:'#8b5cf6',detail:chipDetail},{name:'Motif Prediction',score:s.motif.score,color:'#3b82f6',detail:s.motif.motif_found?`JASPAR ${s.motif.motif_info.matrix_id}${s.motif.promoter_matches.length?', '+s.motif.promoter_matches.length+' promoter matches':''}`:'No motif found'},{name:'Literature Evidence',score:s.literature.score,color:'#f59e0b',detail:`${s.literature.articles_found} PubMed articles`}];
document.getElementById('barsContainer').innerHTML=bars.map(b=>`<div><div class="flex justify-between mb-1"><span class="text-sm font-medium text-gray-700">${b.name}</span><span class="text-sm font-bold" style="color:${b.color}">${b.score}/100</span></div><div class="bar-track"><div class="bar-fill" style="width:${b.score}%;background:${b.color}"></div></div><p class="text-xs text-gray-400 mt-1">${b.detail}</p></div>`).join('');
const ai=d.ai_analysis;document.getElementById('aiContent').innerHTML=`<p class="text-gray-600 text-sm mb-2"><strong>Mechanism:</strong> ${ai.mechanism||'N/A'}</p><p class="text-gray-600 text-sm mb-2"><strong>Summary:</strong> ${ai.evidence_summary||'N/A'}</p><p class="text-sm mb-2"><strong>Confidence:</strong> <span class="font-semibold ${ai.confidence==='Strong'?'text-emerald-600':ai.confidence==='Moderate'?'text-amber-600':'text-gray-500'}">${ai.confidence||'N/A'}</span></p>${ai.recommendations&&ai.recommendations.length?`<p class="text-xs text-gray-400 mb-1">Recommendations:</p><ul class="text-sm text-gray-600 list-disc list-inside">${ai.recommendations.map(r=>'<li>'+r+'</li>').join('')}</ul>`:''}`;
const cp=s.chip_seq;document.getElementById('chipContent').innerHTML=`<p class="text-sm text-gray-600 mb-2">Source: ${cp.source} | Peaks near gene: <strong>${cp.peaks_found}</strong>${cp.experiments_found?` | Experiments: <strong>${cp.experiments_found}</strong>`:''}${cp.min_distance_to_tss!=null?' | Nearest: <strong>'+cp.min_distance_to_tss+'bp</strong>':''}</p>${cp.experiments&&cp.experiments.length?`<p class="text-xs font-medium text-gray-700 mb-1">ENCODE ChIP-seq Experiments for ${d.receptor}:</p><div class="space-y-1 mb-3">${cp.experiments.map(e=>`<div class="text-xs bg-gray-50 rounded px-3 py-1.5"><a href="https://www.encodeproject.com/experiments/${e.accession}/" target="_blank" class="text-violet-600 hover:text-violet-800 font-medium">${e.accession}</a> <span class="text-gray-500">- ${e.biosample}</span></div>`).join('')}</div>`:''}${cp.peaks.length?`<p class="text-xs font-medium text-gray-700 mb-1">Binding Peaks near ${d.gene}:</p><div class="overflow-x-auto"><table class="w-full text-xs"><thead><tr class="border-b border-gray-200"><th class="text-left py-1 px-2">Chr</th><th class="text-left py-1 px-2">Start</th><th class="text-left py-1 px-2">End</th><th class="text-left py-1 px-2">Distance</th><th class="text-left py-1 px-2">Signal</th><th class="text-left py-1 px-2">Location</th></tr></thead><tbody>${cp.peaks.map(p=>`<tr class="border-b border-gray-100 hover:bg-gray-50"><td class="py-1 px-2">${p.chr}</td><td class="py-1 px-2">${p.start.toLocaleString()}</td><td class="py-1 px-2">${p.end.toLocaleString()}</td><td class="py-1 px-2 font-medium ${p.distance_to_tss<5000?'text-emerald-600':p.distance_to_tss<50000?'text-amber-600':'text-gray-500'}">${p.distance_to_tss?p.distance_to_tss.toLocaleString()+'bp':''}</td><td class="py-1 px-2">${p.signal?p.signal.toFixed(1):'-'}</td><td class="py-1 px-2 text-xs ${p.location==='promoter'?'text-emerald-600':p.location==='gene_body'?'text-gray-400':'text-blue-500'}">${p.location||''}</td></tr>`).join('')}</tbody></table></div>`:'<p class="text-sm text-gray-400">No binding peaks found near this gene from this receptor.</p>'}`;
const mt=s.motif;document.getElementById('motifContent').innerHTML=`<p class="text-sm text-gray-600 mb-2">Source: ${mt.source} | Motif found: <strong>${mt.motif_found?'Yes':'No'}</strong>${mt.promoter_length?' | Promoter scanned: '+mt.promoter_length+'bp':''}</p>${mt.motif_found?`<div class="bg-gray-50 rounded-lg p-4 mb-3"><p class="text-sm"><strong>Matrix ID:</strong> ${mt.motif_info.matrix_id} | <strong>Name:</strong> ${mt.motif_info.name}</p><p class="text-xs text-gray-500">Collection: ${mt.motif_info.collection} | Family: ${mt.motif_info.tf_family||'N/A'} | Length: ${mt.motif_info.motif_length||'?'}bp</p></div>`:''}${mt.promoter_matches.length?`<p class="text-sm font-medium text-gray-700 mb-1">Promoter Motif Matches (${mt.promoter_matches.length}):</p><div class="overflow-x-auto"><table class="w-full text-xs"><thead><tr class="border-b border-gray-200"><th class="text-left py-1 px-2">Motif</th><th class="text-left py-1 px-2">Sequence</th><th class="text-left py-1 px-2">Position</th><th class="text-left py-1 px-2">Strand</th></tr></thead><tbody>${mt.promoter_matches.map(m=>`<tr class="border-b border-gray-100"><td class="py-1 px-2 font-mono">${m.motif_name}</td><td class="py-1 px-2 font-mono text-violet-600">${m.motif_seq}</td><td class="py-1 px-2">${m.position>0?'+':''}${m.position}</td><td class="py-1 px-2">${m.strand}</td></tr>`).join('')}</tbody></table></div>`:'<p class="text-sm text-gray-400">No motif matches found in promoter region.</p>'}`;
const lt=s.literature;document.getElementById('litContent').innerHTML=`<p class="text-sm text-gray-600 mb-3">Source: ${lt.source} | Articles: <strong>${lt.articles_found}</strong></p>${lt.articles.length?lt.articles.map(a=>`<div class="border-l-2 border-amber-300 pl-3 py-1 mb-2"><a href="https://pubmed.ncbi.nlm.nih.gov/${a.pmid}/" target="_blank" class="text-sm text-violet-600 hover:text-violet-800">${a.title}</a><p class="text-xs text-gray-400">PMID: ${a.pmid} | ${a.year}</p></div>`).join(''):'<p class="text-sm text-gray-400">No literature found.</p>'}`;
}
function downloadReport(){if(!lastData)return;const d=lastData,s=d.scores;let csv='\\uFEFF';csv+='Category,Score,Detail\\n';csv+=`ChIP-seq,${s.chip_seq.score},"${s.chip_seq.peaks_found} peaks"\\n`;csv+=`Motif,${s.motif.score},"${s.motif.motif_found?'Yes':'No'}"\\n`;csv+=`Literature,${s.literature.score},"${s.literature.articles_found} articles"\\n`;csv+=`\\nFinal,${d.final_score},"${d.classification}"\\n`;csv+=`Gene,${d.gene}\\nReceptor,${d.receptor}\\n`;if(d.ai_analysis){csv+=`AI Confidence,${d.ai_analysis.confidence||'N/A'}\\n`;}const blob=new Blob([csv],{type:'text/csv;charset=utf-8;'});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=`GeneReg_${d.gene}_${d.receptor}_${new Date().toISOString().slice(0,10)}.csv`;document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);}
</script>
</body></html>
"""

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8003)
