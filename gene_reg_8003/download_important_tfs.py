#!/usr/bin/env python3
"""
Download missing important TF ChIP-seq peak files from ENCODE
Focuses on nuclear receptors and key TFs for gene regulation prediction
"""
import os, json, time, requests, gzip, shutil, re
from collections import defaultdict

ENCODE_API = "https://www.encodeproject.org"
OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
HEADERS = {"Accept": "application/json"}

# Important TFs/receptors for gene regulation prediction
IMPORTANT_TFS = [
    # Nuclear receptors
    "ESR1", "ESR2", "AR", "PGR", "NR3C1", "NR3C2", "NR3C3",
    "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG",
    "PPARA", "PPARD", "PPARG",
    "VDR", "FXR", "NR1H3", "NR1H2", "NR1I2",
    "THRA", "THRB",
    "HNF4A", "HNF4G",
    # Key TFs
    "FOXA1", "FOXA2", "FOXA3",
    "GATA1", "GATA2", "GATA3", "GATA4", "GATA6",
    "MYC", "MAX", "MXD1",
    "SP1", "SP2",
    "NFYA", "NFYB",
    "E2F1", "E2F4",
    "JUN", "JUND", "FOS", "FOSL1", "FOSL2",
    "MYB", "REL", "RELA", "NFKB1",
    "TP53", "TP63", "TP73",
    "STAT1", "STAT3", "STAT5A", "STAT5B",
    "BRD4", "EP300", "CREBBP",
    "TBP", "TAF1",
    "CEBPB", "CEBPD",
    "YY1", "YY2",
    "KLF1", "KLF4", "KLF5",
    "RUNX1", "RUNX2", "RUNX3",
    "SMAD1", "SMAD2", "SMAD3", "SMAD4",
    "TCF7L2", "LEF1",
    "NOTCH1", "RBPJ",
    "HOXA9", "HOXB4",
    "SOX2", "SOX9", "SOX10", "SOX17",
    "OCT4", "POU5F1", "NANOG",
    "NEUROD1", "ASCL1",
    "IRF1", "IRF4",
    "BATF", "BATF3",
    "BCL11A", "BCL11B",
    "REST", "CTCF",
    "RAD21", "SMC3", "STAG1",
    "POLR2A", "POLR2B",
    "EP300", "CREBBP", "KAT2B",
    # Additional important
    "EGR1", "EGR2",
    "ATF2", "ATF3", "ATF4", "ATF7",
    "NRF1", "NRF2", "NFE2L2",
    "USF1", "USF2",
    "MAFK", "BACH1", "BACH2",
    "MAZ", "SP3",
]

def download_and_decompress(url, filepath):
    gz_path = filepath + ".gz"
    try:
        r = requests.get(url, timeout=180, stream=True, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        try:
            with gzip.open(gz_path, "rb") as f_in:
                with open(filepath, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
        except (gzip.BadGzipFile, OSError):
            os.rename(gz_path, filepath)
        return True
    except Exception as e:
        for p in [gz_path, filepath]:
            if os.path.exists(p):
                os.remove(p)
        print(f"    Download failed: {e}")
        return False

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Check which TFs we already have
    existing = set()
    for d in os.listdir(OUT_DIR):
        dp = os.path.join(OUT_DIR, d)
        if os.path.isdir(dp) and not d.startswith("_"):
            beds = [f for f in os.listdir(dp) if f.endswith(".bed")]
            if beds:
                existing.add(d)
    
    print(f"Already have {len(existing)} TFs with data")
    
    total_dl = 0
    total_fail = 0
    
    for tf in IMPORTANT_TFS:
        tf_safe = tf.replace("/", "_").replace(" ", "_")
        tf_dir = os.path.join(OUT_DIR, tf_safe)
        
        if tf_safe in existing:
            # Check how many files
            beds = [f for f in os.listdir(tf_dir) if f.endswith(".bed")]
            print(f"  {tf}: already have {len(beds)} files, skipping")
            continue
        
        print(f"\n{tf}: searching ENCODE...", end=" ", flush=True)
        
        # Search for narrowPeak IDR files for this TF
        try:
            r = requests.get(f"{ENCODE_API}/search/", params={
                "type": "File",
                "file_type": "bed narrowPeak",
                "output_type": "IDR thresholded peaks",
                "assay_title": "TF ChIP-seq",
                "target.label": tf,
                "status": "released",
                "format": "json",
                "limit": 100,
                "frame": "object"
            }, headers=HEADERS, timeout=30)
            
            if r.status_code != 200:
                print(f"API error {r.status_code}")
                continue
            
            files = r.json().get("@graph", [])
            if not files:
                print("no files found")
                continue
            
            print(f"{len(files)} files", end="", flush=True)
            os.makedirs(tf_dir, exist_ok=True)
            
            ok = 0
            fail = 0
            for f in files:
                acc = f.get("accession", "")
                if not acc:
                    continue
                fp = os.path.join(tf_dir, f"{acc}.bed")
                if os.path.exists(fp):
                    ok += 1
                    continue
                
                # Try CDN URL first, then ENCODE href
                cm = f.get("cloud_metadata", {})
                url = ""
                if isinstance(cm, dict):
                    url = cm.get("url", "")
                if not url:
                    href = f.get("href", "")
                    if href:
                        url = f"{ENCODE_API}{href}"
                
                if not url:
                    continue
                
                if download_and_decompress(url, fp):
                    ok += 1
                    total_dl += 1
                else:
                    fail += 1
                    total_fail += 1
                time.sleep(0.15)
            
            print(f" ({ok} ok, {fail} fail)")
            
        except Exception as e:
            print(f"Error: {e}")
        
        time.sleep(0.2)
    
    print(f"\nDone! Downloaded: {total_dl}, Failed: {total_fail}")

if __name__ == "__main__":
    main()
