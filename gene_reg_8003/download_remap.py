#!/usr/bin/env python3
"""
Download AR and other key TF ChIP-seq peak data from ReMap database
ReMap has much more TF binding data than ENCODE
URL: https://remap.univ-amu.fr/
"""
import os, gzip, shutil, requests, time
from collections import defaultdict

OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
REMAP_BASE = "https://remap.univ-amu.fr/storage/v2/consensus"

# Key TFs that need more data
TFS_TO_DOWNLOAD = [
    "AR",      # Androgen receptor - critical, only 3 files in ENCODE
    "PGR",     # Progesterone receptor
    "ESR2",    # Estrogen receptor beta
    "STAT1",   # Signal transducer
    "STAT3",   # Signal transducer
    "SP1",     # General TF
    "NFYA",    # General TF
    "NFYB",    # General TF
    "RELA",    # NF-kB
]

def download_remap_tf(tf_name):
    """Download ReMap consensus peaks for a TF"""
    # ReMap URL format: https://remap.univ-amu.fr/storage/v2/consensus/{tf_name}.bed.gz
    url = f"{REMAP_BASE}/{tf_name}.bed.gz"
    tf_safe = tf_name.replace("/", "_").replace(" ", "_")
    tf_dir = os.path.join(OUT_DIR, tf_safe)
    os.makedirs(tf_dir, exist_ok=True)
    
    # Check if already exists
    out_file = os.path.join(tf_dir, f"remap_{tf_safe}.bed")
    if os.path.exists(out_file):
        print(f"  {tf_name}: already have ReMap file")
        return True
    
    print(f"  Downloading ReMap {tf_name}...", end=" ", flush=True)
    try:
        r = requests.get(url, timeout=120, stream=True, headers={"User-Agent": "Mozilla/5.0"})
        if r.status_code != 200:
            print(f"not available ({r.status_code})")
            return False
        
        gz_path = out_file + ".gz"
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        
        # Decompress
        try:
            with gzip.open(gz_path, "rb") as f_in:
                with open(out_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
        except (gzip.BadGzipFile, OSError):
            os.rename(gz_path, out_file)
        
        # Count peaks
        count = 0
        with open(out_file) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    count += 1
        print(f"{count} peaks")
        return True
    except Exception as e:
        print(f"error: {e}")
        for p in [out_file, out_file + ".gz"]:
            if os.path.exists(p): os.remove(p)
        return False

def main():
    print("Downloading ReMap consensus peaks for key TFs...")
    ok = 0
    fail = 0
    for tf in TFS_TO_DOWNLOAD:
        if download_remap_tf(tf):
            ok += 1
        else:
            fail += 1
        time.sleep(0.5)
    
    print(f"\nReMap download complete: {ok} ok, {fail} fail")

if __name__ == "__main__":
    main()
