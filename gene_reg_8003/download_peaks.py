#!/usr/bin/env python3
"""
Download ENCODE TF ChIP-seq peak files - per-TF strategy
For each TF, search experiments -> get peak file URLs -> download
"""
import os, json, time, requests, gzip, shutil, pickle
from collections import defaultdict

ENCODE_API = "https://www.encodeproject.org"
OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
INDEX_PATH = os.path.join(OUT_DIR, "_peak_index.pkl")

def download_and_decompress(url, filepath):
    gz_path = filepath + ".gz"
    try:
        r = requests.get(url, timeout=120, stream=True)
        r.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        with gzip.open(gz_path, "rb") as f_in:
            with open(filepath, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)
        return True
    except Exception as e:
        for p in [gz_path, filepath]:
            if os.path.exists(p): os.remove(p)
        return False

def get_tf_list():
    """Get all TFs with ChIP-seq data from ENCODE facets"""
    r = requests.get(f"{ENCODE_API}/search/", params={
        "type": "Experiment",
        "assay_title": "TF ChIP-seq",
        "status": "released",
        "format": "json",
        "limit": 0
    }, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    data = r.json()
    tfs = []
    for fc in data.get("facets", []):
        if fc.get("field") == "target.label":
            for term in fc.get("terms", []):
                tfs.append(term["key"])
            break
    return tfs

def get_peak_files_for_tf(tf_name):
    """Get peak file URLs for a specific TF"""
    peak_files = []
    limit = 25
    offset = 0
    
    while True:
        r = requests.get(f"{ENCODE_API}/search/", params={
            "type": "Experiment",
            "assay_title": "TF ChIP-seq",
            "target.label": tf_name,
            "status": "released",
            "format": "json",
            "limit": limit,
            "from": offset
        }, headers={"Accept": "application/json"}, timeout=30)
        if r.status_code != 200:
            break
        
        exps = r.json().get("@graph", [])
        if not exps:
            break
        
        for exp in exps:
            exp_acc = exp.get("accession", "")
            bio = exp.get("biosample_ontology", {})
            biosample = bio.get("term_name", "unknown") if isinstance(bio, dict) else "unknown"
            
            # Get files from experiment
            try:
                fr = requests.get(f"{ENCODE_API}/experiments/{exp_acc}/", params={
                    "format": "json"
                }, headers={"Accept": "application/json"}, timeout=20)
                if fr.status_code != 200:
                    continue
                detail = fr.json()
                for f in detail.get("files", []):
                    if not isinstance(f, dict):
                        continue
                    ft = f.get("file_type", "")
                    ot = f.get("output_type", "")
                    href = f.get("href", "")
                    acc = f.get("accession", "")
                    sz = f.get("file_size", 0)
                    if "narrowPeak" in ft and "IDR" in ot and href:
                        peak_files.append({
                            "url": f"{ENCODE_API}{href}",
                            "accession": acc,
                            "file_size": sz,
                            "biosample": biosample
                        })
            except:
                pass
            time.sleep(0.1)
        
        offset += limit
        if len(exps) < limit:
            break
        time.sleep(0.3)
    
    return peak_files

def build_index():
    print("\nBuilding peak index...")
    index = {}
    total_peaks = 0
    
    for tf_dir_name in os.listdir(OUT_DIR):
        tf_dir = os.path.join(OUT_DIR, tf_dir_name)
        if not os.path.isdir(tf_dir) or tf_dir_name.startswith("_"):
            continue
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
                        if len(parts) < 3: continue
                        chrom = parts[0].replace("chr", "")
                        try:
                            start, end = int(parts[1]), int(parts[2])
                        except ValueError:
                            continue
                        if chrom not in index:
                            index[chrom] = []
                        index[chrom].append((start, end, tf_dir_name))
                        total_peaks += 1
            except:
                pass
    
    for chrom in index:
        index[chrom].sort(key=lambda x: x[0])
    
    with open(INDEX_PATH, "wb") as f:
        pickle.dump(index, f)
    
    print(f"Index: {total_peaks} peaks across {len(index)} chromosomes")
    print(f"Saved to: {INDEX_PATH}")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Get all TF names
    print("Getting TF list from ENCODE...")
    tfs = get_tf_list()
    print(f"Found {len(tfs)} TFs with ChIP-seq data")
    
    # Check which TFs already downloaded
    done_tfs = set()
    for d in os.listdir(OUT_DIR):
        if os.path.isdir(os.path.join(OUT_DIR, d)) and not d.startswith("_"):
            bed_count = len([f for f in os.listdir(os.path.join(OUT_DIR, d)) if f.endswith(".bed")])
            if bed_count > 0:
                done_tfs.add(d)
    print(f"Already downloaded: {len(done_tfs)} TFs")
    
    total_dl = 0
    total_fail = 0
    start = time.time()
    
    for i, tf in enumerate(tfs):
        tf_safe = tf.replace("/", "_").replace(" ", "_")
        if tf_safe in done_tfs:
            continue
        
        print(f"\n[{i+1}/{len(tfs)}] {tf}", end=" ... ", flush=True)
        
        peak_files = get_peak_files_for_tf(tf)
        if not peak_files:
            print("no peak files")
            continue
        
        print(f"{len(peak_files)} files", end="", flush=True)
        tf_dir = os.path.join(OUT_DIR, tf_safe)
        os.makedirs(tf_dir, exist_ok=True)
        
        ok = 0
        fail = 0
        for pf in peak_files:
            fp = os.path.join(tf_dir, f"{pf['accession']}.bed")
            if os.path.exists(fp):
                ok += 1
                continue
            if download_and_decompress(pf["url"], fp):
                ok += 1
                total_dl += 1
            else:
                fail += 1
                total_fail += 1
            time.sleep(0.2)
        
        print(f" ({ok} ok, {fail} fail)")
        done_tfs.add(tf_safe)
        
        if (i + 1) % 20 == 0:
            elapsed = time.time() - start
            rate = (i + 1) / elapsed * 60
            print(f"  === {i+1}/{len(tfs)} TFs done | {total_dl} files downloaded | {rate:.0f} TFs/min ===")
    
    elapsed = time.time() - start
    print(f"\nDone! {total_dl} downloaded, {total_fail} failed in {elapsed/60:.0f} minutes")
    
    build_index()

if __name__ == "__main__":
    main()
