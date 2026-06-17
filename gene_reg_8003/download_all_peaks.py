#!/usr/bin/env python3
"""
Download ALL narrowPeak files for TF ChIP-seq from ENCODE
Including: IDR thresholded, conservative IDR, optimal IDR, peaks, pseudoreplicated IDR
This dramatically increases coverage for TFs like AR
"""
import os, json, time, requests, gzip, shutil, re
from collections import defaultdict

ENCODE_API = "https://www.encodeproject.org"
OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
HEADERS = {"Accept": "application/json"}
DELAY = 0.15

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
            if os.path.exists(p): os.remove(p)
        return False

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Check existing files
    existing = set()
    tf_file_counts = defaultdict(int)
    for d in os.listdir(OUT_DIR):
        dp = os.path.join(OUT_DIR, d)
        if os.path.isdir(dp) and not d.startswith("_"):
            for f in os.listdir(dp):
                if f.endswith(".bed"):
                    existing.add(f.replace(".bed", ""))
                    tf_file_counts[d] += 1
    
    print(f"Already have {len(existing)} files across {len(tf_file_counts)} TFs")
    
    # Search for ALL narrowPeak files (not just IDR)
    output_types = [
        "IDR thresholded peaks",
        "conservative IDR thresholded peaks",
        "optimal IDR thresholded peaks",
        "peaks",
        "pseudoreplicated IDR thresholded peaks",
    ]
    
    all_files = []
    for ot in output_types:
        offset = 0
        limit = 500
        ot_count = 0
        while True:
            try:
                r = requests.get(f"{ENCODE_API}/search/", params={
                    "type": "File",
                    "file_type": "bed narrowPeak",
                    "output_type": ot,
                    "assay_title": "TF ChIP-seq",
                    "status": "released",
                    "format": "json",
                    "limit": limit,
                    "from": offset,
                    "frame": "object"
                }, headers=HEADERS, timeout=60)
                r.raise_for_status()
                data = r.json()
            except Exception as e:
                print(f"  API error: {e}")
                time.sleep(5)
                continue
            
            files = data.get("@graph", [])
            if not files:
                break
            
            for f in files:
                acc = f.get("accession", "")
                href = f.get("href", "")
                if not acc or not href:
                    continue
                
                # Extract TF from target
                target = f.get("target", "")
                tf_name = ""
                if isinstance(target, str):
                    m = re.search(r'/targets/([^/]+)', target)
                    if m:
                        tf_name = m.group(1).replace("-human", "").replace("-mouse", "").replace("-rat", "")
                elif isinstance(target, dict):
                    tf_name = target.get("label", "")
                
                if not tf_name:
                    continue
                
                # Only GRCh38
                assembly = f.get("assembly", "")
                if assembly != "GRCh38":
                    continue
                
                cm = f.get("cloud_metadata", {})
                url = ""
                if isinstance(cm, dict):
                    url = cm.get("url", "")
                if not url:
                    url = f"{ENCODE_API}{href}"
                
                all_files.append({
                    "url": url,
                    "accession": acc,
                    "tf": tf_name,
                    "output_type": ot,
                    "file_size": f.get("file_size", 0),
                })
                ot_count += 1
            
            total = data.get("total", 0)
            offset += limit
            if offset >= total or len(files) < limit:
                break
            time.sleep(DELAY)
        
        print(f"  {ot}: found {ot_count} files")
    
    print(f"\nTotal files to consider: {len(all_files)}")
    
    # Deduplicate by accession
    seen = set()
    unique_files = []
    for f in all_files:
        if f["accession"] not in seen:
            seen.add(f["accession"])
            unique_files.append(f)
    print(f"Unique accessions: {len(unique_files)}")
    
    # Group by TF
    tf_groups = defaultdict(list)
    for f in unique_files:
        tf_groups[f["tf"]].append(f)
    
    print(f"Covering {len(tf_groups)} TFs")
    top_tfs = sorted(tf_groups.items(), key=lambda x: len(x[1]), reverse=True)[:10]
    for tf, files in top_tfs:
        print(f"  {tf}: {len(files)} files")
    
    # Download
    total_dl = 0
    total_skip = 0
    total_fail = 0
    start = time.time()
    
    for tf_name in sorted(tf_groups.keys()):
        files = tf_groups[tf_name]
        tf_safe = tf_name.replace("/", "_").replace(" ", "_")
        tf_dir = os.path.join(OUT_DIR, tf_safe)
        os.makedirs(tf_dir, exist_ok=True)
        
        ok = 0
        skip = 0
        fail = 0
        for f in files:
            if f["accession"] in existing:
                skip += 1
                total_skip += 1
                continue
            fp = os.path.join(tf_dir, f"{f['accession']}.bed")
            if os.path.exists(fp):
                skip += 1
                total_skip += 1
                continue
            if download_and_decompress(f["url"], fp):
                ok += 1
                total_dl += 1
                existing.add(f["accession"])
            else:
                fail += 1
                total_fail += 1
            time.sleep(DELAY)
        
        if ok > 0:
            print(f"  {tf_name}: +{ok} ({skip} skip, {fail} fail)")
    
    elapsed = time.time() - start
    print(f"\nDone! Downloaded: {total_dl}, Skipped: {total_skip}, Failed: {total_fail}")
    print(f"Time: {elapsed/60:.1f} min")

if __name__ == "__main__":
    main()
