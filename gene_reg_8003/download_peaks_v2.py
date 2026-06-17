#!/usr/bin/env python3
"""
Download ENCODE TF ChIP-seq peak files - OPTIMIZED v2
Directly query ENCODE Files API for narrowPeak+IDR files
Skip per-experiment detail requests (10x faster)
"""
import os, json, time, requests, gzip, shutil, pickle, sys
from collections import defaultdict

ENCODE_API = "https://www.encodeproject.org"
OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
INDEX_PATH = os.path.join(OUT_DIR, "_peak_index.pkl")

HEADERS = {"Accept": "application/json"}
# Be polite to ENCODE API
REQUEST_DELAY = 0.15

def log(msg):
    print(msg, flush=True)

def download_and_decompress(url, filepath):
    gz_path = filepath + ".gz"
    try:
        r = requests.get(url, timeout=120, stream=True)
        r.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        # Try gzip decompress first
        try:
            with gzip.open(gz_path, "rb") as f_in:
                with open(filepath, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
        except (gzip.BadGzipFile, OSError):
            # Not actually gzipped, just rename
            os.rename(gz_path, filepath)
        return True
    except Exception as e:
        for p in [gz_path, filepath]:
            if os.path.exists(p):
                os.remove(p)
        return False

def get_all_peak_files():
    """Get ALL narrowPeak IDR files from ENCODE in bulk - no per-experiment requests!"""
    all_files = []
    limit = 500
    offset = 0
    
    while True:
        log(f"  Fetching files batch (offset={offset})...")
        try:
            r = requests.get(f"{ENCODE_API}/search/", params={
                "type": "File",
                "file_type": "bed narrowPeak",
                "output_type": "IDR thresholded peaks",
                "status": "released",
                "format": "json",
                "limit": limit,
                "from": offset,
                "frame": "object"
            }, headers=HEADERS, timeout=60)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            log(f"  API error: {e}, retrying in 5s...")
            time.sleep(5)
            continue
        
        files = data.get("@graph", [])
        if not files:
            break
        
        for f in files:
            if not isinstance(f, dict):
                continue
            href = f.get("href", "")
            acc = f.get("accession", "")
            if not href or not acc:
                continue
            
            # Extract TF name from dataset target
            dataset = f.get("dataset", "")
            target = f.get("target", {})
            tf_name = ""
            if isinstance(target, dict):
                tf_name = target.get("label", "")
            elif isinstance(target, str):
                tf_name = target.replace("/targets/", "").replace("/", "")
            
            # Extract biosample
            biosample = ""
            bio = f.get("biosample_ontology", {})
            if isinstance(bio, dict):
                biosample = bio.get("term_name", "")
            
            assembly = f.get("assembly", "unknown")
            file_size = f.get("file_size", 0)
            
            all_files.append({
                "url": f"{ENCODE_API}{href}",
                "accession": acc,
                "tf": tf_name or "unknown",
                "biosample": biosample,
                "assembly": assembly,
                "file_size": file_size
            })
        
        total = data.get("total", 0)
        offset += limit
        log(f"  Got {len(all_files)}/{total} files so far")
        
        if offset >= total or len(files) < limit:
            break
        time.sleep(REQUEST_DELAY)
    
    return all_files

def build_index():
    log("\nBuilding peak index...")
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
                        if len(parts) < 3:
                            continue
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
    
    log(f"Index: {total_peaks} peaks across {len(index)} chromosomes")
    log(f"Saved to: {INDEX_PATH}")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Step 1: Bulk fetch all peak file metadata
    log("=" * 60)
    log("ENCODE Peak Downloader v2 - Bulk Files API")
    log("=" * 60)
    log("\nStep 1: Fetching all narrowPeak+IDR file metadata from ENCODE...")
    peak_files = get_all_peak_files()
    log(f"\nFound {len(peak_files)} peak files total")
    
    # Group by TF
    tf_groups = defaultdict(list)
    for pf in peak_files:
        tf_groups[pf["tf"]].append(pf)
    log(f"Covering {len(tf_groups)} TFs")
    
    # Check what's already downloaded
    done_files = set()
    for tf_dir_name in os.listdir(OUT_DIR):
        tf_dir = os.path.join(OUT_DIR, tf_dir_name)
        if os.path.isdir(tf_dir) and not tf_dir_name.startswith("_"):
            for fname in os.listdir(tf_dir):
                if fname.endswith(".bed"):
                    done_files.add(fname.replace(".bed", ""))
    log(f"Already downloaded: {len(done_files)} files")
    
    # Step 2: Download
    log(f"\nStep 2: Downloading files...")
    total_dl = 0
    total_skip = 0
    total_fail = 0
    start = time.time()
    
    for tf_name, files in sorted(tf_groups.items()):
        tf_safe = tf_name.replace("/", "_").replace(" ", "_")
        tf_dir = os.path.join(OUT_DIR, tf_safe)
        os.makedirs(tf_dir, exist_ok=True)
        
        ok = 0
        skip = 0
        fail = 0
        for pf in files:
            fp = os.path.join(tf_dir, f"{pf['accession']}.bed")
            if os.path.exists(fp) or pf["accession"] in done_files:
                skip += 1
                total_skip += 1
                continue
            if download_and_decompress(pf["url"], fp):
                ok += 1
                total_dl += 1
            else:
                fail += 1
                total_fail += 1
            time.sleep(REQUEST_DELAY)
        
        if ok > 0 or skip > 0:
            log(f"  {tf_name}: {ok} downloaded, {skip} skipped, {fail} failed")
        
        # Progress every 50 TFs
        if total_dl > 0 and total_dl % 100 == 0:
            elapsed = time.time() - start
            rate = total_dl / elapsed * 60
            log(f"  === Progress: {total_dl} downloaded | {rate:.0f} files/min ===")
    
    elapsed = time.time() - start
    log(f"\n{'=' * 60}")
    log(f"Download complete!")
    log(f"  Downloaded: {total_dl}")
    log(f"  Skipped (existing): {total_skip}")
    log(f"  Failed: {total_fail}")
    log(f"  Time: {elapsed/60:.1f} minutes")
    log(f"{'=' * 60}")
    
    # Step 3: Build index
    log("\nStep 3: Building peak index for fast lookup...")
    build_index()
    log("\nAll done!")

if __name__ == "__main__":
    main()
