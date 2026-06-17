#!/usr/bin/env python3
"""
Download ENCODE TF ChIP-seq peak files - OPTIMIZED v3
Strategy: 
1. Search Experiments API for TF ChIP-seq with narrowPeak files (bulk, no per-experiment requests)
2. Use file search to get download URLs directly
3. Match experiments to files via accession
"""
import os, json, time, requests, gzip, shutil, pickle, re, sys
from collections import defaultdict

ENCODE_API = "https://www.encodeproject.org"
OUT_DIR = "/www/wwwroot/gene_reg/encode_peaks"
INDEX_PATH = os.path.join(OUT_DIR, "_peak_index.pkl")

HEADERS = {"Accept": "application/json"}
REQUEST_DELAY = 0.15

def log(msg):
    print(msg, flush=True)

def download_and_decompress(url, filepath):
    gz_path = filepath + ".gz"
    try:
        r = requests.get(url, timeout=180, stream=True, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        # Try gzip decompress
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
        return False

def step1_get_experiment_file_map():
    """
    Step 1: Search Files API for narrowPeak+IDR files that belong to TF ChIP-seq experiments.
    Use assay_title filter to get only TF ChIP-seq files.
    """
    log("Step 1: Fetching TF ChIP-seq narrowPeak files from ENCODE Files API...")
    
    all_files = []
    limit = 500
    offset = 0
    
    while True:
        try:
            r = requests.get(f"{ENCODE_API}/search/", params={
                "type": "File",
                "file_type": "bed narrowPeak",
                "output_type": "IDR thresholded peaks",
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
            log(f"  API error at offset {offset}: {e}, retrying in 5s...")
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
            
            # Extract dataset accession from dataset field (e.g. "/experiments/ENCSR000AAA/")
            dataset = f.get("dataset", "")
            dataset_acc = ""
            if isinstance(dataset, str):
                m = re.search(r'(ENCSR\w+)', dataset)
                if m:
                    dataset_acc = m.group(1)
            
            assembly = f.get("assembly", "unknown")
            file_size = f.get("file_size", 0)
            
            # Use cloud_metadata URL if available (faster CDN)
            cloud = f.get("cloud_metadata", {})
            cdn_url = ""
            if isinstance(cloud, dict):
                cdn_url = cloud.get("url", "")
            
            all_files.append({
                "url": cdn_url or f"{ENCODE_API}{href}",
                "accession": acc,
                "dataset_acc": dataset_acc,
                "assembly": assembly,
                "file_size": file_size
            })
        
        total = data.get("total", 0)
        offset += limit
        log(f"  Got {len(all_files)}/{total} files")
        
        if offset >= total or len(files) < limit:
            break
        time.sleep(REQUEST_DELAY)
    
    log(f"Total files found: {len(all_files)}")
    return all_files

def step2_get_tf_mapping(file_list):
    """
    Step 2: For each unique experiment accession, look up the TF name.
    Use batch experiment search instead of per-experiment requests.
    """
    log("\nStep 2: Mapping experiments to TF names...")
    
    # Get unique experiment accessions
    exp_accs = set(f["dataset_acc"] for f in file_list if f["dataset_acc"])
    log(f"  {len(exp_accs)} unique experiments to look up")
    
    tf_map = {}  # exp_acc -> tf_name
    
    # Search experiments in batches by TF name using facets
    # Alternative: search all TF ChIP-seq experiments and extract TF from target.label
    offset = 0
    limit = 500
    batch_size = 500
    
    while True:
        try:
            r = requests.get(f"{ENCODE_API}/search/", params={
                "type": "Experiment",
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
            log(f"  API error: {e}, retrying in 5s...")
            time.sleep(5)
            continue
        
        exps = data.get("@graph", [])
        if not exps:
            break
        
        for exp in exps:
            exp_acc = exp.get("accession", "")
            target = exp.get("target", {})
            tf_name = ""
            if isinstance(target, dict):
                tf_name = target.get("label", "")
            elif isinstance(target, str):
                # target is a URL like "/targets/CTCF-human/"
                m = re.search(r'/targets/([^/]+)', target)
                if m:
                    tf_name = m.group(1).replace("-human", "").replace("-mouse", "")
            
            if exp_acc and tf_name:
                tf_map[exp_acc] = tf_name
        
        total = data.get("total", 0)
        offset += limit
        log(f"  Mapped {len(tf_map)} experiments (scanned {offset}/{total})")
        
        if offset >= total or len(exps) < limit:
            break
        time.sleep(REQUEST_DELAY)
    
    log(f"  TF mapping complete: {len(tf_map)} experiments -> TF names")
    return tf_map

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
    
    log("=" * 60)
    log("ENCODE Peak Downloader v3 - Optimized Bulk API")
    log("=" * 60)
    
    # Step 1: Get all narrowPeak files for TF ChIP-seq
    peak_files = step1_get_experiment_file_map()
    
    if not peak_files:
        log("No files found! Check API filters.")
        return
    
    # Step 2: Map experiments to TF names
    tf_map = step2_get_tf_mapping(peak_files)
    
    # Merge TF names into file list
    for pf in peak_files:
        pf["tf"] = tf_map.get(pf["dataset_acc"], "unknown")
    
    # Group by TF
    tf_groups = defaultdict(list)
    for pf in peak_files:
        tf_groups[pf["tf"]].append(pf)
    
    log(f"\n  {len(peak_files)} files across {len(tf_groups)} TFs")
    for tf in sorted(tf_groups.keys())[:10]:
        log(f"    {tf}: {len(tf_groups[tf])} files")
    if len(tf_groups) > 10:
        log(f"    ... and {len(tf_groups) - 10} more TFs")
    
    # Step 3: Download
    log(f"\nStep 3: Downloading {len(peak_files)} files...")
    
    # Check existing
    done_files = set()
    for tf_dir_name in os.listdir(OUT_DIR):
        tf_dir = os.path.join(OUT_DIR, tf_dir_name)
        if os.path.isdir(tf_dir) and not tf_dir_name.startswith("_"):
            for fname in os.listdir(tf_dir):
                if fname.endswith(".bed"):
                    done_files.add(fname.replace(".bed", ""))
    
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
        
        if ok > 0:
            elapsed = time.time() - start
            rate = total_dl / elapsed * 60 if elapsed > 0 else 0
            log(f"  {tf_name}: +{ok} downloaded ({fail} fail) [{rate:.0f} files/min]")
        elif skip > 0:
            pass  # Don't log skipped TFs
        elif fail > 0:
            log(f"  {tf_name}: {fail} failed!")
        
        # Progress
        if total_dl > 0 and total_dl % 200 == 0:
            elapsed = time.time() - start
            eta = (len(peak_files) - total_dl - total_skip) / (total_dl / elapsed) / 60
            log(f"  === {total_dl} downloaded, {total_fail} failed | ETA: {eta:.0f} min ===")
    
    elapsed = time.time() - start
    log(f"\n{'=' * 60}")
    log(f"Download complete!")
    log(f"  Downloaded: {total_dl}")
    log(f"  Skipped: {total_skip}")
    log(f"  Failed: {total_fail}")
    log(f"  Time: {elapsed/60:.1f} minutes")
    log(f"{'=' * 60}")
    
    # Step 4: Build index
    log("\nStep 4: Building peak index for fast lookup...")
    build_index()
    log("\nAll done! 🎉")

if __name__ == "__main__":
    main()
