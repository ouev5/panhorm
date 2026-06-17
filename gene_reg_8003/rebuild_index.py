#!/usr/bin/env python3
"""Rebuild peak index from all downloaded bed files"""
import os, pickle
from collections import defaultdict

PEAK_DIR = "/www/wwwroot/gene_reg/encode_peaks"
INDEX_PATH = os.path.join(PEAK_DIR, "_peak_index.pkl")

def build_index():
    print("Building peak index...")
    index = {}
    total_peaks = 0
    total_files = 0
    
    for tf_dir_name in os.listdir(PEAK_DIR):
        tf_dir = os.path.join(PEAK_DIR, tf_dir_name)
        if not os.path.isdir(tf_dir) or tf_dir_name.startswith("_"):
            continue
        for fname in os.listdir(tf_dir):
            if not fname.endswith(".bed"):
                continue
            filepath = os.path.join(tf_dir, fname)
            total_files += 1
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
            except Exception as e:
                print(f"  Error reading {filepath}: {e}")
    
    for chrom in index:
        index[chrom].sort(key=lambda x: x[0])
    
    with open(INDEX_PATH, "wb") as f:
        pickle.dump(index, f)
    
    print(f"Done! {total_peaks} peaks across {len(index)} chromosomes from {total_files} files")
    print(f"Saved to: {INDEX_PATH}")
    
    # Show top chromosomes
    top = sorted(index.items(), key=lambda x: len(x[1]), reverse=True)[:5]
    for chrom, peaks in top:
        print(f"  chr{chrom}: {len(peaks)} peaks")

if __name__ == "__main__":
    build_index()
