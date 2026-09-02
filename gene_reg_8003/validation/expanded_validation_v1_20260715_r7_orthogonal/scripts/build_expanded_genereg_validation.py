#!/usr/bin/env python3
"""Build a provenance-preserving independent validation set for GeneReg.

This script deliberately does not import or modify GeneReg's scoring code.  It
downloads the external source files into a separate validation directory,
constructs a frozen 178-positive/192-negative set, and records every rule used
to select a row.  No DeepSeek call is made.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


RUN_ID = "expanded_validation_v1_20260715_r7_orthogonal"
GTRD_ROOT = "http://gtrd.biouml.org:8888/downloads/current"
CHIP_ATLAS_APP = "https://chip-atlas.org/target_genes?type=tsv"
PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_REST = "https://rest.ensembl.org"
JASPAR_API = "https://jaspar.elixir.no/api/v1"
SPECIES = "Homo sapiens"
ASSEMBLY = "hg38"
LOCAL_PEAK_DIR = Path("/www/wwwroot/gene_reg/encode_peaks")
POSITIVE_PER_TF = 20
NEGATIVE_PER_TF = 24
PROMOTER_UPSTREAM = 2000
PROMOTER_DOWNSTREAM = 500
# A PWM hit is called only at >=95% of the theoretical score range.  This
# operationalizes "no JASPAR motif" as no high-confidence promoter motif and
# avoids treating ubiquitous low-information nuclear-receptor half-sites as
# functional motif support.
JASPAR_RELATIVE_THRESHOLD = 0.95


# The eight TFs are hormone/nuclear-receptor relevant, have local GeneReg
# peak directories, published GTRD promoter target files, and ChIP-Atlas
# Target Genes output.  Exact JASPAR CORE matrix IDs are configured explicitly
# to avoid false substring matches (for example, AR must not resolve to ALX3).
# The matrix ID and PWM threshold are saved for every screening result.
TF_CONFIG: dict[str, dict[str, str]] = {
    "AR": {"uniprot": "P10275", "jaspar_matrix_id": "MA0007.2"},
    "ESR1": {"uniprot": "P03372", "jaspar_matrix_id": "MA0112.4"},
    "PGR": {"uniprot": "P06401", "jaspar_matrix_id": "MA2327.1"},
    "NR3C1": {"uniprot": "P04150", "jaspar_matrix_id": "MA0113.4"},
    "PPARG": {"uniprot": "P37231", "jaspar_matrix_id": "MA0066.2"},
    "HNF4A": {"uniprot": "P41235", "jaspar_matrix_id": "MA0114.5"},
    "RARA": {"uniprot": "P10276", "jaspar_matrix_id": "MA0729.1"},
    "GATA3": {"uniprot": "P23771", "jaspar_matrix_id": "MA0037.3"},
}


# Manually curated recent (2024–2025) high-confidence pairs.  Each source
# reports at least two independent experimental axes for the relevant
# regulation (binding/occupancy and expression or functional perturbation).
# The raw PubMed XML is retained under raw/pubmed/ for audit.
RECENT_LITERATURE_PAIRS: list[dict[str, str]] = [
    *[
        {
            "tf": "HNF4A",
            "gene": gene,
            "pmid": "38909044",
            "year": "2024",
            "evidence_summary": (
                "Manual curation: HNF4A ChIP-seq binding plus tissue-specific "
                "target-gene regulation/functional follow-up reported in beta "
                "cells and hepatocytes."
            ),
        }
        for gene in (
            "ACY3", "HAAO", "HNF1A", "MAP3K11", "ABCD3", "CDKN2AIP",
            "USH1C", "VIL1", "GAD2", "HOPX", "AKAP1",
        )
    ],
    *[
        {
            "tf": "AR",
            "gene": gene,
            "pmid": "40596315",
            "year": "2025",
            "evidence_summary": (
                "Manual curation: AR target-gene regulation supported by RNA-seq, "
                "qPCR/Western validation, and promoter-binding/transactivation assays "
                "in prostate cancer cells."
            ),
        }
        for gene in ("ZFHX3", "KLK3", "FKBP5", "TMPRSS2", "TNK1", "ADAM7", "MAPRE2")
    ],
]


EXCLUDE_GENE_PREFIXES = (
    "MTRNR", "MT-", "RPL", "RPS", "MRPL", "MRPS", "HIST", "LOC", "LINC",
    "MIR", "SNOR", "SCARNA", "OR", "TR", "IG", "HLA-",
)
EXCLUDE_GENE_EXACT = {"", "NA", "N/A"}


class SourceError(RuntimeError):
    """Raised when an external source cannot be retrieved reproducibly."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def http_get(session: requests.Session, url: str, *, timeout: int = 90, params: dict[str, str] | None = None) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            response = session.get(url, params=params, timeout=timeout, headers={"User-Agent": "AHormoneDB-GeneReg-validation/1.0"})
            response.raise_for_status()
            return response
        except Exception as exc:  # Requests raises different transient errors.
            last_error = exc
            if attempt < 3:
                time.sleep(attempt * 2)
    raise SourceError(f"Failed to retrieve {url}: {last_error}")


def http_post_json(session: requests.Session, url: str, payload: dict[str, Any], *, timeout: int = 90) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            response = session.post(
                url,
                json=payload,
                timeout=timeout,
                headers={"User-Agent": "AHormoneDB-GeneReg-validation/1.0", "Accept": "application/json"},
            )
            response.raise_for_status()
            return response
        except Exception as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(attempt * 2)
    raise SourceError(f"Failed to POST to {url}: {last_error}")


def download(session: requests.Session, url: str, destination: Path, manifest: list[dict[str, Any]]) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        last_error: Exception | None = None
        for attempt in range(1, 4):
            try:
                with session.get(url, stream=True, timeout=(20, 300), headers={"User-Agent": "AHormoneDB-GeneReg-validation/1.0"}) as response:
                    response.raise_for_status()
                    with destination.open("wb") as output:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                output.write(chunk)
                break
            except Exception as exc:
                last_error = exc
                destination.unlink(missing_ok=True)
                if attempt == 3:
                    raise SourceError(f"Failed to download {url}: {last_error}")
                time.sleep(attempt * 3)
    manifest.append({
        "retrieved_at_utc": utc_now(),
        "url": url,
        "path": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
    })
    return destination


def gtrd_target_url(uniprot: str) -> str:
    return (
        f"{GTRD_ROOT}/intervals/target_genes/Homo%20sapiens/"
        f"genes%20promoter%5B-1000,+100%5D/{uniprot}.txt"
    )


def chip_atlas_url(session: requests.Session, tf: str) -> str:
    payload = {"condition": {"genome": ASSEMBLY, "antigen": tf, "distance": "1"}}
    response = session.post(CHIP_ATLAS_APP, json=payload, timeout=60, headers={"User-Agent": "AHormoneDB-GeneReg-validation/1.0"})
    response.raise_for_status()
    url = response.json().get("url", "")
    if not url.startswith("https://"):
        raise SourceError(f"ChIP-Atlas did not return a data URL for {tf}: {response.text[:500]}")
    return url


def load_gtrd_targets(path: Path) -> dict[str, dict[str, Any]]:
    values: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = (row.get("Gene symbol") or "").strip().upper()
            ensembl_id = (row.get("ID") or "").strip().split(".")[0]
            try:
                site_count = int(row.get("SiteCount") or 0)
            except ValueError:
                site_count = 0
            if gene and ensembl_id:
                values[gene] = {"ensembl_id": ensembl_id, "gtrd_site_count": site_count}
    return values


def local_chip_atlas_experiment_ids(tf: str) -> set[str]:
    """Return ChIP-Atlas SRX datasets already present in GeneReg's peak corpus.

    The live GeneReg project stores these files as ``chipatlas_SRX....bed``.
    They must not provide the orthogonal ChIP-Atlas confirmation used for a
    GTRD-labelled validation positive.
    """
    tf_dir = LOCAL_PEAK_DIR / tf
    if not tf_dir.is_dir():
        return set()
    return {
        item.stem.removeprefix("chipatlas_")
        for item in tf_dir.glob("chipatlas_SRX*.bed")
    }


def load_chip_atlas_targets(path: Path, local_experiment_ids: set[str]) -> dict[str, dict[str, Any]]:
    values: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 2 or header[0] != "Target_genes":
            raise SourceError(f"Unexpected ChIP-Atlas target-gene format: {path}")
        experiment_columns = [
            (index, value.split("|", 1)[0])
            for index, value in enumerate(header[2:-1], 2)
            if value.split("|", 1)[0].startswith("SRX")
        ]
        independent_experiment_columns = [
            (index, accession)
            for index, accession in experiment_columns
            if accession not in local_experiment_ids
        ]
        for row in reader:
            if len(row) < 2:
                continue
            gene = row[0].strip().upper()
            try:
                average = float(row[1])
            except ValueError:
                average = 0.0
            experiment_scores: list[float] = []
            independent_experiment_scores: list[float] = []
            for index, _ in experiment_columns:
                value = row[index] if index < len(row) else ""
                try:
                    experiment_scores.append(float(value))
                except ValueError:
                    pass
            for index, _ in independent_experiment_columns:
                value = row[index] if index < len(row) else ""
                try:
                    independent_experiment_scores.append(float(value))
                except ValueError:
                    pass
            values[gene] = {
                "chip_atlas_average": round(average, 6),
                "chip_atlas_nonzero_experiments": sum(value > 0 for value in experiment_scores),
                "chip_atlas_independent_nonzero_experiments": sum(value > 0 for value in independent_experiment_scores),
                "chip_atlas_independent_experiment_universe": len(independent_experiment_columns),
            }
    return values


def valid_gene_symbol(gene: str) -> bool:
    if gene in EXCLUDE_GENE_EXACT or any(gene.startswith(prefix) for prefix in EXCLUDE_GENE_PREFIXES):
        return False
    return bool(re.fullmatch(r"[A-Z][A-Z0-9-]{1,19}", gene))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def fetch_jaspar_matrix(session: requests.Session, tf: str, matrix_id: str, cache_dir: Path) -> dict[str, Any]:
    cache_path = cache_dir / f"{tf}_{matrix_id}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    matrix = http_get(session, f"{JASPAR_API}/matrix/{matrix_id}").json()
    if str(matrix.get("name") or "").upper() != tf.upper():
        raise SourceError(f"JASPAR matrix {matrix_id} does not exactly correspond to TF {tf}")
    payload = {"tf": tf, "matrix_id": matrix_id, "matrix": matrix, "source_url": f"{JASPAR_API}/matrix/{matrix_id}"}
    write_json(cache_path, payload)
    return payload


def promoter_sequence(session: requests.Session, ensembl_id: str, cache_dir: Path) -> str:
    cache_path = cache_dir / f"{ensembl_id}.fa"
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8").strip()
    lookup = http_get(session, f"{ENSEMBL_REST}/lookup/id/{ensembl_id}", params={"expand": "0"}).json()
    chromosome = str(lookup.get("seq_region_name") or "")
    start, end, strand = int(lookup.get("start") or 0), int(lookup.get("end") or 0), int(lookup.get("strand") or 1)
    if not chromosome or not start or not end:
        raise SourceError(f"Incomplete Ensembl location for {ensembl_id}")
    if strand == 1:
        left, right = max(1, start - PROMOTER_UPSTREAM), start + PROMOTER_DOWNSTREAM
    else:
        left, right = max(1, end - PROMOTER_DOWNSTREAM), end + PROMOTER_UPSTREAM
    sequence = http_get(
        session,
        f"{ENSEMBL_REST}/sequence/region/human/{chromosome}:{left}..{right}:1",
        params={"content-type": "text/plain"},
    ).text.upper().replace("\n", "").strip()
    if not sequence or not re.fullmatch(r"[ACGTN]+", sequence):
        raise SourceError(f"Invalid promoter sequence for {ensembl_id}")
    cache_path.write_text(sequence + "\n", encoding="utf-8")
    return sequence


def pwm_hits(sequence: str, matrix_payload: dict[str, Any]) -> tuple[int, float]:
    """Return high-stringency PWM hits on both strands and the raw threshold."""
    pfm = matrix_payload["matrix"].get("pfm") or {}
    columns = list(zip(pfm.get("A", []), pfm.get("C", []), pfm.get("G", []), pfm.get("T", [])))
    if not columns:
        raise SourceError(f"Missing PFM in JASPAR matrix {matrix_payload['matrix_id']}")
    weights: list[dict[str, float]] = []
    max_score = 0.0
    min_score = 0.0
    for a, c, g, t in columns:
        counts = {"A": float(a), "C": float(c), "G": float(g), "T": float(t)}
        total = sum(counts.values()) + 0.4
        column = {base: __import__("math").log2((count + 0.1) / total / 0.25) for base, count in counts.items()}
        weights.append(column)
        max_score += max(column.values())
        min_score += min(column.values())
    threshold = min_score + JASPAR_RELATIVE_THRESHOLD * (max_score - min_score)
    motif_len = len(weights)
    def count_hits(seq: str) -> int:
        hits = 0
        for i in range(0, len(seq) - motif_len + 1):
            window = seq[i:i + motif_len]
            if "N" in window:
                continue
            score = sum(weights[pos][base] for pos, base in enumerate(window))
            if score >= threshold:
                hits += 1
        return hits
    return count_hits(sequence) + count_hits(reverse_complement(sequence)), round(threshold, 6)


def pubmed_count(session: requests.Session, tf: str, gene: str) -> tuple[int, str]:
    query = f'"{tf}"[Title/Abstract] AND "{gene}"[Title/Abstract]'
    result = http_get(session, f"{PUBMED_EUTILS}/esearch.fcgi", params={"db": "pubmed", "term": query, "retmax": "0", "retmode": "json"}).json()
    try:
        return int(result["esearchresult"]["count"]), query
    except (KeyError, TypeError, ValueError) as exc:
        raise SourceError(f"Unexpected PubMed result for {query}: {exc}")


def fetch_pubmed_xml(session: requests.Session, pmid: str, raw_dir: Path) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    destination = raw_dir / f"PMID_{pmid}.xml"
    if destination.exists():
        return
    response = http_get(session, f"{PUBMED_EUTILS}/efetch.fcgi", params={"db": "pubmed", "id": pmid, "retmode": "xml"})
    destination.write_bytes(response.content)
    root = ET.fromstring(response.content)
    date = root.findtext(".//JournalIssue/PubDate/Year") or root.findtext(".//JournalIssue/PubDate/MedlineDate") or ""
    if not re.match(r"20(23|24|25)", date):
        raise SourceError(f"PMID {pmid} is outside the 2023–2025 curation interval: {date}")


def select_gtrd_positives(
    gtrd_by_tf: dict[str, dict[str, dict[str, Any]]],
    chip_by_tf: dict[str, dict[str, dict[str, Any]]],
    excluded_pairs: set[tuple[str, str]],
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for tf in TF_CONFIG:
        candidates: list[dict[str, Any]] = []
        for gene, source in gtrd_by_tf[tf].items():
            chip = chip_by_tf[tf].get(gene, {
                "chip_atlas_average": 0.0,
                "chip_atlas_nonzero_experiments": 0,
                "chip_atlas_independent_nonzero_experiments": 0,
                "chip_atlas_independent_experiment_universe": 0,
            })
            if (
                (tf, gene) in excluded_pairs
                or not valid_gene_symbol(gene)
                or source["gtrd_site_count"] < 2
                or chip["chip_atlas_average"] <= 0
                # At least two ChIP-Atlas SRX experiments that are not stored
                # in the current GeneReg peak corpus.  Together with the GTRD
                # >=2-meta-cluster condition, this gives a conservative
                # multi-evidence positive label.
                or chip["chip_atlas_independent_nonzero_experiments"] < 2
            ):
                continue
            candidates.append({"gene": gene, **source, **chip})
        candidates.sort(key=lambda item: (-item["gtrd_site_count"], -item["chip_atlas_average"], item["gene"]))
        if len(candidates) < POSITIVE_PER_TF:
            raise SourceError(f"Only {len(candidates)} eligible GTRD/ChIP-Atlas positives for {tf}")
        for item in candidates[:POSITIVE_PER_TF]:
            selected.append({
                "label": 1,
                "subset": "GTRD_independent",
                "tf": tf,
                "gene": item["gene"],
                "ensembl_id": item["ensembl_id"],
                "gtrd_site_count": item["gtrd_site_count"],
                "chip_atlas_average": item["chip_atlas_average"],
                "chip_atlas_nonzero_experiments": item["chip_atlas_nonzero_experiments"],
                "chip_atlas_independent_nonzero_experiments": item["chip_atlas_independent_nonzero_experiments"],
                "chip_atlas_independent_experiment_universe": item["chip_atlas_independent_experiment_universe"],
                "evidence_summary": "GTRD promoter target with at least two meta-clusters, confirmed by at least two non-zero ChIP-Atlas SRX experiments absent from the current GeneReg peak corpus.",
                "selection_rule": "gtrd_site_count>=2 AND chip_atlas_average>0 AND chip_atlas_independent_nonzero_experiments>=2",
                "pmid": "",
                "publication_year": "",
            })
    return selected


def chunked(items: list[Any], size: int) -> list[list[Any]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def fetch_promoters_bulk(session: requests.Session, ensembl_ids: list[str], cache_dir: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Cache promoter regions with Ensembl bulk lookup and bulk region APIs."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    sequences: dict[str, str] = {}
    errors: dict[str, str] = {}
    missing: list[str] = []
    for ensembl_id in ensembl_ids:
        cache_path = cache_dir / f"{ensembl_id}.fa"
        if cache_path.exists():
            sequences[ensembl_id] = cache_path.read_text(encoding="utf-8").strip()
        else:
            missing.append(ensembl_id)
    locations: dict[str, dict[str, Any]] = {}
    for batch in chunked(missing, 100):
        try:
            result = http_post_json(session, f"{ENSEMBL_REST}/lookup/id", {"ids": batch}).json()
        except SourceError as exc:
            for ensembl_id in batch:
                errors[ensembl_id] = str(exc)
            continue
        for ensembl_id in batch:
            item = result.get(ensembl_id)
            if not isinstance(item, dict):
                errors[ensembl_id] = "Ensembl bulk lookup returned no record"
                continue
            try:
                chromosome = str(item["seq_region_name"])
                start, end, strand = int(item["start"]), int(item["end"]), int(item["strand"])
                if strand == 1:
                    left, right = max(1, start - PROMOTER_UPSTREAM), start + PROMOTER_DOWNSTREAM
                else:
                    left, right = max(1, end - PROMOTER_DOWNSTREAM), end + PROMOTER_UPSTREAM
                locations[ensembl_id] = {"region": f"{chromosome}:{left}..{right}:1"}
            except (KeyError, TypeError, ValueError) as exc:
                errors[ensembl_id] = f"Invalid Ensembl location: {exc}"
    region_to_ids: dict[str, list[str]] = defaultdict(list)
    for ensembl_id, item in locations.items():
        region_to_ids[item["region"]].append(ensembl_id)
    for regions in chunked(list(region_to_ids), 50):
        try:
            result = http_post_json(session, f"{ENSEMBL_REST}/sequence/region/human", {"regions": regions}).json()
        except SourceError as exc:
            for region in regions:
                for ensembl_id in region_to_ids[region]:
                    errors[ensembl_id] = str(exc)
            continue
        returned: set[str] = set()
        for item in result:
            region = item.get("query")
            sequence = str(item.get("seq") or "").upper().replace("\n", "")
            if region not in region_to_ids or not sequence or not re.fullmatch(r"[ACGTN]+", sequence):
                continue
            for ensembl_id in region_to_ids[region]:
                sequences[ensembl_id] = sequence
                (cache_dir / f"{ensembl_id}.fa").write_text(sequence + "\n", encoding="utf-8")
                returned.add(ensembl_id)
        for region in regions:
            for ensembl_id in region_to_ids[region]:
                if ensembl_id not in returned:
                    errors[ensembl_id] = "Ensembl bulk sequence returned no promoter"
    return sequences, errors


def select_negatives(
    session: requests.Session,
    gtrd_by_tf: dict[str, dict[str, dict[str, Any]]],
    chip_by_tf: dict[str, dict[str, dict[str, Any]]],
    all_genes: dict[str, str],
    raw_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    promoter_cache = raw_dir / "ensembl_promoters"
    jaspar_cache = raw_dir / "jaspar"
    promoter_cache.mkdir(parents=True, exist_ok=True)
    jaspar_cache.mkdir(parents=True, exist_ok=True)
    negatives: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    manual_pairs = {(item["tf"], item["gene"]) for item in RECENT_LITERATURE_PAIRS}
    for tf, config in TF_CONFIG.items():
        matrix = fetch_jaspar_matrix(session, tf, TF_CONFIG[tf]["jaspar_matrix_id"], jaspar_cache)
        candidates = [
            (gene, ensembl_id) for gene, ensembl_id in all_genes.items()
            if gene not in gtrd_by_tf[tf]
            and chip_by_tf[tf].get(gene, {}).get("chip_atlas_average", 0.0) <= 0
            and valid_gene_symbol(gene)
            and (tf, gene) not in manual_pairs
        ]
        # Stable pseudo-random ordering avoids cherry-picking while remaining reproducible.
        rng = random.Random(f"{RUN_ID}:{tf}")
        rng.shuffle(candidates)
        picked = 0
        candidate_index = 0
        # Use Ensembl bulk lookup/region endpoints for the promoter/PWM phase,
        # then retain the PubMed rate limit only for the motif-free subset.
        while picked < NEGATIVE_PER_TF and candidate_index < len(candidates):
            batch = candidates[candidate_index:candidate_index + 96]
            candidate_index += len(batch)
            sequence_by_id, errors_by_id = fetch_promoters_bulk(session, [ensembl_id for _, ensembl_id in batch], promoter_cache)
            for gene, ensembl_id in sorted(batch, key=lambda item: (item[0], item[1])):
                if picked >= NEGATIVE_PER_TF:
                    break
                row: dict[str, Any] = {
                    "tf": tf,
                    "gene": gene,
                    "ensembl_id": ensembl_id,
                    "gtrd_absent": True,
                    "chip_atlas_average": 0.0,
                    "chip_atlas_nonzero_experiments": 0,
                    "chip_atlas_independent_nonzero_experiments": 0,
                    "chip_atlas_independent_experiment_universe": matrix.get("chip_atlas_independent_experiment_universe", ""),
                    "jaspar_matrix_id": matrix["matrix_id"],
                    "jaspar_source_url": matrix["source_url"],
                    "jaspar_relative_threshold": JASPAR_RELATIVE_THRESHOLD,
                    "pubmed_checked_at_utc": utc_now(),
                }
                if ensembl_id in errors_by_id or ensembl_id not in sequence_by_id:
                    row["screen_status"] = "rejected_source_error"
                    row["error"] = errors_by_id.get(ensembl_id, "No promoter sequence returned")
                    audit.append(row)
                    continue
                hits, raw_threshold = pwm_hits(sequence_by_id[ensembl_id], matrix)
                row["jaspar_pwm_raw_threshold"] = raw_threshold
                row["jaspar_pwm_hits"] = hits
                if row.get("jaspar_pwm_hits", 0) != 0:
                    row["screen_status"] = "rejected_jaspar_hit"
                    audit.append(row)
                    continue
                try:
                    count, query = pubmed_count(session, tf, row["gene"])
                except SourceError as exc:
                    row["screen_status"] = "rejected_pubmed_error"
                    row["error"] = str(exc)
                    audit.append(row)
                    continue
                row["pubmed_query"] = query
                row["pubmed_count"] = count
                if count != 0:
                    row["screen_status"] = "rejected_pubmed_cooccurrence"
                    audit.append(row)
                else:
                    row["screen_status"] = "accepted_negative"
                    audit.append(row)
                    negatives.append({
                        "label": 0,
                        "subset": "strict_negative",
                        "tf": tf,
                        "gene": row["gene"],
                        "ensembl_id": row["ensembl_id"],
                        "gtrd_site_count": 0,
                        "chip_atlas_average": 0.0,
                        "chip_atlas_nonzero_experiments": 0,
                        "chip_atlas_independent_nonzero_experiments": 0,
                        "chip_atlas_independent_experiment_universe": 0,
                        "jaspar_matrix_id": matrix["matrix_id"],
                        "jaspar_pwm_hits": 0,
                        "jaspar_relative_threshold": JASPAR_RELATIVE_THRESHOLD,
                        "pubmed_query": query,
                        "pubmed_count": 0,
                        "pmid": "",
                        "publication_year": "",
                        "evidence_summary": "Strict negative: absent from the TF-specific GTRD promoter target list, ChIP-Atlas score=0, no high-stringency JASPAR PWM hit in the promoter window, and zero PubMed title/abstract co-occurrence.",
                        "selection_rule": "gtrd_absent AND chip_atlas_average=0 AND jaspar_pwm_hits=0 AND pubmed_count=0",
                    })
                    picked += 1
                # NCBI asks that unauthenticated E-utilities users stay below 3 calls/s.
                time.sleep(0.36)
            print(f"{tf}: strict negatives accepted {picked}/{NEGATIVE_PER_TF}; promoter candidates screened {candidate_index}", flush=True)
        if picked != NEGATIVE_PER_TF:
            raise SourceError(f"Only {picked}/{NEGATIVE_PER_TF} strict negatives were found for {tf}")
    return negatives, audit


def build(output_dir: Path, reuse_raw_from: Path | None = None) -> None:
    if output_dir.exists():
        raise SourceError(f"Output directory already exists; refusing to overwrite preserved data: {output_dir}")
    raw_dir = output_dir / "raw"
    derived_dir = output_dir / "derived"
    logs_dir = output_dir / "logs"
    manifest: list[dict[str, Any]] = []
    session = requests.Session()
    session.headers.update({"Accept": "application/json, text/plain, */*"})
    try:
        if reuse_raw_from:
            prior_raw = reuse_raw_from / "raw"
            if not prior_raw.is_dir():
                raise SourceError(f"Reusable raw-data directory is missing: {prior_raw}")
            output_dir.mkdir(parents=True, exist_ok=False)
            shutil.copytree(prior_raw, raw_dir, copy_function=os.link)
        else:
            raw_dir.mkdir(parents=True, exist_ok=False)
        # Keep source documentation and the two GTRD metadata tables needed to
        # audit release provenance and identify ENCODE-labelled experiments.
        for name in ("README",):
            download(session, f"{GTRD_ROOT}/intervals/target_genes/{name}", raw_dir / "gtrd" / name, manifest)
        for name in ("ChIP-seq.metadata.txt", "external_references.metadata.txt"):
            download(session, f"{GTRD_ROOT}/metadata/{name}", raw_dir / "gtrd" / "metadata" / name, manifest)

        gtrd_by_tf: dict[str, dict[str, dict[str, Any]]] = {}
        chip_by_tf: dict[str, dict[str, dict[str, Any]]] = {}
        chip_urls: dict[str, str] = {}
        for tf, config in TF_CONFIG.items():
            gtrd_file = download(session, gtrd_target_url(config["uniprot"]), raw_dir / "gtrd" / "target_genes" / f"{tf}_{config['uniprot']}.tsv", manifest)
            gtrd_by_tf[tf] = load_gtrd_targets(gtrd_file)
            source_url = chip_atlas_url(session, tf)
            chip_urls[tf] = source_url
            chip_file = download(session, source_url, raw_dir / "chip_atlas" / f"{tf}_hg38_tss1kb.tsv", manifest)
            chip_by_tf[tf] = load_chip_atlas_targets(chip_file, local_chip_atlas_experiment_ids(tf))

        # Source XML for each manually curated PMID is retained verbatim.
        for pmid in sorted({item["pmid"] for item in RECENT_LITERATURE_PAIRS}):
            fetch_pubmed_xml(session, pmid, raw_dir / "pubmed")
        manual_pairs = [dict(item) for item in RECENT_LITERATURE_PAIRS]
        write_csv(derived_dir / "manual_recent_literature_pairs.csv", manual_pairs, list(manual_pairs[0]))

        manual_pair_keys = {(item["tf"], item["gene"]) for item in manual_pairs}
        gtrd_positives = select_gtrd_positives(gtrd_by_tf, chip_by_tf, manual_pair_keys)
        if len(gtrd_positives) != len(TF_CONFIG) * POSITIVE_PER_TF:
            raise SourceError(f"Expected 160 GTRD positives, got {len(gtrd_positives)}")

        # Create canonical entries for recent manual positives.  Add GTRD and
        # ChIP-Atlas values where available, but do not use them as selection
        # criteria for this literature subset.
        literature_positives: list[dict[str, Any]] = []
        for item in manual_pairs:
            gtrd = gtrd_by_tf.get(item["tf"], {}).get(item["gene"], {})
            chip = chip_by_tf.get(item["tf"], {}).get(item["gene"], {})
            literature_positives.append({
                "label": 1,
                "subset": "recent_literature_manual",
                "tf": item["tf"],
                "gene": item["gene"],
                "ensembl_id": gtrd.get("ensembl_id", ""),
                "gtrd_site_count": gtrd.get("gtrd_site_count", ""),
                "chip_atlas_average": chip.get("chip_atlas_average", ""),
                "chip_atlas_nonzero_experiments": chip.get("chip_atlas_nonzero_experiments", ""),
                "chip_atlas_independent_nonzero_experiments": chip.get("chip_atlas_independent_nonzero_experiments", ""),
                "chip_atlas_independent_experiment_universe": chip.get("chip_atlas_independent_experiment_universe", ""),
                "jaspar_matrix_id": "",
                "jaspar_pwm_hits": "",
                "jaspar_relative_threshold": "",
                "pubmed_query": f"PMID:{item['pmid']}",
                "pubmed_count": 1,
                "pmid": item["pmid"],
                "publication_year": item["year"],
                "evidence_summary": item["evidence_summary"],
                "selection_rule": "Manual 2023–2025 publication curation; source reports binding/occupancy plus expression, perturbation, or functional evidence.",
            })
        if len(literature_positives) != 18:
            raise SourceError(f"Expected 18 manually curated positives, got {len(literature_positives)}")

        all_genes: dict[str, str] = {}
        for targets in gtrd_by_tf.values():
            for gene, source in targets.items():
                if valid_gene_symbol(gene):
                    all_genes.setdefault(gene, source["ensembl_id"])
        negatives, negative_audit = select_negatives(session, gtrd_by_tf, chip_by_tf, all_genes, raw_dir)
        if len(negatives) != len(TF_CONFIG) * NEGATIVE_PER_TF:
            raise SourceError(f"Expected 192 negatives, got {len(negatives)}")

        positives = gtrd_positives + literature_positives
        if len(positives) != 178:
            raise SourceError(f"Expected 178 positives, got {len(positives)}")
        combined = []
        for index, row in enumerate(positives + negatives, 1):
            combined.append({"record_id": f"{RUN_ID}_{index:04d}", "species": SPECIES, "used_in_training": "false", **row})
        # Explicit checks prevent accidental class imbalance or pair duplication.
        if sum(int(row["label"]) for row in combined) != 178 or sum(int(row["label"]) == 0 for row in combined) != 192:
            raise SourceError("Final class counts do not equal 178 positive and 192 negative records")
        pair_keys = [(row["tf"], row["gene"]) for row in combined]
        if len(pair_keys) != len(set(pair_keys)):
            raise SourceError("Duplicate TF–gene pair found in final validation set")

        fields = [
            "record_id", "label", "subset", "species", "tf", "gene", "ensembl_id", "used_in_training",
            "gtrd_site_count", "chip_atlas_average", "chip_atlas_nonzero_experiments", "jaspar_matrix_id",
            "chip_atlas_independent_nonzero_experiments", "chip_atlas_independent_experiment_universe",
            "jaspar_pwm_hits", "jaspar_relative_threshold", "pubmed_query", "pubmed_count", "pmid",
            "publication_year", "evidence_summary", "selection_rule",
        ]
        combined.sort(key=lambda row: (int(row["label"]) * -1, row["subset"], row["tf"], row["gene"]))
        write_csv(derived_dir / "combined_test_set_178_positive_192_negative.csv", combined, fields)
        write_csv(derived_dir / "negative_screening_audit.csv", negative_audit, sorted({key for row in negative_audit for key in row}))
        write_json(derived_dir / "chip_atlas_urls.json", chip_urls)
        write_json(
            derived_dir / "local_chipatlas_accessions_excluded.json",
            {tf: sorted(local_chip_atlas_experiment_ids(tf)) for tf in TF_CONFIG},
        )

        summary = {
            "run_id": RUN_ID,
            "created_at_utc": utc_now(),
            "model_data_isolation": "The GeneReg rule-based score is not trained on this validation directory. No labels from this directory are read by app.py.",
            "positive_count": len(positives),
            "positive_gtrd_independent_count": len(gtrd_positives),
            "positive_recent_literature_manual_count": len(literature_positives),
            "negative_count": len(negatives),
            "total_count": len(combined),
            "positive_rule": "GTRD subset: >=2 promoter meta-clusters plus >=2 non-zero ChIP-Atlas SRX experiments whose accessions are absent from the live GeneReg peak corpus. Recent-literature subset: manual 2024–2025 curation with reported binding/occupancy and expression, perturbation, or functional evidence.",
            "negative_rule": "Absent from TF-specific GTRD promoter targets, ChIP-Atlas target score=0, zero high-stringency JASPAR PWM promoter hits, and zero PubMed title/abstract co-occurrence.",
            "jaspar_relative_pwm_threshold": JASPAR_RELATIVE_THRESHOLD,
            "tf_distribution": {
                tf: {
                    "gtrd_positive": sum(row["subset"] == "GTRD_independent" and row["tf"] == tf for row in combined),
                    "negative": sum(row["subset"] == "strict_negative" and row["tf"] == tf for row in combined),
                }
                for tf in TF_CONFIG
            },
        }
        write_json(derived_dir / "validation_summary.json", summary)
        write_json(derived_dir / "source_manifest.json", manifest)
        (output_dir / "README.md").write_text(
            "# Expanded GeneReg independent validation set\n\n"
            "This directory is immutable input/output for the expanded validation set. It is kept outside the active GeneReg scoring corpus (`encode_peaks/`).\n\n"
            "- `raw/`: downloaded GTRD, ChIP-Atlas, JASPAR, Ensembl, and PubMed source data.\n"
            "- `derived/combined_test_set_178_positive_192_negative.csv`: frozen joint test set.\n"
            "- `derived/negative_screening_audit.csv`: all examined negative candidates and rejection reasons.\n"
            "- `derived/source_manifest.json`: retrieval URL, timestamp, byte count and SHA-256 for every downloaded source file.\n"
            "- `derived/local_chipatlas_accessions_excluded.json`: SRX accessions already in GeneReg's peak corpus and excluded from the GTRD positive confirmation.\n\n"
            "The current GeneReg application is rule-based and does not read this directory; the test set is therefore held out from its scoring corpus.\n",
            encoding="utf-8",
        )
        print(json.dumps(summary, indent=2))
    except Exception:
        # Preserve downloaded inputs for audit, while making failed builds easy
        # to distinguish from a successful frozen release.
        (output_dir / "BUILD_FAILED").write_text("Build failed; inspect logs and raw/ before rerunning with a new output directory.\n", encoding="utf-8")
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=f"/www/wwwroot/gene_reg/validation/{RUN_ID}")
    parser.add_argument("--reuse-raw-from", help="Existing failed/partial run whose raw directory will be hard-linked into the new run")
    args = parser.parse_args()
    build(Path(args.output_dir), Path(args.reuse_raw_from) if args.reuse_raw_from else None)


if __name__ == "__main__":
    main()
