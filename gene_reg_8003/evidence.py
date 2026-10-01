"""Sequence helpers matching the manuscript promoter window and PWM cutoff."""
import math

UPSTREAM = 2000
DOWNSTREAM = 500
RELATIVE_THRESHOLD = 0.95


def promoter_region(gene_info):
    """Ensembl 1-based inclusive region, oriented in the gene's direction."""
    strand = int(gene_info.get("strand", 1))
    # MyGene genomic_pos coordinates are 0-based, half-open.
    tss = int(gene_info["start"]) + 1 if strand == 1 else int(gene_info["end"])
    if strand == 1:
        left, right = max(1, tss - UPSTREAM), tss + DOWNSTREAM - 1
    else:
        left, right = max(1, tss - DOWNSTREAM + 1), tss + UPSTREAM
    return f"{gene_info['chr']}:{left}..{right}:{strand}", (tss - left if strand == 1 else right - tss)


def scan_pwm(sequence, pfm, matrix_id, tss_offset=UPSTREAM):
    """Scan both strands with 0.1 pseudocount and a uniform background."""
    columns = list(zip(*(pfm.get(base, []) for base in "ACGT")))
    if not columns or len({len(pfm.get(base, [])) for base in "ACGT"}) != 1:
        raise ValueError("JASPAR matrix has an incomplete position-frequency matrix.")
    weights = []
    for counts in columns:
        total = sum(float(count) for count in counts) + 0.4
        weights.append({base: math.log2((float(count) + 0.1) / total / 0.25)
                        for base, count in zip("ACGT", counts)})
    low = sum(min(column.values()) for column in weights)
    high = sum(max(column.values()) for column in weights)
    if high <= low:
        return []
    threshold = low + RELATIVE_THRESHOLD * (high - low)
    size = len(weights)
    matches = []
    sequence = sequence.upper()
    for index in range(len(sequence) - size + 1):
        window = sequence[index:index + size]
        if any(base not in "ACGT" for base in window):
            continue
        reverse = window.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for strand, candidate in (("+", window), ("-", reverse)):
            score = sum(column[base] for column, base in zip(weights, candidate))
            if score >= threshold:
                matches.append({"motif_name": matrix_id, "motif_seq": candidate,
                                "position": index - tss_offset, "strand": strand,
                                "pwm_score": round(score, 6),
                                "relative_score": round((score - low) / (high - low), 6)})
    return matches
