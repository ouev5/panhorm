"""Frozen manuscript scoring rule; never fit normalization on a user query."""
import json
import math
from pathlib import Path

MODEL_CONFIG = json.loads(Path(__file__).with_name("model_config.json").read_text(encoding="utf-8"))
WEIGHTS = MODEL_CONFIG["weights"]
NORMALIZATION = MODEL_CONFIG["normalization"]
THRESHOLD = MODEL_CONFIG["threshold_normalized"]


def score_evidence(chip: float, motif: float, literature: float) -> dict:
    """Return the weighted score and label using development-frozen parameters.

    Normalized values outside [0, 1] are retained: new observations may fall
    outside the development range. Decisions use full precision, not rounded
    display values. The score is an evidence priority, not a probability.
    """
    components = {"chip": float(chip), "motif": float(motif), "literature": float(literature)}
    if any(not math.isfinite(value) or not 0 <= value <= 100 for value in components.values()):
        raise ValueError("Evidence component scores must be finite values between 0 and 100.")
    raw = math.fsum(WEIGHTS[name] * value for name, value in components.items())
    normalized = (raw - NORMALIZATION["min"]) / (NORMALIZATION["max"] - NORMALIZATION["min"])
    return {
        "weighted_score": raw,
        "normalized_score": normalized,
        # Algebraically equivalent raw comparison avoids subtract/divide
        # roundoff rejecting a value exactly at the inclusive boundary.
        "prediction": int(raw >= MODEL_CONFIG["threshold_raw"]),
        "threshold_normalized": THRESHOLD,
        "threshold_raw": MODEL_CONFIG["threshold_raw"],
        "weights": dict(WEIGHTS),
        "model_id": MODEL_CONFIG["model_id"],
    }
