"""Regenerate checksums for the GeneReg manuscript release artifacts."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "release_manifest.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    previous = json.loads(MANIFEST.read_text())
    names = {entry["path"] for entry in previous["files"]}
    names.update({"README.md", "README_v2.6.md", "model_config.json", "scoring.py", "evidence.py",
                  "requirements.txt", "requirements-analysis.txt",
                  "analysis_v2.6/manuscript_benchmark_415.csv",
                  "analysis_v2.6/manuscript_metrics.json"})
    names.update(str(path.relative_to(ROOT)) for path in (ROOT / "tests").glob("*.py"))
    names.update(str(path.relative_to(ROOT)) for path in (ROOT / "analysis_v2.6/scripts").glob("*.py"))
    names.update(str(path.relative_to(ROOT)) for path in (ROOT / "analysis_v2.6/results").glob("*.xlsx"))
    manifest = {
        "release_id": "AHormoneDB-GeneReg-v2.6-weighted-evidence",
        "model": json.loads((ROOT / "model_config.json").read_text()),
        "app_sha256": sha256(ROOT / "app.py"),
        "benchmark_sha256": sha256(ROOT / "analysis_v2.6/manuscript_benchmark_415.csv"),
        "validation_output_sha256": sha256(ROOT / "analysis_v2.6/manuscript_metrics.json"),
        "manuscript_sha256": "b1814954c502ec34d31d8f8ce3e56d32ff01ae23c789e4740baf98321122abd2",
        "files": [{"path": name, "bytes": (ROOT / name).stat().st_size, "sha256": sha256(ROOT / name)} for name in sorted(names)],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Updated {MANIFEST} with {len(manifest['files'])} artifact checksums.")


if __name__ == "__main__":
    main()
