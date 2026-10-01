import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

import app
from evidence import promoter_region, scan_pwm
from reproduce_manuscript import build_summary, read_corpus
from scoring import score_evidence
from rerun_415_leakfree import merge_415

ROOT = Path(__file__).resolve().parents[1]


def test_all_published_predictions_and_test_metrics():
    summary = build_summary(read_corpus())
    assert summary["metrics"]["f1"] == pytest.approx(0.8941176470588236)
    assert summary["metrics"]["roc_auc"] == pytest.approx(0.9229651162790697)
    assert summary["metrics"]["pr_auc"] == pytest.approx(0.9240180870185092)
    assert summary == json.loads((ROOT / "analysis_v2.6/manuscript_metrics.json").read_text())


def test_published_rows_match_frozen_source_component_scores():
    source = {row["record_id"]: row for row in merge_415()}
    rows = read_corpus()
    assert len(source) == len(rows) == 415
    for row in rows:
        for key in ("tf", "gene", "label", "chip_score", "motif_score", "literature_score"):
            assert row[key] == source[row["record_id"]][key]


def test_decision_uses_unrounded_normalization():
    # Both display as 41.925 to three decimals, but lie on opposite sides.
    assert score_evidence(0, (41.925 - 0.0001) / 0.7, 0)["prediction"] == 0
    assert score_evidence(0, (41.925 + 0.0001) / 0.7, 0)["prediction"] == 1
    assert score_evidence(41.925, 41.925, 41.925)["prediction"] == 1
    assert score_evidence(0, 0, 0)["normalized_score"] < 0
    assert score_evidence(100, 100, 100)["normalized_score"] > 1


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, 101])
def test_invalid_evidence_is_rejected(value):
    with pytest.raises(ValueError):
        score_evidence(value, 35, 0)


def stub_evidence(monkeypatch, chip=100, motif=60, literature=75):
    monkeypatch.setattr(app.gf, "fetch", lambda _: {"symbol": "UBB", "chr": "1", "start": 10000, "end": 20000, "tss": 10000})
    monkeypatch.setattr(app.ca, "analyze", lambda *args: {"score": chip, "peaks": [], "peaks_found": 0})
    monkeypatch.setattr(app.ma, "analyze", lambda *args: {"score": motif, "promoter_matches": []})
    monkeypatch.setattr(app.la, "analyze", lambda *args: {"score": literature, "articles": []})
    monkeypatch.setattr(app, "_chip_atlas_stats", lambda *args: {"avg": 10000, "max": 10000, "nz": 100, "tot": 100})


def test_api_matches_weighted_rule_and_ai_cannot_change_it(monkeypatch):
    stub_evidence(monkeypatch)
    calls = []
    monkeypatch.setattr(app.ai_a, "analyze", lambda *args: calls.append(args) or {"confidence": "Strong", "score_adjustment": 100})
    client = TestClient(app.app)
    default = client.get("/api/analyze", params={"gene": "UBB", "receptor": "NR3C1"})
    assert default.status_code == 200
    assert not calls
    data = default.json()
    # No extra total-score penalty for a housekeeping target or ChIP-Atlas bonus.
    assert data["raw_score"] == 67
    assert data["normalized_score"] == 1
    assert data["prediction"] == 1
    explained = client.post("/api/analyze", json={"gene": "UBB", "receptor": "NR3C1", "include_ai": True}).json()
    assert calls
    assert explained["ai_bonus"] == 0
    for key in ("raw_score", "final_score", "normalized_score", "prediction", "weights"):
        assert explained[key] == data[key]


def test_validation_api_presents_83_pair_internal_test():
    data = TestClient(app.app).get("/api/validation").json()
    assert data["corpus"]["records"] == 415
    assert data["corpus"]["hash_matches_manifest"]
    assert data["development_set"]["records"] == 332
    assert data["test_set"] == {"records": 83, "positives": 40, "evidence_absent": 43}
    assert data["model"]["threshold_normalized"] == 0.41


def test_negative_prediction_is_not_an_evidence_level_cutoff(monkeypatch):
    stub_evidence(monkeypatch, chip=0, motif=59, literature=0)
    data = TestClient(app.app).get("/api/analyze", params={"gene": "GENE", "receptor": "NR3C1"}).json()
    assert data["raw_score"] == 41.3  # Above 40, but below 41.925.
    assert data["prediction"] == 0
    assert data["classification"] == "Below prioritization threshold"


@pytest.mark.parametrize("strand,expected", [(1, "1:8001..10500:1"), (-1, "1:19501..22000:-1")])
def test_promoter_window_coordinates_and_orientation(strand, expected):
    region, offset = promoter_region({"chr": "1", "start": 10000, "end": 20000, "strand": strand})
    assert region == expected
    assert offset == 2000


def test_pwm_matches_both_strands_and_skips_ambiguous_bases():
    # A deliberately asymmetric ACG motif makes strand orientation testable.
    pfm = {"A": [100, 0, 0], "C": [0, 100, 0], "G": [0, 0, 100], "T": [0, 0, 0]}
    matches = scan_pwm("ACGNNCGT", pfm, "test", tss_offset=3)
    assert [(match["position"], match["strand"]) for match in matches] == [(-3, "+"), (2, "-")]
    assert all(match["relative_score"] >= 0.95 for match in matches)
    assert scan_pwm("AAAAAA", pfm, "test") == []


def test_literature_uses_total_count_not_retmax(monkeypatch):
    calls = []
    xml = '<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>1</PMID><Article><ArticleTitle>TF and GENE</ArticleTitle><Journal><JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue></Journal></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'
    def request(url, **kwargs):
        calls.append((url, kwargs))
        if "esearch" in url:
            return SimpleNamespace(raise_for_status=lambda: None, json=lambda: {"esearchresult": {"count": "35", "idlist": ["1"]}})
        return SimpleNamespace(text=xml)
    monkeypatch.setattr(app.requests, "get", request)
    result = app.LiteratureAnalyzer().analyze("GENE", "TF")
    assert result["articles_found"] == 35
    assert result["score"] == 90
    assert result["articles"][0]["year"] == "2020"
    assert '"GENE"[Title/Abstract]' in calls[0][1]["params"]["term"]
    assert "2010/01/01" in result["query"] and "2025/12/31" in result["query"]


def test_motif_analyzer_uses_pwm_and_promoter_endpoint(monkeypatch, tmp_path):
    calls = []
    pfm = {"A": [100, 0, 0], "C": [0, 100, 0], "G": [0, 0, 100], "T": [0, 0, 0]}
    def request(url, **kwargs):
        calls.append((url, kwargs))
        if "/sequence/region/" in url:
            payload = {"seq": "N" * 2000 + "ACG" + "N" * 497}
        elif url.endswith("/matrix/"):
            payload = {"results": [{"matrix_id": "test.1", "name": "TF", "collection": "CORE"}]}
        else:
            payload = {"pfm": pfm}
        return SimpleNamespace(status_code=200, raise_for_status=lambda: None, json=lambda: payload)
    monkeypatch.setattr(app.requests, "get", request)
    monkeypatch.setattr(app, "PROMOTER_CACHE_DIR", tmp_path)
    result = app.MotifAnalyzer().analyze("GENE", "TF", {"chr": "1", "start": 10000, "end": 20000, "strand": 1, "ensembl_id": "ENSGtest"})
    assert result["promoter_length"] == 2500
    assert result["promoter_matches"][0]["position"] == 0
    assert result["score"] == 41
    assert any("/sequence/region/human/1:8001..10500:1" in url for url, _ in calls)
    assert calls[0][1]["params"]["release"] == "2024"
