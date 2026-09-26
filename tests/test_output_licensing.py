from biocore.providers.base import Finding, Tier, Category
from geneask.annotators import alphagenome_vep as ag, alphamissense as am, alphagenome_atlas as atlas


def finding():
    return Finding("1-100-A-G", "variant_lookup", "d", Tier.SPECULATIVE, [Category.CLINICAL],
                   detail={"research_candidate": True})


def test_commercial_ag_am_gates_before_cache_or_mirror(monkeypatch):
    monkeypatch.setenv("DNAREPORT_OUTPUT_MODE", "commercial")
    monkeypatch.setenv("ALPHAGENOME_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.setattr(ag, "_cache_con", lambda *a: (_ for _ in ()).throw(AssertionError("cache accessed")))
    assert ag.score_variant("1-100-A-G", "unused") is None
    p = ag.Pacing()
    assert ag.annotate_findings([finding()], pacing=p) == 0
    assert p.as_dict()["status"] == "license_blocked" and p.as_dict()["reason"]
    assert am.lookup("1-100-A-G") is None
    assert am.findings_for("1-100-A-G") == []
    status = {}
    assert am.annotate_findings([finding()], status=status) == 0
    assert status["status"] == "license_blocked"


def test_commercial_atlas_ignores_api_cache_even_avi_only(monkeypatch, tmp_path):
    monkeypatch.setenv("DNAREPORT_OUTPUT_MODE", "commercial")
    monkeypatch.setattr(atlas, "_cache_read", lambda *a: (_ for _ in ()).throw(AssertionError("cache accessed")))
    monkeypatch.setattr(atlas, "_run_remote", lambda *a: (_ for _ in ()).throw(AssertionError("API accessed")))
    monkeypatch.setattr(atlas, "lookup_local_avi", lambda *a: {"status": "unavailable"})
    result = atlas.query_variant("1-100-A-G", requested_scorers=["AVI_SCORE"])
    assert result["status"] == "license_blocked" and result["license"]["commercial_allowed"] is False


def test_commercial_atlas_retains_only_local_avi(monkeypatch):
    monkeypatch.setenv("DNAREPORT_OUTPUT_MODE", "commercial")
    monkeypatch.setattr(atlas, "_run_remote", lambda *a: (_ for _ in ()).throw(AssertionError("API accessed")))
    monkeypatch.setattr(atlas, "lookup_local_avi", lambda vid: {"variant_id": vid, "status": "complete", "avi_score": .7,
        "provenance": "alphagenome_atlas_local_avi", "tracks": [{"scorer": "AVI_SCORE", "raw_score": .7}]})
    result = atlas.query_variant("1-100-A-G")
    assert result["avi_score"] == .7 and result["remote_status"] == "license_blocked"
    assert result["missing_scorers"] == ["AVI_SCORE_FEATURE_IMPORTANCE"]
