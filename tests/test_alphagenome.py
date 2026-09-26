"""AlphaGenome enrichment: no-op-when-absent + cache behavior (no key, no network)."""
import os, sqlite3
from geneask.annotators.alphagenome_vep import (annotate_findings, _enabled, _parse_vid,
                                                 score_variant, _is_uncertain)
from biocore.providers.base import Finding, Tier, Category


def test_disabled_without_flag(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_ENABLED", raising=False)
    monkeypatch.setenv("ALPHA_GENOME_KEY", "x")
    assert _enabled() is False           # key alone isn't enough; needs opt-in flag
    fs = [Finding("1-1-A-G", "clinvar", "d", Tier.SPECULATIVE, [Category.CLINICAL],
                  detail={"clinical_significance": "Uncertain significance"})]
    assert annotate_findings(fs) == 0    # no-op


def test_disabled_without_key(monkeypatch):
    monkeypatch.setenv("ALPHAGENOME_ENABLED", "1")
    monkeypatch.delenv("ALPHA_GENOME_KEY", raising=False)
    assert _enabled() is False


def test_parse_and_uncertain():
    assert _parse_vid("13-32316419-C-G") == ("chr13", 32316419, "C", "G")
    assert _parse_vid("cg123") is None
    unc = Finding("1-1-A-G", "s", "d", Tier.SPECULATIVE, [Category.CLINICAL],
                  detail={"clinical_significance": "Uncertain significance"})
    cert = Finding("1-1-A-G", "s", "d", Tier.ROBUST, [Category.CLINICAL],
                   detail={"clinical_significance": "Pathogenic"})
    assert _is_uncertain(unc) and not _is_uncertain(cert)


def test_cache_hit_no_network(tmp_path):
    # a pre-populated cache row is returned without any client/network call
    db = str(tmp_path / "ag.db")
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE ag(variant_id TEXT PRIMARY KEY, summary TEXT)")
    import json
    con.execute("INSERT INTO ag VALUES (?,?)", ("1-1-A-G", json.dumps({"variant_id": "1-1-A-G",
                "top_modality": "RNA_SEQ", "quantile_score": 0.8, "direction": "increase", "n_tracks": 5})))
    con.commit(); con.close()
    s = score_variant("1-1-A-G", api_key="unused", cache_db=db)
    assert s["top_modality"] == "RNA_SEQ" and s["quantile_score"] == 0.8


def test_resource_exhausted_is_recognised_structurally_and_by_message():
    """RESOURCE_EXHAUSTED is the ONLY quota signal AlphaGenome gives — they publish
    no numbers on purpose. The client wraps errors differently across versions, so
    missing it would mean silently hammering an API that is already refusing."""
    from geneask.annotators.alphagenome_vep import _is_resource_exhausted

    class _Code:
        name = "RESOURCE_EXHAUSTED"
    class _GrpcError(Exception):
        def code(self): return _Code()

    assert _is_resource_exhausted(_GrpcError("quota"))
    assert _is_resource_exhausted(RuntimeError("StatusCode.RESOURCE_EXHAUSTED: quota"))
    assert not _is_resource_exhausted(RuntimeError("connection reset"))


def test_quota_exhaustion_stops_the_run_and_caches_nothing(tmp_path, monkeypatch):
    """One refusal must end the run, not cost one variant. And the non-answer must
    not be cached: it says nothing about the variant, and caching it would make an
    upstream busy minute permanent."""
    from geneask.annotators.alphagenome_vep import Pacing
    monkeypatch.setenv("ALPHAGENOME_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "x")
    # The vendor client isn't installed in the test env, and an ImportError would
    # be indistinguishable from the refusal we're actually testing for. Stub the
    # modules score_variant imports so the failure under test is the API's.
    import sys, types
    for name in ("alphagenome", "alphagenome.data", "alphagenome.models"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    sys.modules["alphagenome.data"].genome = types.SimpleNamespace(
        Variant=lambda **k: None, Interval=lambda **k: None)
    sys.modules["alphagenome.models"].variant_scorers = types.SimpleNamespace(
        RECOMMENDED_VARIANT_SCORERS={}, tidy_scores=lambda s: None)
    # client_available() imports this name; without it the stub reads as "package
    # not installed" and the run would stop for the wrong reason.
    sys.modules["alphagenome.models"].dna_client = types.SimpleNamespace(create=lambda k: None)
    calls = []
    def boom(api_key):
        calls.append(api_key)
        raise RuntimeError("StatusCode.RESOURCE_EXHAUSTED: out of quota")
    monkeypatch.setattr("geneask.annotators.alphagenome_vep._client", boom)
    db = str(tmp_path / "ag.db")
    fs = [Finding(f"1-{i}-A-G", "clinvar", "d", Tier.SPECULATIVE, [Category.CLINICAL],
                  detail={"clinical_significance": "Uncertain significance"})
          for i in range(6)]
    pacing = Pacing()
    assert annotate_findings(fs, cache_db=db, pacing=pacing) == 0
    assert len(calls) == 1 and pacing.halted
    con = sqlite3.connect(db)
    assert con.execute("SELECT count(*) FROM ag").fetchone()[0] == 0
    con.close()


def test_no_variant_cap_by_default(monkeypatch):
    """The everyday control is the rate and the time budget, not a variant count —
    a flat cap truncates a big report for reasons unrelated to the upstream, and
    does it identically at two reports a day and at two thousand."""
    from geneask.annotators.alphagenome_vep import default_pacing
    monkeypatch.delenv("ALPHAGENOME_MAX_VARIANTS", raising=False)
    assert default_pacing(":memory:").max_variants is None
    monkeypatch.setenv("ALPHAGENOME_MAX_VARIANTS", "3")
    assert default_pacing(":memory:").max_variants == 3


def test_missing_client_is_surfaced_not_swallowed(monkeypatch):
    """score_variant catches every exception so one bad variant can't fail a
    report — which also swallows ImportError. With the feature switched ON and the
    vendor package absent, that made a deployment fault look identical to "no
    regulatory effects found". It must be reported instead."""
    import geneask.annotators.alphagenome_vep as ag
    monkeypatch.setenv("ALPHAGENOME_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "x")
    monkeypatch.setattr(ag, "client_available", lambda: False)
    p = ag.Pacing()
    fs = [Finding("1-1-A-G", "clinvar", "d", Tier.SPECULATIVE, [Category.CLINICAL],
                  detail={"clinical_significance": "Uncertain significance"})]
    assert ag.annotate_findings(fs, pacing=p) == 0
    assert p.client_missing and p.spent == 0    # flagged, and no slot wasted


def test_novel_candidates_and_coverage(monkeypatch):
    from geneask.annotators import alphagenome_vep as ag
    monkeypatch.setenv("ALPHAGENOME_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.setattr(ag, "client_available", lambda: True)
    monkeypatch.setattr(ag, "score_variant", lambda marker, *a, **kw:
                        {"top_modality": "SPLICE_SITES", "quantile_score": .98, "direction": "unknown"})
    def finding(marker, **detail):
        return Finding(marker, "novel_variant", "d", Tier.SPECULATIVE, [Category.CLINICAL], detail=detail)
    fs = [finding("1-100-A-G", novel_candidate=True),
          finding("1-100-A-G", novel_candidate=True),
          finding("1-200-A-G", novel_candidate=True),
          finding("1-300-A-G", novel_candidate=True, clinical_significance="Benign"),
          finding("1-400-A-G", novel_candidate=True, clinical_significance="Pathogenic")]
    p = ag.Pacing(max_variants=1)
    assert ag.annotate_findings(fs, pacing=p) == 2
    assert p.as_dict()["eligible"] == 2
    assert (p.scored, p.skipped, p.status) == (1, 1, "partial")
    assert "increase" not in fs[0].description
    assert "alphagenome" not in fs[-1].detail


def test_summary_unsigned_metadata_and_nonfinite():
    from geneask.annotators.alphagenome_vep import summarize_scores
    class Frame:
        columns = ["quantile_score"]
        def __len__(self): return 2
        def to_dict(self, orientation):
            return [{"quantile_score": .99, "raw_score": .4, "output_type": "SPLICE_SITES",
                     "gene_name": "TAL1", "biosample_name": "blood"},
                    {"quantile_score": float("nan")}]
    s = summarize_scores(Frame(), "1-100-A-G")
    assert s["direction"] == "unknown"
    assert s["gene_name"] == "TAL1" and s["biosample_name"] == "blood"
    assert s["n_tracks"] == 1 and len(s["top_tracks"]) == 1
    assert s["provenance"] == "alphagenome_api"


def test_legacy_cache_direction_is_untrusted(tmp_path):
    import json
    from geneask.annotators.alphagenome_vep import Pacing
    db = str(tmp_path / "cache.db")
    with sqlite3.connect(db) as con:
        con.execute("CREATE TABLE ag(variant_id TEXT PRIMARY KEY, summary TEXT)")
        con.execute("INSERT INTO ag VALUES (?,?)", ("1-100-A-G", json.dumps({"direction": "increase"})))
    p = Pacing()
    assert score_variant("1-100-A-G", "unused", db, p)["direction"] == "unknown"
    assert p.cache_hits == 1 and p.spent == 0
