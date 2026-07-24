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
