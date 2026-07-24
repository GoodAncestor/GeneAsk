"""gnomAD AF annotator: cache behavior + in-place finding annotation (no network)."""
import sqlite3
from geneask.annotators.gnomad_freq import allele_frequency, annotate_findings, _cache_con
from biocore.providers.base import Finding, Tier, Category


def test_cache_hit_no_network(tmp_path):
    db = str(tmp_path / "af.db")
    con = _cache_con(db)
    con.execute("INSERT INTO af VALUES (?,?)", ("1-100-A-G", 0.35)); con.commit(); con.close()
    # cached value returned without any API call
    assert allele_frequency("1-100-A-G", cache_db=db) == 0.35
    # cached miss (-1) returns None, also no API call
    con = _cache_con(db); con.execute("INSERT INTO af VALUES (?,?)", ("1-200-C-T", -1.0)); con.commit(); con.close()
    assert allele_frequency("1-200-C-T", cache_db=db) is None


def test_annotate_only_variant_markers(tmp_path):
    db = str(tmp_path / "af.db")
    con = _cache_con(db)
    con.execute("INSERT INTO af VALUES (?,?)", ("1-100-A-G", 0.35)); con.commit(); con.close()
    fs = [Finding("1-100-A-G", "clinvar", "gene X: pathogenic", Tier.MODERATE, [Category.CLINICAL], detail={}),
          Finding("cg999", "ewas", "methylation", Tier.ROBUST, [Category.AGING], detail={})]
    n = annotate_findings(fs, cache_db=db)
    assert n == 1
    assert fs[0].detail["gnomad_af"] == 0.35
    assert "35.0% of people" in fs[0].description
    assert "gnomad_af" not in fs[1].detail    # CpG marker untouched
