"""Targeted ClinVar mirror lookups.

Screening one person's callset used to pull all ~4.2M mirror rows into a dict to
answer a few thousand questions — seconds of CPU and gigabytes of memory on every
report. These pin that the targeted path returns the same answers.
"""
import sqlite3
import pytest
from geneask.annotators import clinvar_mirror as cm


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    db = tmp_path / "clinvar.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO meta VALUES ('schema_version', '2')")
    con.execute("""CREATE TABLE variants(
        variant_id TEXT PRIMARY KEY, clinvar_variation_id TEXT, gene TEXT,
        clinical_significance TEXT, review_status TEXT, gold_stars INTEGER,
        conditions TEXT, condition_ids TEXT, molecular_consequence TEXT,
        origin TEXT, allele_id TEXT)""")
    con.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)", [
        ("1-100-A-G", "v1", "BRCA1", "Pathogenic", "reviewed by expert panel", 3,
         "[]", "[]", None, "[]", None),
        ("2-200-C-T", "v2", "MTHFR", "Benign", "criteria provided", 1,
         "[]", "[]", None, "[]", None),
        ("3-300-G-A", "v3", "TP53", "Likely pathogenic", "criteria provided", 2,
         "[]", "[]", None, "[]", None),
    ])
    con.commit()
    con.close()
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(db))
    return str(db)


def test_targeted_lookup_agrees_with_loading_everything(mirror):
    full = cm.load_panel_from_mirror()
    hit = cm.lookup_from_mirror(["1-100-A-G", "3-300-G-A"])
    assert hit["1-100-A-G"] == full["1-100-A-G"]
    assert hit["3-300-G-A"] == full["3-300-G-A"]
    assert "2-200-C-T" not in hit


def test_no_mirror_is_None_not_empty(tmp_path, monkeypatch):
    """The caller falls back to the bundled panel on None. Returning {} instead
    would read as 'the mirror says this person has nothing', which is the
    opposite of 'there is no mirror' — and would silently drop the panel."""
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(tmp_path / "absent.db"))
    assert cm.lookup_from_mirror(["1-100-A-G"]) is None


def test_mirror_present_but_nothing_matched_is_empty_not_None(mirror):
    assert cm.lookup_from_mirror(["9-999-T-C"]) == {}


def test_more_variants_than_sqlites_parameter_limit(mirror):
    ids = [f"9-{i}-A-G" for i in range(3000)] + ["1-100-A-G"]
    got = cm.lookup_from_mirror(ids)
    assert set(got) == {"1-100-A-G"}


def test_screening_uses_the_mirror_without_loading_the_bundled_panel(mirror, monkeypatch):
    """The bundled 393k-variant panel was gunzipped and JSON-parsed on every
    report and then discarded, because the mirror superseded it."""
    from geneask.interpret import clinvar_screen as cs
    called = []
    monkeypatch.setattr(cs, "load_panel", lambda *a, **k: called.append(1) or {})
    out = cs.screen_findings([{"variant_id": "1-100-A-G", "genotype": "A/G",
                               "platform": "WGS"}])
    assert not called, "bundled panel was loaded despite the mirror answering"
    assert len(out) == 1
    assert "BRCA1" in out[0].description
