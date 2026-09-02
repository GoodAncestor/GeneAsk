"""A ClinVar finding must link to its own ClinVar record.

Both the bundled panel and the full mirror already carry
`clinvar_variation_id` — the screen simply never put it on the Finding, so
bio-core's renderer fell through to `f.link or SOURCES['clinvar'].url` and every
clinical finding in a genome report pointed at the ClinVar homepage. A reader
following a pathogenic call had no way back to the submitters, the review status
or the assertion they were being shown.
"""
from geneask.interpret.clinvar_screen import screen_findings

_PANEL = {
    "BRCA2": {"variants": [
        {"variant_id": "13-32340301-A-G", "clinvar_variation_id": "51062",
         "clinical_significance": "Pathogenic",
         "review_status": "reviewed by expert panel", "gold_stars": 3},
        # a record the panel carries without a variation id
        {"variant_id": "13-32340999-C-T", "clinical_significance": "Pathogenic",
         "review_status": "criteria provided, single submitter", "gold_stars": 1},
    ]},
}


def _screen(variant_id):
    carried = [{"variant_id": variant_id, "genotype": "A/G", "platform": "WGS"}]
    found = screen_findings(carried, panel=_PANEL)
    assert len(found) == 1
    return found[0]


def test_link_points_at_the_variant_record(monkeypatch, tmp_path):
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(tmp_path / "absent.db"))
    f = _screen("13-32340301-A-G")
    assert f.link == "https://www.ncbi.nlm.nih.gov/clinvar/variation/51062/"


def test_variation_id_is_carried_in_detail(monkeypatch, tmp_path):
    """So the JSON export and MCP consumers can resolve the record themselves
    rather than re-deriving it by scraping the link."""
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(tmp_path / "absent.db"))
    assert _screen("13-32340301-A-G").detail["clinvar_variation_id"] == "51062"


def test_mirror_backed_findings_link_to_the_record_too(tmp_path, monkeypatch):
    """The deployed workers screen against the full mirror, not the bundled
    panel, so the panel-only test above would not have caught a regression on
    the path every real report actually takes."""
    import sqlite3
    db = tmp_path / "clinvar.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO meta VALUES ('schema_version', '2')")
    con.execute("""CREATE TABLE variants(
        variant_id TEXT PRIMARY KEY, clinvar_variation_id TEXT, gene TEXT,
        clinical_significance TEXT, review_status TEXT, gold_stars INTEGER,
        conditions TEXT, condition_ids TEXT, molecular_consequence TEXT,
        origin TEXT, allele_id TEXT)""")
    con.execute("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                ("13-32340301-A-G", "51062", "BRCA2", "Pathogenic",
                 "reviewed by expert panel", 3, "[]", "[]", None, "[]", None))
    con.commit()
    con.close()
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(db))
    found = screen_findings([{"variant_id": "13-32340301-A-G",
                              "genotype": "A/G", "platform": "WGS"}])
    assert found[0].link == "https://www.ncbi.nlm.nih.gov/clinvar/variation/51062/"


def test_no_variation_id_falls_back_to_the_gene_search(monkeypatch, tmp_path):
    """Not the homepage: a gene-scoped ClinVar search is still somewhere the
    reader can look. It must not pretend to be the variant's own record."""
    monkeypatch.setenv("CLINVAR_MIRROR_DB", str(tmp_path / "absent.db"))
    f = _screen("13-32340999-C-T")
    assert f.link == "https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA2%5Bgene%5D"
    assert "/variation/" not in f.link
