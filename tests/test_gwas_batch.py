"""Batched GWAS mirror lookups.

The per-rsID entry point opened a fresh SQLite connection for every SNP on the
array. On a 650k-variant consumer export that was 31 of a 48 second report — two
thirds of the runtime — so the batch path is what keeps a large array inside the
CDN's time ceiling. These pin the behaviour that made it safe to switch.
"""
import sqlite3
import pytest
from geneask.annotators import gwas_catalog as gw


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    db = tmp_path / "gwas.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO meta VALUES ('schema_version', '2')")
    con.execute("""CREATE TABLE assoc(rsid TEXT, chrom TEXT, pos INTEGER,
        risk_allele TEXT, trait TEXT, mapped_trait TEXT, efo_uri TEXT,
        pvalue REAL, raf REAL, pmid TEXT, effect_type TEXT, effect REAL,
        ci_text TEXT, accession TEXT, initial_n TEXT, replication_n TEXT,
        mapped_gene TEXT, pvalue_text TEXT)""")
    empty = (None,) * 8
    con.executemany("INSERT INTO assoc VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
        ("rs1", "1", 100, "A", "trait one", "trait one", "http://e/EFO_1", 1e-9, 0.3, "111") + empty,
        ("rs1", "1", 100, "G", "trait two", "trait two", "http://e/EFO_2", 1e-6, 0.2, "222") + empty,
        ("rs2", "2", 200, "T", "trait three", "trait three", "http://e/EFO_3", 1e-12, 0.1, "333") + empty,
    ])
    con.execute("CREATE INDEX idx_rsid ON assoc(rsid)")
    con.commit()
    con.close()
    monkeypatch.setenv("GWAS_MIRROR_DB", str(db))
    gw._conns.by_path = {}          # the connection cache is per-thread and sticky
    return str(db)


def test_batch_returns_the_same_rows_the_single_lookup_does(mirror):
    """The batch path replaced the per-rsID one in the report; if they disagreed,
    switching would have silently changed what people are told."""
    batch = gw.mirror_lookup_many(["rs1", "rs2"])
    assert batch["rs1"] == gw.mirror_lookup("rs1")
    assert batch["rs2"] == gw.mirror_lookup("rs2")


def test_missing_rsids_are_absent_rather_than_empty(mirror):
    got = gw.mirror_lookup_many(["rs1", "rs_nope"])
    assert "rs_nope" not in got
    assert len(got["rs1"]) == 2


def test_duplicate_rsids_are_queried_once(mirror):
    """A callset repeats rsIDs; querying each occurrence is wasted work."""
    got = gw.mirror_lookup_many(["rs1", "rs1", "rs1"])
    assert len(got["rs1"]) == 2          # not 6


def test_more_rsids_than_sqlites_parameter_limit(mirror):
    """SQLite rejects a statement with more than 999 parameters, so the batch has
    to chunk. A whole array is three orders of magnitude past that."""
    rsids = [f"rs{i}" for i in range(5000)] + ["rs1"]
    got = gw.mirror_lookup_many(rsids)
    assert len(got["rs1"]) == 2


def test_no_mirror_is_empty_not_an_error(mirror, monkeypatch, tmp_path):
    monkeypatch.setenv("GWAS_MIRROR_DB", str(tmp_path / "absent.db"))
    assert gw.mirror_lookup_many(["rs1"]) == {}


def test_findings_from_rows_matches_findings_for(mirror):
    rows = gw.mirror_lookup("rs1")
    a = gw.findings_from_rows("rs1", rows, carried_alleles={"A"})
    b = gw.findings_for("rs1", carried_alleles={"A"})
    assert [f.description for f in a] == [f.description for f in b]
    assert any("you carry it: yes" in f.description for f in a)
