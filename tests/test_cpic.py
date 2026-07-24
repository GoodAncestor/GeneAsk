"""CPIC PGx annotator: mirror lookup + tiering (mini fixture, no network)."""
import sqlite3
from geneask.annotators.cpic_pgx import recommendations_for_gene, _LEVEL_TIER
from biocore.providers.base import Tier


def _mini(tmp_path):
    db = tmp_path / "pgx.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE drug(drugid TEXT PRIMARY KEY, name TEXT)")
    con.execute("INSERT INTO drug VALUES ('RxNorm:1','clopidogrel')")
    con.execute("CREATE TABLE pair(genesymbol TEXT, drugid TEXT, cpiclevel TEXT, guidelineid INTEGER)")
    con.execute("INSERT INTO pair VALUES ('CYP2C19','RxNorm:1','A',100)")
    con.execute("""CREATE TABLE recommendation(guidelineid INTEGER, drugid TEXT, phenotypes TEXT,
                   drugrecommendation TEXT, classification TEXT, implications TEXT, lookupkey TEXT)""")
    con.execute("INSERT INTO recommendation VALUES (100,'RxNorm:1','{\"CYP2C19\":\"Poor Metabolizer\"}',"
                "'Use alternative antiplatelet','Strong','{}','{}')")
    con.execute("CREATE TABLE allele(genesymbol TEXT, name TEXT, functionalstatus TEXT)")
    con.commit(); con.close()
    return str(db)


def test_recommendations_for_gene(tmp_path):
    db = _mini(tmp_path)
    fs = recommendations_for_gene("CYP2C19", db_path=db)
    assert len(fs) == 1
    f = fs[0]
    assert f.detail["drug"] == "clopidogrel"
    assert f.detail["cpic_level"] == "A"
    assert f.tier == Tier.ROBUST          # level A -> robust
    assert f.detail["topic"] == "pharmacogenomic"
    assert "clopidogrel" in f.description


def test_phenotype_filter(tmp_path):
    db = _mini(tmp_path)
    assert len(recommendations_for_gene("CYP2C19", phenotype="Poor Metabolizer", db_path=db)) == 1
    assert len(recommendations_for_gene("CYP2C19", phenotype="Ultrarapid", db_path=db)) == 0


def test_absent_gene_and_no_mirror(tmp_path):
    assert recommendations_for_gene("NADK", db_path=_mini(tmp_path)) == []
    assert recommendations_for_gene("CYP2C19", db_path=str(tmp_path/"nope.db")) == []
