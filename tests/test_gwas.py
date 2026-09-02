"""GWAS Catalog annotator: mirror lookup + tiering + risk-allele note (mini fixture)."""
import sqlite3
from pathlib import Path
from geneask.annotators.gwas_catalog import mirror_lookup, findings_for, _tier
from biocore.providers.base import Tier


def _mini(tmp_path):
    db = tmp_path / "gwas.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
    con.execute("INSERT INTO meta VALUES ('schema_version', '2')")
    con.execute("""CREATE TABLE assoc(rsid TEXT, chrom TEXT, pos INTEGER,
        risk_allele TEXT, trait TEXT, mapped_trait TEXT, efo_uri TEXT,
        pvalue REAL, raf TEXT, pmid TEXT, effect_type TEXT, effect REAL,
        ci_text TEXT, accession TEXT, initial_n TEXT, replication_n TEXT,
        mapped_gene TEXT, pvalue_text TEXT)""")
    con.executemany("INSERT INTO assoc VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
        ("rs704", "2", 100, "G", "LDL cholesterol", "LDL cholesterol",
         "http://www.ebi.ac.uk/efo/EFO_0004611", 4e-8, "0.4", "12345678",
         None, None, None, None, None, None, None, None),
        ("rs704", "2", 100, "G", "heart disease", "coronary artery disease",
         "http://www.ebi.ac.uk/efo/EFO_0000378", 2e-3, "0.4", "22222222",
         None, None, None, None, None, None, None, None)])
    con.execute("CREATE INDEX idx_rsid ON assoc(rsid)")
    con.commit(); con.close()
    return str(db)


def test_lookup_and_findings(tmp_path):
    db = _mini(tmp_path)
    assert len(mirror_lookup("rs704", db)) == 2
    fs = findings_for("rs704", carried_alleles={"G", "A"}, db_path=db)
    assert len(fs) == 2
    # genome-wide significant p=4e-8 -> robust; suggestive p=2e-3 -> speculative
    tiers = sorted(f.tier.value for f in fs)
    assert "robust" in tiers and "speculative" in tiers
    # risk-allele carriage surfaced
    assert any("you carry it: yes" in f.description for f in fs)
    assert all(f.detail["modality"] == "genome" for f in fs)


def test_tier_thresholds():
    assert _tier(1e-9) == Tier.ROBUST         # genome-wide significant
    assert _tier(1e-6) == Tier.MODERATE       # suggestive
    assert _tier(0.01) == Tier.SPECULATIVE
    assert _tier(None) == Tier.SPECULATIVE


def test_absent_rsid(tmp_path):
    assert findings_for("rs999", db_path=_mini(tmp_path)) == []
