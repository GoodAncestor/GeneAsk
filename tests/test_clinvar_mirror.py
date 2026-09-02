"""Full ClinVar mirror: VCF parse -> SQLite, index shape, screen upgrade."""
import gzip, sqlite3, os
from pathlib import Path
from geneask.annotators.clinvar_mirror import (
    build_mirror, load_panel_from_mirror, lookup_from_mirror, mirror_status,
    _stars, _gene,
)


_VCF = """##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t930200\t12345\tG\tA\t.\t.\tALLELEID=5001;GENEINFO=SAMD11:148398;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Hereditary_breast_ovarian_cancer_syndrome|not_provided;CLNDISDB=MedGen:C0677776,OMIM:604370|MedGen:CN517202;MC=SO:0001574|splice_acceptor_variant;ORIGIN=1
1\t930245\t12346\tG\tA\t.\t.\tGENEINFO=SAMD11:148398;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter
2\t500\t99\tA\t.\t.\t.\tCLNSIG=Pathogenic
"""


def _write_vcf(tmp_path):
    p = tmp_path / "clinvar.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_VCF)
    return p


def test_build_and_index(tmp_path, monkeypatch):
    _write_vcf(tmp_path)   # pre-place so build_mirror skips download
    db = str(tmp_path / "cv.db")
    s = build_mirror(db_path=db, workdir=str(tmp_path))
    assert s["variants"] == 2          # the ALT='.' row is skipped
    idx = load_panel_from_mirror(db)
    assert "1-930200-G-A" in idx
    rec = idx["1-930200-G-A"]
    assert rec["gene"] == "SAMD11"
    assert rec["gold_stars"] == 3      # reviewed by expert panel
    assert "Pathogenic" in rec["clinical_significance"]


def test_v2_columns_are_kept(tmp_path):
    _write_vcf(tmp_path)
    db = str(tmp_path / "cv.db")
    build_mirror(db_path=db, workdir=str(tmp_path))
    rec = load_panel_from_mirror(db)["1-930200-G-A"]
    assert rec["conditions"] == ["Hereditary breast ovarian cancer syndrome"]
    assert rec["condition_ids"] == ["MedGen:C0677776", "OMIM:604370"]
    assert rec["molecular_consequence"] == "splice_acceptor_variant"
    assert rec["origin"] == ["germline"]
    assert rec["allele_id"] == "5001"
    rec2 = load_panel_from_mirror(db)["1-930245-G-A"]
    assert rec2["conditions"] == [] and rec2["molecular_consequence"] is None


def test_v1_database_is_refused(tmp_path):
    db = str(tmp_path / "old.db")
    con = sqlite3.connect(db)
    con.execute("""CREATE TABLE variants(variant_id TEXT PRIMARY KEY,
        clinvar_variation_id TEXT, gene TEXT, clinical_significance TEXT,
        review_status TEXT, gold_stars INTEGER)""")
    con.execute("INSERT INTO variants VALUES ('1-1-A-G','1','G','Pathogenic','x',1)")
    con.commit()
    con.close()
    assert lookup_from_mirror(["1-1-A-G"], db_path=db) is None
    assert load_panel_from_mirror(db) is None
    status = mirror_status(db)
    assert status.health.value == "unavailable"
    assert "schema" in (status.note or "").lower()


def test_stars_and_gene():
    assert _stars("reviewed by expert panel") == 3
    assert _stars("practice guideline") == 4
    assert _stars("no assertion criteria provided") == 0
    assert _gene("BRCA2:675|ZAR1L:100") == "BRCA2"


def test_screen_prefers_mirror(tmp_path, monkeypatch):
    # when CLINVAR_MIRROR_DB points at a built mirror, index_by_variant_id uses it
    _write_vcf(tmp_path)
    db = str(tmp_path / "cv.db")
    build_mirror(db_path=db, workdir=str(tmp_path))
    monkeypatch.setenv("CLINVAR_MIRROR_DB", db)
    from geneask.interpret.clinvar_screen import index_by_variant_id
    idx = index_by_variant_id({})     # empty bundled panel; mirror should win
    assert "1-930200-G-A" in idx and idx["1-930200-G-A"]["gene"] == "SAMD11"


def test_repeated_mondo_prefix_is_collapsed():
    from geneask.annotators.clinvar_mirror import _condition_ids, _norm_id
    assert _norm_id("MONDO:MONDO:0012933") == "MONDO:0012933"
    assert _norm_id("MedGen:C0677776") == "MedGen:C0677776"
    assert _condition_ids("MONDO:MONDO:0012933,MedGen:C2675520|MedGen:CN517202",
                          "Breast-ovarian_cancer|not_provided") == ["MONDO:0012933", "MedGen:C2675520"]
