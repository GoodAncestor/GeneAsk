"""AlphaMissense annotator: mirror lookup + tiering + finding annotation (fixture)."""
import gzip, sqlite3
from geneask.annotators.alphamissense import build_mirror, lookup, findings_for, annotate_findings, _CLASS_TIER
from biocore.providers.base import Finding, Tier, Category

_TSV = ("#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n"
        "chr1\t100\tG\tT\thg38\tQ1\tENST1\tR567M\t0.98\tlikely_pathogenic\n"
        "chr1\t200\tA\tC\thg38\tQ2\tENST2\tV2L\t0.10\tlikely_benign\n")


def _mini(tmp_path):
    with gzip.open(tmp_path / "AlphaMissense_hg38.tsv.gz", "wt") as fh:
        fh.write(_TSV)
    db = str(tmp_path / "am.db")
    build_mirror(db_path=db, workdir=str(tmp_path))
    return db


def test_build_and_lookup(tmp_path):
    db = _mini(tmp_path)
    r = lookup("1-100-G-T", db)
    assert r["am_class"] == "likely_pathogenic" and r["pathogenicity"] == 0.98
    assert r["protein_variant"] == "R567M"


def test_findings_for_tier(tmp_path):
    db = _mini(tmp_path)
    fs = findings_for("1-100-G-T", db_path=db)
    assert fs and fs[0].tier == Tier.MODERATE      # likely_pathogenic -> moderate
    assert "likely pathogenic" in fs[0].description


def test_annotate_in_place(tmp_path):
    db = _mini(tmp_path)
    fs = [Finding("1-100-G-T", "clinvar", "gene X uncertain", Tier.SPECULATIVE, [Category.CLINICAL], detail={}),
          Finding("cg42", "ewas", "methylation", Tier.ROBUST, [Category.AGING], detail={})]
    n = annotate_findings(fs, db_path=db)
    assert n == 1
    assert fs[0].detail["alphamissense"]["class"] == "likely_pathogenic"
    assert "alphamissense" not in fs[1].detail        # CpG untouched


def test_no_mirror_noop(tmp_path):
    assert lookup("1-100-G-T", str(tmp_path / "nope.db")) is None
    assert annotate_findings([], db_path=str(tmp_path / "nope.db")) == 0


def test_mirror_status_and_missing_is_not_benign(tmp_path):
    from geneask.annotators.alphamissense import mirror_status
    missing = str(tmp_path / "missing.db")
    assert mirror_status(missing)["status"] == "unavailable"
    empty = tmp_path / "empty.db"
    empty.touch()
    assert mirror_status(str(empty))["status"] == "invalid"
    db = _mini(tmp_path)
    status = {}
    fs = [Finding("1-999-A-G", "novel_variant", "d", Tier.SPECULATIVE, [Category.CLINICAL], detail={})]
    assert annotate_findings(fs, db, status=status) == 0
    assert status["not_found"] == 1 and status["scored"] == 0
    assert "alphamissense" not in fs[0].detail
    status = {}
    annotate_findings(fs, missing, status=status)
    assert status["skipped"] == 1 and status["not_found"] == 0
