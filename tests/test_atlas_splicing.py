import pytest
from biocore.providers.base import Finding, Tier, Category
from geneask.annotators import atlas_splicing, atlas_avi

REAL_HEADER = "#CHROM\tPOS\tREF\tALT\talphagenome_splicing"


def _fixture(tmp_path, header=REAL_HEADER, name="splice.tsv"):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / name
    path.write_text(header + "\nchr1\t65409\tA\tC\t0.003052\nchr1\t65409\tA\tG\t0.8\nchr1\t70000\tC\tT\tnan\n")
    bgz = str(path) + ".gz"
    pysam.tabix_compress(str(path), bgz, force=True)
    pysam.tabix_index(bgz, seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)
    return bgz


def test_exact_alleles_with_verified_header(tmp_path, monkeypatch):
    monkeypatch.delenv("DNAREPORT_OUTPUT_MODE", raising=False)
    path = _fixture(tmp_path)
    assert atlas_splicing.local_status(path)["status"] == "ready"
    status = {}
    rows = atlas_splicing.lookup_many(["1-65409-A-G", "chr1-65409-A-C", "1-65409-A-T", "1-70000-C-T"], path, status=status)
    assert set(rows) == {"1-65409-A-G", "1-65409-A-C"}
    assert rows["1-65409-A-G"]["splicing_score"] == .8
    assert rows["1-65409-A-G"]["license"]["commercial_allowed"] is False
    assert status["scored"] == 2 and status["not_found"] == 2


def test_other_headers_fail_closed(tmp_path, monkeypatch):
    monkeypatch.delenv("DNAREPORT_OUTPUT_MODE", raising=False)
    path = _fixture(tmp_path, "#CHROM\tPOS\tREF\tALT\tscore")
    assert atlas_splicing.local_status(path)["status"] == "unsupported_schema"
    assert atlas_splicing.lookup_many(["1-65409-A-G"], path) == {}


def test_commercial_mode_blocks_before_reading(tmp_path, monkeypatch):
    monkeypatch.setenv("DNAREPORT_OUTPUT_MODE", "commercial")
    monkeypatch.setattr(atlas_splicing, "_file", lambda *a: pytest.fail("file opened in commercial mode"))
    status = {}
    assert atlas_splicing.lookup_many(["1-65409-A-G"], status=status) == {}
    assert status["status"] == "license_blocked"


def test_annotate_findings_attaches_record(tmp_path, monkeypatch):
    monkeypatch.delenv("DNAREPORT_OUTPUT_MODE", raising=False)
    monkeypatch.setenv("ALPHAGENOME_ATLAS_SPLICING_FILE", _fixture(tmp_path))
    f = Finding("1-65409-A-G", "variant_lookup", "x", Tier.SPECULATIVE, [Category.CLINICAL], detail={})
    assert atlas_splicing.annotate_findings([f]) == 1
    assert f.detail["alphagenome_atlas_splicing"]["splicing_score"] == .8


def test_avi_reader_accepts_official_avi_header(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "avi.tsv"
    path.write_text("#CHROM\tPOS\tREF\tALT\traw_score\tPHRED\nchr1\t10001\tT\tA\t-0.03868\t1.06466\n")
    bgz = str(path) + ".gz"
    pysam.tabix_compress(str(path), bgz, force=True)
    pysam.tabix_index(bgz, seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)
    row = atlas_avi.lookup("1-10001-T-A", bgz)
    assert row["avi_score"] == -0.03868 and row["avi_phred"] == 1.06466
