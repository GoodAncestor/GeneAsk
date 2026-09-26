import pytest
from geneask.annotators import atlas_avi


def _fixture(tmp_path, header="#CHROM\tPOS\tREF\tALT\tAVI_SCORE\tPHRED"):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "avi.tsv"
    path.write_text(header + "\nchr1\t100\tA\tG\t0.75\t12\nchr1\t200\tC\tT\t0.12\t2\n")
    bgz = str(path) + ".gz"
    pysam.tabix_compress(str(path), bgz, force=True)
    pysam.tabix_index(bgz, seq_col=0, start_col=1, end_col=1, zerobased=False, force=True)
    return bgz


def test_exact_local_alleles_bulk_and_phred(tmp_path):
    path = _fixture(tmp_path)
    assert atlas_avi.local_status(path)["status"] == "ready"
    status = {}
    rows = atlas_avi.lookup_many(["1-100-A-G", "chr1-200-C-T", "1-100-A-C", "2-100-A-G"], path, status=status)
    assert set(rows) == {"1-100-A-G", "1-200-C-T"}
    assert rows["1-100-A-G"]["avi_phred"] == 12
    assert rows["1-100-A-G"]["avi_score"] == .75
    assert status["not_found"] == 2
    assert atlas_avi.lookup("1-100-A-C", path)["status"] == "not_found"


def test_unknown_schema_is_explicit_and_configurable(tmp_path, monkeypatch):
    path = _fixture(tmp_path, "#a\tb\tc\td\te\tf")
    assert atlas_avi.local_status(path)["status"] == "unsupported_schema"
    monkeypatch.setenv("ALPHAGENOME_ATLAS_AVI_COLUMNS", "CHROM,POS,REF,ALT,AVI_SCORE,PHRED")
    assert atlas_avi.lookup("1-100-A-G", path)["avi_score"] == .75


def test_offline_query_uses_local_without_network(tmp_path, monkeypatch):
    from geneask.annotators import alphagenome_atlas as atlas
    path = _fixture(tmp_path)
    monkeypatch.setenv("ALPHAGENOME_ATLAS_AVI_FILE", path)
    monkeypatch.setenv("ALPHAGENOME_ATLAS_CACHE_DB", str(tmp_path / "cache.db"))
    monkeypatch.setattr(atlas, "_run_remote", lambda *a: pytest.fail("offline network attempt"))
    result = atlas.query_variant("1-100-A-G", offline=True)
    assert result["avi_score"] == .75
    assert result["missing_scorers"] == ["AVI_SCORE_FEATURE_IMPORTANCE"]
    assert result["remote_status"] == "offline"


def test_official_raw_score_header_alias(tmp_path):
    path = _fixture(tmp_path, "#chrom\tpos\tref\talt\traw_score\tphred")
    assert atlas_avi.lookup("1-100-A-G", path)["avi_score"] == .75
