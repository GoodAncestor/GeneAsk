"""Position-based genome-build detection (overrides lying headers)."""
from geneask.parsers.base import GenotypeRecord
from geneask.parsers.build_detect import detect_build, _anchors


def _recs(build_key):
    return [GenotypeRecord(a["rsid"], a["chrom"], a[build_key], "A", "G")
            for a in _anchors().values()]


def test_detect_38_and_37():
    assert detect_build(_recs("pos38")) == "38"
    assert detect_build(_recs("pos37")) == "37"


def test_undecidable_without_markers():
    assert detect_build([GenotypeRecord("rsNONE", "1", 12345, "A", "G")]) is None
