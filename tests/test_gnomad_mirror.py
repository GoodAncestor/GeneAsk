"""gnomAD AF mirror: key shape, resumability, and mirror-first reader behavior
(no network — all VCFs are pre-placed synthetic fixtures)."""
import gzip, sqlite3
from pathlib import Path
from geneask.annotators.gnomad_mirror import build_mirror
from geneask.annotators import gnomad_freq
from biocore.providers.base import Finding, Tier, Category
from biocore.variants.carried import carried_variants


_CHR21 = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t100\t.\tA\tG\t.\t.\tAC=10;AN=100;AF=0.1
21\t200\t.\tC\tT\t.\t.\tAC=1;AN=1000;AF=0.001
21\t300\t.\tG\tA,C\t.\t.\tAC=5,3;AN=200;AF=0.025,0.015
"""

_CHR22 = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
22\t400\t.\tT\tA\t.\t.\tAC=2;AN=200;AF=0.01
"""


def _place(wd, chrom, text):
    p = Path(wd) / f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
    with gzip.open(p, "wt") as fh:
        fh.write(text)
    return p


# The real gnomAD v4 files say "chr21", not "21" — verified against the published
# chrY sites file. A fixture without the prefix cannot tell a working strip from a
# missing one, so one case carries the real shape.
_CHR21_PREFIXED = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr21\t100\t.\tA\tG\t.\t.\tAC=10;AN=100;AF=0.1
"""


def test_the_chr_prefix_upstream_actually_uses_is_stripped(tmp_path):
    _place(tmp_path, "21", _CHR21_PREFIXED)
    db = str(tmp_path / "prefixed.db")
    build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21"])
    import sqlite3
    keys = [r[0] for r in sqlite3.connect(db).execute("select variant_id from af")]
    assert keys == ["21-100-A-G"], keys


def test_keys_match_carried_variants_shape(tmp_path):
    """The mirror's key format must be exactly what carried_variants() emits —
    'chrom-pos-ref-alt', no 'chr' prefix — or every mirror lookup misses."""
    _place(tmp_path, "21", _CHR21)
    db = str(tmp_path / "gnomad.db")
    s = build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21"])
    assert s["variants_added"] == 4     # 3 sites, one biallelic split into 2

    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr21\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"
    )
    carried = carried_variants(str(vcf))
    assert carried[0]["variant_id"] == "21-100-A-G"

    con = sqlite3.connect(db)
    row = con.execute("SELECT af FROM af WHERE variant_id=?", (carried[0]["variant_id"],)).fetchone()
    con.close()
    assert row is not None and abs(row[0] - 0.1) < 1e-9


def test_multiallelic_split(tmp_path):
    _place(tmp_path, "21", _CHR21)
    db = str(tmp_path / "gnomad.db")
    build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21"])
    con = sqlite3.connect(db)
    assert con.execute("SELECT af FROM af WHERE variant_id='21-300-G-A'").fetchone()[0] == 0.025
    assert con.execute("SELECT af FROM af WHERE variant_id='21-300-G-C'").fetchone()[0] == 0.015
    con.close()


def test_resumable_across_runs(tmp_path):
    """A chromosome already marked done must not be re-downloaded or re-parsed —
    resuming a killed build must not restart it from zero."""
    _place(tmp_path, "21", _CHR21)
    _place(tmp_path, "22", _CHR22)
    db = str(tmp_path / "gnomad.db")

    s1 = build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21"])
    assert s1["chroms_built"] == ["21"] and s1["chroms_skipped"] == []

    # The build deletes each source file once ingested — chr21 alone is 7.76 GB in
    # the real release, so keeping them would leave hundreds of GB beside a mirror
    # of a few dozen. Its absence here is therefore the assertion: chr21 cannot be
    # re-read on the next run, so the run below must skip it from `progress` alone.
    assert not (tmp_path / "gnomad.genomes.v4.1.sites.chr21.vcf.bgz").exists()

    s2 = build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21", "22"])
    assert s2["chroms_skipped"] == ["21"]
    assert s2["chroms_built"] == ["22"]
    assert s2["variants_added"] == 1     # only chr22's one site

    con = sqlite3.connect(db)
    n = con.execute("SELECT count(*) FROM af").fetchone()[0]
    con.close()
    assert n == 5   # chr21's 4 + chr22's 1, both present after the resumed run


def test_bounded_build_caps_variants_and_leaves_chrom_unfinished(tmp_path):
    _place(tmp_path, "21", _CHR21)
    db = str(tmp_path / "gnomad.db")
    s = build_mirror(db_path=db, workdir=str(tmp_path), chroms=["21"], max_variants=2)
    assert s["variants_added"] == 2
    assert s["chroms_built"] == []     # cut short: not recorded as complete

    con = sqlite3.connect(db)
    done = con.execute("SELECT chrom FROM progress WHERE done=1").fetchall()
    con.close()
    assert done == []   # a resumed run will redo chr21 from scratch, not skip it


def test_mirror_hit_needs_no_network(tmp_path, monkeypatch):
    """Once a mirror exists, allele_frequency() must never touch the live API —
    that is the entire point of building one."""
    def _boom(*a, **k):
        raise AssertionError("live gnomAD API called despite a mirror being present")
    monkeypatch.setattr("urllib.request.urlopen", _boom)

    _place(tmp_path, "21", _CHR21)
    mirror = str(tmp_path / "gnomad.db")
    build_mirror(db_path=mirror, workdir=str(tmp_path), chroms=["21"])

    assert gnomad_freq.allele_frequency("21-100-A-G", mirror_db=mirror) == 0.1
    # a variant absent from the mirror is a real miss, still no live call
    assert gnomad_freq.allele_frequency("21-999-A-G", mirror_db=mirror) is None


def test_annotate_findings_mirror_first_bulk_no_budget(tmp_path, monkeypatch):
    def _boom(*a, **k):
        raise AssertionError("live gnomAD API called despite a mirror being present")
    monkeypatch.setattr("urllib.request.urlopen", _boom)

    _place(tmp_path, "21", _CHR21)
    mirror = str(tmp_path / "gnomad.db")
    build_mirror(db_path=mirror, workdir=str(tmp_path), chroms=["21"])

    fs = [Finding("21-100-A-G", "clinvar", "gene X: pathogenic", Tier.MODERATE, [Category.CLINICAL], detail={}),
          Finding("21-999-A-G", "clinvar", "unmirrored", Tier.MODERATE, [Category.CLINICAL], detail={}),
          Finding("cg999", "ewas", "methylation", Tier.ROBUST, [Category.AGING], detail={})]
    n = gnomad_freq.annotate_findings(fs, mirror_db=mirror)
    assert n == 1
    assert fs[0].detail["gnomad_af"] == 0.1
    assert "10.0% of people" in fs[0].description
    assert "gnomad_af" not in fs[1].detail
    assert "gnomad_af" not in fs[2].detail


def test_live_fallback_unchanged_without_a_mirror(tmp_path, monkeypatch):
    """No mirror on disk: the original cached/live path must still run exactly as
    before — nothing here is allowed to regress a box that hasn't built one yet."""
    calls = []
    class _R:
        def __init__(self, body): self._b = body
        def read(self): return self._b
        def __enter__(self): return self
        def __exit__(self, *a): return False
    def fake(req, timeout=None):
        calls.append(req)
        return _R(b'{"data":{"variant":{"genome":{"af":0.42},"exome":null}}}')
    monkeypatch.setattr("urllib.request.urlopen", fake)

    cache = str(tmp_path / "cache.db")
    missing_mirror = str(tmp_path / "no_such_mirror.db")
    af = gnomad_freq.allele_frequency("1-1-A-G", cache_db=cache, mirror_db=missing_mirror)
    assert af == 0.42
    assert len(calls) == 1
