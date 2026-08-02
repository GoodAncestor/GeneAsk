"""Annotated-VCF output: header provenance/escaping, Number=A multi-allelic
correctness, and round-trip safety for records GeneAsk can't annotate."""
import gzip
import re

from geneask.interpret.annotated_vcf import (
    annotate_vcf, _meta_field, _escape_meta, _safe_token, _append_info,
)

_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n'
)

# a tiny fake ClinVar panel, gene-keyed like the bundled one, with one entry
# whose clinical_significance has a comma (the classic real-world ClinVar
# string) to exercise the value-sanitizer, and one clean pathogenic hit.
_PANEL = {"BRCA1": {"release": "2099-01-01", "variants": [
    {"variant_id": "1-100-A-G", "gene": "BRCA1",
     "clinical_significance": "Pathogenic, low penetrance",
     "review_status": "criteria provided, single submitter", "gold_stars": 1},
]}, "TP53": {"release": "2099-01-01", "variants": [
    {"variant_id": "1-200-C-T", "gene": "TP53",
     "clinical_significance": "Benign", "review_status": "reviewed by expert panel",
     "gold_stars": 3},
]}}


def _write(tmp_path, body_lines, name="in.vcf"):
    p = tmp_path / name
    p.write_text(_HEADER + "\n".join(body_lines) + "\n")
    return str(p)


def _info_dict(info_field):
    return dict(kv.split("=", 1) for kv in info_field.split(";") if "=" in kv)


def _read_out(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return fh.read().splitlines()


# ---------------------------------------------------------------- header ---

def test_header_declares_every_info_field_with_correct_number_type(tmp_path):
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    lines = _read_out(out)
    info_lines = {re.match(r"##INFO=<ID=([^,]+),", l).group(1): l
                  for l in lines if l.startswith("##INFO=<ID=GA_")}
    expected = {"GA_CLNSIG": ("A", "String"), "GA_CLNGENE": ("A", "String"),
                "GA_CLNREVSTAT": ("A", "String"), "GA_CLNSTARS": ("A", "Integer"),
                "GA_AM_CLASS": ("A", "String"), "GA_AM_SCORE": ("A", "Float"),
                "GA_AM_PROT": ("A", "String")}
    assert set(info_lines) == set(expected)
    for fid, (num, typ) in expected.items():
        assert f"Number={num}" in info_lines[fid]
        assert f"Type={typ}" in info_lines[fid]
        assert 'Description="' in info_lines[fid]


def test_header_has_provenance_line_before_chrom(tmp_path):
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    lines = _read_out(out)
    prov = [l for l in lines if l.startswith("##geneaskAnnotateVCF=")]
    assert len(prov) == 1
    assert "Tool=geneask" in prov[0]
    assert "Version=" in prov[0]
    assert re.search(r"RunTimestamp=\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", prov[0])
    assert "clinvar:explicit_panel" in prov[0]
    chrom_idx = next(i for i, l in enumerate(lines) if l.startswith("#CHROM"))
    prov_idx = next(i for i, l in enumerate(lines) if l.startswith("##geneaskAnnotateVCF="))
    assert prov_idx < chrom_idx    # meta lines precede the column header, always


def test_original_meta_lines_preserved_and_not_duplicated(tmp_path):
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    lines = _read_out(out)
    assert '##fileformat=VCFv4.2' in lines
    assert '##FILTER=<ID=PASS,Description="All filters passed">' in lines
    assert sum(1 for l in lines if l.startswith("##INFO=<ID=GA_CLNSIG,")) == 1


# ------------------------------------------------------------ escaping ----

def test_meta_field_quotes_only_when_needed():
    assert _meta_field("plain") == "plain"
    assert _meta_field("has,comma") == '"has,comma"'
    assert _meta_field('has"quote') == '"has\\"quote"'
    assert _meta_field("has=equals") == '"has=equals"'
    assert _meta_field("has space") == '"has space"'


def test_escape_meta_backslash_and_quote():
    assert _escape_meta('a "quoted" \\ value') == 'a \\"quoted\\" \\\\ value'


def test_commas_equals_quotes_in_description_dont_corrupt_header(tmp_path):
    """A malicious/weird panel entry with comma+quote+equals in its
    clinical_significance must not break header parsing of the *value line*
    it lands in (the INFO/Description declarations are fixed strings we
    control; this proves the per-record VALUE sanitizer handles the same
    hazards the header-escaping requirement calls out)."""
    panel = {"WEIRD": {"release": "2099-01-01", "variants": [
        {"variant_id": "1-100-A-G", "gene": "WEIRD",
         "clinical_significance": 'odd"value=with,commas', "review_status": "x", "gold_stars": 1},
    ]}}
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=panel)
    lines = _read_out(out)
    rec = next(l for l in lines if l.startswith("1\t100\t"))
    info = _info_dict(rec.split("\t")[7])
    # sanitized: no raw comma/quote/equals survives inside the value slot
    val = info["GA_CLNSIG"]
    assert '"' not in val and "=" not in val and "," not in val
    # and the record line still splits into exactly 10 tab fields (untouched structure)
    assert len(rec.split("\t")) == 10


def test_safe_token_and_append_info():
    assert _safe_token(None) == "."
    assert _safe_token("") == "."
    assert _safe_token("Likely pathogenic, low penetrance") == "Likely_pathogenic_low_penetrance"
    assert _append_info(".", {"X": "1"}) == "X=1"
    assert _append_info("DP=10", {"X": "1", "Y": "2"}) == "DP=10;X=1;Y=2"


# --------------------------------------------------------- multi-allelic --

def test_multiallelic_number_a_alignment(tmp_path):
    """ALT=G,T where only G hits the panel: GA_CLNSIG must be 'Pathogenic...,.'."""
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG,T\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    rec = next(l for l in _read_out(out) if l.startswith("1\t100\t"))
    info = _info_dict(rec.split("\t")[7])
    sig_vals = info["GA_CLNSIG"].split(",")
    gene_vals = info["GA_CLNGENE"].split(",")
    assert len(sig_vals) == 2 == len(gene_vals)     # one slot per ALT allele
    assert sig_vals[0].startswith("Pathogenic") and sig_vals[1] == "."
    assert gene_vals[0] == "BRCA1" and gene_vals[1] == "."


def test_clean_hit_full_field_set(tmp_path):
    vcf = _write(tmp_path, ["1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1"])
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    rec = next(l for l in _read_out(out) if l.startswith("1\t200\t"))
    info = _info_dict(rec.split("\t")[7])
    assert info["GA_CLNSIG"] == "Benign"
    assert info["GA_CLNGENE"] == "TP53"
    assert info["GA_CLNSTARS"] == "3"
    assert stats["clinvar_hits"] == 1
    assert stats["annotated"] == 1
    assert stats["records"] == 1
    assert stats["alt_alleles"] == 1


# ------------------------------------------------------- round-trip -------

def test_unannotatable_records_pass_through_byte_identical_info(tmp_path):
    """A record whose (chrom,pos,ref,alt) matches nothing gets no GA_* keys
    at all — original INFO field untouched, not even a placeholder '.'"""
    vcf = _write(tmp_path, ["1\t999\trs9\tA\tG\t.\tPASS\tDP=30;AF=0.5\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    rec = next(l for l in _read_out(out) if l.startswith("1\t999\t"))
    assert rec.split("\t")[7] == "DP=30;AF=0.5"     # byte-identical, no GA_* added


def test_symbolic_and_breakend_alt_pass_through_unchanged(tmp_path):
    vcf = _write(tmp_path, [
        "1\t500\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1",
        "1\t600\t.\tG\tG]17:1584563]\t.\tPASS\t.\tGT\t0/1",
        "1\t700\t.\tA\t*\t.\tPASS\t.\tGT\t0/1",
        "1\t800\t.\tA\t.\t.\tPASS\t.\tGT\t0/0",
    ])
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    lines = [l for l in _read_out(out) if l.startswith("1\t")]
    assert len(lines) == 4
    for orig, got in zip(["1\t500\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1",
                          "1\t600\t.\tG\tG]17:1584563]\t.\tPASS\t.\tGT\t0/1",
                          "1\t700\t.\tA\t*\t.\tPASS\t.\tGT\t0/1",
                          "1\t800\t.\tA\t.\t.\tPASS\t.\tGT\t0/0"], lines):
        assert orig == got
    assert stats["annotated"] == 0


def test_records_never_dropped_or_reordered(tmp_path):
    body = [f"1\t{100+i}\trs{i}\tA\tG\t.\tPASS\t.\tGT\t0/1" for i in range(20)]
    vcf = _write(tmp_path, body)
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    out_recs = [l for l in _read_out(out) if l.startswith("1\t")]
    assert len(out_recs) == 20
    assert [l.split("\t")[1] for l in out_recs] == [str(100 + i) for i in range(20)]


def test_sample_and_format_columns_untouched(tmp_path):
    vcf = _write(tmp_path, ["1\t200\trs2\tC\tT\t99\tPASS\tDP=5\tGT:AD:DP\t1/1:0,12:12"])
    out = str(tmp_path / "out.vcf")
    annotate_vcf(vcf, out, clinvar_panel=_PANEL)
    rec = next(l for l in _read_out(out) if l.startswith("1\t200\t"))
    fields = rec.split("\t")
    assert fields[0:4] == ["1", "200", "rs2", "C"]
    assert fields[5:7] == ["99", "PASS"]      # QUAL, FILTER untouched
    assert fields[8] == "GT:AD:DP"
    assert fields[9] == "1/1:0,12:12"


def test_gz_input_and_output_roundtrip(tmp_path):
    p = tmp_path / "in.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_HEADER + "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\n")
    out = str(tmp_path / "out.vcf.gz")
    stats = annotate_vcf(str(p), out, clinvar_panel=_PANEL)
    assert stats["clinvar_hits"] == 1
    lines = _read_out(out)
    assert any(l.startswith("1\t200\t") for l in lines)
    assert any(l.startswith("##geneaskAnnotateVCF=") for l in lines)


# ------------------------------------------------- alphamissense layer ----

def test_alphamissense_layer_noop_when_no_mirror(tmp_path):
    vcf = _write(tmp_path, ["1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(vcf, out, clinvar_panel=_PANEL,
                         alphamissense_db=str(tmp_path / "nope.db"))
    assert stats["alphamissense_hits"] == 0
    lines = _read_out(out)
    prov = next(l for l in lines if l.startswith("##geneaskAnnotateVCF="))
    assert "alphamissense:not_available" in prov
    rec = next(l for l in lines if l.startswith("1\t200\t"))
    info = _info_dict(rec.split("\t")[7])
    assert info["GA_AM_CLASS"] == "."


def test_alphamissense_layer_used_when_mirror_built(tmp_path):
    from geneask.annotators.alphamissense import build_mirror
    tsv = ("#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\t"
           "am_pathogenicity\tam_class\n"
           "chr1\t200\tC\tT\thg38\tQ1\tENST1\tR567M\t0.91\tlikely_pathogenic\n")
    with gzip.open(tmp_path / "AlphaMissense_hg38.tsv.gz", "wt") as fh:
        fh.write(tsv)
    db = str(tmp_path / "am.db")
    build_mirror(db_path=db, workdir=str(tmp_path))

    vcf = _write(tmp_path, ["1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(vcf, out, clinvar_panel=_PANEL, alphamissense_db=db)
    assert stats["alphamissense_hits"] == 1
    rec = next(l for l in _read_out(out) if l.startswith("1\t200\t"))
    info = _info_dict(rec.split("\t")[7])
    assert info["GA_AM_CLASS"] == "likely_pathogenic"
    assert info["GA_AM_SCORE"] == "0.9100"
    assert info["GA_AM_PROT"] == "R567M"
    prov = next(l for l in _read_out(out) if l.startswith("##geneaskAnnotateVCF="))
    assert "alphamissense:mirror" in prov


# ------------------------------------------------ chrom-prefix handling ---

def test_chr_prefixed_input_still_matches_bare_panel(tmp_path):
    """UCSC-style 'chr1' input must still hit a panel keyed bare '1-...' —
    and the output CHROM column must still read 'chr1' verbatim (never rewritten)."""
    vcf = ('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
          'chr1\t200\trs2\tC\tT\t.\tPASS\t.\n')
    p = tmp_path / "in.vcf"
    p.write_text(vcf)
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(str(p), out, clinvar_panel=_PANEL)
    assert stats["clinvar_hits"] == 1
    rec = next(l for l in _read_out(out) if l.startswith("chr1\t"))
    assert rec.split("\t")[0] == "chr1"


def test_bundled_panel_default_used_when_no_override(tmp_path):
    """Without an explicit panel, the real bundled 157-gene panel loads (or
    the full mirror if one happens to be built) — smoke test the default path."""
    vcf = _write(tmp_path, ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"])
    out = str(tmp_path / "out.vcf")
    stats = annotate_vcf(vcf, out)
    assert stats["records"] == 1
    prov = next(l for l in _read_out(out) if l.startswith("##geneaskAnnotateVCF="))
    assert "clinvar:bundled_157gene_panel" in prov or "clinvar:full_mirror" in prov
