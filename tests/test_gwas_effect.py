from geneask.annotators.gwas_catalog import _effect, build_mirror, findings_for


def test_effect_parsing_follows_catalog_rule():
    assert _effect("1.32", "[1.21-1.44]") == ("or", 1.32)
    assert _effect("0.045", "[0.03-0.06] unit increase") == ("beta", 0.045)
    assert _effect("", "") == (None, None)
    assert _effect("NR", "[NR]") == (None, None)


def _tsv(tmp_path, cols):
    header = [f"c{i}" for i in range(38)]
    row = [""] * 38
    for index, value in cols.items():
        row[index] = value
    path = tmp_path / "gwas_assoc.tsv"
    path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")


def test_mirror_keeps_effect_and_finding_carries_it(tmp_path):
    _tsv(tmp_path, {
        1: "12345", 7: "Type 2 diabetes", 8: "10,000 European",
        9: "5,000 European", 11: "10", 12: "114758349", 14: "TCF7L2",
        20: "rs7903146-T", 21: "rs7903146", 26: "0.3", 27: "1e-40",
        29: "1E-40", 30: "1.37", 31: "[1.31-1.43]",
        34: "type II diabetes mellitus",
        35: "http://www.ebi.ac.uk/efo/EFO_0001360", 36: "GCST000001",
    })
    db = str(tmp_path / "g.db")
    build_mirror(db_path=db, workdir=str(tmp_path))
    finding, = findings_for("rs7903146", carried_alleles={"T"}, db_path=db)
    detail = finding.detail
    assert detail["effect_type"] == "or" and detail["effect"] == 1.37
    assert detail["ci_text"] == "[1.31-1.43]"
    assert detail["accession"] == "GCST000001"
    assert detail["initial_n"] == "10,000 European"
    assert detail["mapped_gene"] == "TCF7L2"
    assert detail["chrom"] == "10" and detail["pos"] == 114758349
    assert detail["risk_allele_carried"] is True
    finding2, = findings_for("rs7903146", carried_alleles={"C"}, db_path=db)
    assert finding2.detail["risk_allele_carried"] is False
    finding3, = findings_for("rs7903146", carried_alleles=None, db_path=db)
    assert finding3.detail["risk_allele_carried"] is None
