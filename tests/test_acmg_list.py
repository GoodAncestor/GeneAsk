"""The ACMG SF v3.2 list travels as data and answers by gene."""
from geneask.interpret.lists import acmg_sf, ACMG_SF_VERSION, all_acmg_sf


def test_list_is_the_81_genes_of_sf_v3_2():
    assert ACMG_SF_VERSION == "SF v3.2"
    assert len(all_acmg_sf()) == 81


def test_brca2_and_hfe_entries():
    b = acmg_sf("BRCA2")
    assert b["condition"] == "Hereditary breast and ovarian cancer syndrome"
    assert b["inheritance"] == "dominant"
    h = acmg_sf("hfe")
    assert h["inheritance"] == "recessive"
    assert acmg_sf("MTHFR") is None


def test_ttn_carries_its_truncating_only_note():
    assert "truncating" in acmg_sf("TTN")["note"]
