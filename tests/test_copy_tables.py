"""The copy tables: lookup rules, fallbacks, and the house voice."""
import re
from pathlib import Path
from geneask.interpret.copy import (classify, platform_class, next_step, gene_function,
                                    condition_phrase, copy_meta, _TABLES)

BANNED = re.compile(r"honest|it's not |not a [A-Z]|the point is|worth (noting|being explicit)"
                    r"|in other words|crucially|importantly,|notably,|by design|you have|you will",
                    re.I)


def test_classify_precedence_matches_render_direction():
    assert classify("Conflicting classifications of pathogenicity; risk factor") == "conflicting"
    assert classify("Uncertain significance") == "vus"
    assert classify("Pathogenic; Affects") == "plp"
    assert classify("Likely pathogenic") == "plp"
    assert classify("Benign; drug response") == "drug_response"
    assert classify("risk factor") == "risk_factor"
    assert classify("Benign") == "benign"
    assert classify("association") == "other"
    assert classify("") == "other"


def test_platform_and_next_step():
    assert platform_class("WGS") == "wgs_wes" and platform_class("ARRAY") == "array"
    a = next_step("plp", "array", "het", stars=2)
    assert "clinical test" in a["next_step"].lower()
    assert "2 of 4" in a["how_sure"]
    w = next_step("plp", "wgs_wes", "het", stars=2)
    assert "genetic counsel" in w["next_step"].lower()
    assert "whole-genome" in w["how_sure"]
    b = next_step("benign", "wgs_wes", "hom", stars=1)
    assert "no risk" not in b["next_step"].lower()
    for cls in ("vus", "conflicting", "drug_response", "risk_factor"):
        assert next_step(cls, "array", None)["next_step"]
    assert next_step("other", "array", None)["next_step"] == ""
    assert "unknown number" in next_step("plp", "wgs_wes", None)["how_sure"]


def test_gene_function_and_condition():
    g = gene_function("BRCA2")
    assert g["url"].startswith("https://www.ncbi.nlm.nih.gov/gene/") and len(g["sentence"]) < 160
    assert gene_function("brca2") == g
    assert gene_function("NOTAGENE") is None
    c = condition_phrase(["MedGen:C0677776"], "Hereditary breast ovarian cancer syndrome", "BRCA2")
    assert c["name"] == "Hereditary breast and ovarian cancer syndrome"
    assert c["inheritance"] == "dominant"
    c2 = condition_phrase([], "Some_rare_condition", "ZZZ9")
    assert c2["name"] == "Some rare condition" and c2["inheritance"] is None


def test_every_table_sentence_passes_the_voice_check():
    for name, path in _TABLES.items():
        text = Path(path).read_text()
        hits = [m.group(0) for m in BANNED.finditer(text)]
        assert not hits, f"{name}: {hits}"
        for s in re.findall(r'"([^"]{40,})"', text):
            for sentence in re.split(r"(?<=[.!?])\s+", s):
                assert len(sentence.split()) <= 25, f"{name}: too long: {sentence}"


def test_gene_table_covers_the_panel_and_the_acmg_list():
    import gzip, json
    from geneask.interpret.lists import all_acmg_sf
    panel = json.load(gzip.open(Path(__file__).resolve().parents[1]
                                / "geneask/data/reference/clinvar_panel_157genes.json.gz", "rt"))
    for g in set(panel) | set(all_acmg_sf()):
        assert gene_function(g), g


def test_meta_says_who_reviewed():
    m = copy_meta("clinical_next_step")
    assert m["version"] and isinstance(m["reviewed_by"], list)
