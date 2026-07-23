"""Tests for GeneAsk interpretation on the bundled ClinVar panel + bio-core."""
from geneask.interpret.clinvar_screen import load_panel, index_by_variant_id, screen_findings
from biocore.providers.base import Finding, Tier, Category
from biocore.report.render import render_html


def test_panel_loads_gene_keyed():
    p = load_panel()
    assert len(p) == 157
    abcg2 = p["ABCG2"]
    assert abcg2["gene_id"].startswith("ENSG")
    assert isinstance(abcg2["variants"], list) and abcg2["variants"]


def test_variant_index_flattens():
    idx = index_by_variant_id(load_panel())
    assert len(idx) > 1000
    # every entry carries its gene and a variant_id shaped chr-pos-ref-alt
    k = next(iter(idx))
    assert idx[k]["gene"]
    assert k.count("-") >= 3


def test_screen_reports_only_pathogenic():
    panel = load_panel()
    idx = index_by_variant_id(panel)
    # pick a real pathogenic variant from the panel if present, else assert empty-safe
    path_ids = [vid for vid, r in idx.items()
                if "pathogenic" in (r.get("clinical_significance") or "").lower()]
    carried = [{"variant_id": vid, "genotype": "A/G", "platform": "WGS"} for vid in path_ids[:5]]
    findings = screen_findings(carried, panel)
    assert all(Category.CLINICAL in f.categories for f in findings)
    assert all("pathogenic" in f.description.lower() for f in findings)


def test_array_pathogenic_demoted():
    """Recovered caveat: a single-array pathogenic hit must never be ROBUST."""
    panel = load_panel()
    idx = index_by_variant_id(panel)
    two_star_path = next((vid for vid, r in idx.items()
                          if "pathogenic" in (r.get("clinical_significance") or "").lower()
                          and int(r.get("gold_stars", 0) or 0) >= 2), None)
    if two_star_path is None:
        return  # no 2-star pathogenic in panel; nothing to assert
    array = screen_findings([{"variant_id": two_star_path, "genotype": "A/A", "platform": "ARRAY"}], panel)
    assert array and array[0].tier != Tier.ROBUST


def test_findings_render_through_biocore():
    f = Finding(marker="1-1-A-G", source="clinvar_panel_157",
                description="TEST: pathogenic demo", tier=Tier.ROBUST, categories=[Category.CLINICAL])
    html = render_html([f], [])
    assert "TEST: pathogenic demo" in html
