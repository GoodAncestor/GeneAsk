"""screen_findings() puts record and call data into detail fields."""
from geneask.interpret.clinvar_screen import screen_findings


_PANEL = {"BRCA2": {"variants": [{
    "variant_id": "13-32316419-CAG-C",
    "clinical_significance": "Likely pathogenic",
    "review_status": "criteria provided, multiple submitters, no conflicts",
    "gold_stars": 2,
    "clinvar_variation_id": "51063",
    "conditions": ["Hereditary breast ovarian cancer syndrome"],
    "condition_ids": ["MedGen:C0677776"],
    "molecular_consequence": "splice_acceptor_variant",
    "origin": ["germline"],
    "allele_id": "5001",
}]}}

_CARRIED = [{
    "variant_id": "13-32316419-CAG-C", "genotype": "C/CAG", "platform": "WGS",
    "zygosity": "het", "filter": "PASS", "qual": 812.5, "gq": 99, "dp": 41,
}]


def test_detail_carries_condition_and_call_fields():
    finding, = screen_findings(_CARRIED, _PANEL)
    detail = finding.detail
    assert detail["conditions"] == ["Hereditary breast ovarian cancer syndrome"]
    assert detail["condition_ids"] == ["MedGen:C0677776"]
    assert detail["molecular_consequence"] == "splice_acceptor_variant"
    assert detail["origin"] == ["germline"]
    assert detail["zygosity"] == "het" and detail["genotype"] == "C/CAG"
    assert detail["filter"] == "PASS" and detail["qual"] == 812.5
    assert detail["gq"] == 99 and detail["dp"] == 41
    assert finding.source == "clinvar_panel_157"


def test_missing_optional_fields_read_as_empty():
    panel = {"G": {"variants": [{
        "variant_id": "1-1-A-G", "clinical_significance": "Pathogenic",
        "gold_stars": 1,
    }]}}
    finding, = screen_findings(
        [{"variant_id": "1-1-A-G", "genotype": "A/G", "platform": "ARRAY"}],
        panel,
    )
    assert finding.detail["conditions"] == []
    assert finding.detail["zygosity"] is None
