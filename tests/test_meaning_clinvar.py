"""A ClinVar finding's four parts come from its fields and the copy tables."""
from biocore.providers.base import Finding, Tier, Category
from geneask.interpret.meaning import interpret, interpret_clinvar


def _brca2(**over):
    d = {"gene": "BRCA2", "topic": "cancer", "modality": "genome",
         "clinical_significance": "Likely pathogenic", "gold_stars": 2,
         "review_status": "criteria provided, multiple submitters, no conflicts",
         "clinvar_variation_id": "51063", "platform": "WGS", "zygosity": "het",
         "conditions": ["Hereditary breast ovarian cancer syndrome"],
         "condition_ids": ["MedGen:C0677776"], "molecular_consequence": "splice_acceptor_variant",
         "gnomad": {"af": 3.2e-6, "ac": 5, "an": 1560000, "version": "v4.1"}}
    d.update(over)
    return Finding(marker="13-32316419-CAG-C", source="clinvar_mirror", description="x",
                   tier=Tier.ROBUST, categories=[Category.CLINICAL], detail=d,
                   link="https://www.ncbi.nlm.nih.gov/clinvar/variation/51063/", pmids=[])


def test_four_parts_come_from_fields_and_tables():
    ip = interpret_clinvar(_brca2())
    assert ip.condition == "Hereditary breast and ovarian cancer syndrome"
    assert ip.zygosity == "het"
    assert "BRCA2 makes a protein" in ip.found
    assert "likely pathogenic" in ip.found.lower()
    assert "splice acceptor variant" in ip.found
    assert ip.can_mean.startswith("One altered copy was read.")
    assert "one altered copy is enough" in ip.can_mean.lower()
    assert "2 of 4 review stars" in ip.how_sure
    assert "5 of 1,560,000" in ip.how_sure
    assert "genetic counsel" in ip.next_step.lower()
    assert ip.copy_version == "1" and ip.reviewed_by == []


def test_array_call_asks_for_confirmation_and_recessive_het_is_a_carrier():
    ip = interpret_clinvar(_brca2(platform="ARRAY", gene="HFE",
                                  conditions=["Hemochromatosis type 1"], condition_ids=[], gnomad=None))
    assert "clinical test" in ip.next_step.lower()
    assert "carrier" in ip.can_mean.lower()
    assert "sampled chromosomes" not in ip.how_sure


def test_unknown_gene_and_condition_fall_back_without_inventing():
    ip = interpret_clinvar(_brca2(gene="ZZZ9", conditions=["Some_rare_condition"], condition_ids=[]))
    assert ip.found.startswith("A change in ZZZ9.")
    assert ip.condition == "Some rare condition"
    assert "how the condition is inherited" in ip.can_mean


def test_vus_and_conflicting_say_no_action():
    v = interpret_clinvar(_brca2(clinical_significance="Uncertain significance"))
    assert "No action follows" in v.next_step and "but not this change" in v.can_mean
    c = interpret_clinvar(_brca2(clinical_significance="Conflicting classifications of pathogenicity"))
    assert "disagree" in c.how_sure and "Others do not" in c.can_mean


def test_interpret_sets_chain_and_skips_unknown_sources():
    f = _brca2()
    other = Finding(marker="rs1", source="somethingelse", description="x", tier=Tier.UNKNOWN,
                    categories=[Category.TRAIT])
    assert interpret([f, other]) == 1
    kinds = [c.kind for c in f.evidence_chain]
    assert kinds[:3] == ["variant", "gene", "condition"]
    assert any(c.kind == "assertion" and "51063" in (c.url or "") for c in f.evidence_chain)
    assert any(c.kind == "paper" and "ACMG" in c.label for c in f.evidence_chain)
    assert other.interpretation is None
