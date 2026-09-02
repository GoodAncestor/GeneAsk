"""GWAS and CPIC findings are explained in group terms, never as personal risk."""
from biocore.providers.base import Finding, Tier, Category
from geneask.interpret.meaning import interpret


def _gwas(**over):
    d = {"p": 1e-40, "trait": "Type 2 diabetes", "efo": "http://www.ebi.ac.uk/efo/EFO_0001360",
         "risk_allele": "T", "raf": "0.3", "risk_allele_carried": True,
         "effect_type": "or", "effect": 1.37, "ci_text": "[1.31-1.43]", "accession": "GCST000001",
         "initial_n": "10,000 European ancestry individuals", "replication_n": "5,000 European ancestry individuals",
         "mapped_gene": "TCF7L2", "topic": "other", "modality": "genome"}
    d.update(over)
    return Finding(marker="rs7903146", source="gwas_catalog", description="x", tier=Tier.ROBUST,
                   categories=[Category.TRAIT], detail=d,
                   link="https://www.ebi.ac.uk/gwas/variants/rs7903146", pmids=["12345"])


def test_gwas_or_is_said_in_words_and_never_as_personal_risk():
    f = _gwas(); interpret([f]); ip = f.interpretation
    assert "you carry the allele" in ip.found.lower()
    assert "1.37 times" in ip.can_mean and "odds" in ip.can_mean
    assert "not a personal risk" in ip.can_mean.lower()
    assert "10,000" in ip.how_sure and "Replication" in ip.how_sure
    assert ip.next_step == ""
    assert [c.kind for c in f.evidence_chain][:2] == ["variant", "trait"]
    assert any(c.id == "GCST000001" for c in f.evidence_chain)


def test_gwas_beta_uncarried_and_missing_effect():
    f = _gwas(effect_type="beta", effect=0.045, ci_text="[0.03-0.06] unit increase"); interpret([f])
    assert "0.045" in f.interpretation.can_mean and "unit increase" in f.interpretation.can_mean
    g = _gwas(risk_allele_carried=False); interpret([g])
    assert "you do not carry" in g.interpretation.found.lower()
    h = _gwas(effect_type=None, effect=None, replication_n=""); interpret([h])
    assert "effect size is not recorded" in h.interpretation.can_mean.lower()
    assert "no replication sample" in h.interpretation.how_sure


def test_cpic_names_the_drug_and_the_missing_diplotype():
    f = Finding(marker="CYP2C19", source="cpic", description="x", tier=Tier.ROBUST,
                categories=[Category.CLINICAL],
                detail={"gene": "CYP2C19", "drug": "clopidogrel", "cpic_level": "A",
                        "topic": "pharmacogenomic", "modality": "genome",
                        "phenotypes": ["Poor metabolizer"], "classification": "Strong"})
    interpret([f]); ip = f.interpretation
    assert "clopidogrel" in ip.found and "CYP2C19" in ip.found
    assert "did not determine" in ip.how_sure.lower()
    assert "pharmacist" in ip.next_step.lower()
    assert f.evidence_chain[0].kind == "gene"


def test_cpic_names_the_called_diplotype_and_specific_recommendation():
    f = Finding(
        marker="CYP2C19",
        source="cpic",
        description="x",
        tier=Tier.ROBUST,
        categories=[Category.CLINICAL],
        detail={
            "gene": "CYP2C19",
            "drug": "clopidogrel",
            "cpic_level": "A",
            "classification": "Strong",
            "recommendation": "Consider another antiplatelet",
            "diplotype": "*1/*2",
            "phenotype": "Intermediate Metabolizer",
            "activity_score": 1.0,
            "diplotype_source": "PharmCAT 3.4.0",
        },
    )

    interpret([f])
    interpretation = f.interpretation

    assert interpretation.found == (
        "Your CYP2C19 type is *1/*2 (Intermediate Metabolizer), called by PharmCAT."
    )
    assert "PharmCAT 3.4.0" in interpretation.how_sure
    assert "level" not in interpretation.how_sure          # level and score stay in detail
    assert "prescriber" in interpretation.next_step
    assert "clopidogrel" in interpretation.can_mean
