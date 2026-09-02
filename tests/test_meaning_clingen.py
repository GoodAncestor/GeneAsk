from biocore.providers.base import Category, Finding, Tier

from geneask.annotators import clingen
from geneask.interpret.meaning import interpret


def test_clingen_inheritance_and_evidence_reach_the_meaning(monkeypatch):
    monkeypatch.setattr(
        clingen,
        "inheritance_for",
        lambda gene, condition_ids: "recessive" if gene == "RARE1" else None,
    )
    monkeypatch.setattr(
        clingen,
        "lookup",
        lambda gene: {
            "validity": [
                {
                    "disease": "Rare recessive condition",
                    "mondo": "MONDO:0000001",
                    "classification": "Definitive",
                    "report_url": "https://clinicalgenome.org/validity/rare1",
                }
            ],
            "actionability": [],
        },
    )
    monkeypatch.setattr(
        clingen,
        "actionability_for",
        lambda gene: {
            "gene": gene,
            "context": "adult",
            "score": 10,
            "report_url": "https://clinicalgenome.org/actionability/rare1",
        },
    )
    finding = Finding(
        marker="1-100-A-G",
        source="clinvar_mirror",
        description="x",
        tier=Tier.ROBUST,
        categories=[Category.CLINICAL],
        detail={
            "gene": "RARE1",
            "clinical_significance": "Pathogenic",
            "conditions": ["Rare recessive condition"],
            "condition_ids": ["MONDO:0000001"],
            "zygosity": "het",
            "gold_stars": 2,
        },
    )

    interpret([finding])

    assert "carrier" in finding.interpretation.can_mean
    labels = [link.label for link in finding.evidence_chain]
    assert "ClinGen: Definitive for Rare recessive condition" in labels
    assert "ClinGen actionability score 10" in labels
