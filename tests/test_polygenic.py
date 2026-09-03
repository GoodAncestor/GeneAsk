import math
import sqlite3

import pytest

from biocore.providers.base import Category, Finding, Tier
from geneask.annotators.gnomad_mirror import lookup_by_position
from geneask.interpret.polygenic import CAVEAT, trait_scores


def _finding(marker, effect_type, effect, carried, zygosity, gene):
    return Finding(
        marker=marker,
        source="gwas_catalog",
        description="association",
        tier=Tier.ROBUST,
        categories=[Category.TRAIT],
        detail={
            "trait": "Type 2 diabetes",
            "efo": "EFO_0001360",
            "effect_type": effect_type,
            "effect": effect,
            "risk_allele_carried": carried,
            "zygosity": zygosity,
            "mapped_gene": gene,
        },
    )


def test_three_variant_score_matches_the_hand_calculation():
    findings = [
        _finding("1-100-A-G", "or", 1.5, True, "het", "GENE1"),
        _finding("1-200-C-T", "or", 1.2, True, "hom", "GENE2"),
        _finding("1-300-G-A", "beta", 0.3, False, "hom", "GENE3"),
    ]
    frequencies = {"1-100-A-G": 0.3, "1-200-C-T": 0.5, "1-300-G-A": 0.1}

    score, = trait_scores(findings, af_lookup=lambda ids: frequencies)

    weights = [math.log(1.5), math.log(1.2), 0.3]
    expected_mean = sum(2 * p * weight for p, weight in zip((0.3, 0.5, 0.1), weights))
    expected_var = sum(
        2 * p * (1 - p) * weight * weight
        for p, weight in zip((0.3, 0.5, 0.1), weights)
    )
    expected_score = weights[0] + 2 * weights[1]
    expected_z = (expected_score - expected_mean) / math.sqrt(expected_var)

    assert score.n_variants == 3
    assert score.n_with_af == 3
    assert score.score == pytest.approx(expected_score, abs=0.001)
    assert score.mean == pytest.approx(expected_mean, abs=0.001)
    assert score.sd == pytest.approx(math.sqrt(expected_var), abs=0.001)
    assert score.z == pytest.approx(expected_z, abs=0.001)
    assert score.percentile == 80
    assert score.direction_word == "higher"
    assert score.top[0][0:2] == ("1-100-A-G", "GENE1")
    assert score.caveat == CAVEAT


def test_score_is_omitted_with_fewer_than_three_frequency_rows():
    findings = [
        _finding("1-100-A-G", "or", 1.5, True, "het", "GENE1"),
        _finding("1-200-C-T", "or", 1.2, True, "hom", "GENE2"),
        _finding("1-300-G-A", "beta", 0.3, False, "hom", "GENE3"),
    ]

    assert trait_scores(
        findings,
        af_lookup=lambda ids: {"1-100-A-G": 0.3, "1-200-C-T": 0.5},
    ) == []


def test_unknown_carriage_and_invalid_odds_ratios_do_not_count():
    findings = [
        _finding("1-100-A-G", "or", 0, True, "het", "GENE1"),
        _finding("1-200-C-T", "beta", 0.2, None, "het", "GENE2"),
    ]

    assert trait_scores(findings, af_lookup=lambda ids: {}) == []


def test_position_lookup_returns_every_alt_at_the_site(tmp_path):
    database = tmp_path / "gnomad.db"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE af(variant_id TEXT PRIMARY KEY, af REAL)")
        connection.execute("INSERT INTO af VALUES ('1-100-A-G', 0.3)")
        connection.execute("INSERT INTO af VALUES ('1-100-A-T', 0.1)")
        connection.execute("INSERT INTO af VALUES ('1-101-C-T', 0.2)")

    assert lookup_by_position("chr1", 100, db_path=str(database)) == {
        "1-100-A-G": 0.3,
        "1-100-A-T": 0.1,
    }
