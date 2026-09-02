# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Published gene lists the report leans on for triage and phrasing.

ACMG SF: the genes the American College of Medical Genetics tells clinical
laboratories to report when they find a pathogenic variant, because something
can be done about it. The list is data, so a new version is a file change.
Gene and phenotype come from ClinVar's ACMG page; inheritance from the SF v3.2
paper (PMID 37347242).
"""
from __future__ import annotations
import json, functools
from pathlib import Path

_ACMG = Path(__file__).resolve().parents[1] / "data" / "reference" / "acmg_sf_v3_2.json"

ACMG_SF_VERSION = "SF v3.2"
ACMG_SF_URL = "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"


@functools.lru_cache(maxsize=1)
def _acmg() -> dict:
    return json.loads(_ACMG.read_text())


def all_acmg_sf() -> dict:
    """{gene: {condition, inheritance, note}} for every gene on the list."""
    return dict(_acmg()["rows"])


def acmg_sf(gene: str) -> dict | None:
    """The list entry for one gene, or None when the gene is not on the list."""
    key = (gene or "").strip().upper()
    row = _acmg()["rows"].get(key)
    return {"gene": key, **row} if row else None
