# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Reviewed copy for clinical findings. The tables are data; this module only
looks them up and fills the templates.

Each table carries `version`, `reviewed_by` and `reviewed_on`. A report copies
those onto every Interpretation it renders, so a reader can see whether a
person has read the wording that reached them.
"""
from __future__ import annotations
import json, functools
from pathlib import Path
from .lists import acmg_sf

_DIR = Path(__file__).resolve().parents[1] / "data" / "copy"
_TABLES = {"clinical_next_step": _DIR / "clinical_next_step.json",
           "gene_function": _DIR / "gene_function.json",
           "condition_phrasing": _DIR / "condition_phrasing.json"}


@functools.lru_cache(maxsize=None)
def _load(name: str) -> dict:
    return json.loads(_TABLES[name].read_text())


def copy_meta(name: str) -> dict:
    t = _load(name)
    return {"version": str(t.get("version", "")),
            "reviewed_by": list(t.get("reviewed_by") or [])}


def classify(sig: str) -> str:
    """Classification class for a ClinVar significance string.

    Same precedence as biocore.report.render.direction(): uncertainty first,
    then disease, then drug response, then benign. "Conflicting classifications
    of pathogenicity" contains "pathogenic", which is why uncertainty is tested
    first.
    """
    s = (sig or "").lower()
    if not s:
        return "other"
    if "conflicting" in s:
        return "conflicting"
    if "uncertain" in s or "not provided" in s or "no classification" in s:
        return "vus"
    if "pathogenic" in s:
        return "plp"
    if "risk factor" in s or "risk allele" in s:
        return "risk_factor"
    if "drug response" in s:
        return "drug_response"
    if "benign" in s:
        return "benign"
    return "other"


def platform_class(platform: str | None) -> str:
    return "array" if (platform or "").upper() == "ARRAY" else "wgs_wes"


_PLATFORM_WORDS = {"wgs_wes": "whole-genome or exome sequencing",
                   "array": "a consumer genotyping array"}


def next_step(cls: str, plat: str, zygosity: str | None, *,
              stars: int | None = None) -> dict:
    """{how_sure, next_step} for a classification class, platform and zygosity.
    Falls back from the most specific key to the least."""
    rows = _load("clinical_next_step")["rows"]
    fill = {"stars": "an unknown number" if stars is None else str(stars),
            "platform_words": _PLATFORM_WORDS.get(plat, _PLATFORM_WORDS["wgs_wes"])}
    for key in (f"{cls}|{plat}|{zygosity or 'any'}", f"{cls}|{plat}|any",
                f"{cls}|any|{zygosity or 'any'}", f"{cls}|any|any", "other|any|any"):
        if key in rows:
            return {k: v.format(**fill) for k, v in rows[key].items()}
    return {"how_sure": "", "next_step": ""}


def gene_function(gene: str) -> dict | None:
    """{sentence, source, url} for a gene symbol, or None when there is no row."""
    return _load("gene_function")["rows"].get((gene or "").strip().upper())


def condition_phrase(condition_ids: list[str], fallback_name: str, gene: str) -> dict:
    """Return the condition name, ClinGen inheritance, and reviewed URL."""
    rows = _load("condition_phrasing")["rows"]
    row = None
    for cid in condition_ids or []:
        if cid in rows:
            row = dict(rows[cid])
            break
    if row is None:
        by_gene = rows.get((gene or "").strip().upper())
        row = dict(by_gene) if by_gene else None

    from geneask.annotators import clingen
    inheritance = clingen.inheritance_for(gene, condition_ids)
    if row:
        if inheritance:
            row["inheritance"] = inheritance
        return row
    entry = acmg_sf(gene)
    if entry:
        return {
            "name": entry["condition"],
            "inheritance": inheritance or entry["inheritance"],
            "url": None,
        }
    return {"name": (fallback_name or "").replace("_", " ").strip(),
            "inheritance": inheritance, "url": None}
