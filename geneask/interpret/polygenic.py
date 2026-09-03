# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Offline polygenic position summaries from GWAS effects and gnomAD AF."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path


CAVEAT = (
    "Population reference is the whole gnomAD set; effect sizes and frequencies "
    "vary by ancestry."
)


@dataclass
class TraitScore:
    trait: str
    efo: str | None
    n_variants: int
    n_with_af: int
    score: float
    mean: float
    sd: float
    z: float | None
    percentile: int | None
    direction_word: str
    top: list = field(default_factory=list)
    caveat: str = CAVEAT


def _normal_trait(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value or "").lower()).strip()


def _weight(detail: dict) -> float | None:
    try:
        effect = float(detail.get("effect"))
    except (TypeError, ValueError):
        return None
    effect_type = str(detail.get("effect_type") or "").lower()
    if effect_type == "or":
        return math.log(effect) if effect > 0 else None
    if effect_type == "beta":
        return effect
    return None


def _dosage(detail: dict) -> int | None:
    carried = detail.get("risk_allele_carried")
    if carried is None:
        return None
    if carried is False:
        return 0
    return 2 if str(detail.get("zygosity") or "").lower() == "hom" else 1


def _is_variant_id(value: str) -> bool:
    parts = str(value or "").split("-", 3)
    return len(parts) == 4 and parts[1].isdigit()


def _entry(finding) -> dict | None:
    detail = finding.detail or {}
    weight = _weight(detail)
    dosage = _dosage(detail)
    trait = str(detail.get("trait") or "").strip()
    if weight is None or dosage is None or not trait:
        return None
    marker = str(finding.marker or "")
    variant_id = str(detail.get("variant_id") or "")
    if not _is_variant_id(variant_id) and _is_variant_id(marker):
        variant_id = marker
    chrom = detail.get("chrom")
    pos = detail.get("pos")
    risk_allele = str(detail.get("risk_allele") or "").upper()
    if not variant_id and chrom and pos and detail.get("ref"):
        alt = detail.get("alt") or risk_allele
        if alt:
            variant_id = f"{chrom}-{pos}-{detail['ref']}-{alt}"
    lookup_key = variant_id or marker
    efo = str(detail.get("efo") or "").strip() or None
    gene = str(detail.get("mapped_gene") or "").split(" - ")[0].split(",")[0].strip()
    return {
        "trait": trait,
        "efo": efo,
        "marker": marker,
        "lookup_key": lookup_key,
        "variant_id": variant_id,
        "chrom": chrom,
        "pos": pos,
        "risk_allele": risk_allele,
        "gene": gene,
        "weight": weight,
        "dosage": dosage,
    }


def _af_value(value) -> float | None:
    if isinstance(value, dict):
        value = value.get("af")
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and 0 <= value <= 1 else None


def _cache_frequencies(entries: list[dict], frequencies: dict[str, float]) -> None:
    from geneask.annotators.gnomad_freq import _cache_con, _cache_path

    path = _cache_path()
    if not Path(path).exists():
        return
    con = _cache_con(path)
    try:
        for entry in entries:
            if entry["lookup_key"] in frequencies:
                continue
            if entry["variant_id"]:
                row = con.execute(
                    "SELECT af FROM af WHERE variant_id=?", (entry["variant_id"],)
                ).fetchone()
                if row:
                    af = _af_value(row[0])
                    if af is not None:
                        frequencies[entry["lookup_key"]] = af
                continue
            if not entry["chrom"] or not entry["pos"]:
                continue
            chrom = str(entry["chrom"])
            chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
            rows = con.execute(
                "SELECT variant_id, af FROM af WHERE variant_id LIKE ?",
                (f"{chrom}-{entry['pos']}-%",),
            )
            for variant_id, raw_af in rows:
                if str(variant_id).rsplit("-", 1)[-1].upper() == entry["risk_allele"]:
                    af = _af_value(raw_af)
                    if af is not None:
                        frequencies[entry["lookup_key"]] = af
                        break
    finally:
        con.close()


def _default_frequencies(entries: list[dict]) -> dict[str, float]:
    from geneask.annotators.gnomad_mirror import lookup_by_position, lookup_many

    exact = [entry["variant_id"] for entry in entries if entry["variant_id"]]
    mirror_rows = lookup_many(exact)
    frequencies = {}
    for entry in entries:
        af = _af_value(mirror_rows.get(entry["variant_id"]))
        if af is not None:
            frequencies[entry["lookup_key"]] = af
    for entry in entries:
        if entry["lookup_key"] in frequencies or not entry["chrom"] or not entry["pos"]:
            continue
        rows = lookup_by_position(entry["chrom"], entry["pos"])
        for variant_id, raw_af in rows.items():
            if variant_id.rsplit("-", 1)[-1].upper() == entry["risk_allele"]:
                af = _af_value(raw_af)
                if af is not None:
                    frequencies[entry["lookup_key"]] = af
                    break
    _cache_frequencies(entries, frequencies)
    return frequencies


def trait_scores(findings, *, af_lookup=None) -> list[TraitScore]:
    """Build one population-position score for each sufficiently covered trait."""
    groups: dict[str, list[dict]] = {}
    for finding in findings or []:
        if getattr(finding, "source", "") != "gwas_catalog":
            continue
        entry = _entry(finding)
        if entry is None:
            continue
        key = f"efo:{entry['efo']}" if entry["efo"] else f"trait:{_normal_trait(entry['trait'])}"
        groups.setdefault(key, []).append(entry)

    entries = [entry for group in groups.values() for entry in group]
    if af_lookup is None:
        frequencies = _default_frequencies(entries)
    else:
        raw = af_lookup([entry["lookup_key"] for entry in entries]) or {}
        frequencies = {
            key: af for key, value in raw.items() if (af := _af_value(value)) is not None
        }

    scores = []
    for group in groups.values():
        with_af = [(entry, frequencies.get(entry["lookup_key"])) for entry in group]
        with_af = [(entry, af) for entry, af in with_af if af is not None]
        if len(with_af) < 3:
            continue
        score = sum(entry["weight"] * entry["dosage"] for entry, _af in with_af)
        mean = sum(2 * af * entry["weight"] for entry, af in with_af)
        variance = sum(
            2 * af * (1 - af) * entry["weight"] ** 2 for entry, af in with_af
        )
        sd = math.sqrt(variance)
        z = (score - mean) / sd if sd > 0 else None
        if z is None:
            percentile = None
            direction = "about average"
        else:
            raw_percentile = 100 * 0.5 * (1 + math.erf(z / math.sqrt(2)))
            percentile = max(10, min(90, int(round(raw_percentile / 10) * 10)))
            direction = "higher" if z >= 0.5 else "lower" if z <= -0.5 else "about average"
        top = sorted(
            [
                (entry["marker"], entry["gene"], entry["weight"] * entry["dosage"])
                for entry, _af in with_af
            ],
            key=lambda item: abs(item[2]),
            reverse=True,
        )[:5]
        first = group[0]
        scores.append(
            TraitScore(
                trait=first["trait"],
                efo=first["efo"],
                n_variants=len(group),
                n_with_af=len(with_af),
                score=score,
                mean=mean,
                sd=sd,
                z=z,
                percentile=percentile,
                direction_word=direction,
                top=top,
                caveat=CAVEAT,
            )
        )
    return scores
