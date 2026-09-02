# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""ClinVar clinical screen against the bundled 157-gene panel.

Loads the cached ClinVar panel (ACMG SF v3.2 secondary-findings genes + carrier
+ pharmacogenes, 393,806 variants) and reports which pathogenic/likely-pathogenic
variants the individual carries, as bio-core Findings tagged CLINICAL.

The panel is bundled as reference data (data/reference/clinvar_panel_157genes.json.gz)
so the screen runs offline. Evidence tier follows ClinVar review status:
multi-star assertions -> ROBUST, single/none -> MODERATE, conflicting -> SPECULATIVE.

IMPORTANT (the recovered headline result): consumer-array 'pathogenic' calls are
frequently false positives — 10 of 11 in the source project were refuted by WGS.
This screen therefore reports the calling platform and defers to the unified
callset's multi-source confidence; a single-array pathogenic hit is never
presented as robust.
"""
from __future__ import annotations
import gzip, json
from pathlib import Path
from urllib.parse import quote
from biocore.providers.base import Finding, Tier, Category

_PANEL = Path(__file__).resolve().parents[1] / "data" / "reference" / "clinvar_panel_157genes.json.gz"

# Hereditary cancer-predisposition genes (ACMG SF v3.2 cancer set + common
# syndrome genes). A pathogenic hit here is tagged topic='cancer' so it groups
# and filters under the Cancer subject. Not exhaustive — extended as the panel grows.
_CANCER_GENES = {
    "BRCA1", "BRCA2", "PALB2", "TP53", "PTEN", "STK11", "CDH1", "APC", "MUTYH",
    "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",            # Lynch / mismatch-repair
    "RB1", "VHL", "MEN1", "RET", "NF1", "NF2", "TSC1", "TSC2", "SMAD4", "BMPR1A",
    "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX", "TMEM127", "WT1", "BAP1",
    "CHEK2", "ATM", "NBN", "BARD1", "BRIP1", "RAD51C", "RAD51D",  # breast/ovarian
    "CDKN2A", "CDK4", "MITF",                            # melanoma
    "FH", "FLCN", "MET", "HOXB13", "POLD1", "POLE", "GREM1", "NTHL1",
}


def load_panel(path: str | None = None) -> dict:
    """Load the bundled ClinVar panel, gene-keyed:
    {gene: {gene_id, release, n, variants: [{variant_id, clinical_significance,
             gold_stars, review_status, major_consequence, pos, ...}, ...]}}.
    """
    p = Path(path) if path else _PANEL
    with gzip.open(p, "rt") as fh:
        return json.load(fh)


def index_by_variant_id(panel: dict, wanted=None) -> dict:
    """Flatten to {variant_id: {gene, **variant_record}} keyed on 'chrom-pos-ref-alt'
    (GRCh38, matching the callset key).

    Prefers the FULL ClinVar mirror when a worker has built it (CLINVAR_MIRROR_DB):
    the complete ~4.2M-variant set replaces the bundled 157-gene panel as a drop-in.
    Falls back to flattening the bundled panel when no mirror exists.

    `wanted` is the set of variant_ids the caller will actually look up. Given
    one, only those rows are read: screening a callset against the mirror used to
    pull all ~4.2M rows into memory to answer a few thousand questions, which
    cost seconds of CPU and gigabytes of RSS on every single report.
    """
    try:
        if wanted is not None:
            from ..annotators.clinvar_mirror import lookup_from_mirror
            hit = lookup_from_mirror(wanted)
            if hit is not None:       # None means no mirror; {} means no matches
                return hit
        else:
            from ..annotators.clinvar_mirror import load_panel_from_mirror
            full = load_panel_from_mirror()
            if full:
                return full
    except Exception:
        pass
    idx = {}
    for gene, entry in panel.items():
        for v in entry.get("variants", []):
            rec = dict(v)
            rec["gene"] = gene
            idx[v["variant_id"]] = rec
    return idx


def _abbrev_allele(geno: str, keep: int = 8) -> str:
    """Shorten long indel alleles for display: a 331 bp deletion becomes
    'A/AGGAGG…(331bp)' instead of dumping the full sequence. Splits on '/' so a
    diploid genotype abbreviates each allele; short SNP alleles pass through."""
    def short(a: str) -> str:
        a = a.strip()
        return a if len(a) <= keep + 3 else f"{a[:keep]}…({len(a)}bp)"
    if not geno or geno == "?":
        return geno
    return "/".join(short(a) for a in geno.split("/"))


def _record_link(variation_id, gene: str) -> str | None:
    """Deep link to the variant's own ClinVar record.

    Both the bundled panel and the full mirror carry `clinvar_variation_id`, so
    a finding almost always resolves to one page. Without it, fall back to a
    gene-scoped search rather than the ClinVar homepage — a reader shown a
    pathogenic call needs somewhere to check the submitters and the review
    status, and the homepage is not that. The fallback deliberately avoids a
    /variation/ URL so a gene-level guess can never be mistaken for the variant's
    own record.
    """
    vid = str(variation_id or "").strip()
    if vid:
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{quote(vid)}/"
    g = (gene or "").strip()
    if g and g != "?":
        return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(g)}%5Bgene%5D"
    return None


def _tier_from_review(sig: str, gold_stars: int) -> Tier:
    sig = (sig or "").lower()
    if "conflicting" in sig:
        return Tier.SPECULATIVE
    if gold_stars >= 2:
        return Tier.ROBUST
    if gold_stars == 1:
        return Tier.MODERATE
    return Tier.SPECULATIVE


def screen_findings(carried_variant_ids: list[dict], panel: dict | None = None) -> list[Finding]:
    """Screen the individual's carried variants against the panel.

    carried_variant_ids entries: {variant_id: 'chr-pos-ref-alt', genotype, platform}.
    Returns CLINICAL Findings for pathogenic/likely-pathogenic panel hits only.
    """
    wanted = {v.get("variant_id") for v in carried_variant_ids}
    # Ask the mirror first and only fall back to the bundled panel if there
    # isn't one. Loading the panel up front meant every report on a mirrored
    # box gunzipped and JSON-parsed 393,806 variants and then discarded them,
    # because the mirror superseded the result.
    idx = None
    source = "clinvar_panel_157"
    if panel is None:
        try:
            from ..annotators.clinvar_mirror import lookup_from_mirror
            idx = lookup_from_mirror(wanted)      # None when no mirror exists
            if idx is not None:
                source = "clinvar_mirror"
        except Exception:
            idx = None
    if idx is None:
        idx = index_by_variant_id(panel if panel is not None else load_panel(),
                                  wanted=wanted)
    findings = []
    for v in carried_variant_ids:
        rec = idx.get(v.get("variant_id"))
        if not rec:
            continue
        sig = (rec.get("clinical_significance") or "").lower()
        if "pathogenic" not in sig:
            continue  # report only P/LP hits
        stars = int(rec.get("gold_stars", 0) or 0)
        tier = _tier_from_review(sig, stars)
        # a single-array call is demoted regardless of ClinVar stars (recovered caveat:
        # 10 of 11 array 'pathogenic' calls were WGS-refuted in the source project)
        if v.get("platform") == "ARRAY" and tier == Tier.ROBUST:
            tier = Tier.MODERATE
        gene = rec.get("gene", "?")
        # a cancer-predisposition gene -> topic cancer, else generic clinical
        topic = "cancer" if gene in _CANCER_GENES else "clinical"
        geno = _abbrev_allele(v.get("genotype", "?"))
        vid = rec.get("clinvar_variation_id")
        findings.append(Finding(
            marker=v.get("variant_id", "?"), source=source,
            description=f"{gene}: {rec.get('clinical_significance', sig)}"
                        f" (genotype {geno}, {v.get('platform','?')})",
            tier=tier, categories=[Category.CLINICAL],
            # the record itself, not the ClinVar homepage: bio-core's renderer
            # uses `f.link or SOURCES['clinvar'].url`, so leaving this unset sent
            # every clinical finding in a genome report to the same generic page.
            link=_record_link(vid, gene),
            detail={"gene": gene, "topic": topic, "modality": "genome",
                    "clinical_significance": rec.get("clinical_significance", sig),
                    "review_status": rec.get("review_status"),
                    # carried so JSON/MCP consumers resolve the record directly
                    # instead of parsing it back out of the link
                    "clinvar_variation_id": vid,
                    "gold_stars": stars, "platform": v.get("platform"),
                    "conditions": list(rec.get("conditions") or []),
                    "condition_ids": list(rec.get("condition_ids") or []),
                    "molecular_consequence": rec.get("molecular_consequence"),
                    "origin": list(rec.get("origin") or []),
                    "allele_id": rec.get("allele_id"),
                    "genotype": v.get("genotype"),
                    "zygosity": v.get("zygosity"),
                    "filter": v.get("filter"), "qual": v.get("qual"),
                    "gq": v.get("gq"), "dp": v.get("dp")}))
    return findings
