# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Trait interpretation — genotype lookup with multi-source confidence.

Given a trait table (rsid, effect_allele, trait, gene) and the unified
multi-source callset, report the individual's genotype at each SNP and map the
recovered multi-source confidence onto bio-core evidence tiers so GeneAsk
reports render through the same renderer as MethylAsk.

Refactored from the recovered trait_report.py. The bcftools query and the
plus-strand base-pair reporting (23andMe/Promethease convention) are preserved;
the addition is emitting bio-core Finding objects.
"""
from __future__ import annotations
import subprocess, csv
from pathlib import Path
from biocore.providers.base import Finding, Tier, Category

# a small bundled table of well-established, non-alarming consumer-genetics traits
# (caffeine/alcohol metabolism, lactase persistence, earwax type, ...). Used when
# a caller doesn't supply its own table.
DEFAULT_TRAIT_TABLE = str(Path(__file__).resolve().parents[1]
                          / "data" / "reference" / "trait_table.csv")


def _run(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


def load_callset_by_rsid(vcf: str, rsids) -> dict:
    """Single bcftools pass -> {rsid: {chrom,pos,ref,alt,gt_bases,platf,nsrc,conc,srcs}}.

    Genotypes are plus-strand base pairs (consumer raw-data convention).
    """
    want = set(rsids)
    out = {}
    q = _run(
        "bcftools query -f "
        "'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%GT]\\t%INFO/PLATF\\t%INFO/NSRC\\t%INFO/CONC\\t%INFO/SRCS\\n' "
        + vcf)
    for line in q.split("\n"):
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        chrom, pos, rsid, ref, alt, gt, platf, nsrc, conc, srcs = f[:10]
        if rsid not in want:
            continue
        alts = alt.split(",")
        def base(i):
            if i in (".", ""):
                return "."
            if i == "0":
                return ref
            try:
                return alts[int(i) - 1]
            except (ValueError, IndexError):
                return "."
        g = gt.replace("|", "/").split("/")
        if len(g) == 1:
            bases = base(g[0])
        else:
            bs = [base(x) for x in g]
            called = [b for b in bs if b != "."]
            bases = "/".join(sorted(bs)) if "." not in bs else (
                "/".join(called + ["."]) if called else "./.")
        out[rsid] = dict(chrom=chrom, pos=pos, ref=ref, alt=alt, gt_bases=bases,
                         platf=platf, nsrc=int(nsrc), conc=int(conc), srcs=srcs)
    return out


def _tier_from_confidence(rec) -> Tier:
    """Map recovered multi-source confidence onto a bio-core evidence tier.

    array+WGS agree, or >=3 sources agree -> ROBUST; single-platform -> MODERATE;
    array-vs-WGS conflict -> SPECULATIVE (flag the disagreement, don't hide it).
    """
    if rec is None:
        return Tier.UNKNOWN
    if rec["platf"] == "BOTH" and rec["conc"] == 1:
        return Tier.ROBUST
    if rec["platf"] == "BOTH" and rec["conc"] == 0:
        return Tier.SPECULATIVE  # array vs WGS disagree
    if rec["nsrc"] >= 3 and rec["conc"] == 1:
        return Tier.ROBUST
    if rec["platf"] in ("WGS", "ARRAY"):
        return Tier.MODERATE
    return Tier.SPECULATIVE


def _findings_from_calls(traits: list, calls: dict, source: str) -> list[Finding]:
    """Build trait Findings from a {rsid: rec} callset. rec needs 'gt_bases'
    ('A/G') and is passed to _tier_from_confidence for the tier."""
    findings = []
    for t in traits:
        rec = calls.get(t.get("rsid"))
        if rec is None:
            continue  # the person's data doesn't cover this SNP -> don't report it
        geno = rec["gt_bases"]
        ea = t.get("effect_allele", "").strip()
        carries = (ea in geno.split("/")) if ea else None
        desc = f"{t.get('trait','trait')} ({t.get('gene','')}): genotype {geno}"
        if carries is not None:
            desc += f"; carries effect allele {ea}: {'yes' if carries else 'no'}"
        findings.append(Finding(
            marker=t.get("rsid", "?"), source=source,
            description=desc, tier=_tier_from_confidence(rec),
            categories=[Category.TRAIT]))
    return findings


def trait_findings(vcf: str, trait_table: str) -> list[Finding]:
    """Interpret a trait table against a unified VCF callset -> bio-core Findings."""
    traits = list(csv.DictReader(open(trait_table)))
    rsids = [t["rsid"] for t in traits if t.get("rsid")]
    calls = load_callset_by_rsid(vcf, rsids)
    return _findings_from_calls(traits, calls, source="unified_callset")


def trait_findings_from_parse(parsed, trait_table: str) -> list[Finding]:
    """Interpret a trait table against a parsed array callset (no bcftools/VCF).
    Matches by rsID directly from the vendor parser's GenotypeRecords."""
    traits = list(csv.DictReader(open(trait_table)))
    by_rsid = {r.rsid: r for r in parsed.records if not r.is_nocall}
    calls = {}
    for t in traits:
        r = by_rsid.get(t.get("rsid"))
        if r is None:
            continue
        gt = "/".join(sorted(a for a in (r.allele1, r.allele2) if a))
        # array calls carry no per-source confidence; mark platform ARRAY so the
        # tier heuristic treats them as single-source array (not WGS-confirmed).
        calls[r.rsid] = {"gt_bases": gt, "platf": "ARRAY", "nsrc": 1,
                         "conc": 1, "srcs": "array"}
    return _findings_from_calls(traits, calls, source="array_callset")
