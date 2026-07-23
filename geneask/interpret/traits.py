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
from biocore.providers.base import Finding, Tier, Category


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


def trait_findings(vcf: str, trait_table: str) -> list[Finding]:
    """Interpret a trait table against the unified callset -> bio-core Findings."""
    traits = list(csv.DictReader(open(trait_table)))
    rsids = [t["rsid"] for t in traits if t.get("rsid")]
    calls = load_callset_by_rsid(vcf, rsids)
    findings = []
    for t in traits:
        rec = calls.get(t.get("rsid"))
        geno = rec["gt_bases"] if rec else "not_in_callset"
        ea = t.get("effect_allele", "").strip()
        carries = (ea in geno.split("/")) if (rec and ea) else None
        desc = f"{t.get('trait','trait')} ({t.get('gene','')}): genotype {geno}"
        if carries is not None:
            desc += f"; carries effect allele {ea}: {'yes' if carries else 'no'}"
        findings.append(Finding(
            marker=t.get("rsid", "?"), source="unified_callset",
            description=desc, tier=_tier_from_confidence(rec),
            categories=[Category.TRAIT]))
    return findings
