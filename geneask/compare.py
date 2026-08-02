# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Compare a person's genome across multiple tests/platforms.

The "do my tests agree?" feature: given a multi-sample VCF where each sample is
one test of the same person (e.g. five 23andMe chips + a whole-genome sequence),
this reads per-sample base-pair genotypes and drives bio-core's organism-
agnostic comparison primitives — pairwise concordance, discordance typing, and
KING-robust relatedness (the sanity check that every "sample" really is the
same person and not a mix-up).

Genotype comparison mechanism lives in bio-core (biocore.compare.genotype_calls);
this module is the GeneAsk-specific glue: read genotypes from a VCF into the
{sample: {marker: "A/G"}} form bio-core expects, and phrase the result for a
consumer report. Recovered from the "six genomes" concordance analysis.
"""
from __future__ import annotations
import subprocess

from biocore.compare.genotype_calls import (
    concordance_matrix, concordance_pair, discordance_breakdown, king_relatedness,
)
from biocore.providers.base import Finding, Tier, Category


def _run(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True,
                          text=True).stdout


def load_calls_by_sample(vcf: str) -> dict:
    """Read a multi-sample VCF -> {sample: {marker_id: "A/G"}}.

    Genotypes are emitted as plus-strand base-pair strings (sorted alleles),
    matching the consumer raw-data convention; missing is "./." or ".". marker_id
    is the VCF ID column (rsID) when present, else chrom:pos.
    """
    hdr = _run(f"bcftools query -l {vcf}").split()
    calls = {s: {} for s in hdr}
    fmt = "'%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n'"
    q = _run(f"bcftools query -f {fmt} {vcf}")
    for line in q.split("\n"):
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 5 + len(hdr):
            continue
        chrom, pos, rid, ref, alt = f[:5]
        gts = f[5:5 + len(hdr)]
        alts = alt.split(",")
        marker = rid if rid not in (".", "") else f"{chrom}:{pos}"

        def base(i):
            if i in (".", ""):
                return "."
            if i == "0":
                return ref
            try:
                return alts[int(i) - 1]
            except (ValueError, IndexError):
                return "."

        for sample, gt in zip(hdr, gts):
            g = gt.replace("|", "/").split("/")
            if len(g) == 1:
                bases = base(g[0])
            else:
                bs = [base(x) for x in g]
                called = [b for b in bs if b != "."]
                bases = "/".join(sorted(bs)) if "." not in bs else (
                    "/".join(called + ["."]) if called else "./.")
            calls[sample][marker] = bases
    return calls


def compare_samples(vcf: str, *, autosomal_only: bool = True) -> dict:
    """Full multi-test comparison from a VCF. Returns a dict with:
      samples, concordance (matrix + labels + overlap), discordance (per pair),
      relatedness (KING per pair). Autosomal restriction for KING drops chrX/Y/M
      markers (identified by marker id prefix) to keep kinship unbiased.
    """
    calls = load_calls_by_sample(vcf)
    labels, mat, overlap = concordance_matrix(calls)

    auto = None
    if autosomal_only:
        auto = set()
        for m in {k for c in calls.values() for k in c}:
            head = m.split(":")[0].lower().lstrip("chr")
            # keep autosomes (numeric) and unprefixed rsIDs (assume autosomal
            # unless clearly a sex/mito contig)
            if head.isdigit() or m.lower().startswith("rs"):
                auto.add(m)

    from itertools import combinations
    disc, king = {}, {}
    for a, b in combinations(labels, 2):
        cats, _ex, nsh = discordance_breakdown(calls[a], calls[b])
        disc[(a, b)] = dict(categories=cats, n_shared=nsh,
                            n_discordant=sum(cats.values()))
        king[(a, b)] = king_relatedness(calls[a], calls[b], restrict_markers=auto)
    return dict(samples=labels, concordance=dict(labels=labels, matrix=mat,
                overlap=overlap), discordance=disc, relatedness=king)


def compare_findings(vcf: str) -> list:
    """Phrase the comparison as bio-core Findings for a report: a self-identity
    check (KING) and an overall concordance statement. Descriptive, not clinical.
    """
    res = compare_samples(vcf)
    findings = []
    labels = res["samples"]
    # identity check: min kinship across pairs; all should be ~0.5 (same person)
    if res["relatedness"]:
        kins = [v["kinship"] for v in res["relatedness"].values()
                if v["kinship"] == v["kinship"]]  # drop NaN
        if kins:
            lo = min(kins)
            same = lo > 0.35   # comfortably above parent-child (~0.25)
            findings.append(Finding(
                marker="identity_check", source="geneask.compare",
                description=(f"All {len(labels)} tests are the same individual "
                             f"(min pairwise KING kinship {lo:.3f}, expected ~0.5 for self)"
                             if same else
                             f"Warning: a test pair has kinship {lo:.3f}, below the self "
                             f"threshold — possible sample mix-up"),
                tier=(Tier.ROBUST if same else Tier.MODERATE),
                categories=[Category.TRAIT]))
    # overall concordance (mean of off-diagonal)
    mat, labs = res["concordance"]["matrix"], res["concordance"]["labels"]
    offdiag = [mat[i][j] for i in range(len(labs)) for j in range(len(labs))
               if i < j and mat[i][j] == mat[i][j]]
    if offdiag:
        mean_c = sum(offdiag) / len(offdiag)
        findings.append(Finding(
            marker="cross_test_concordance", source="geneask.compare",
            description=f"Tests agree {100*mean_c:.2f}% of the time on shared, "
                        f"fully-called markers (mean over {len(offdiag)} pairs); "
                        f"disagreements are typically technical, not biological",
            tier=Tier.MODERATE, categories=[Category.TRAIT]))
    return findings
