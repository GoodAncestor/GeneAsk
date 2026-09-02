# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Compose what a genome finding means, from its fields and the reviewed copy.

Nothing here invents a fact. Each sentence comes from a field on the finding
or from a copy row that names its source. The composer fills
``Finding.interpretation`` (four parts) and ``Finding.evidence_chain`` (the
path from the variant to its evidence), and leaves ``description`` alone for
callers that only know the old sentence.
"""
from __future__ import annotations
from biocore.providers.base import Finding, Interpretation, ChainLink
from .copy import (classify, platform_class, next_step, gene_function,
                   condition_phrase, copy_meta)
from .lists import acmg_sf, ACMG_SF_URL, ACMG_SF_VERSION

_CLINVAR_SOURCES = {"clinvar", "clinvar_mirror", "clinvar_panel_157"}

_ZYG = {"het": "One altered copy", "hom": "Two altered copies",
        "hemi": "One copy, on the X or Y chromosome",
        "unknown": "An unknown number of copies"}


def _sig_words(sig: str) -> str:
    return (sig or "").replace(";", ",").replace("_", " ").strip().lower()


def _can_mean_clinvar(cls: str, zyg: str | None, cond: dict) -> str:
    name = cond["name"] or "the named condition"
    inh = cond.get("inheritance")
    if cls == "plp":
        if inh == "dominant":
            if zyg == "hom":
                return (f"Two altered copies were read. For a dominant condition one is enough "
                        f"to raise the chance of {name}. It does not mean the condition is present.")
            return (f"One altered copy is enough to raise the chance of {name}. "
                    f"It does not mean the condition is present or certain.")
        if inh == "recessive":
            if zyg == "hom":
                return (f"Two altered copies were read. For {name} that is the pattern that can "
                        f"cause the condition. A clinician can say whether it applies.")
            return (f"One altered copy makes a person a carrier of {name}. A carrier does not "
                    f"usually have the condition. It matters for children.")
        if inh == "x_linked":
            return (f"This change sits on the X chromosome. Its effect on {name} depends on sex "
                    f"and on the number of copies. A clinician can say whether it applies.")
        return (f"Laboratories link this change to {name}. Whether it applies depends on how "
                f"the condition is inherited and on the number of copies.")
    if cls == "vus":
        return f"Laboratories have linked this gene to {name}, but not this change."
    if cls == "conflicting":
        return f"Some laboratories link this change to {name}. Others do not."
    if cls == "benign":
        return f"This change does not cause {name}."
    if cls == "drug_response":
        return f"This change can alter how some medicines work. ClinVar names {name}."
    if cls == "risk_factor":
        return f"This change shifts the odds of {name}. It does not decide the outcome."
    return f"ClinVar records an association with {name} without a classification."


def interpret_clinvar(f: Finding) -> Interpretation:
    d = f.detail or {}
    gene = str(d.get("gene") or "").strip()
    sig = d.get("clinical_significance", "")
    cls = classify(sig)
    plat = platform_class(d.get("platform"))
    zyg = d.get("zygosity")
    stars = d.get("gold_stars")
    cond = condition_phrase(list(d.get("condition_ids") or []),
                            (d.get("conditions") or [""])[0], gene)
    gf = gene_function(gene)

    found = [gf["sentence"] if gf else f"A change in {gene or 'a gene'}."]
    found.append(f"ClinVar classifies this change as {_sig_words(sig)} for {cond['name']}.")
    if d.get("molecular_consequence"):
        found.append(f"The change is a {str(d['molecular_consequence']).replace('_', ' ')}.")

    ns = next_step(cls, plat, zyg, stars=None if stars is None else int(stars))
    how = ns["how_sure"]
    g = d.get("gnomad") or {}
    if g.get("ac") is not None and g.get("an"):
        ver = f" {g['version']}" if g.get("version") else ""
        how += (f" gnomAD{ver} saw this change in {int(g['ac']):,} of "
                f"{int(g['an']):,} sampled chromosomes.")

    can = _can_mean_clinvar(cls, zyg, cond)
    if cls == "plp" and zyg in _ZYG:
        can = f"{_ZYG[zyg]} was read. " + can

    cites = []
    if gf:
        cites.append(ChainLink(kind="gene", label=gene, url=gf["url"]))
    if cond.get("url"):
        cites.append(ChainLink(kind="condition", label=cond["name"], url=cond["url"]))
    if f.link:
        vid = d.get("clinvar_variation_id")
        cites.append(ChainLink(kind="assertion", label="ClinVar record",
                               id=f"ClinVar:{vid}" if vid else None, url=f.link))
    meta = copy_meta("clinical_next_step")
    return Interpretation(found=" ".join(found), can_mean=can, how_sure=how,
                          next_step=ns["next_step"], condition=cond["name"] or None,
                          condition_ids=list(d.get("condition_ids") or []),
                          zygosity=zyg, citations=cites,
                          copy_version=meta["version"], reviewed_by=meta["reviewed_by"])


def _chain_clinvar(f: Finding, ip: Interpretation) -> list[ChainLink]:
    d = f.detail or {}
    gene = str(d.get("gene") or "")
    gf = gene_function(gene)
    vid = d.get("clinvar_variation_id")
    chain = [
        ChainLink(kind="variant", label=f.marker,
                  id=f"ClinVar:{vid}" if vid else None, url=f.link),
        ChainLink(kind="gene", label=gene or "gene",
                  url=(gf or {}).get("url")
                  or f"https://www.ncbi.nlm.nih.gov/gene/?term={gene}%5Bsym%5D"),
        ChainLink(kind="condition", label=ip.condition or "condition not named",
                  id=(ip.condition_ids or [None])[0],
                  url=next((c.url for c in ip.citations if c.kind == "condition"), None)),
    ]
    if f.link:
        chain.append(ChainLink(kind="assertion",
                               label=str(d.get("review_status") or "ClinVar assertion"),
                               id=f"stars:{d.get('gold_stars')}" if d.get("gold_stars") is not None else None,
                               url=f.link))
    if acmg_sf(gene):
        chain.append(ChainLink(kind="paper", label=f"ACMG {ACMG_SF_VERSION} reportable gene",
                               url=ACMG_SF_URL))
    for p in f.pmids or []:
        chain.append(ChainLink(kind="paper", label=f"PMID {p}", id=f"PMID:{p}",
                               url=f"https://pubmed.ncbi.nlm.nih.gov/{p}/"))
    return chain


def interpret_gwas(f: Finding) -> Interpretation:
    d = f.detail or {}
    trait = d.get("trait") or "this trait"
    ra = d.get("risk_allele") or ""
    carried = d.get("risk_allele_carried")
    allele = f"the {ra} allele" if ra else "an allele"
    if carried is True:
        found = f"A study links {allele} at {f.marker} to {trait}. You carry the allele."
    elif carried is False:
        found = f"A study links {allele} at {f.marker} to {trait}. You do not carry the allele."
    else:
        found = (f"A study links {f.marker} to {trait}. Whether you carry the studied allele "
                 f"was not determined.")
    et, ev = d.get("effect_type"), d.get("effect")
    if et == "or" and ev:
        can = (f"People with the allele had {float(ev):g} times the odds of {trait} in the study. "
               f"This is a group average. It is not a personal risk.")
    elif et == "beta" and ev is not None:
        unit = (str(d.get("ci_text") or "").split("]")[-1].strip()) or "units"
        can = (f"Each copy of the allele moved {trait} by {float(ev):g} {unit} on average in "
               f"the study. This is a group average. It is not a personal prediction.")
    else:
        can = (f"The study reports an association with {trait}. The effect size is not recorded "
               f"in the catalogue, so its size cannot be said here.")
    p = d.get("p")
    try:
        how = f"The association reached p = {float(p):.0e}." if p not in (None, "") else "The p-value is not recorded."
    except (TypeError, ValueError):
        how = "The p-value is not recorded."
    if d.get("initial_n"):
        how += f" The study read {d['initial_n']}."
    if d.get("replication_n"):
        how += f" Replication: {d['replication_n']}."
    else:
        how += " The catalogue records no replication sample."
    cites = [ChainLink(kind="paper", label=f"PMID {pm}", id=f"PMID:{pm}",
                       url=f"https://pubmed.ncbi.nlm.nih.gov/{pm}/") for pm in (f.pmids or [])]
    return Interpretation(found=found, can_mean=can, how_sure=how, next_step="",
                          condition=None, condition_ids=[], zygosity=None, citations=cites,
                          copy_version="inline", reviewed_by=[])


def _chain_gwas(f: Finding, ip: Interpretation) -> list[ChainLink]:
    d = f.detail or {}
    chain = [ChainLink(kind="variant", label=f.marker, url=f.link),
             ChainLink(kind="trait", label=str(d.get("trait") or "trait"),
                       id=d.get("efo"), url=d.get("efo"))]
    if d.get("mapped_gene"):
        g = str(d["mapped_gene"]).split(" - ")[0].split(",")[0].strip()
        chain.append(ChainLink(kind="gene", label=g,
                               url=f"https://www.ncbi.nlm.nih.gov/gene/?term={g}%5Bsym%5D"))
    if d.get("accession"):
        acc = str(d["accession"])
        chain.append(ChainLink(kind="assertion", label=f"GWAS Catalog study {acc}", id=acc,
                               url=f"https://www.ebi.ac.uk/gwas/studies/{acc}"))
    chain.extend(ip.citations)
    return chain


def interpret_cpic(f: Finding) -> Interpretation:
    d = f.detail or {}
    gene, drug = str(d.get("gene") or f.marker), str(d.get("drug") or "a medicine")
    found = (f"You carry a change in {gene}. CPIC publishes prescribing guidance for {drug} "
             f"by {gene} status.")
    can = (f"For some {gene} types, {drug} works differently or needs a different dose. "
           f"The guidance depends on the type, which this report did not determine.")
    how = (f"CPIC evidence level {d.get('cpic_level') or '?'}. This report did not determine "
           f"your {gene} type (diplotype), so the specific recommendation cannot be selected here.")
    nxt = f"Tell a prescriber or pharmacist that you carry a {gene} change before starting {drug}."
    return Interpretation(found=found, can_mean=can, how_sure=how, next_step=nxt,
                          condition=None, condition_ids=[], zygosity=None,
                          citations=[ChainLink(kind="assertion", label=f"CPIC {gene}–{drug}",
                                               url="https://cpicpgx.org/genes-drugs/")],
                          copy_version="inline", reviewed_by=[])


def interpret(findings: list[Finding]) -> int:
    """Fill interpretation and evidence_chain in place. Returns how many were filled."""
    n = 0
    for f in findings:
        d = f.detail or {}
        if f.source in _CLINVAR_SOURCES and d.get("clinical_significance"):
            f.interpretation = interpret_clinvar(f)
            f.evidence_chain = _chain_clinvar(f, f.interpretation)
        elif f.source == "gwas_catalog":
            f.interpretation = interpret_gwas(f)
            f.evidence_chain = _chain_gwas(f, f.interpretation)
        elif f.source == "cpic":
            f.interpretation = interpret_cpic(f)
            f.evidence_chain = [ChainLink(kind="gene", label=str(d.get("gene") or f.marker)),
                                *f.interpretation.citations]
        else:
            continue
        n += 1
    return n
