"""Turn parsed array genotypes into carried-variant records — panel-anchored.

Array data gives (chrom, pos, allele1, allele2) but no REF/ALT. The ClinVar
panel already knows REF/ALT for every variant it screens (variant_id =
chrom-pos-ref-alt, GRCh38). So rather than anchor every array site to a whole
reference FASTA, we anchor only at panel positions: for each panel variant, if
the person's call at that (chrom,pos) carries the panel's ALT allele, emit that
variant_id as carried. This produces exactly the {variant_id, genotype, platform}
records GeneAsk's clinvar screen consumes — no FASTA, no bcftools, no VCF step.

Coordinate note: array builds vary (GRCh37/38); the panel is GRCh38. A build
mismatch means positions won't line up and nothing matches — the caller should
lift or detect build first. This converter assumes GRCh38-aligned input.
"""
from __future__ import annotations
from .base import ParseResult


def _index_panel_by_pos(panel: dict) -> dict:
    """{(chrom, pos): [(ref, alt, variant_id), ...]} from the gene-keyed panel."""
    idx: dict = {}
    for gene_rec in panel.values():
        for v in gene_rec.get("variants", []):
            vid = v.get("variant_id", "")
            parts = vid.split("-")
            if len(parts) != 4:
                continue
            chrom, pos, ref, alt = parts
            try:
                key = (chrom, int(pos))
            except ValueError:
                continue
            idx.setdefault(key, []).append((ref, alt, vid))
    return idx


def carried_from_parse(parsed: ParseResult, panel: dict,
                       platform: str = "ARRAY") -> list[dict]:
    """Match parsed array records against the panel; return carried-variant dicts
    {variant_id, genotype, platform} for panel variants whose ALT the person carries."""
    pos_idx = _index_panel_by_pos(panel)
    out: list[dict] = []
    for r in parsed.records:
        if r.is_nocall:
            continue
        hits = pos_idx.get((r.chrom, r.pos))
        if not hits:
            continue
        called = {a for a in (r.allele1, r.allele2) if a}   # bases the person carries
        for ref, alt, vid in hits:
            # indel-anchor protection (cf. allelix ADR-0011): an array readout is a
            # single base, so it must only match SNV substitutions (single-base
            # ref AND alt). A single-base call must never match an anchor-base
            # indel like 'TC->T', where 'T' is just the anchor, not a real allele.
            if len(ref) != 1 or len(alt) != 1:
                continue
            if alt in called:
                geno = "/".join(sorted(a for a in (r.allele1, r.allele2) if a)) or alt
                out.append({"variant_id": vid, "genotype": geno, "platform": platform})
    return out
