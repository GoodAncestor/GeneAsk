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


def _match_at(pos_idx: dict, chrom: str, pos: int, called: set,
              geno: str, platform: str, out: list) -> int:
    """Emit carried-variant dicts for panel SNVs at (chrom,pos) whose ALT is
    carried. Returns the number of hits emitted."""
    hits = pos_idx.get((chrom, pos))
    if not hits:
        return 0
    n = 0
    for ref, alt, vid in hits:
        # indel-anchor protection (cf. allelix ADR-0011): an array readout is a
        # single base, so it must only match SNV substitutions (single-base ref
        # AND alt). A single-base call must never match an anchor-base indel like
        # 'TC->T', where 'T' is just the anchor, not a real allele.
        if len(ref) != 1 or len(alt) != 1:
            continue
        if alt in called:
            out.append({"variant_id": vid, "genotype": geno, "platform": platform})
            n += 1
    return n


def carried_from_parse(parsed: ParseResult, panel: dict,
                       platform: str = "ARRAY") -> list[dict]:
    """Match parsed array records against the GRCh38 ClinVar panel; return
    carried-variant dicts {variant_id, genotype, platform}.

    Build handling (the panel is GRCh38):
      - build '38'      -> match coordinates as-is
      - build '37'      -> lift each coordinate 37->38, then match
      - build 'unknown' -> try BOTH as-is and lifted, keep whichever the record
                           hits (handles MyHappyGenes, whose header lies about build)
    A build-37 upload matched against GRCh38 without liftover silently returns
    nothing; that failure mode is why this is coordinate-aware, not naive.
    """
    from .lift import lift_37_to_38
    from .build_detect import detect_build
    pos_idx = _index_panel_by_pos(panel)
    out: list[dict] = []
    build = (parsed.build or "unknown").replace("GRCh", "").replace("hg", "")

    # header build claims are sometimes wrong (MyHappyGenes says 37.1, ships 38),
    # so when unknown — or always, as a cross-check — vote from marker-SNP
    # positions. A confident position-based call overrides an unknown header.
    detected = detect_build(parsed.records)
    if build not in ("37", "38") and detected:
        build = detected
        parsed.notes.append(f"genome build auto-detected as GRCh{detected} from marker positions")
    elif detected and detected != build:
        parsed.notes.append(f"header build GRCh{build} disagrees with marker positions "
                            f"(GRCh{detected}); using detected GRCh{detected}")
        build = detected
    lifted_used = unliftable = 0

    for r in parsed.records:
        if r.is_nocall:
            continue
        called = {a for a in (r.allele1, r.allele2) if a}
        geno = "/".join(sorted(called)) or (next(iter(called), ""))

        if build == "38":
            _match_at(pos_idx, r.chrom, r.pos, called, geno, platform, out)
        elif build == "37":
            lp = lift_37_to_38(r.chrom, r.pos)
            if lp is None:
                unliftable += 1
                continue
            lifted_used += 1
            _match_at(pos_idx, lp[0], lp[1], called, geno, platform, out)
        else:  # unknown: try as-is (maybe already 38) then lifted (maybe 37)
            hit = _match_at(pos_idx, r.chrom, r.pos, called, geno, platform, out)
            if not hit:
                lp = lift_37_to_38(r.chrom, r.pos)
                if lp is not None:
                    _match_at(pos_idx, lp[0], lp[1], called, geno, platform, out)

    # surface what we did, so the report can note it rather than silently 0-match
    if build == "37":
        note = f"input build GRCh37: lifted to GRCh38 ({lifted_used} sites"
        if unliftable:
            note += f", {unliftable} unliftable"
        note += ") before ClinVar matching"
        parsed.notes.append(note)
    elif build == "unknown":
        parsed.notes.append("input build unknown: matched at GRCh38 and, where that "
                            "missed, at lifted GRCh37->38 coordinates")
    return out
