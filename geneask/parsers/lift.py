# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""GRCh37 -> GRCh38 coordinate liftover for consumer-genotype uploads.

The ClinVar panel is GRCh38, but most consumer arrays (AncestryDNA, FTDNA,
MyHeritage, Living DNA, and many 23andMe exports) are GRCh37. A GRCh37 position
does not line up with a GRCh38 panel coordinate, so without liftover a build-37
upload silently matches nothing. This lifts each (chrom, pos) 37->38 using the
UCSC hg19ToHg38 chain, bundled offline as package data.
"""
from __future__ import annotations
from pathlib import Path
from functools import lru_cache

_CHAIN = Path(__file__).resolve().parents[1] / "data" / "reference" / "chains" / "hg19ToHg38.over.chain.gz"


@lru_cache(maxsize=1)
def _lo():
    """Lazily construct the LiftOver from the bundled chain. Returns None if
    pyliftover isn't installed or the chain is missing (caller degrades)."""
    try:
        from pyliftover import LiftOver
    except ImportError:
        return None
    if not _CHAIN.exists():
        return None
    return LiftOver(str(_CHAIN))


def lift_37_to_38(chrom: str, pos: int):
    """Lift a 1-based GRCh37 (chrom, pos) to GRCh38. Returns (chrom, pos) with no
    'chr' prefix, or None if unliftable / liftover unavailable. pyliftover uses
    0-based coordinates and 'chr'-prefixed contigs."""
    lo = _lo()
    if lo is None:
        return None
    c = chrom if chrom.startswith("chr") else "chr" + chrom
    res = lo.convert_coordinate(c, pos - 1)   # 0-based in
    if not res:
        return None
    new_chrom, new_pos0 = res[0][0], res[0][1]
    nc = new_chrom[3:] if new_chrom.startswith("chr") else new_chrom
    return (nc, new_pos0 + 1)                  # back to 1-based
