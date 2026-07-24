"""Genome-build auto-detection from position data (cf. allelix ADR-0021).

Consumer arrays report an rsID and a position, but the header's build claim is
sometimes wrong (MyHappyGenes labels itself 37.1 but ships GRCh38). Rather than
trust the header, we vote: for a small set of universally-genotyped marker SNPs
whose GRCh37 and GRCh38 positions differ (bundled, fetched from Ensembl), check
which build's coordinate the file reports for each. The majority wins.

Returns '37', '38', or None (too few markers seen to decide) — the caller keeps
its header value or falls back to try-both when None.
"""
from __future__ import annotations
import json
from pathlib import Path
from functools import lru_cache

_ANCHORS = Path(__file__).resolve().parents[1] / "data" / "reference" / "build_anchors.json"


@lru_cache(maxsize=1)
def _anchors() -> dict:
    """{rsid: {'chrom', 'pos37', 'pos38'}} — marker SNPs with build-distinct pos."""
    try:
        return {a["rsid"]: a for a in json.load(open(_ANCHORS))}
    except Exception:
        return {}


def detect_build(records, min_votes: int = 2) -> str | None:
    """Vote build from parsed GenotypeRecords by matching marker-SNP positions.
    records: iterable of objects with .rsid and .pos. Needs >= min_votes markers
    seen to return a call; otherwise None (undecidable)."""
    anch = _anchors()
    if not anch:
        return None
    v37 = v38 = 0
    for r in records:
        a = anch.get(r.rsid)
        if a is None:
            continue
        if r.pos == a["pos38"]:
            v38 += 1
        elif r.pos == a["pos37"]:
            v37 += 1
    if v37 + v38 < min_votes:
        return None
    if v38 > v37:
        return "38"
    if v37 > v38:
        return "37"
    return None   # tie -> undecidable
