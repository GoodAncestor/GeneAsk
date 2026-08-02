# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""gnomAD allele-frequency enrichment — per-variant, cached (no bulk mirror).

gnomAD is terabytes; mirroring it is the wrong design for a consumer report that
only touches the handful of variants a person actually carries. Instead we query
the gnomAD GraphQL API per variant and cache the answer on disk, so the second
lookup of any variant is instant and offline. AF reframes a scary ClinVar hit
("35% of people carry this") — an enrichment layered onto existing findings,
not a standalone finding source.

Data: gnomAD (Broad Institute), https://gnomad.broadinstitute.org/api
"""
from __future__ import annotations
import os, json, sqlite3, urllib.request, urllib.error
from pathlib import Path

_API = "https://gnomad.broadinstitute.org/api"
_CACHE_ENV = "GNOMAD_CACHE_DB"
_DEFAULT_CACHE = os.path.expanduser("~/.cache/geneask/gnomad_af.db")

_QUERY = """query($vid:String!,$ds:DatasetId!){
  variant(variantId:$vid, dataset:$ds){ genome{af} exome{af} }
}"""


def _cache_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_CACHE_ENV) or _DEFAULT_CACHE


def _cache_con(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE IF NOT EXISTS af(variant_id TEXT PRIMARY KEY, af REAL)")
    return con


def _to_gnomad_id(variant_id: str) -> str:
    """'chrom-pos-ref-alt' is already gnomAD's variantId shape (chr prefix stripped)."""
    return variant_id[3:] if variant_id.lower().startswith("chr") else variant_id


def allele_frequency(variant_id: str, dataset: str = "gnomad_r4",
                     cache_db: str | None = None, timeout: int = 20) -> float | None:
    """Population allele frequency for 'chrom-pos-ref-alt' (GRCh38). Cached on disk;
    returns the max of genome/exome AF, or None if unknown / API unreachable.
    A cached None (miss) is stored as -1.0 so we don't re-hit the API for it."""
    cache = _cache_path(cache_db)
    con = _cache_con(cache)
    try:
        row = con.execute("SELECT af FROM af WHERE variant_id=?", (variant_id,)).fetchone()
        if row is not None:
            return None if row[0] < 0 else row[0]
        af = None
        try:
            body = json.dumps({"query": _QUERY,
                               "variables": {"vid": _to_gnomad_id(variant_id), "ds": dataset}}).encode()
            req = urllib.request.Request(_API, data=body,
                                         headers={"Content-Type": "application/json",
                                                  "User-Agent": "curl/8"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                v = (json.loads(r.read()).get("data") or {}).get("variant")
            if v:
                afs = [x["af"] for x in (v.get("genome"), v.get("exome"))
                       if x and x.get("af") is not None]
                af = max(afs) if afs else None
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError):
            return None    # API down: don't cache, let a later run try again
        con.execute("INSERT OR REPLACE INTO af VALUES (?,?)",
                    (variant_id, af if af is not None else -1.0))
        con.commit()
        return af
    finally:
        con.close()


def annotate_findings(findings, cache_db: str | None = None):
    """Attach population AF to each finding whose marker is a 'chrom-pos-ref-alt'
    variant id, in place: sets f.detail['gnomad_af'] and appends a plain-language
    frequency note to the description. Findings whose marker isn't a variant id
    (CpG probes, rsIDs) are left untouched. Returns the count annotated."""
    n = 0
    for f in findings:
        m = f.marker or ""
        parts = m.split("-")
        if len(parts) != 4 or not parts[1].isdigit():
            continue    # not a chrom-pos-ref-alt variant id
        af = allele_frequency(m, cache_db=cache_db)
        if af is None:
            continue
        if f.detail is None:
            f.detail = {}
        f.detail["gnomad_af"] = af
        pct = af * 100
        freq = (f"{pct:.1f}% of people" if pct >= 0.1
                else f"~{pct:.3f}% of people (rare)")
        f.description = f"{f.description} — carried by {freq} (gnomAD)"
        n += 1
    return n
