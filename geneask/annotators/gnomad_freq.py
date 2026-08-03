# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""gnomAD allele-frequency enrichment — per-variant, cached (no bulk mirror).

gnomAD is terabytes; mirroring it is the wrong design for a consumer report that
only touches the handful of variants a person actually carries. Instead we query
the gnomAD GraphQL API per variant and cache the answer on disk, so the second
lookup of any variant is instant and offline. AF reframes a scary ClinVar hit
("35% of people carry this") — an enrichment layered onto existing findings,
not a standalone finding source.

**This is somebody else's rate-limited API, called from a public endpoint.** One
report annotates every variant finding it has, and a genome full of rare variants
misses the cache on nearly all of them, so the live-call count per report is
bounded by the report, not by us. Two limits therefore apply, both per report:

  - GNOMAD_MAX_LOOKUPS (default 50) caps live calls. Cache hits are free and
    don't count — the cap is on what we ask Broad for, not on what we answer.
  - A 429 stops the rest of the run outright rather than continuing into a
    limiter that is already saying no. Before this a 429 looked exactly like
    "variant unknown", was not cached (correctly), and so every later report
    retried the same variants into the same wall.

Findings past the cap keep their cached AF where there is one and simply go
un-annotated otherwise; nothing fails, and status is reported by the caller.

Data: gnomAD (Broad Institute), https://gnomad.broadinstitute.org/api
"""
from __future__ import annotations
import os, json, sqlite3, urllib.request, urllib.error
from pathlib import Path

_API = "https://gnomad.broadinstitute.org/api"
_CACHE_ENV = "GNOMAD_CACHE_DB"
_MAX_ENV = "GNOMAD_MAX_LOOKUPS"
_DEFAULT_MAX = 50
_DEFAULT_CACHE = os.path.expanduser("~/.cache/geneask/gnomad_af.db")


class LiveBudget:
    """How many live API calls this report may still make. Shared by every lookup
    in one run; `halted` is the 429 circuit-breaker, which is separate from simple
    exhaustion so the caller can tell 'we stopped asking' from 'they stopped
    answering'."""

    def __init__(self, limit: int):
        self.remaining = max(0, limit)
        self.halted = False
        self.spent = 0

    def allows(self) -> bool:
        return self.remaining > 0 and not self.halted

    def spend(self) -> None:
        self.remaining -= 1
        self.spent += 1

    def halt(self) -> None:
        self.halted = True
        self.remaining = 0

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
                     cache_db: str | None = None, timeout: int = 20,
                     budget: "LiveBudget | None" = None) -> float | None:
    """Population allele frequency for 'chrom-pos-ref-alt' (GRCh38). Cached on disk;
    returns the max of genome/exome AF, or None if unknown / API unreachable.
    A cached None (miss) is stored as -1.0 so we don't re-hit the API for it.

    With a budget, a cache miss that has no allowance left returns None WITHOUT
    calling out and without caching — the answer is unknown to us, not known to be
    absent, and caching it would poison the variant permanently."""
    cache = _cache_path(cache_db)
    con = _cache_con(cache)
    try:
        row = con.execute("SELECT af FROM af WHERE variant_id=?", (variant_id,)).fetchone()
        if row is not None:
            return None if row[0] < 0 else row[0]
        if budget is not None and not budget.allows():
            return None
        af = None
        try:
            body = json.dumps({"query": _QUERY,
                               "variables": {"vid": _to_gnomad_id(variant_id), "ds": dataset}}).encode()
            req = urllib.request.Request(_API, data=body,
                                         headers={"Content-Type": "application/json",
                                                  "User-Agent": "curl/8"})
            if budget is not None:
                budget.spend()      # count the attempt, not the success: a call
                                    # that times out still cost them a request
            with urllib.request.urlopen(req, timeout=timeout) as r:
                v = (json.loads(r.read()).get("data") or {}).get("variant")
            if v:
                afs = [x["af"] for x in (v.get("genome"), v.get("exome"))
                       if x and x.get("af") is not None]
                af = max(afs) if afs else None
        except urllib.error.HTTPError as e:
            if e.code == 429 and budget is not None:
                budget.halt()       # they are rate-limiting us; stop asking
            return None    # don't cache: a limiter's "no" is not "variant unknown"
        except (urllib.error.URLError, TimeoutError, ValueError):
            return None    # API down: don't cache, let a later run try again
        con.execute("INSERT OR REPLACE INTO af VALUES (?,?)",
                    (variant_id, af if af is not None else -1.0))
        con.commit()
        return af
    finally:
        con.close()


def default_budget() -> LiveBudget:
    """One report's allowance of live gnomAD calls, from GNOMAD_MAX_LOOKUPS."""
    try:
        limit = int(os.environ.get(_MAX_ENV, _DEFAULT_MAX))
    except ValueError:
        limit = _DEFAULT_MAX
    return LiveBudget(limit)


def annotate_findings(findings, cache_db: str | None = None,
                      budget: "LiveBudget | None" = None):
    """Attach population AF to each finding whose marker is a 'chrom-pos-ref-alt'
    variant id, in place: sets f.detail['gnomad_af'] and appends a plain-language
    frequency note to the description. Findings whose marker isn't a variant id
    (CpG probes, rsIDs) are left untouched. Returns the count annotated.

    Pass a budget to inspect afterwards how much of the allowance a report used,
    or whether gnomAD rate-limited it; omit it for the default per-report cap."""
    if budget is None:
        budget = default_budget()
    n = 0
    for f in findings:
        m = f.marker or ""
        parts = m.split("-")
        if len(parts) != 4 or not parts[1].isdigit():
            continue    # not a chrom-pos-ref-alt variant id
        af = allele_frequency(m, cache_db=cache_db, budget=budget)
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
