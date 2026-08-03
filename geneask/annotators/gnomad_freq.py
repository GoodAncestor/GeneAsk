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
bounded by the report, not by us.

We pace instead of capping. gnomAD publishes its limit in its own server code
(broadinstitute/gnomad-browser): 30 requests per calendar minute per client IP,
and the `variant` field we query costs 1 of the 300-per-minute cost budget. There
is no daily quota. So a report may make as many lookups as it likes provided it
makes them at that rate — a flat per-report cap would truncate a big report at a
number unrelated to anything upstream, and would do it just as hard at two reports
a day as at two thousand.

  - GNOMAD_MAX_PER_MINUTE (default 25, against their 30) is the rate, held in
    SQLite so the app's two uvicorn workers share one allowance — the limit is per
    IP, and an in-process counter would let each worker spend all of it.
  - GNOMAD_TIME_BUDGET_S (default 45) bounds how long ONE report will wait for
    that rate. This is a page-latency decision, not a quota one.
  - A 429 stops the run outright rather than continuing into a limiter that is
    already saying no. Before this a 429 looked exactly like "variant unknown",
    was not cached (correctly), and so every later report retried the same
    variants into the same wall.

Findings past the deadline keep their cached AF where there is one and simply go
un-annotated otherwise; nothing fails, and status is reported by the caller.

Data: gnomAD (Broad Institute), https://gnomad.broadinstitute.org/api
"""
from __future__ import annotations
import os, json, sqlite3, urllib.request, urllib.error
from pathlib import Path

_API = "https://gnomad.broadinstitute.org/api"
_CACHE_ENV = "GNOMAD_CACHE_DB"
_RATE_ENV = "GNOMAD_MAX_PER_MINUTE"
_TIME_ENV = "GNOMAD_TIME_BUDGET_S"
# Their enforced limit is 30/min per IP; 25 leaves room for anything else on this
# address (the sibling uvicorn worker's own SQLite view, a retry in flight).
_DEFAULT_RATE = 25
_DEFAULT_TIME_BUDGET = 45.0
_DEFAULT_CACHE = os.path.expanduser("~/.cache/geneask/gnomad_af.db")


class LiveBudget:
    """One report's allowance of live gnomAD calls — a rate and a waiting time,
    not a count.

    `halted` is the 429 circuit-breaker and is deliberately distinct from running
    out of deadline: the caller needs to tell "we stopped asking" from "they
    stopped answering", because only the second is a fact about gnomAD.

    A limit of 0 (or a None limiter) means 'never call out', which is what the
    tests and any offline caller want."""

    def __init__(self, limit: int | None = None, limiter=None, deadline=None):
        # `limit` is the legacy count form, kept because it is the clearest way to
        # express "no live calls at all" and to pin behaviour in tests.
        self.remaining = limit if limit is not None else -1   # -1 = uncounted
        self.limiter = limiter
        self.deadline = deadline
        self.halted = False
        self.spent = 0

    def allows(self) -> bool:
        """Cheap pre-check: is a live call permitted at all? Does not reserve the
        slot — acquire() does that, because the wait belongs next to the call."""
        if self.halted or self.remaining == 0:
            return False
        return not (self.deadline is not None and self.deadline.expired())

    def acquire(self) -> bool:
        """Reserve one live call, waiting for the rate if needed."""
        if not self.allows():
            return False
        if self.limiter is not None and not self.limiter.acquire(self.deadline):
            return False      # the wait would outlast this report's patience
        return True

    def spend(self) -> None:
        if self.remaining > 0:
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
        if budget is not None and not budget.acquire():
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


def default_budget(cache_db: str | None = None) -> LiveBudget:
    """One report's allowance: gnomAD's published rate, and this report's patience.

    The rate counter lives in the AF cache database because that file is already
    mounted into every process that makes these calls — which is what makes the
    allowance shared rather than per-process."""
    from biocore.net.pace import RateLimiter, Deadline
    def _num(env, default, cast):
        try:
            return cast(os.environ.get(env, default))
        except (TypeError, ValueError):
            return default
    rate = _num(_RATE_ENV, _DEFAULT_RATE, int)
    seconds = _num(_TIME_ENV, _DEFAULT_TIME_BUDGET, float)
    return LiveBudget(limiter=RateLimiter(_cache_path(cache_db), "gnomad", rate),
                      deadline=Deadline(seconds))


def annotate_findings(findings, cache_db: str | None = None,
                      budget: "LiveBudget | None" = None):
    """Attach population AF to each finding whose marker is a 'chrom-pos-ref-alt'
    variant id, in place: sets f.detail['gnomad_af'] and appends a plain-language
    frequency note to the description. Findings whose marker isn't a variant id
    (CpG probes, rsIDs) are left untouched. Returns the count annotated.

    Pass a budget to inspect afterwards how much of the allowance a report used,
    or whether gnomAD rate-limited it; omit it for the default per-report cap."""
    if budget is None:
        budget = default_budget(cache_db)
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
