# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""gnomAD allele-frequency enrichment — mirror-first, live API as a fallback.

Originally per-variant live GraphQL, paced and time-boxed (see the budget
machinery below), on the theory that gnomAD is terabytes and mirroring the
whole thing is the wrong design for a report that only touches the handful of
variants one person carries. Measured consequence: a report with 3,947 variant
findings got AF for 82 of them — whichever the loop reached before the 45s
budget ran out, a silently different 82 every run. That is an undisclosed
sample on a health report, not an engineering tradeoff worth keeping.

gnomad_mirror.py now builds a local SQLite of the same shape as the cache table
below (variant_id -> AF, -1.0 = known absent) from gnomAD's public sites VCFs.
When that mirror exists, allele_frequency() answers from it exclusively: no
rate limiter, no time budget, no live call, hit or miss. The live path below —
budget, rate limiter, 429 circuit breaker, disk cache — stays intact and is
still exactly what runs on a box where the mirror hasn't been built yet
(refresh:gnomad hasn't run, or hasn't finished its first chromosome), so nothing
regresses before the mirror exists. AF reframes a scary ClinVar hit ("35% of
people carry this") — an enrichment layered onto existing findings, not a
standalone finding source.

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
_MIRROR_ENV = "GNOMAD_MIRROR_DB"
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


def _mirror_path(explicit: str | None = None) -> str:
    if explicit:
        return explicit
    if os.environ.get(_MIRROR_ENV):
        return os.environ[_MIRROR_ENV]
    from geneask.annotators.gnomad_mirror import _DEFAULT_DB
    return _DEFAULT_DB


# SQLite's default parameter ceiling is 999; stay under it with room to spare
# (same limit clinvar_mirror.lookup_from_mirror works around, for the same reason).
_MIRROR_CHUNK = 900


def _mirror_lookup_many(variant_ids, mirror_db: str) -> dict:
    """{variant_id: af} for every id the mirror has an AF row for, af<0 excluded
    (that means 'looked up during the build and confirmed absent from gnomAD',
    which is indistinguishable from 'not asked about' to a caller that only
    wants a number to annotate with)."""
    wanted = sorted({v for v in variant_ids if v})
    if not wanted:
        return {}
    con = sqlite3.connect(mirror_db)
    out = {}
    try:
        for i in range(0, len(wanted), _MIRROR_CHUNK):
            chunk = wanted[i:i + _MIRROR_CHUNK]
            q = ("SELECT variant_id, af FROM af WHERE variant_id IN (%s)"
                 % ",".join("?" * len(chunk)))
            for vid, af in con.execute(q, chunk):
                if af is not None and af >= 0:
                    out[vid] = af
    except sqlite3.OperationalError:
        return {}
    finally:
        con.close()
    return out


def _to_gnomad_id(variant_id: str) -> str:
    """'chrom-pos-ref-alt' is already gnomAD's variantId shape (chr prefix stripped)."""
    return variant_id[3:] if variant_id.lower().startswith("chr") else variant_id


def allele_frequency(variant_id: str, dataset: str = "gnomad_r4",
                     cache_db: str | None = None, timeout: int = 20,
                     budget: "LiveBudget | None" = None,
                     mirror_db: str | None = None) -> float | None:
    """Population allele frequency for 'chrom-pos-ref-alt' (GRCh38).

    Mirror-first: if GNOMAD_MIRROR_DB (or mirror_db) exists on disk, the answer
    comes from it exclusively — hit or miss, no live call, no rate limiter, no
    budget spent. A mirror is built to cover the whole genome, so a miss against
    a mirror that exists means "gnomAD doesn't have this variant", not "we ran
    out of time to ask" — the two rate/budget knobs below exist only for a box
    that has no mirror yet, where they still mean exactly what they used to.

    Without a mirror, falls back to the original per-variant live GraphQL call,
    cached on disk so the second lookup of any variant is instant and offline.
    A cached None (miss) is stored as -1.0 so we don't re-hit the API for it.

    With a budget, a cache miss that has no allowance left returns None WITHOUT
    calling out and without caching — the answer is unknown to us, not known to be
    absent, and caching it would poison the variant permanently."""
    mirror = _mirror_path(mirror_db)
    if Path(mirror).exists():
        hits = _mirror_lookup_many([variant_id], mirror)
        return hits.get(variant_id)
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
                      budget: "LiveBudget | None" = None,
                      mirror_db: str | None = None):
    """Attach population AF to each finding whose marker is a 'chrom-pos-ref-alt'
    variant id, in place: sets f.detail['gnomad_af'] and appends a plain-language
    frequency note to the description. Findings whose marker isn't a variant id
    (CpG probes, rsIDs) are left untouched. Returns the count annotated.

    Mirror-first, same rule as allele_frequency(): when the mirror exists this
    is ONE chunked SQLite query for every variant marker in `findings`, and no
    rate limiter or budget applies — the number of findings a report can afford
    to annotate stops being bounded by a 45-second clock. Pass a budget to
    inspect afterwards how much of the live allowance a report used (only
    consulted on a box with no mirror); omit it for the default per-report cap."""
    mirror = _mirror_path(mirror_db)
    markers = [f.marker or "" for f in findings]
    variant_markers = [m for m in markers if len(m.split("-")) == 4 and m.split("-")[1].isdigit()]

    if Path(mirror).exists():
        hits = _mirror_lookup_many(variant_markers, mirror)
        n = 0
        for f, m in zip(findings, markers):
            af = hits.get(m)
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

    if budget is None:
        budget = default_budget(cache_db)
    variant_marker_set = set(variant_markers)
    n = 0
    for f, m in zip(findings, markers):
        if m not in variant_marker_set:
            continue    # not a chrom-pos-ref-alt variant id
        af = allele_frequency(m, cache_db=cache_db, budget=budget, mirror_db=mirror_db)
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
