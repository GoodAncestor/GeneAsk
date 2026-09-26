# SPDX-License-Identifier: Apache-2.0
"""Bounded AlphaGenome Atlas lookups with persistent cache and local AVI support."""
from __future__ import annotations
import json
import math
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from .alphagenome_vep import _parse_vid, _is_resource_exhausted

DEFAULT_SCORERS = ("AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE")
CACHE_VERSION = "atlas-summary-v2"


def _number(name, default, minimum, maximum):
    try:
        value = float(os.getenv(name, default))
        return max(minimum, min(maximum, value)) if math.isfinite(value) else default
    except (ValueError, TypeError):
        return default


def _base(variant_id):
    return {"model": "AlphaGenome Atlas", "variant_id": variant_id,
            "assembly": "GRCh38", "status": "disabled", "tracks": [], "cache_hit": False}


def _normalize(variant_id):
    parsed = _parse_vid(variant_id)
    if not parsed:
        return None
    chrom, pos, ref, alt = parsed
    if chrom not in {"chr" + str(i) for i in range(1, 23)} | {"chrX", "chrY"}:
        return None
    if pos < 1 or len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return None
    return f"{chrom[3:]}-{pos}-{ref}-{alt}"


def _cache_path(explicit=None):
    return explicit or os.getenv("ALPHAGENOME_ATLAS_CACHE_DB", os.path.expanduser("~/.cache/geneask/alphagenome_atlas.db"))


def _cache_key(variant_id, scorers):
    return json.dumps([CACHE_VERSION, os.getenv("ALPHAGENOME_ATLAS_DATA_VERSION", "2026-09"), variant_id, sorted(scorers)])


def _cache_read(path, key):
    if not Path(path).is_file():
        return None
    try:
        with sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True, timeout=2) as con:
            row = con.execute("SELECT summary, stored FROM atlas WHERE cache_key=?", (key,)).fetchone()
        ttl = _number("ALPHAGENOME_ATLAS_CACHE_TTL_S", 30 * 86400, 0, 365 * 86400)
        if row and time.time() - row[1] <= ttl:
            result = json.loads(row[0])
            if result.get("status") == "complete":
                return {**result, "cache_hit": True, "remote_attempted": False}
    except (sqlite3.Error, ValueError, OSError):
        pass
    return None


def _cache_write(path, key, result):
    if result.get("status") != "complete":
        return
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(path, timeout=2) as con:
            con.execute("CREATE TABLE IF NOT EXISTS atlas(cache_key TEXT PRIMARY KEY, summary TEXT, stored REAL)")
            con.execute("INSERT OR REPLACE INTO atlas VALUES (?,?,?)", (key, json.dumps(result, allow_nan=False), time.time()))
    except (sqlite3.Error, ValueError, OSError):
        # Prediction remains useful when its optional cache cannot be written.
        result["cache_write_failed"] = True


def _summarize(scores, variant_id, scorers):
    bounded = []
    total = 0
    avi = None
    covered = set()
    for name, data in scores.items():
        top = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                raw = float(data.X[i, j])
                if not math.isfinite(raw):
                    continue
                total += 1
                covered.add(str(name))
                track = {"scorer": str(name), "raw_score": raw, "direction": "unknown"}
                if str(name) == "AVI_SCORE":
                    avi = raw
                quantiles = data.layers.get("quantiles")
                if quantiles is not None and math.isfinite(float(quantiles[i, j])):
                    track["quantile_score"] = float(quantiles[i, j])
                for metadata in (data.obs.iloc[i], data.var.iloc[j]):
                    for key in ("gene_id", "gene_name", "biosample_name", "ontology_curie", "name", "feature", "feature_name", "variant_scorer"):
                        value = metadata.get(key)
                        if value is not None and str(value) not in ("nan", "<NA>"):
                            track[key] = str(value)
                # Some SDK scorers encode attribution names in a labeled index.
                # Numeric row counters are not biological feature labels.
                if str(name) == "AVI_SCORE_FEATURE_IMPORTANCE" and not any(track.get(k) for k in ("name", "feature", "feature_name")):
                    for metadata, index in ((data.var, j), (data.obs, i)):
                        labels = getattr(metadata, "index", None)
                        if labels is not None:
                            label = labels[index]
                            if isinstance(label, str) and label.strip() and not label.replace(".", "", 1).isdigit():
                                track["feature_name"] = label
                                break
                top.append(track)
                # Bounded memory even with tens of thousands of tracks. Raw
                # scores are ranked only within their own scorer, never pooled.
                top.sort(key=lambda t: abs(t.get("quantile_score", t["raw_score"])), reverse=True)
                del top[5:]
        bounded.extend(top)
    missing = [s for s in scorers if s not in covered]
    result = {**_base(variant_id), "status": ("partial" if missing else "complete") if total else "not_found",
              "tracks": bounded, "n_tracks": total, "requested_scorers": scorers,
              "missing_scorers": missing,
              "provenance": "alphagenome_atlas_api", "data_version": os.getenv("ALPHAGENOME_ATLAS_DATA_VERSION", "2026-09"),
              "queried_at": datetime.now(timezone.utc).isoformat(),
              "score_explanation": "Research variant-impact score, not a personal disease probability."}
    if avi is not None:
        result["avi_score"] = avi
    return result


def _remote_query(variant_id, scorers):
    """Child-only SDK execution. Parent terminates the process at its deadline."""
    result = _base(variant_id)
    try:
        from alphagenome.atlas import atlas
        from alphagenome.data import genome
    except ImportError:
        return {**result, "status": "client_missing"}
    try:
        client = atlas.create(os.environ["ALPHA_GENOME_KEY"], timeout=10)
        chrom, pos, ref, alt = _parse_vid(variant_id)
        scores = client.query_variant(genome.Variant(chromosome=chrom, position=pos,
                                      reference_bases=ref, alternate_bases=alt), requested_scorers=scorers)
        return _summarize(scores, variant_id, scorers)
    except Exception as exc:
        return {**result, "status": "rate_limited" if _is_resource_exhausted(exc) else "failed"}


def _run_remote(variant_id, scorers, timeout_s):
    try:
        process = subprocess.run([sys.executable, "-m", "geneask.annotators.alphagenome_atlas"],
            input=json.dumps({"variant_id": variant_id, "scorers": scorers}), text=True,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout_s, check=False)
        if process.returncode:
            return {**_base(variant_id), "status": "failed"}
        result = json.loads(process.stdout)
        if not isinstance(result, dict):
            raise ValueError("Invalid worker response")
        return result
    except subprocess.TimeoutExpired:
        return {**_base(variant_id), "status": "timeout"}
    except (OSError, ValueError):
        return {**_base(variant_id), "status": "failed"}


def query_variant(variant_id: str, *, requested_scorers=None, offline=False,
                  cache_db=None, timeout_s=None, _local_result=None) -> dict:
    normalized = _normalize(variant_id)
    if normalized is None:
        return {**_base(variant_id), "status": "not_applicable"}
    variant_id = normalized
    scorers = requested_scorers if requested_scorers is not None else [s.strip() for s in
        os.getenv("ALPHAGENOME_ATLAS_SCORERS", ",".join(DEFAULT_SCORERS)).split(",") if s.strip()]
    if not scorers or len(scorers) > 20 or any(not isinstance(s, str) or len(s) > 200 for s in scorers):
        return {**_base(variant_id), "status": "scorers_required"}
    scorers = sorted(set(scorers))
    path, key = _cache_path(cache_db), _cache_key(variant_id, scorers)
    cached = _cache_read(path, key)
    if cached:
        return cached
    local = lookup_local_avi(variant_id) if _local_result is None else _local_result
    local_success = local.get("status") == "complete" and "AVI_SCORE" in scorers
    if local_success:
        local = {**local, "requested_scorers": scorers,
                 "missing_scorers": [s for s in scorers if s != "AVI_SCORE"]}
        if not local["missing_scorers"]:
            return local
    def unavailable(state):
        return ({**local, "remote_status": state} if local_success else
                {**_base(variant_id), "status": state, "local_status": local["status"]})
    if offline:
        return unavailable("offline")
    if os.getenv("ALPHAGENOME_ATLAS_ENABLED", "").lower() not in ("1", "true", "yes", "on"):
        return unavailable("disabled")
    if not os.getenv("ALPHA_GENOME_KEY"):
        return unavailable("missing_key")
    timeout = _number("ALPHAGENOME_ATLAS_TIMEOUT_S", 20, .1, 120) if timeout_s is None else max(.01, min(float(timeout_s), 120))
    result = _run_remote(variant_id, scorers, timeout)
    if result["status"] not in ("complete", "partial") and local_success:
        return {**local, "remote_status": result["status"], "remote_attempted": True}
    if result["status"] == "partial" and local_success and "AVI_SCORE" in result.get("missing_scorers", []):
        missing = [s for s in result["missing_scorers"] if s != "AVI_SCORE"]
        result = {**result, "avi_score": local["avi_score"],
                  "tracks": result.get("tracks", []) + local.get("tracks", []),
                  "n_tracks": result.get("n_tracks", 0) + local.get("n_tracks", 1),
                  "missing_scorers": missing, "status": "partial" if missing else "complete",
                  "local_avi": True, "local_source_url": local.get("source_url"),
                  "provenance": "alphagenome_atlas_api_and_local_avi"}
    result["remote_attempted"] = True
    if not result.get("local_avi"):
        _cache_write(path, key, result)
    return result


def annotate_findings(findings, *, offline=False, status=None, max_variants=None, cache_db=None):
    groups = {}
    for finding in findings:
        variant_id = _normalize(finding.marker or "")
        if variant_id:
            groups.setdefault(variant_id, []).append(finding)
    limit = int(_number("ALPHAGENOME_ATLAS_MAX_VARIANTS", 10, 1, 50)) if max_variants is None else max(0, min(50, int(max_variants)))
    deadline = time.monotonic() + _number("ALPHAGENOME_ATLAS_TIME_BUDGET_S", 25, .1, 300)
    from .atlas_avi import lookup_many
    local_status = {}
    local_records = lookup_many(groups, status=local_status)
    coverage = dict(status="complete", eligible=len(groups), scored=0, failed=0,
                    skipped=0, cache_hits=0, local_hits=0, not_found=0, partial=0,
                    remote_attempted=0, local_status=local_status.get("status"))
    reasons = {}
    enriched = 0
    halted = False
    def priority(item):
        variant_id, targets = item
        details = [finding.detail or {} for finding in targets]
        if any(d.get("research_candidate") is True for d in details):
            rank = 0
        elif any(d.get("novel_candidate") is True or
                 "uncertain" in str(d.get("clinical_significance", "")).lower() or
                 "conflicting" in str(d.get("clinical_significance", "")).lower()
                 for d in details):
            rank = 1
        else:
            rank = 2
        avi = local_records.get(variant_id, {}).get("avi_score", -math.inf)
        return rank, -avi, variant_id
    for variant_id, targets in sorted(groups.items(), key=priority):
        remaining = deadline - time.monotonic()
        permit_remote = not offline and remaining > 0 and coverage["remote_attempted"] < limit and not halted
        local = local_records.get(variant_id) or {"status": "not_found" if local_status.get("available") else local_status.get("status", "unavailable")}
        result = query_variant(variant_id, offline=not permit_remote, cache_db=cache_db,
            timeout_s=max(.01, min(remaining, _number("ALPHAGENOME_ATLAS_TIMEOUT_S", 20, .1, 120))), _local_result=local)
        state = result["status"]
        reasons[state] = reasons.get(state, 0) + 1
        attempted = bool(result.get("remote_attempted")) and not result.get("cache_hit")
        coverage["remote_attempted"] += int(attempted)
        remote_state = result.get("remote_status", state)
        if remote_state == "rate_limited":
            halted = True
        if state in ("complete", "partial"):
            coverage["scored"] += 1
            coverage["cache_hits"] += int(result.get("cache_hit", False))
            coverage["local_hits"] += int(result.get("provenance") == "alphagenome_atlas_local_avi" or result.get("local_avi", False))
            coverage["partial"] += int(state == "partial" or bool(result.get("missing_scorers")))
            if attempted and remote_state in ("failed", "timeout", "rate_limited"):
                coverage["failed"] += 1
            for finding in targets:
                if finding.detail is None:
                    finding.detail = {}
                finding.detail["alphagenome_atlas"] = result
                enriched += 1
        elif state == "not_found":
            coverage["not_found"] += 1
        elif state in ("failed", "timeout", "rate_limited"):
            coverage["failed"] += 1
        else:
            coverage["skipped"] += 1
    if halted:
        coverage["status"] = "rate_limited"
    elif not coverage["scored"] and len(reasons) == 1:
        coverage["status"] = next(iter(reasons))
    elif coverage["failed"] or coverage["skipped"] or coverage["partial"]:
        coverage["status"] = "partial"
    coverage["reasons"] = reasons
    if status is not None:
        status.update(coverage)
    return enriched


def lookup_local_avi(variant_id, path=None):
    """Local-only indexed score lookup. Format is specified in the companion module."""
    from .atlas_avi import lookup
    return lookup(variant_id, path=path)


if __name__ == "__main__":
    # SDK logging must not corrupt the one JSON result returned to the parent.
    import contextlib
    payload = json.load(sys.stdin)
    with contextlib.redirect_stdout(sys.stderr):
        output = _remote_query(payload["variant_id"], payload["scorers"])
    print(json.dumps(output, allow_nan=False))
