"""AlphaGenome regulatory variant-effect enrichment — key-gated, cached, capped.

Every other source we integrate is a LOOKUP (ClinVar/GWAS/gnomAD) — it can only
speak to variants someone already catalogued. AlphaGenome PREDICTS the regulatory
effect of a variant from DNA sequence, so it can say something about the
non-coding / uncertain-significance variants the catalogues can't. That's where
its cancer relevance sits (regulatory variants near oncogenes, e.g. TAL1).

Design (deliberately constrained — the key is a rate-limited non-commercial key
on a possibly-public site):
  - key read SERVER-SIDE only from ALPHA_GENOME_KEY; never sent to a client,
    never bundled. No key -> every function is a no-op. This is also the license
    boundary: the non-commercial constraint travels with the KEY, not this code.
  - opt-in: does nothing unless ALPHAGENOME_ENABLED is truthy AND a key is set.
  - per-report cap (ALPHAGENOME_MAX_VARIANTS, default 10): only the top-N
    uncertain variants are ever scored, so one report can't drain the quota.
  - disk cache keyed by variant: a variant is scored at most once, ever.
Enrichment only — it annotates existing findings, it is not a finding source.
"""
from __future__ import annotations
import os, json, sqlite3
from pathlib import Path

_KEY_ENV = "ALPHA_GENOME_KEY"
_ENABLED_ENV = "ALPHAGENOME_ENABLED"
_CACHE_ENV = "ALPHAGENOME_CACHE_DB"
_MAX_ENV = "ALPHAGENOME_MAX_VARIANTS"
_DEFAULT_CACHE = os.path.expanduser("~/.cache/geneask/alphagenome.db")
_DEFAULT_MAX = 10
# interval width centered on the variant; 100KB is enough for local regulatory
# context and far cheaper than the 1MB max.
_SEQ_LEN = 131072


def _enabled() -> bool:
    return bool(os.environ.get(_KEY_ENV)) and \
        os.environ.get(_ENABLED_ENV, "").lower() in ("1", "true", "yes", "on")


def _cache_con(path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE IF NOT EXISTS ag(variant_id TEXT PRIMARY KEY, summary TEXT)")
    return con


def _cache_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_CACHE_ENV) or _DEFAULT_CACHE


def _parse_vid(variant_id: str):
    """'chrom-pos-ref-alt' -> (chrom, pos, ref, alt) or None for non-SNV/indel ids."""
    p = variant_id.split("-")
    if len(p) != 4 or not p[1].isdigit():
        return None
    chrom, pos, ref, alt = p
    if not (ref and alt):
        return None
    return (chrom if chrom.startswith("chr") else "chr" + chrom, int(pos), ref, alt)


def _client(api_key: str):
    from alphagenome.models import dna_client
    return dna_client.create(api_key)


def _recommended_scorers():
    from alphagenome.models import variant_scorers
    # RECOMMENDED_VARIANT_SCORERS is a dict of named scorers; take the set the
    # library ships (bounded by MAX_VARIANT_SCORERS_PER_REQUEST=20).
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    return scorers[:20]


def score_variant(variant_id: str, api_key: str, cache_db: str | None = None) -> dict | None:
    """Score one variant's regulatory effect. Returns a small summary dict
    {variant_id, top_modality, max_abs_score, n_tracks} or None. Cached on disk."""
    parsed = _parse_vid(variant_id)
    if parsed is None:
        return None
    cache = _cache_path(cache_db)
    con = _cache_con(cache)
    try:
        row = con.execute("SELECT summary FROM ag WHERE variant_id=?", (variant_id,)).fetchone()
        if row is not None:
            return json.loads(row[0]) if row[0] else None
        chrom, pos, ref, alt = parsed
        summary = None
        try:
            from alphagenome.data import genome
            from alphagenome.models import variant_scorers as vs
            client = _client(api_key)
            variant = genome.Variant(chromosome=chrom, position=pos,
                                     reference_bases=ref, alternate_bases=alt)
            half = _SEQ_LEN // 2
            interval = genome.Interval(chromosome=chrom, start=max(0, pos - half), end=pos + half)
            scores = client.score_variant(interval=interval, variant=variant,
                                          variant_scorers=_recommended_scorers())
            df = vs.tidy_scores(scores)
            if df is not None and len(df):
                # summarize using quantile_score (normalized -1..+1: the effect's
                # percentile vs genome background — interpretable), NOT raw_score
                # (unnormalized model output, up to ~1e5, meaningless to a reader).
                col = "quantile_score" if "quantile_score" in df.columns else None
                if col:
                    df = df.assign(_abs=df[col].abs())
                    top = df.loc[df["_abs"].idxmax()]
                    modcol = "output_type" if "output_type" in df.columns else "variant_scorer"
                    q = float(top[col])
                    summary = {"variant_id": variant_id,
                               "top_modality": str(top.get(modcol, "regulatory")),
                               "quantile_score": round(q, 3),
                               "direction": "increase" if q > 0 else "decrease",
                               "n_tracks": int(len(df))}
        except Exception:
            return None    # API/library error: don't cache, allow a later retry
        con.execute("INSERT OR REPLACE INTO ag VALUES (?,?)",
                    (variant_id, json.dumps(summary) if summary else ""))
        con.commit()
        return summary
    finally:
        con.close()


def _is_uncertain(f) -> bool:
    """A variant finding the catalogues couldn't resolve — the case AlphaGenome
    adds value for: ClinVar 'uncertain significance' or 'conflicting'."""
    sig = str((f.detail or {}).get("clinical_significance", "")).lower()
    return ("uncertain" in sig) or ("conflicting" in sig)


def annotate_findings(findings, cache_db: str | None = None) -> int:
    """Score the top-N uncertain variant findings and attach the predicted
    regulatory effect in place. No-op unless enabled + key present. Returns count."""
    if not _enabled():
        return 0
    api_key = os.environ.get(_KEY_ENV)
    try:
        max_n = int(os.environ.get(_MAX_ENV, _DEFAULT_MAX))
    except ValueError:
        max_n = _DEFAULT_MAX
    # candidates: variant-id markers that are uncertain/non-catalogued
    cands = [f for f in findings
             if _parse_vid(f.marker or "") is not None and _is_uncertain(f)][:max_n]
    n = 0
    for f in cands:
        s = score_variant(f.marker, api_key, cache_db=cache_db)
        if not s:
            continue
        if f.detail is None:
            f.detail = {}
        f.detail["alphagenome"] = s
        q = abs(s.get("quantile_score", 0))
        strength = "strong" if q >= 0.9 else ("moderate" if q >= 0.5 else "weak")
        f.description = (f"{f.description} — AlphaGenome predicts a {strength} "
                         f"regulatory effect (predicted {s['direction']} in "
                         f"{s['top_modality']}, quantile {s['quantile_score']:+.2f})")
        n += 1
    return n
