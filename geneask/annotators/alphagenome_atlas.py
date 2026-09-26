# SPDX-License-Identifier: Apache-2.0
"""Optional, explicit single-SNV lookup using the official AlphaGenome Atlas API.

Disabled by default. This adapter never falls back to fresh model inference and
never treats missing coverage or an upstream failure as a benign prediction.
"""
from __future__ import annotations
import math
import os
from datetime import datetime, timezone
from .alphagenome_vep import _parse_vid, _is_resource_exhausted


def query_variant(variant_id: str, *, requested_scorers: list[str] | None = None) -> dict:
    result = {"model": "AlphaGenome Atlas", "variant_id": variant_id,
              "assembly": "GRCh38", "status": "disabled", "tracks": []}
    if os.getenv("ALPHAGENOME_ATLAS_ENABLED", "").lower() not in ("1", "true", "yes", "on"):
        return result
    if not os.getenv("ALPHA_GENOME_KEY"):
        return {**result, "status": "missing_key"}
    parsed = _parse_vid(variant_id)
    if not parsed or len(parsed[2]) != 1 or len(parsed[3]) != 1 or any(b not in "ACGT" for b in parsed[2:]) or parsed[1] < 1:
        return {**result, "status": "not_applicable"}
    scorers = requested_scorers or [s.strip() for s in os.getenv("ALPHAGENOME_ATLAS_SCORERS", "").split(",") if s.strip()]
    if not scorers or len(scorers) > 20:
        return {**result, "status": "scorers_required"}
    try:
        from alphagenome.atlas import atlas
        from alphagenome.data import genome
    except ImportError:
        return {**result, "status": "client_missing"}
    try:
        client = atlas.create(os.environ["ALPHA_GENOME_KEY"], timeout=30)
        chrom, pos, ref, alt = parsed
        scores = client.query_variant(genome.Variant(chromosome=chrom, position=pos,
                                      reference_bases=ref, alternate_bases=alt),
                                      requested_scorers=scorers)
        tracks = []
        total = 0
        for name, data in scores.items():
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    raw = float(data.X[i, j])
                    if not math.isfinite(raw):
                        continue
                    total += 1
                    track = {"scorer": str(name), "raw_score": raw, "direction": "unknown"}
                    quantiles = data.layers.get("quantiles")
                    if quantiles is not None and math.isfinite(float(quantiles[i, j])):
                        track["quantile_score"] = float(quantiles[i, j])
                    for metadata in (data.obs.iloc[i], data.var.iloc[j]):
                        for key in ("gene_id", "gene_name", "biosample_name", "ontology_curie"):
                            value = metadata.get(key)
                            if value is not None and str(value) not in ("nan", "<NA>"):
                                track[key] = str(value)
                    tracks.append(track)
        # Scores from different scorers need not share units. Preserve bounded
        # per-scorer examples instead of ranking raw AVI vs molecular effects.
        bounded = []
        for name in scores:
            group = [t for t in tracks if t["scorer"] == name]
            group.sort(key=lambda t: abs(t.get("quantile_score", 0)), reverse=True)
            bounded.extend(group[:5])
        return {**result, "status": "complete" if total else "not_found", "tracks": bounded,
                "n_tracks": total, "requested_scorers": scorers,
                "provenance": "alphagenome_atlas_api", "queried_at": datetime.now(timezone.utc).isoformat()}
    except Exception as exc:
        return {**result, "status": "rate_limited" if _is_resource_exhausted(exc) else "failed"}
