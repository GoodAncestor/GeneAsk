# SPDX-License-Identifier: Apache-2.0
"""Local, indexed AVI score access; no download or unindexed whole-file scans.

The official download is tabix-indexed but its column schema must be inspected
at installation. Headers (or ALPHAGENOME_ATLAS_AVI_COLUMNS) must unambiguously
identify chromosome, 1-based position, reference, alternate and AVI score.
Unknown formats fail explicitly, never guessing column positions or assembly.
"""
from __future__ import annotations
from collections import defaultdict
import math
import os
from pathlib import Path
import re
from biocore.licensing import prediction_license

SOURCE_URL = "https://deepmind.google.com/science/alphagenome/_/download/atlas/avi_scores_snvs_tabix.zip"
ALIASES = {
    "chrom": {"chrom", "chromosome", "chr"},
    "pos": {"pos", "position", "position1", "position1based"},
    "ref": {"ref", "reference", "referenceallele", "referencebases"},
    "alt": {"alt", "alternate", "alternative", "alternateallele", "alternatebases"},
    "avi": {"avi", "aviscore", "score", "rawscore"},
    "phred": {"phred", "phredscore", "aviphred", "aviscorephred", "phredscaledscore"},
}


def _files(path=None):
    raw = path or os.getenv("ALPHAGENOME_ATLAS_AVI_FILE")
    if not raw:
        return []
    p = Path(raw)
    if p.is_dir():
        return sorted(f for f in p.iterdir() if f.is_file() and f.suffix in (".gz", ".bgz", ".bgzf"))
    return [p] if p.is_file() else []


def _columns(header):
    configured = os.getenv("ALPHAGENOME_ATLAS_AVI_COLUMNS", "")
    candidates = [configured.split(",")] if configured else [h.lstrip("#").strip().split("\t") for h in header]
    for names in reversed(candidates):
        indices = {}
        for index, name in enumerate(names):
            cleaned = re.sub(r"[^a-z0-9]", "", name.lower())
            for canonical, aliases in ALIASES.items():
                if cleaned in aliases:
                    if canonical in indices:
                        return None
                    indices[canonical] = index
        if {"chrom", "pos", "ref", "alt", "avi"}.issubset(indices):
            return indices
    return None


def local_status(path=None):
    files = _files(path)
    result = {"status": "unavailable", "available": False, "assembly": "GRCh38", "source_url": SOURCE_URL}
    if not files:
        return result
    try:
        import pysam
    except ImportError:
        return {**result, "status": "client_missing"}
    try:
        for file in files:
            if not (Path(str(file) + ".tbi").is_file() or Path(str(file) + ".csi").is_file()):
                return {**result, "status": "index_missing"}
            with pysam.TabixFile(str(file)) as tabix:
                if not _columns(tabix.header):
                    return {**result, "status": "unsupported_schema"}
        return {**result, "status": "ready", "available": True, "files": len(files)}
    except (OSError, ValueError):
        return {**result, "status": "invalid"}


def lookup_many(variant_ids, path=None, *, status=None):
    """One handle per shard, grouped 64kb seeks; return only exact allele matches.

    Accepts a single BGZF or directory of indexed chromosome shards. Callers can
    pass a mutable status to distinguish missing installation from true misses.
    """
    from .alphagenome_atlas import _normalize
    coverage = local_status(path)
    wanted = {_normalize(v) for v in variant_ids}
    wanted.discard(None)
    coverage.update(eligible=len(wanted), scored=0, not_found=0)
    if not coverage["available"]:
        if status is not None:
            status.update(coverage)
        return {}
    import pysam
    windows = defaultdict(dict)
    for vid in wanted:
        chrom, pos, ref, alt = vid.split("-")
        pos = int(pos)
        windows[(chrom, (pos - 1) // 65536)][vid] = pos
    out = {}
    try:
        for file in _files(path):
            with pysam.TabixFile(str(file)) as tabix:
                columns = _columns(tabix.header)
                contigs = set(tabix.contigs)
                for (chrom, window), entries in sorted(windows.items()):
                    contig = chrom if chrom in contigs else "chr" + chrom
                    if contig not in contigs:
                        continue
                    for line in tabix.fetch(contig, min(entries.values()) - 1, max(entries.values())):
                        row = line.split("\t")
                        if len(row) <= max(columns.values()):
                            raise ValueError("Unsupported short AVI row")
                        c = row[columns["chrom"]].removeprefix("chr")
                        p = int(row[columns["pos"]])
                        vid = f'{c}-{p}-{row[columns["ref"]]}-{row[columns["alt"]]}'
                        if vid not in entries:
                            continue
                        score = float(row[columns["avi"]])
                        if not math.isfinite(score):
                            continue
                        value = {"model": "AlphaGenome Atlas", "variant_id": vid, "assembly": "GRCh38",
                                 "status": "complete", "avi_score": score, "cache_hit": False,
                                 "provenance": "alphagenome_atlas_local_avi", "source_url": SOURCE_URL,
                                 "license": prediction_license("alphagenome_atlas_local_avi"),
                                 "data_version": os.getenv("ALPHAGENOME_ATLAS_DATA_VERSION", "2026-09"),
                                 "tracks": [{"scorer": "AVI_SCORE", "raw_score": score, "direction": "unknown"}],
                                 "n_tracks": 1, "score_explanation": "Research variant-impact score, not a personal disease probability."}
                        if "phred" in columns:
                            phred = float(row[columns["phred"]])
                            if math.isfinite(phred):
                                value["avi_phred"] = phred
                        out[vid] = value
    except (OSError, ValueError, IndexError):
        coverage.update(status="invalid", available=False)
        # Never report partial/corrupt reads as authoritative misses.
        out = {}
    coverage["scored"] = len(out)
    coverage["not_found"] = len(wanted) - len(out) if coverage["available"] else 0
    if status is not None:
        status.update(coverage)
    return out


def lookup(variant_id, path=None):
    from .alphagenome_atlas import _normalize, _base
    normalized = _normalize(variant_id)
    if normalized is None:
        return {**_base(variant_id), "status": "not_applicable"}
    status = {}
    records = lookup_many([normalized], path, status=status)
    return records.get(normalized) or {**_base(normalized),
        "status": "not_found" if status["available"] else status["status"]}
