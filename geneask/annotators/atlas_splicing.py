# SPDX-License-Identifier: Apache-2.0
"""Local, indexed AlphaGenome merged-splicing scores (non-commercial only).

Schema verified from the official archive on 2026-09-27
(`combined_alphagenome_splicing_snvs.tsv.gz` in `combined_splicing_snvs_tabix.zip`):
`#CHROM POS REF ALT alphagenome_splicing`, GRCh38, `chr`-prefixed contigs.
The score is the Atlas paper's merged splicing score (sites + usage +
junctions/5, each a maximum over coordinates and genes): 0 means no predicted
change; there is no published clinical cutoff. Any other header fails closed.

Commercial output mode returns nothing before any file is opened.
"""
from __future__ import annotations
from collections import defaultdict
import math
import os
from pathlib import Path
from biocore.licensing import commercial_mode, prediction_license

SOURCE_URL = "https://deepmind.google.com/science/alphagenome/_/download/atlas/combined_splicing_snvs_tabix.zip"
HEADER = ["CHROM", "POS", "REF", "ALT", "alphagenome_splicing"]
EXPLANATION = ("Predicted change in how this gene's RNA is spliced (AlphaGenome merged splicing score). "
               "0 means no predicted change; higher means a larger predicted change. "
               "No clinical cutoff is published; this is research evidence, not a diagnosis.")


def _file(path=None):
    raw = path or os.getenv("ALPHAGENOME_ATLAS_SPLICING_FILE")
    return Path(raw) if raw and Path(raw).is_file() else None


def _header_ok(header) -> bool:
    return bool(header) and header[-1].lstrip("#").strip().split("\t") == HEADER


def local_status(path=None):
    result = {"status": "unavailable", "available": False, "assembly": "GRCh38", "source_url": SOURCE_URL}
    if commercial_mode():
        return {**result, "status": "license_blocked", "output_mode": "commercial"}
    file = _file(path)
    if file is None:
        return result
    try:
        import pysam
    except ImportError:
        return {**result, "status": "client_missing"}
    if not (Path(str(file) + ".tbi").is_file() or Path(str(file) + ".csi").is_file()):
        return {**result, "status": "index_missing"}
    try:
        with pysam.TabixFile(str(file)) as tabix:
            if not _header_ok(list(tabix.header)):
                return {**result, "status": "unsupported_schema"}
    except (OSError, ValueError):
        return {**result, "status": "invalid"}
    return {**result, "status": "ready", "available": True}


def lookup_many(variant_ids, path=None, *, status=None):
    """Exact allele matches only, grouped into 64 kb windows like the AVI reader."""
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
        chrom, pos, _, _ = vid.split("-")
        windows[(chrom, (int(pos) - 1) // 65536)][vid] = int(pos)
    out = {}
    try:
        with pysam.TabixFile(str(_file(path))) as tabix:
            contigs = set(tabix.contigs)
            for (chrom, _), entries in sorted(windows.items()):
                contig = chrom if chrom in contigs else "chr" + chrom
                if contig not in contigs:
                    continue
                for line in tabix.fetch(contig, min(entries.values()) - 1, max(entries.values())):
                    row = line.split("\t")
                    if len(row) < len(HEADER):
                        raise ValueError("Unsupported short splicing row")
                    vid = "-".join((row[0].removeprefix("chr"), row[1], row[2].upper(), row[3].upper()))
                    if vid not in entries:
                        continue
                    score = float(row[4])
                    if not math.isfinite(score):
                        continue
                    out[vid] = {"model": "AlphaGenome Atlas merged splicing", "variant_id": vid,
                                "assembly": "GRCh38", "status": "complete", "splicing_score": score,
                                "provenance": "alphagenome_atlas_local_splicing", "source_url": SOURCE_URL,
                                "license": prediction_license("alphagenome_atlas_splicing"),
                                "data_version": os.getenv("ALPHAGENOME_ATLAS_DATA_VERSION", "2026-09"),
                                "score_explanation": EXPLANATION}
    except (OSError, ValueError, IndexError):
        coverage.update(status="invalid", available=False)
        out = {}
    coverage["scored"] = len(out)
    coverage["not_found"] = len(wanted) - len(out) if coverage["available"] else 0
    if status is not None:
        status.update(coverage)
    return out


def annotate_findings(findings, *, status=None):
    """Attach `detail['alphagenome_atlas_splicing']`; the output policy strips it in commercial mode."""
    from .alphagenome_atlas import _normalize
    targets = defaultdict(list)
    for finding in findings:
        vid = _normalize(finding.marker or "")
        if vid:
            targets[vid].append(finding)
    scores = lookup_many(list(targets), status=status)
    for vid, record in scores.items():
        for finding in targets[vid]:
            if finding.detail is None:
                finding.detail = {}
            finding.detail["alphagenome_atlas_splicing"] = record
    return len(scores)
