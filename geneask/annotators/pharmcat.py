# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2026 GoodAncestor
"""Call PharmCAT for sequencing VCFs and parse its phenotype JSON.

PharmCAT runs inside report scratch. The caller keeps only diplotype facts in
the report and lets the worker delete every PharmCAT intermediate with scratch.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


PHARMCAT_VERSION = "3.4.0"
_SEQUENCING_PLATFORMS = {"WGS", "WES"}

# Shape of PharmCAT 3.4.0's `<sample>.phenotype.json`, read from a real run on
# alien02 on 2026-09-02 (`matcherMetadata.namedAlleleMatcherVersion` 2.0.0):
#   {"matcherMetadata": {...},
#    "geneReports": {"CYP2C19": {"geneSymbol": "CYP2C19", "callSource": "MATCHER",
#                                "recommendationDiplotypes": [{"label": "*1/*2",
#                                    "phenotypes": ["Intermediate Metabolizer"],
#                                    "activityScore": 1.0, ...}], ...}, ...},
#    "unannotatedGeneCalls": [...]}
# A gene the matcher could not call carries label "Unknown/Unknown" and
# phenotypes ["No Result"]; that is "not called", never a result.
_UNCALLED_LABELS = ("unknown",)
_UNCALLED_PHENOTYPES = ("no result", "indeterminate", "")


class DiplotypeCalls(dict):
    """A gene-keyed call mapping with a non-personal pipeline note."""

    def __init__(self, *args, note: str = "", **kwargs):
        super().__init__(*args, **kwargs)
        self.note = note


def available() -> bool:
    return shutil.which("pharmcat_pipeline") is not None


def platform_ok(platform) -> bool:
    return str(platform or "").upper() in _SEQUENCING_PLATFORMS


def _activity(value):
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def _called(diplotype: dict) -> dict | None:
    """One recommendation diplotype as a call, or None when it is not a call."""
    label = str(diplotype.get("label") or "").strip()
    phenotypes = [str(x).strip() for x in (diplotype.get("phenotypes") or []) if str(x).strip()]
    phenotype = phenotypes[0] if phenotypes else ""
    if not label or any(u in label.lower() for u in _UNCALLED_LABELS):
        return None
    if phenotype.lower() in _UNCALLED_PHENOTYPES:
        return None
    return {"diplotype": label, "phenotype": phenotype,
            "activity_score": _activity(diplotype.get("activityScore"))}


def _parse(document) -> DiplotypeCalls:
    calls = DiplotypeCalls()
    reports = (document or {}).get("geneReports") if isinstance(document, dict) else None
    if not isinstance(reports, dict):
        return calls
    for gene, entry in reports.items():
        if not isinstance(entry, dict):
            continue
        symbol = str(entry.get("geneSymbol") or gene).upper()
        for diplotype in entry.get("recommendationDiplotypes") or []:
            if not isinstance(diplotype, dict):
                continue
            call = _called(diplotype)
            if call:
                calls[symbol] = {**call, "source": f"PharmCAT {PHARMCAT_VERSION}"}
                break
    return calls


def _json_candidates(scratch: str) -> list[Path]:
    paths = list(Path(scratch).rglob("*.json"))
    return sorted(
        paths,
        key=lambda path: (
            "phenotype" not in path.name.lower(),
            "report" not in path.name.lower(),
            str(path),
        ),
    )


def call_diplotypes(
    vcf_path: str, scratch: str, *, timeout: float = 600
) -> DiplotypeCalls:
    """Run PharmCAT once and return complete gene diplotype and phenotype calls."""
    command = [
        "pharmcat_pipeline",
        vcf_path,
        "-o",
        scratch,
        "-reporterJson",
    ]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        seconds = f"{timeout:g}"
        return DiplotypeCalls(
            note=(
                "PharmCAT does not produce diplotypes before the "
                f"{seconds}-second limit."
            )
        )
    except OSError as error:
        return DiplotypeCalls(
            note=f"PharmCAT does not run on this worker ({type(error).__name__})."
        )

    if completed.returncode != 0:
        return DiplotypeCalls(
            note=(
                "PharmCAT does not produce diplotypes because the pipeline returns "
                f"status {completed.returncode}."
            )
        )

    found_json = False
    for path in _json_candidates(scratch):
        found_json = True
        try:
            calls = _parse(json.loads(path.read_text()))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if calls:
            return calls
    if found_json:
        return DiplotypeCalls(
            note="PharmCAT records no complete diplotype and phenotype pair."
        )
    return DiplotypeCalls(note="PharmCAT does not produce a phenotype JSON file.")
