# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2026 GoodAncestor
"""Call PharmCAT for sequencing VCFs and parse its phenotype JSON.

PharmCAT runs inside report scratch. The caller keeps only diplotype facts in
the report and lets the worker delete every PharmCAT intermediate with scratch.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path


PHARMCAT_VERSION = "3.4.0"
_SEQUENCING_PLATFORMS = {"WGS", "WES"}

# PharmCAT documentation was unreachable from studio8t on 2026-09-02. These
# aliases include the fields named in the approved plan and common nested forms.
# The Task C4 real run must confirm them against an unmodified phenotype JSON.
_GENE_FIELDS = ("gene", "geneSymbol", "genesymbol")
_DIPLOTYPE_FIELDS = ("diplotype", "diplotypeName", "diplotypeString")
_PHENOTYPE_FIELDS = ("phenotype", "phenotypeName", "metabolizerStatus")
_ACTIVITY_FIELDS = ("activityScore", "activity_score")


class DiplotypeCalls(dict):
    """A gene-keyed call mapping with a non-personal pipeline note."""

    def __init__(self, *args, note: str = "", **kwargs):
        super().__init__(*args, **kwargs)
        self.note = note


def available() -> bool:
    return shutil.which("pharmcat_pipeline") is not None


def platform_ok(platform) -> bool:
    return str(platform or "").upper() in _SEQUENCING_PLATFORMS


def _normal(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


def _field(record: dict, names: tuple[str, ...]):
    wanted = {_normal(name) for name in names}
    for key, value in record.items():
        if _normal(key) in wanted:
            return value
    return None


def _text(value) -> str:
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, dict):
        for key in ("name", "value", "label", "term"):
            text = _text(value.get(key))
            if text:
                return text
    if isinstance(value, list):
        values = [_text(item) for item in value]
        values = [item for item in values if item]
        if len(values) == 2 and all(item.startswith("*") for item in values):
            return "/".join(values)
        return values[0] if values else ""
    return ""


def _activity(value):
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def _looks_like_gene(value: str) -> bool:
    return bool(re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", str(value or "")))


def _records(value, inherited_gene: str = ""):
    if isinstance(value, list):
        for item in value:
            yield from _records(item, inherited_gene)
        return
    if not isinstance(value, dict):
        return

    gene = _text(_field(value, _GENE_FIELDS)) or inherited_gene
    diplotype = _text(_field(value, _DIPLOTYPE_FIELDS))
    phenotype = _text(_field(value, _PHENOTYPE_FIELDS))
    if gene and diplotype and phenotype:
        yield {
            "gene": gene.upper(),
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": _activity(_field(value, _ACTIVITY_FIELDS)),
        }

    for key, child in value.items():
        child_gene = gene
        if not child_gene and _looks_like_gene(key) and isinstance(child, (dict, list)):
            child_gene = key
        yield from _records(child, child_gene)


def _parse(document) -> DiplotypeCalls:
    calls = DiplotypeCalls()
    for record in _records(document):
        calls[record["gene"]] = {
            "diplotype": record["diplotype"],
            "phenotype": record["phenotype"],
            "activity_score": record["activity_score"],
            "source": f"PharmCAT {PHARMCAT_VERSION}",
        }
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
