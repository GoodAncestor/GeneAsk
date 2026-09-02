# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Offline ClinGen gene-validity and actionability snapshot.

ClinGen publishes these curated downloads under CC0 1.0. The builder detects
each header because the validity CSV starts with a human-readable preamble.

Sources:
  https://search.clinicalgenome.org/kb/gene-validity/download
  https://actionability.clinicalgenome.org/ac/Adult/api/summ?format=tsv
  https://actionability.clinicalgenome.org/ac/Pediatric/api/summ?format=tsv
"""
from __future__ import annotations

import csv
import os
import re
import sqlite3
import urllib.request
from datetime import date
from pathlib import Path

from biocore.providers.base import Health, ProviderStatus


VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
ADULT_URL = "https://actionability.clinicalgenome.org/ac/Adult/api/summ?format=tsv"
PEDIATRIC_URL = (
    "https://actionability.clinicalgenome.org/ac/Pediatric/api/summ?format=tsv"
)
SCHEMA_VERSION = "1"
_DEFAULT_DB = "/data/clingen/clingen_mirror.db"
_DB_ENV = "CLINGEN_MIRROR_DB"

_VALIDITY_FILE = "clingen_gene_validity.csv"
_ADULT_FILE = "clingen_actionability_adult.tsv"
_PEDIATRIC_FILE = "clingen_actionability_pediatric.tsv"

# The vendor hosts were unreachable from studio8t on 2026-09-02. These aliases
# follow the approved plan and fixture. Header detection makes a changed live
# shape fail clearly instead of skipping a fixed number of preamble lines.
_GENE_FIELDS = ("GENE SYMBOL", "GENE", "GENE(S)")
_DISEASE_FIELDS = ("DISEASE LABEL", "DISEASE", "CONDITION")
_SCORE_FIELDS = ("SCORE", "ACTIONABILITY SCORE", "TOTAL SCORE")
_REPORT_FIELDS = ("ONLINE REPORT", "REPORT URL", "REPORT")

_CLASSIFICATION_RANK = {
    "definitive": 7,
    "strong": 6,
    "moderate": 5,
    "limited": 4,
    "disputed": 3,
    "refuted": 2,
    "no known disease relationship": 1,
}


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def _normal(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(value or "").upper())


def _value(row: dict[str, str], *names: str) -> str:
    normalized = {_normal(key): value.strip() for key, value in row.items()}
    for name in names:
        value = normalized.get(_normal(name))
        if value is not None:
            return value
    return ""


def _rows(path: Path, delimiter: str, required: tuple[tuple[str, ...], ...]):
    """Yield named rows after finding a header that contains every field group."""
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = None
        for candidate in reader:
            present = {_normal(cell) for cell in candidate}
            if all(any(_normal(name) in present for name in group) for group in required):
                header = candidate
                break
        if header is None:
            raise ValueError(f"ClinGen header not found in {path.name}")
        for values in reader:
            if not any(str(value).strip() for value in values):
                continue
            padded = values + [""] * max(0, len(header) - len(values))
            yield dict(zip(header, padded))


def _download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size:
        return
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(request, timeout=600) as response, destination.open(
        "wb"
    ) as output:
        while chunk := response.read(1 << 20):
            output.write(chunk)


def _number(value: str) -> float | None:
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def _genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in re.split(r"[,;|]", value) if gene.strip()]


def build_mirror(db_path: str | None = None, workdir: str | None = None) -> dict:
    """Download ClinGen snapshots and build the schema-v1 SQLite mirror."""
    database = _db_path(db_path)
    Path(database).parent.mkdir(parents=True, exist_ok=True)
    source_dir = Path(workdir or Path(database).parent)
    source_dir.mkdir(parents=True, exist_ok=True)
    validity_path = source_dir / _VALIDITY_FILE
    adult_path = source_dir / _ADULT_FILE
    pediatric_path = source_dir / _PEDIATRIC_FILE
    _download(VALIDITY_URL, validity_path)
    _download(ADULT_URL, adult_path)
    _download(PEDIATRIC_URL, pediatric_path)

    validity_rows = []
    for row in _rows(
        validity_path,
        ",",
        ((_GENE_FIELDS), (_DISEASE_FIELDS), (("MOI",)), (("CLASSIFICATION",))),
    ):
        gene = _value(row, *_GENE_FIELDS).upper()
        disease = _value(row, *_DISEASE_FIELDS)
        if not gene or not disease:
            continue
        validity_rows.append(
            (
                gene,
                disease,
                _value(row, "DISEASE ID (MONDO)", "MONDO", "DISEASE ID"),
                _value(row, "MOI"),
                _value(row, "CLASSIFICATION"),
                _value(row, "GCEP"),
                _value(row, "CLASSIFICATION DATE", "DATE"),
                _value(row, *_REPORT_FIELDS),
            )
        )

    actionability_rows = []
    for context, path in (("adult", adult_path), ("pediatric", pediatric_path)):
        for row in _rows(
            path,
            "\t",
            ((_GENE_FIELDS), (_DISEASE_FIELDS), (_SCORE_FIELDS)),
        ):
            disease = _value(row, *_DISEASE_FIELDS)
            score = _number(_value(row, *_SCORE_FIELDS))
            for gene in _genes(_value(row, *_GENE_FIELDS)):
                actionability_rows.append(
                    (
                        gene,
                        disease,
                        context,
                        score,
                        _value(row, *_REPORT_FIELDS),
                    )
                )

    retrieved_on = date.today().isoformat()
    with sqlite3.connect(database) as connection:
        connection.execute("DROP TABLE IF EXISTS validity")
        connection.execute("DROP TABLE IF EXISTS actionability")
        connection.execute("DROP TABLE IF EXISTS meta")
        connection.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
        connection.executemany(
            "INSERT INTO meta VALUES (?, ?)",
            (("schema_version", SCHEMA_VERSION), ("retrieved_on", retrieved_on)),
        )
        connection.execute(
            """CREATE TABLE validity(
                gene TEXT, disease TEXT, mondo TEXT, moi TEXT,
                classification TEXT, gcep TEXT, date TEXT, report_url TEXT,
                UNIQUE(gene, disease, mondo, moi, classification, gcep))"""
        )
        connection.execute(
            """CREATE TABLE actionability(
                gene TEXT, disease TEXT, context TEXT, score REAL, report_url TEXT,
                UNIQUE(gene, disease, context))"""
        )
        connection.executemany(
            "INSERT OR REPLACE INTO validity VALUES (?,?,?,?,?,?,?,?)", validity_rows
        )
        connection.executemany(
            "INSERT OR REPLACE INTO actionability VALUES (?,?,?,?,?)",
            actionability_rows,
        )
        connection.execute("CREATE INDEX validity_gene ON validity(gene)")
        connection.execute(
            "CREATE INDEX actionability_gene ON actionability(gene, context)"
        )
    return {
        "source": "clingen",
        "validity": len(validity_rows),
        "actionability": len(actionability_rows),
        "retrieved_on": retrieved_on,
        "db": database,
    }


def _schema_ok(connection: sqlite3.Connection) -> bool:
    try:
        row = connection.execute(
            "SELECT value FROM meta WHERE key='schema_version'"
        ).fetchone()
    except sqlite3.OperationalError:
        return False
    return bool(row) and str(row[0]) == SCHEMA_VERSION


def lookup(gene: str, db_path: str | None = None) -> dict:
    """Return ClinGen validity and actionability rows for one gene."""
    database = _db_path(db_path)
    if not Path(database).exists():
        return {}
    with sqlite3.connect(database) as connection:
        connection.row_factory = sqlite3.Row
        if not _schema_ok(connection):
            return {}
        try:
            validity = connection.execute(
                "SELECT * FROM validity WHERE gene=?", (str(gene).upper(),)
            ).fetchall()
            actionability = connection.execute(
                "SELECT * FROM actionability WHERE gene=?", (str(gene).upper(),)
            ).fetchall()
        except sqlite3.OperationalError:
            return {}
    return {
        "validity": [dict(row) for row in validity],
        "actionability": [dict(row) for row in actionability],
    }


def _inheritance(moi: str) -> str | None:
    normalized = _normal(moi)
    if normalized in {"AD", "AUTOSOMALDOMINANT", "AUTOSOMALDOMINANTINHERITANCE"}:
        return "dominant"
    if normalized in {"AR", "AUTOSOMALRECESSIVE", "AUTOSOMALRECESSIVEINHERITANCE"}:
        return "recessive"
    if normalized in {"XL", "XLINKED", "XLINKEDINHERITANCE"}:
        return "x_linked"
    return None


def inheritance_for(
    gene: str, condition_ids: list[str], db_path: str | None = None
) -> str | None:
    """Return inheritance from a matching MONDO row, then the strongest gene row."""
    rows = lookup(gene, db_path=db_path).get("validity", [])
    wanted = {str(identifier).upper() for identifier in (condition_ids or [])}
    matching = [row for row in rows if str(row.get("mondo") or "").upper() in wanted]
    candidates = matching or sorted(
        rows,
        key=lambda row: _CLASSIFICATION_RANK.get(
            str(row.get("classification") or "").lower(), 0
        ),
        reverse=True,
    )
    for row in candidates:
        inheritance = _inheritance(row.get("moi") or "")
        if inheritance:
            return inheritance
    return None


def actionability_for(gene: str, db_path: str | None = None) -> dict | None:
    """Return the highest-scored adult actionability row for one gene."""
    rows = [
        row
        for row in lookup(gene, db_path=db_path).get("actionability", [])
        if row.get("context") == "adult"
    ]
    if not rows:
        return None
    return max(
        rows,
        key=lambda row: (
            row.get("score") is not None,
            row.get("score") if row.get("score") is not None else float("-inf"),
        ),
    )


def mirror_status(db_path: str | None = None) -> ProviderStatus:
    """Return mirror health and require a rebuild for every older schema."""
    database = _db_path(db_path)
    if not Path(database).exists():
        return ProviderStatus(
            name="clingen",
            health=Health.UNAVAILABLE,
            note="Build the ClinGen mirror before using it.",
        )
    with sqlite3.connect(database) as connection:
        if not _schema_ok(connection):
            return ProviderStatus(
                name="clingen",
                health=Health.UNAVAILABLE,
                note=f"ClinGen mirror schema v{SCHEMA_VERSION} requires a rebuild.",
            )
        try:
            count = connection.execute("SELECT COUNT(*) FROM validity").fetchone()[0]
            count += connection.execute(
                "SELECT COUNT(*) FROM actionability"
            ).fetchone()[0]
            retrieved = connection.execute(
                "SELECT value FROM meta WHERE key='retrieved_on'"
            ).fetchone()
        except sqlite3.OperationalError:
            return ProviderStatus(
                name="clingen",
                health=Health.UNAVAILABLE,
                note=f"ClinGen mirror schema v{SCHEMA_VERSION} requires a rebuild.",
            )
    return ProviderStatus(
        name="clingen",
        health=Health.OK,
        version=f"schema {SCHEMA_VERSION}",
        fetched_at=retrieved[0] if retrieved else None,
        record_count=count,
    )
