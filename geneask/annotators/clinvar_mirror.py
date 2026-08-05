# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Full ClinVar annotator — local mirror of the weekly GRCh38 VCF.

Replaces the bundled 157-gene panel with the complete ClinVar set (~4.2M
variants, measured 2026-07-24). refresh() downloads clinvar.vcf.gz (~193MB,
GRCh38) and builds a
SQLite keyed by variant_id ('chrom-pos-ref-alt' — the same key the ClinVar
screen and carried-variant extractor already use), so a full-panel screen is a
drop-in for the 157-gene one: index_by_variant_id-shaped lookups, offline.

Data: NCBI ClinVar, public domain.
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
"""
from __future__ import annotations
import os, re, gzip, sqlite3, urllib.request
from pathlib import Path

_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
_DB_ENV = "CLINVAR_MIRROR_DB"
_DEFAULT_DB = "/data/clinvar/clinvar_mirror.db"

# ClinVar review status -> gold stars (the ClinVar star model).
_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def _info(info: str) -> dict:
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def _stars(revstat: str) -> int:
    s = (revstat or "").replace("_", " ").lower()
    return _STARS.get(s, 0)


def _gene(geneinfo: str) -> str:
    # GENEINFO = 'SYMBOL:id|SYMBOL2:id2' -> first symbol
    if not geneinfo:
        return ""
    return geneinfo.split(":")[0].split("|")[0]


def build_mirror(db_path: str | None = None, workdir: str | None = None,
                 max_rows: int | None = None) -> dict:
    """Download the ClinVar GRCh38 VCF and build a SQLite keyed by variant_id.
    max_rows caps ingestion for validation; None ingests the whole file."""
    db = _db_path(db_path)
    Path(db).parent.mkdir(parents=True, exist_ok=True)
    wd = Path(workdir or Path(db).parent)
    wd.mkdir(parents=True, exist_ok=True)
    vcf = wd / "clinvar.vcf.gz"

    if not vcf.exists() or vcf.stat().st_size == 0:
        req = urllib.request.Request(_URL, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req, timeout=600) as r, open(vcf, "wb") as out:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)

    con = sqlite3.connect(db)
    con.execute("DROP TABLE IF EXISTS variants")
    con.execute("""CREATE TABLE variants(
        variant_id TEXT PRIMARY KEY, clinvar_variation_id TEXT, gene TEXT,
        clinical_significance TEXT, review_status TEXT, gold_stars INTEGER)""")
    n = 0
    rows = []
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            if alt == "." or not ref or not alt:
                continue
            info = _info(f[7])
            clnsig = (info.get("CLNSIG", "") or "").replace("_", " ")
            if not clnsig:
                continue
            revstat = info.get("CLNREVSTAT", "")
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
            rows.append((variant_id, vid, _gene(info.get("GENEINFO", "")),
                         clnsig, revstat.replace("_", " "), _stars(revstat)))
            n += 1
            if len(rows) >= 5000:
                con.executemany("INSERT OR REPLACE INTO variants VALUES (?,?,?,?,?,?)", rows)
                rows = []
            if max_rows and n >= max_rows:
                break
        if rows:
            con.executemany("INSERT OR REPLACE INTO variants VALUES (?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return {"source": "clinvar_full", "variants": n, "db": db}


def _as_index(rows) -> dict:
    return {r["variant_id"]: {
        "gene": r["gene"], "clinical_significance": r["clinical_significance"],
        "review_status": r["review_status"], "gold_stars": r["gold_stars"],
        "clinvar_variation_id": r["clinvar_variation_id"]} for r in rows}


# SQLite's default parameter ceiling is 999; stay under it with room to spare.
_CHUNK = 900


def lookup_from_mirror(variant_ids, db_path: str | None = None) -> dict | None:
    """Index only the variants asked for: {variant_id: {...}} for those present.

    variant_id is the table's PRIMARY KEY, so these are index seeks. Returns None
    — not an empty dict — when there is no mirror, so the caller can tell "no
    mirror, fall back to the bundled panel" from "mirror present, nothing
    matched", which are opposite outcomes.
    """
    db = _db_path(db_path)
    if not Path(db).exists():
        return None
    wanted = sorted({v for v in variant_ids if v})
    if not wanted:
        return {}
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    out = {}
    try:
        for i in range(0, len(wanted), _CHUNK):
            chunk = wanted[i:i + _CHUNK]
            q = ("SELECT * FROM variants WHERE variant_id IN (%s)"
                 % ",".join("?" * len(chunk)))
            out.update(_as_index(con.execute(q, chunk)))
    except sqlite3.OperationalError:
        return None
    finally:
        con.close()
    return out


def load_panel_from_mirror(db_path: str | None = None) -> dict | None:
    """The whole mirror as one index, or None if it hasn't been built.

    This materialises ~4.2M rows into a dict — several seconds of CPU and
    multiple GB of memory. Screening one person's callset needs a few thousand of
    those rows, so that path uses lookup_from_mirror() instead; this remains for
    tools that genuinely want the entire set, and for the CLI.
    """
    db = _db_path(db_path)
    if not Path(db).exists():
        return None
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute("SELECT * FROM variants").fetchall()
    except sqlite3.OperationalError:
        return None
    finally:
        con.close()
    return _as_index(rows)
