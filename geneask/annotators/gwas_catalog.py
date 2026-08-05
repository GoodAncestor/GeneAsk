# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""NHGRI-EBI GWAS Catalog annotator — local mirror + per-rsID lookup.

The Catalog is a single ~735MB TSV of curated SNP-trait associations (rsID,
position, risk allele, trait + EFO, p-value, PMID). refresh() downloads the
bulk file (a zip with one TSV inside), builds a SQLite keyed by rsID, so a
person's carried rsIDs get instant trait-association lookups offline.

Data: NHGRI-EBI GWAS Catalog, CC BY 4.0. Bulk associations download:
  https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false
"""
from __future__ import annotations
import os, io, sqlite3, zipfile, urllib.request
import threading as _threading
from pathlib import Path
from biocore.providers.base import Finding, Tier, Category, ProviderStatus, Health

_URL = "https://www.ebi.ac.uk/gwas/api/search/downloads/associations/v1.0.2?split=false"
_DB_ENV = "GWAS_MIRROR_DB"
_DEFAULT_DB = "/data/gwas/gwas_mirror.db"

# column indices in the "alt-full" associations TSV (validated 2026-07-24, 38 cols)
_C = {"pubmedid": 1, "trait": 7, "chr": 11, "pos": 12, "risk_allele": 20,
      "snps": 21, "raf": 26, "pvalue": 27, "mapped_trait": 34, "mapped_uri": 35}


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def _risk_base(strongest: str) -> str:
    """'rs1234-A' -> 'A'; risk allele '?' or multi -> '' (unknown)."""
    if "-" in strongest:
        a = strongest.rsplit("-", 1)[1].strip()
        return a if a in ("A", "C", "G", "T") else ""
    return ""


def _tier(pvalue: float | None) -> Tier:
    """GWAS convention: p < 5e-8 is genome-wide significant (robust); 5e-8..1e-5
    suggestive (moderate); weaker or missing -> speculative."""
    if pvalue is None:
        return Tier.SPECULATIVE
    if pvalue < 5e-8:
        return Tier.ROBUST
    if pvalue < 1e-5:
        return Tier.MODERATE
    return Tier.SPECULATIVE


def build_mirror(db_path: str | None = None, workdir: str | None = None,
                 max_rows: int | None = None) -> dict:
    """Download the associations TSV and build a SQLite keyed by rsID.
    max_rows caps ingestion for validation; None ingests the full file."""
    db = _db_path(db_path)
    Path(db).parent.mkdir(parents=True, exist_ok=True)
    wd = Path(workdir or Path(db).parent)
    wd.mkdir(parents=True, exist_ok=True)
    tsv = wd / "gwas_assoc.tsv"

    if not tsv.exists() or tsv.stat().st_size == 0:
        req = urllib.request.Request(_URL, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req, timeout=600) as r:
            zf = zipfile.ZipFile(io.BytesIO(r.read()))
        member = zf.namelist()[0]
        with open(tsv, "wb") as out:
            out.write(zf.read(member))

    con = sqlite3.connect(db)
    con.execute("DROP TABLE IF EXISTS assoc")
    con.execute("""CREATE TABLE assoc(
        rsid TEXT, chrom TEXT, pos INTEGER, risk_allele TEXT, trait TEXT,
        mapped_trait TEXT, efo_uri TEXT, pvalue REAL, raf TEXT, pmid TEXT)""")
    n = 0
    with open(tsv, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        rows = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= _C["mapped_uri"]:
                continue
            snps = f[_C["snps"]].strip()
            # a row can list multiple SNPs; index each rs-id separately
            for rs in [s.strip() for s in snps.replace(";", ",").replace(" x ", ",").split(",")]:
                if not rs.startswith("rs"):
                    continue
                try:
                    pv = float(f[_C["pvalue"]]) if f[_C["pvalue"]] else None
                except ValueError:
                    pv = None
                try:
                    pos = int(f[_C["pos"]].split(";")[0]) if f[_C["pos"]] else None
                except ValueError:
                    pos = None
                rows.append((rs, f[_C["chr"]].split(";")[0], pos,
                             _risk_base(f[_C["risk_allele"]]), f[_C["trait"]],
                             f[_C["mapped_trait"]], f[_C["mapped_uri"]], pv,
                             f[_C["raf"]], f[_C["pubmedid"]]))
                n += 1
            if len(rows) >= 5000:
                con.executemany("INSERT INTO assoc VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
                rows = []
            if max_rows and n >= max_rows:
                break
        if rows:
            con.executemany("INSERT INTO assoc VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    con.execute("CREATE INDEX idx_rsid ON assoc(rsid)")
    con.commit()
    con.close()
    return {"source": "gwas_catalog", "associations": n, "db": db}


# One connection per thread, reused. Opening a fresh sqlite3.connect() per lookup
# is what made the array path quadratic in practice: a consumer export carries
# ~650k rsIDs, so the caller paid a connect/schema-load/close for each one, and
# the connection churn — not the indexed query — was the cost. Connections are
# not safe to share between threads, and the web front door now runs analysis in
# a threadpool, so the cache is thread-local rather than a module global.
_conns = _threading.local()


def _conn(db: str):
    cache = getattr(_conns, "by_path", None)
    if cache is None:
        cache = _conns.by_path = {}
    con = cache.get(db)
    if con is None:
        con = sqlite3.connect(db)
        con.row_factory = sqlite3.Row
        cache[db] = con
    return con


# SQLite's default parameter ceiling is 999; stay under it with room to spare.
_CHUNK = 900


def mirror_lookup_many(rsids, db_path: str | None = None) -> dict:
    """{rsid: [row, ...]} for many rsIDs in a handful of queries.

    The single-rsid entry point below is kept for callers with one variant in
    hand, but anything iterating a whole callset should come through here: it is
    the difference between one query per SNP and one per 900.
    """
    db = _db_path(db_path)
    if not Path(db).exists():
        return {}
    wanted = sorted({r for r in rsids if r})
    if not wanted:
        return {}
    out: dict[str, list] = {}
    con = _conn(db)
    try:
        for i in range(0, len(wanted), _CHUNK):
            chunk = wanted[i:i + _CHUNK]
            q = ("SELECT * FROM assoc WHERE rsid IN (%s)"
                 % ",".join("?" * len(chunk)))
            for r in con.execute(q, chunk):
                out.setdefault(r["rsid"], []).append(dict(r))
    except sqlite3.OperationalError:
        return {}
    return out


def mirror_lookup(rsid: str, db_path: str | None = None) -> list[dict]:
    db = _db_path(db_path)
    if not Path(db).exists():
        return []
    try:
        rows = _conn(db).execute(
            "SELECT * FROM assoc WHERE rsid=?", (rsid,)).fetchall()
    except sqlite3.OperationalError:
        return []
    return [dict(r) for r in rows]


def findings_from_rows(rsid: str, rows, carried_alleles: set | None = None) -> list:
    """Build Findings from already-fetched mirror rows, so a caller that batched
    its lookups does not have to go back to the database per rsID."""
    out = []
    for row in rows:
        trait = row.get("mapped_trait") or row.get("trait") or "trait"
        ra = row.get("risk_allele") or ""
        carries = (ra in carried_alleles) if (carried_alleles and ra) else None
        desc = f"Associated with {trait}"
        if ra:
            desc += f" (risk allele {ra}"
            if carries is not None:
                desc += f"; you carry it: {'yes' if carries else 'no'}"
            desc += ")"
        efo = (row.get("efo_uri") or "").rsplit("/", 1)[-1]
        out.append(Finding(
            marker=rsid, source="gwas_catalog", description=desc,
            tier=_tier(row.get("pvalue")), categories=[Category.TRAIT],
            detail={"p": row.get("pvalue"), "trait": trait, "efo": efo,
                    "risk_allele": ra, "raf": row.get("raf"),
                    "topic": "other", "modality": "genome"},
            link=f"https://www.ebi.ac.uk/gwas/variants/{rsid}",
            pmids=[str(row["pmid"])] if row.get("pmid") else []))
    return out


def findings_for(rsid: str, carried_alleles: set | None = None,
                 db_path: str | None = None) -> list[Finding]:
    """GWAS trait-association Findings for a carried rsID. If carried_alleles is
    given, note whether the person carries the risk allele.

    Convenience for a caller holding one variant. Anything walking a whole
    callset should use mirror_lookup_many() + findings_from_rows() instead.
    """
    return findings_from_rows(rsid, mirror_lookup(rsid, db_path), carried_alleles)
