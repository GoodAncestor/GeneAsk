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


def mirror_lookup(rsid: str, db_path: str | None = None) -> list[dict]:
    db = _db_path(db_path)
    if not Path(db).exists():
        return []
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    try:
        rows = con.execute("SELECT * FROM assoc WHERE rsid=?", (rsid,)).fetchall()
    except sqlite3.OperationalError:
        return []
    finally:
        con.close()
    return [dict(r) for r in rows]


def findings_for(rsid: str, carried_alleles: set | None = None,
                 db_path: str | None = None) -> list[Finding]:
    """GWAS trait-association Findings for a carried rsID. If carried_alleles is
    given, note whether the person carries the risk allele."""
    out = []
    for row in mirror_lookup(rsid, db_path):
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
