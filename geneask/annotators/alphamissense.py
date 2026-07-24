"""AlphaMissense annotator — local mirror of precomputed missense pathogenicity.

AlphaMissense scores (nearly) every possible missense substitution in the human
proteome for pathogenicity. It fills the gap ClinVar leaves: a missense variant
nobody has clinically curated still gets a principled benign/pathogenic call.

refresh() downloads AlphaMissense_hg38.tsv.gz (~643MB gz, ~71M rows) and builds a
SQLite keyed by variant_id ('chrom-pos-ref-alt', GRCh38 — same key as ClinVar/
gnomAD), so a carried variant gets an instant offline pathogenicity lookup.

The built mirror is multi-GB (bigger than ClinVar) — build on local disk, then
copy to shared storage (see worker README).

Data: AlphaMissense (Google DeepMind), CC BY-NC-SA 4.0 (non-commercial).
  https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz
"""
from __future__ import annotations
import os, gzip, sqlite3, urllib.request
from pathlib import Path
from biocore.providers.base import Finding, Tier, Category

_URL = "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"
_DB_ENV = "ALPHAMISSENSE_MIRROR_DB"
_DEFAULT_DB = "/data/alphamissense/alphamissense_mirror.db"

# am_class -> tier. 'likely_pathogenic' is a strong computational call; 'ambiguous'
# is genuinely uncertain; 'likely_benign' is reported as benign context.
_CLASS_TIER = {"likely_pathogenic": Tier.MODERATE, "pathogenic": Tier.MODERATE,
               "ambiguous": Tier.SPECULATIVE, "likely_benign": Tier.SPECULATIVE,
               "benign": Tier.SPECULATIVE}


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def build_mirror(db_path: str | None = None, workdir: str | None = None,
                 max_rows: int | None = None) -> dict:
    """Download AlphaMissense hg38 TSV and build a SQLite keyed by variant_id.
    max_rows caps ingestion for validation; None ingests the whole file (~71M)."""
    db = _db_path(db_path)
    Path(db).parent.mkdir(parents=True, exist_ok=True)
    wd = Path(workdir or Path(db).parent)
    wd.mkdir(parents=True, exist_ok=True)
    tsv = wd / "AlphaMissense_hg38.tsv.gz"

    if not tsv.exists() or tsv.stat().st_size == 0:
        req = urllib.request.Request(_URL, headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req, timeout=1200) as r, open(tsv, "wb") as out:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)

    con = sqlite3.connect(db)
    con.execute("PRAGMA journal_mode=OFF")       # bulk load, rebuilt-from-scratch
    con.execute("PRAGMA synchronous=OFF")
    con.execute("DROP TABLE IF EXISTS am")
    con.execute("""CREATE TABLE am(variant_id TEXT PRIMARY KEY, uniprot_id TEXT,
                   protein_variant TEXT, pathogenicity REAL, am_class TEXT)""")
    n = 0
    rows = []
    with gzip.open(tsv, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[2], f[3]
            chrom = chrom[3:] if chrom.startswith("chr") else chrom
            try:
                pathog = float(f[8])
            except ValueError:
                continue
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
            rows.append((variant_id, f[5], f[7], pathog, f[9]))
            n += 1
            if len(rows) >= 20000:
                con.executemany("INSERT OR REPLACE INTO am VALUES (?,?,?,?,?)", rows)
                rows = []
            if max_rows and n >= max_rows:
                break
        if rows:
            con.executemany("INSERT OR REPLACE INTO am VALUES (?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return {"source": "alphamissense", "variants": n, "db": db}


def lookup(variant_id: str, db_path: str | None = None) -> dict | None:
    db = _db_path(db_path)
    if not Path(db).exists():
        return None
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    try:
        r = con.execute("SELECT * FROM am WHERE variant_id=?", (variant_id,)).fetchone()
    except sqlite3.OperationalError:
        return None
    finally:
        con.close()
    return dict(r) if r else None


def annotate_findings(findings, db_path: str | None = None) -> int:
    """Attach AlphaMissense pathogenicity to variant findings whose marker is a
    'chrom-pos-ref-alt' id, in place. Mirror-first: no-op if the mirror isn't
    built. Adds a note only when a variant is actually in AlphaMissense (missense)."""
    n = 0
    for f in findings:
        m = f.marker or ""
        parts = m.split("-")
        if len(parts) != 4 or not parts[1].isdigit():
            continue
        rec = lookup(m, db_path)
        if not rec:
            continue
        if f.detail is None:
            f.detail = {}
        f.detail["alphamissense"] = {"pathogenicity": rec["pathogenicity"],
                                     "class": rec["am_class"],
                                     "protein_variant": rec["protein_variant"]}
        f.description = (f"{f.description} — AlphaMissense: {rec['am_class'].replace('_',' ')} "
                         f"({rec['protein_variant']}, score {rec['pathogenicity']:.2f})")
        n += 1
    return n


def findings_for(variant_id: str, db_path: str | None = None) -> list[Finding]:
    """A standalone AlphaMissense Finding for a variant (when it isn't already a
    finding from another source), for missense variants of uncertain catalog status."""
    rec = lookup(variant_id, db_path)
    if not rec:
        return []
    tier = _CLASS_TIER.get((rec["am_class"] or "").lower(), Tier.SPECULATIVE)
    return [Finding(
        marker=variant_id, source="alphamissense",
        description=(f"AlphaMissense predicts {rec['am_class'].replace('_',' ')} "
                     f"for {rec['protein_variant']} (score {rec['pathogenicity']:.2f})"),
        tier=tier, categories=[Category.CLINICAL],
        detail={"pathogenicity": rec["pathogenicity"], "am_class": rec["am_class"],
                "protein_variant": rec["protein_variant"], "uniprot_id": rec["uniprot_id"],
                "topic": "clinical", "modality": "genome"},
        link="https://alphamissense.hegelab.org/")]
