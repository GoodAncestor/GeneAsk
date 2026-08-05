# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""gnomAD allele-frequency mirror — local SQLite built from the public sites VCFs.

gnomad_freq.py's per-variant live GraphQL lookup was paced at 25/min under a
45s per-report budget. Measured consequence: a report with 3,947 variant
findings got AF for 82 of them, and which 82 depended only on loop order — not
disclosed, not reproducible, and silently different every run. That is the
wrong trade for a health report; every other reference source here (ClinVar,
AlphaMissense, GWAS Catalog, CPIC) already solved it the same way: mirror once,
answer forever, offline.

Unlike clinvar_mirror.py/alphamissense.py this build is NOT restartable from
zero — a full gnomAD genomes release is dozens of per-chromosome files running
into the tens of GB compressed, and a build that dies on chr17 must not
re-walk chr1-16. Progress is recorded per chromosome (`progress` table), and a
chromosome already marked done is skipped entirely on the next run — no
re-download, no re-parse. A crash mid-chromosome leaves that chromosome's rows
in `af` (harmless: re-ingestion is INSERT OR REPLACE) but not marked done, so
it simply restarts from scratch on the next run rather than corrupting
anything.

Keys are 'chrom-pos-ref-alt' with no 'chr' prefix — the same shape
biocore.variants.carried.carried_variants() emits and the same shape the `af`
cache table in gnomad_freq.py already uses, so the reader needs no translation
between mirror-hit and live-hit answers.

Data: gnomAD (Broad Institute), public GRCh38 sites VCFs.
  https://gnomad.broadinstitute.org/downloads
"""
from __future__ import annotations
import contextlib, os, gzip, sqlite3, urllib.request
from pathlib import Path

_URL_TMPL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/"
            "4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz")
_DB_ENV = "GNOMAD_MIRROR_DB"
_DEFAULT_DB = "/data/gnomad/gnomad_mirror.db"
_CHROMS_ENV = "GNOMAD_MIRROR_CHROMS"     # comma list, for a bounded/testing build
_MAX_VARIANTS_ENV = "GNOMAD_MIRROR_MAX_VARIANTS"

# Chromosome order for the (default, unbounded) full build.
ALL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def _env_chroms() -> list[str] | None:
    raw = os.environ.get(_CHROMS_ENV)
    if not raw:
        return None
    return [c.strip() for c in raw.split(",") if c.strip()]


def _env_max_variants() -> int | None:
    raw = os.environ.get(_MAX_VARIANTS_ENV)
    if not raw:
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def _ensure_schema(con: sqlite3.Connection) -> None:
    # IF NOT EXISTS, deliberately not DROP: unlike the single-file mirrors, a
    # resumed build must keep whatever earlier chromosomes already wrote.
    con.execute("CREATE TABLE IF NOT EXISTS af(variant_id TEXT PRIMARY KEY, af REAL)")
    con.execute("CREATE TABLE IF NOT EXISTS progress("
                "chrom TEXT PRIMARY KEY, variants INTEGER, done INTEGER)")


def _done_chroms(con: sqlite3.Connection) -> dict[str, int]:
    return {r[0]: r[1] for r in
           con.execute("SELECT chrom, variants FROM progress WHERE done=1")}


def _info(info: str) -> dict:
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def _download(chrom: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return       # a prior run already pulled this chromosome's file
    req = urllib.request.Request(_URL_TMPL.format(chrom=chrom), headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=1800) as r, open(dest, "wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)


def _ingest_chrom(con: sqlite3.Connection, vcf_path: Path, remaining: int | None) -> int:
    """Insert every allele in one chromosome's VCF. `remaining` (None = unbounded)
    caps how many more variants THIS BUILD will accept in total, for the bounded
    test path — a real build always passes None."""
    n = 0
    rows = []
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, ref, alt_field = f[0], f[1], f[3], f[4]
            chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
            if alt_field == "." or not ref or not alt_field:
                continue
            alts = alt_field.split(",")
            info = _info(f[7])
            # AF is Number=A: one value per ALT allele, same order as `alts`. A
            # multiallelic site with a missing/short AF list is a malformed line,
            # not ours to guess at, so those alleles are skipped rather than
            # mis-paired to the wrong AF.
            af_list = info.get("AF", "").split(",") if "AF" in info else []
            for i, alt in enumerate(alts):
                if not alt or alt == "*":
                    continue
                if i >= len(af_list):
                    continue
                try:
                    af = float(af_list[i])
                except ValueError:
                    continue
                rows.append((f"{chrom}-{pos}-{ref}-{alt}", af))
                n += 1
                if remaining is not None and n >= remaining:
                    break
            if remaining is not None and n >= remaining:
                break
            if len(rows) >= 20000:
                con.executemany("INSERT OR REPLACE INTO af VALUES (?,?)", rows)
                rows = []
    if rows:
        con.executemany("INSERT OR REPLACE INTO af VALUES (?,?)", rows)
    return n


def build_mirror(db_path: str | None = None, workdir: str | None = None,
                 chroms: list[str] | None = None,
                 max_variants: int | None = None) -> dict:
    """Stream gnomAD's per-chromosome GRCh38 sites VCFs into a SQLite keyed by
    variant_id, resuming chromosome-by-chromosome across runs.

    chroms restricts which chromosomes to build (default: all 24, or
    GNOMAD_MIRROR_CHROMS if set) — the bounded path a test or a smoke build uses.
    max_variants caps the TOTAL number of new variants this call will ingest
    (default: unbounded, or GNOMAD_MIRROR_MAX_VARIANTS if set); a chromosome that
    hits the cap partway through is NOT marked done, so the next call resumes it
    from scratch rather than serving a truncated chromosome as complete.
    """
    db = _db_path(db_path)
    Path(db).parent.mkdir(parents=True, exist_ok=True)
    wd = Path(workdir or Path(db).parent)
    wd.mkdir(parents=True, exist_ok=True)

    todo = chroms if chroms is not None else (_env_chroms() or ALL_CHROMS)
    cap = max_variants if max_variants is not None else _env_max_variants()

    con = sqlite3.connect(db)
    con.execute("PRAGMA journal_mode=OFF")     # bulk load; a crash restarts the chrom, not the row
    con.execute("PRAGMA synchronous=OFF")
    _ensure_schema(con)

    already = _done_chroms(con)
    built: list[str] = []
    skipped: list[str] = []
    total_new = 0
    for chrom in todo:
        if chrom in already:
            skipped.append(chrom)
            continue
        if cap is not None and total_new >= cap:
            break     # cap reached before this chromosome started: leave it untouched, not partial
        vcf = wd / f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
        _download(chrom, vcf)
        remaining = (cap - total_new) if cap is not None else None
        try:
            n = _ingest_chrom(con, vcf, remaining)
        finally:
            # Drop each chromosome's VCF as soon as it is ingested. They are only
            # ever read once, and keeping them is not a small waste: chr21 alone is
            # 7.76 GB compressed, so a full 24-chromosome run would leave several
            # hundred GB of source files sitting next to a mirror of a few dozen.
            # Deleting bounds the build to one chromosome file plus the database.
            # Safe for resumption either way — a chromosome is re-downloaded unless
            # it was recorded complete.
            with contextlib.suppress(OSError):
                vcf.unlink()
        total_new += n
        if cap is not None and total_new >= cap:
            # this chromosome was cut short by the cap: its rows are in `af`
            # (harmless, idempotent) but it is NOT progress-complete
            con.commit()
            break
        con.execute("INSERT OR REPLACE INTO progress VALUES (?,?,1)", (chrom, n))
        con.commit()
        built.append(chrom)

    con.close()
    return {"source": "gnomad", "variants_added": total_new,
           "chroms_built": built, "chroms_skipped": skipped, "db": db}


def mirror_available(db_path: str | None = None) -> bool:
    return Path(_db_path(db_path)).exists()
