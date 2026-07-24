"""CPIC pharmacogenomics annotator — local mirror of the CPIC knowledge base.

Pharmacogenomics = how a person's genotype affects drug response. The clinical
chain is genotype -> diplotype (gene's two alleles) -> phenotype (e.g. "CYP2C19
poor metabolizer") -> drug recommendation. CPIC publishes this as curated tables.

This mirror covers the gene->drug RECOMMENDATION layer keyed by phenotype:
  - drug            RxNorm drug id -> name
  - pair            gene-drug pairs + CPIC evidence level (A/B/...)
  - recommendation  phenotype-keyed drug guidance (dosing, alternates)
  - allele          per-allele functional status (for phenotype interpretation)

Calling a diplotype from raw genotype (phasing the allele-definition/location
tables) is a documented follow-on; this layer answers "for gene X phenotype Y,
what does CPIC say about drug Z" — the consumer-facing PGx result.

Data: CPIC (cpicpgx.org), CC0. API: https://api.cpicpgx.org/v1/
"""
from __future__ import annotations
import os, json, sqlite3, urllib.request
from pathlib import Path
from biocore.providers.base import Finding, Tier, Category

_API = "https://api.cpicpgx.org/v1/"
_DB_ENV = "PGX_MIRROR_DB"
_DEFAULT_DB = "/data/pgx/pgx_mirror.db"

# CPIC evidence level -> tier (A = strong, prescribing actionable).
_LEVEL_TIER = {"A": Tier.ROBUST, "A/B": Tier.ROBUST, "B": Tier.MODERATE,
               "B/C": Tier.MODERATE, "C": Tier.SPECULATIVE, "D": Tier.SPECULATIVE}


def _db_path(explicit: str | None = None) -> str:
    return explicit or os.environ.get(_DB_ENV) or _DEFAULT_DB


def _fetch_all(table: str, select: str, page: int = 1000) -> list[dict]:
    """Page a CPIC PostgREST table fully via Range headers."""
    out, start = [], 0
    while True:
        url = f"{_API}{table}?select={select}"
        req = urllib.request.Request(url, headers={
            "User-Agent": "curl/8", "Accept": "application/json",
            "Range-Unit": "items", "Range": f"{start}-{start + page - 1}"})
        with urllib.request.urlopen(req, timeout=60) as r:
            batch = json.loads(r.read())
        out.extend(batch)
        if len(batch) < page:
            break
        start += page
    return out


def build_mirror(db_path: str | None = None, workdir: str | None = None) -> dict:
    """Mirror the CPIC drug/pair/recommendation/allele tables into SQLite."""
    db = _db_path(db_path)
    Path(db).parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db)
    c = con.cursor()
    c.execute("DROP TABLE IF EXISTS drug")
    c.execute("CREATE TABLE drug(drugid TEXT PRIMARY KEY, name TEXT)")
    c.executemany("INSERT OR REPLACE INTO drug VALUES (?,?)",
                  [(d["drugid"], d.get("name")) for d in _fetch_all("drug", "drugid,name")])
    c.execute("DROP TABLE IF EXISTS pair")
    c.execute("CREATE TABLE pair(genesymbol TEXT, drugid TEXT, cpiclevel TEXT, guidelineid INTEGER)")
    c.executemany("INSERT INTO pair VALUES (?,?,?,?)",
                  [(p["genesymbol"], p["drugid"], p.get("cpiclevel"), p.get("guidelineid"))
                   for p in _fetch_all("pair", "genesymbol,drugid,cpiclevel,guidelineid")])
    c.execute("DROP TABLE IF EXISTS recommendation")
    c.execute("""CREATE TABLE recommendation(guidelineid INTEGER, drugid TEXT,
                 phenotypes TEXT, drugrecommendation TEXT, classification TEXT,
                 implications TEXT, lookupkey TEXT)""")
    c.executemany("INSERT INTO recommendation VALUES (?,?,?,?,?,?,?)",
                  [(r.get("guidelineid"), r.get("drugid"), json.dumps(r.get("phenotypes")),
                    r.get("drugrecommendation"), r.get("classification"),
                    json.dumps(r.get("implications")), json.dumps(r.get("lookupkey")))
                   for r in _fetch_all("recommendation",
                     "guidelineid,drugid,phenotypes,drugrecommendation,classification,implications,lookupkey")])
    c.execute("DROP TABLE IF EXISTS allele")
    c.execute("CREATE TABLE allele(genesymbol TEXT, name TEXT, functionalstatus TEXT)")
    c.executemany("INSERT INTO allele VALUES (?,?,?)",
                  [(a["genesymbol"], a.get("name"), a.get("clinicalfunctionalstatus"))
                   for a in _fetch_all("allele", "genesymbol,name,clinicalfunctionalstatus")])
    c.execute("CREATE INDEX idx_pair_gene ON pair(genesymbol)")
    c.execute("CREATE INDEX idx_rec_gid ON recommendation(guidelineid)")
    con.commit()
    counts = {t: con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
              for t in ("drug", "pair", "recommendation", "allele")}
    con.close()
    return {"source": "cpic", "db": db, **counts}


def recommendations_for_gene(genesymbol: str, phenotype: str | None = None,
                             db_path: str | None = None) -> list[Finding]:
    """CPIC drug-guidance Findings for a gene (optionally filtered to a phenotype).
    Each finding names the drug + CPIC recommendation, tiered by evidence level."""
    db = _db_path(db_path)
    if not Path(db).exists():
        return []
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    out = []
    try:
        pairs = con.execute("SELECT * FROM pair WHERE genesymbol=?", (genesymbol,)).fetchall()
        for p in pairs:
            drug = con.execute("SELECT name FROM drug WHERE drugid=?", (p["drugid"],)).fetchone()
            drug_name = drug["name"] if drug else p["drugid"]
            recs = con.execute("SELECT * FROM recommendation WHERE guidelineid=? AND drugid=?",
                               (p["guidelineid"], p["drugid"])).fetchall()
            for r in recs:
                phen = json.loads(r["phenotypes"] or "{}")
                if phenotype and phenotype not in json.dumps(phen):
                    continue
                rec_txt = r["drugrecommendation"] or ""
                tier = _LEVEL_TIER.get((p["cpiclevel"] or "").upper(), Tier.SPECULATIVE)
                out.append(Finding(
                    marker=genesymbol, source="cpic",
                    description=f"{drug_name}: {rec_txt[:160]}",
                    tier=tier, categories=[Category.CLINICAL],
                    detail={"gene": genesymbol, "drug": drug_name,
                            "cpic_level": p["cpiclevel"], "topic": "pharmacogenomic",
                            "modality": "genome", "phenotypes": phen,
                            "classification": r["classification"]}))
    finally:
        con.close()
    return out
