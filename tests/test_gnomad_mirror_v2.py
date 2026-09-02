import gzip
import sqlite3

from geneask.annotators.gnomad_mirror import (
    _ensure_schema, _ingest_chrom, lookup_many, mirror_status,
)


_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr13\t32316419\t.\tCAG\tC\t.\tPASS\tAF=0.0000032;AC=5;AN=1560000;nhomalt=0;AF_nfe=0.0000041;AF_afr=0
chr13\t32316500\t.\tA\tG,T\t.\tPASS\tAF=0.1,0.2;AC=10,20;AN=100;nhomalt=1,2;AF_nfe=0.11,0.21
"""


def test_v2_rows_keep_counts_and_populations(tmp_path):
    path = tmp_path / "c.vcf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(_VCF)
    db = str(tmp_path / "g.db")
    con = sqlite3.connect(db)
    _ensure_schema(con)
    count = _ingest_chrom(con, path, None)
    con.commit()
    con.close()
    assert count == 3
    rows = lookup_many(
        ["13-32316419-CAG-C", "13-32316500-A-T"], db_path=db
    )
    row = rows["13-32316419-CAG-C"]
    assert row["af"] == 0.0000032 and row["ac"] == 5
    assert row["an"] == 1560000 and row["nhomalt"] == 0
    assert row["populations"] == {"nfe": 0.0000041, "afr": 0.0}
    row2 = rows["13-32316500-A-T"]
    assert row2["ac"] == 20 and row2["nhomalt"] == 2
    assert row2["populations"]["nfe"] == 0.21


def test_v1_database_is_refused(tmp_path):
    db = str(tmp_path / "old.db")
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE af(variant_id TEXT PRIMARY KEY, af REAL)")
    con.execute("INSERT INTO af VALUES ('13-1-A-G', 0.1)")
    con.commit()
    con.close()
    assert lookup_many(["13-1-A-G"], db_path=db) == {}
    status = mirror_status(db)
    assert status.health.value == "unavailable"
    assert "schema" in (status.note or "").lower()
