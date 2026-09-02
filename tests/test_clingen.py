import sqlite3
from pathlib import Path

from biocore.providers.base import Health

from geneask.annotators import clingen


VALIDITY = """ClinGen Gene-Disease Validity
Downloaded fixture
For testing only
Preamble four
Preamble five
Preamble six
GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE,GCEP
BRCA2,HGNC:1101,Hereditary breast ovarian cancer syndrome,MONDO:0003582,AD,SOP9,Definitive,https://clinicalgenome.org/brca2,2026-01-01,Hereditary Cancer GCEP
MUTYH,HGNC:7527,MUTYH polyposis,MONDO:0017312,AR,SOP9,Definitive,https://clinicalgenome.org/mutyh,2026-01-02,Hereditary Cancer GCEP
F8,HGNC:3546,Hemophilia A,MONDO:0010602,XL,SOP9,Definitive,https://clinicalgenome.org/f8,2026-01-03,Hemostasis GCEP
BRCA2,HGNC:1101,Another BRCA2 condition,MONDO:9999999,AR,SOP9,Limited,https://clinicalgenome.org/brca2-other,2026-01-04,Other GCEP
GENE5,HGNC:5,Condition five,MONDO:0000005,Unknown,SOP9,Moderate,https://clinicalgenome.org/gene5,2026-01-05,Other GCEP
"""

ACTIONABILITY_ADULT = """# docId\ttopicIri\tcurationType\tlatestSearchDate\tlastUpdated\tlastAuthor\tcontext\tcontextIri\trelease\treleaseDate\tgeneOrVariant\tgeneOmim\tdisease\tomim\tstatus-overall\tstatus-stg1\tstatus-stg2\tstatus-scoring\toutcome\toutcomeScoringGroup\tintervention\tinterventionScoringGroup\tseverity\tlikelihood\tnatureOfIntervention\teffectiveness\toverall
AC1\t/AC1\tGene-Condition\td\td\ta\tAdult\thttps://clinicalgenome.org/actionability/brca2\t1.0.0\td\tBRCA2\t600185\tHereditary breast ovarian cancer syndrome\t612555\tReleased\tComplete\tComplete\tComplete\tBreast cancer\tGroupA\tSurveillance\tGroupA\t3\t3\t3\t3\t12AA
AC1\t/AC1\tGene-Condition\td\td\ta\tAdult\thttps://clinicalgenome.org/actionability/brca2\t1.0.0\td\tBRCA2\t600185\tHereditary breast ovarian cancer syndrome\t612555\tReleased\tComplete\tComplete\tComplete\tOvarian cancer\tGroupA\tSurgery\tGroupB\t3\t2\t2\t2\t9CB
AC2\t/AC2\tGene-Condition\td\td\ta\tAdult\thttps://clinicalgenome.org/actionability/mutyh\t1.0.0\td\tMUTYH\t604933\tMUTYH polyposis\t608456\tReleased\tComplete\tComplete\tComplete\tPolyposis\tGroupA\tColonoscopy\tGroupA\t2\t2\t2\t2\t8CB
AC3\t/AC3\tGene-Condition\td\td\ta\tAdult\thttps://clinicalgenome.org/actionability/sdhd\t1.0.0\td\tSDHD\t602690\tHereditary paraganglioma\t168000\tRetracted\tComplete\tComplete\tComplete\tTumour\tGroupA\tImaging\tGroupA\t3\t3\t3\t3\t10CC
"""

ACTIONABILITY_PEDIATRIC = """# docId\ttopicIri\tcurationType\tlatestSearchDate\tlastUpdated\tlastAuthor\tcontext\tcontextIri\trelease\treleaseDate\tgeneOrVariant\tgeneOmim\tdisease\tomim\tstatus-overall\tstatus-stg1\tstatus-stg2\tstatus-scoring\toutcome\toutcomeScoringGroup\tintervention\tinterventionScoringGroup\tseverity\tlikelihood\tnatureOfIntervention\teffectiveness\toverall
AC4\t/AC4\tGene-Condition\td\td\ta\tPediatric\thttps://clinicalgenome.org/actionability/f8\t1.0.0\td\tF8\t300841\tHemophilia A\t306700\tReleased\tComplete\tComplete\tComplete\tBleeding\tGroupA\tFactor replacement\tGroupA\t3\t3\t3\t3\t10AB
"""


def _sources(tmp_path: Path):
    (tmp_path / "clingen_gene_validity.csv").write_text(VALIDITY)
    (tmp_path / "clingen_actionability_adult.tsv").write_text(ACTIONABILITY_ADULT)
    (tmp_path / "clingen_actionability_pediatric.tsv").write_text(
        ACTIONABILITY_PEDIATRIC
    )


def test_build_detects_headers_after_preambles_and_reads_inheritance(tmp_path):
    _sources(tmp_path)
    database = tmp_path / "clingen.db"

    summary = clingen.build_mirror(str(database), str(tmp_path))

    assert summary["validity"] == 5
    assert summary["actionability"] == 3
    assert clingen.inheritance_for(
        "BRCA2", ["MONDO:0003582"], db_path=str(database)
    ) == "dominant"
    assert clingen.inheritance_for(
        "MUTYH", ["MONDO:0017312"], db_path=str(database)
    ) == "recessive"
    assert clingen.inheritance_for(
        "F8", ["MONDO:0010602"], db_path=str(database)
    ) == "x_linked"


def test_lookup_and_adult_actionability_keep_source_fields(tmp_path):
    _sources(tmp_path)
    database = tmp_path / "clingen.db"
    clingen.build_mirror(str(database), str(tmp_path))

    result = clingen.lookup("BRCA2", db_path=str(database))
    actionability = clingen.actionability_for("BRCA2", db_path=str(database))

    assert result["validity"][0]["classification"] == "Definitive"
    assert any(row["mondo"] == "MONDO:0003582" for row in result["validity"])
    assert actionability["score"] == 12
    assert actionability["context"] == "adult"
    assert actionability["report_url"].endswith("/brca2")


def test_highest_classified_gene_row_is_the_fallback(tmp_path):
    _sources(tmp_path)
    database = tmp_path / "clingen.db"
    clingen.build_mirror(str(database), str(tmp_path))

    assert clingen.inheritance_for(
        "BRCA2", ["MONDO:does-not-match"], db_path=str(database)
    ) == "dominant"


def test_old_schema_is_refused_with_a_rebuild_note(tmp_path):
    database = tmp_path / "old.db"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
        connection.execute("INSERT INTO meta VALUES ('schema_version', '0')")

    assert clingen.lookup("BRCA2", db_path=str(database)) == {}
    status = clingen.mirror_status(str(database))
    assert status.health == Health.UNAVAILABLE
    assert "rebuild" in status.note.lower()


def test_actionability_keeps_the_highest_score_and_drops_retracted(tmp_path):
    _sources(tmp_path)
    db = str(tmp_path / "clingen.db")
    clingen.build_mirror(db, str(tmp_path))
    rows = clingen.lookup("BRCA2", db_path=db)["actionability"]
    assert [r["score"] for r in rows if r["context"] == "adult"] == [12.0]
    assert clingen.lookup("SDHD", db_path=db)["actionability"] == []
