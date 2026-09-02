import sqlite3

from geneask.annotators.cpic_pgx import recommendations_for_phenotype


def _mirror(tmp_path):
    database = tmp_path / "pgx.db"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE drug(drugid TEXT PRIMARY KEY, name TEXT)")
        connection.execute("INSERT INTO drug VALUES ('RxNorm:1', 'clopidogrel')")
        connection.execute(
            "CREATE TABLE pair(genesymbol TEXT, drugid TEXT, cpiclevel TEXT, "
            "guidelineid INTEGER)"
        )
        connection.execute(
            "INSERT INTO pair VALUES ('CYP2C19', 'RxNorm:1', 'A', 100)"
        )
        connection.execute(
            """CREATE TABLE recommendation(
                guidelineid INTEGER, drugid TEXT, phenotypes TEXT,
                drugrecommendation TEXT, classification TEXT,
                implications TEXT, lookupkey TEXT)"""
        )
        connection.execute(
            "INSERT INTO recommendation VALUES "
            "(100, 'RxNorm:1', '{\"CYP2C19\":\"Intermediate Metabolizer\"}', "
            "'Consider another antiplatelet', 'Strong', '{}', '{}')"
        )
        connection.execute(
            "INSERT INTO recommendation VALUES "
            "(100, 'RxNorm:1', '{\"CYP2C19\":\"Poor Metabolizer\"}', "
            "'Use another antiplatelet', 'Strong', '{}', '{}')"
        )
    return str(database)


def test_recommendations_select_the_gene_phenotype_and_keep_the_call(tmp_path):
    findings = recommendations_for_phenotype(
        "CYP2C19",
        "Intermediate Metabolizer",
        db_path=_mirror(tmp_path),
        diplotype="*1/*2",
        activity_score=1.0,
        diplotype_source="PharmCAT 3.4.0",
    )

    assert len(findings) == 1
    detail = findings[0].detail
    assert detail["diplotype"] == "*1/*2"
    assert detail["phenotype"] == "Intermediate Metabolizer"
    assert detail["activity_score"] == 1.0
    assert detail["diplotype_source"] == "PharmCAT 3.4.0"
    assert detail["recommendation"] == "Consider another antiplatelet"


def test_recommendations_do_not_guess_when_the_phenotype_is_missing(tmp_path):
    assert recommendations_for_phenotype(
        "CYP2C19", "", db_path=_mirror(tmp_path), diplotype="*1/*2"
    ) == []


def test_recommendations_do_not_guess_the_diplotype_caller(tmp_path):
    findings = recommendations_for_phenotype(
        "CYP2C19",
        "Intermediate Metabolizer",
        db_path=_mirror(tmp_path),
        diplotype="*1/*2",
    )

    assert "diplotype_source" not in findings[0].detail
