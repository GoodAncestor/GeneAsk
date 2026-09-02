"""PharmCAT phenotype parsing against the real 3.4.0 JSON shape.

The fixture is trimmed from a real run on alien02 (2026-09-02); two genes are
edited to called values in the same shape so a call is exercised.
"""
import shutil
import subprocess
from pathlib import Path

from geneask.annotators import pharmcat


FIXTURE = Path(__file__).parent / "fixtures" / "pharmcat.phenotype.json"


def test_call_diplotypes_parses_the_phenotype_json(monkeypatch, tmp_path):
    shutil.copyfile(FIXTURE, tmp_path / "pharmcat.phenotype.json")
    seen = {}

    def run(command, **kwargs):
        seen["command"] = command
        seen.update(kwargs)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(pharmcat.subprocess, "run", run)

    calls = pharmcat.call_diplotypes("sample.vcf", str(tmp_path), timeout=12)

    assert calls["CYP2C19"] == {
        "diplotype": "*1/*2",
        "phenotype": "Intermediate Metabolizer",
        "activity_score": 1.0,
        "source": "PharmCAT 3.4.0",
    }
    assert calls["CYP2D6"]["diplotype"] == "*1/*4"
    assert seen["command"] == [
        "pharmcat_pipeline",
        "sample.vcf",
        "-o",
        str(tmp_path),
        "-reporterJson",
    ]
    assert seen["timeout"] == 12


def test_platform_gate_allows_sequencing_only():
    assert pharmcat.platform_ok("WGS")
    assert pharmcat.platform_ok("wes")
    assert not pharmcat.platform_ok("ARRAY")
    assert not pharmcat.platform_ok(None)


def test_available_uses_the_pipeline_on_path(monkeypatch):
    monkeypatch.setattr(pharmcat.shutil, "which", lambda name: "/bin/pharmcat")
    assert pharmcat.available()
    monkeypatch.setattr(pharmcat.shutil, "which", lambda name: None)
    assert not pharmcat.available()


def test_pipeline_failure_returns_empty_calls_with_a_note(monkeypatch, tmp_path):
    monkeypatch.setattr(
        pharmcat.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 2, "", "bad"),
    )

    calls = pharmcat.call_diplotypes("sample.vcf", str(tmp_path))

    assert calls == {}
    assert "status 2" in calls.note


def test_pipeline_timeout_returns_empty_calls_with_a_note(monkeypatch, tmp_path):
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs["timeout"])

    monkeypatch.setattr(pharmcat.subprocess, "run", timeout)

    calls = pharmcat.call_diplotypes("sample.vcf", str(tmp_path), timeout=1)

    assert calls == {}
    assert "1-second limit" in calls.note


def test_uncalled_genes_are_absent_not_results(monkeypatch, tmp_path):
    shutil.copyfile(FIXTURE, tmp_path / "x.phenotype.json")
    monkeypatch.setattr(pharmcat.subprocess, "run",
                        lambda command, **kw: subprocess.CompletedProcess(command, 0, "", ""))
    calls = pharmcat.call_diplotypes("sample.vcf", str(tmp_path))
    assert "ABCG2" not in calls and "SLCO1B1" not in calls
    assert set(calls) == {"CYP2C19", "CYP2D6"}
