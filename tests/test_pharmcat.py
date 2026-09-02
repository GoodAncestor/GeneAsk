"""PharmCAT phenotype parsing against the plan-named JSON fields.

The official documentation was unreachable from the L6 sandbox. The fixture
uses geneSymbol, diplotype, phenotype, and activityScore from the approved plan.
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
