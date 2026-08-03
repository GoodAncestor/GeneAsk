"""gnomAD AF annotator: cache behavior + in-place finding annotation (no network)."""
import sqlite3
from geneask.annotators.gnomad_freq import (allele_frequency, annotate_findings, _cache_con,
                                            LiveBudget, default_budget)
from biocore.providers.base import Finding, Tier, Category


def test_cache_hit_no_network(tmp_path):
    db = str(tmp_path / "af.db")
    con = _cache_con(db)
    con.execute("INSERT INTO af VALUES (?,?)", ("1-100-A-G", 0.35)); con.commit(); con.close()
    # cached value returned without any API call
    assert allele_frequency("1-100-A-G", cache_db=db) == 0.35
    # cached miss (-1) returns None, also no API call
    con = _cache_con(db); con.execute("INSERT INTO af VALUES (?,?)", ("1-200-C-T", -1.0)); con.commit(); con.close()
    assert allele_frequency("1-200-C-T", cache_db=db) is None


def test_annotate_only_variant_markers(tmp_path):
    db = str(tmp_path / "af.db")
    con = _cache_con(db)
    con.execute("INSERT INTO af VALUES (?,?)", ("1-100-A-G", 0.35)); con.commit(); con.close()
    fs = [Finding("1-100-A-G", "clinvar", "gene X: pathogenic", Tier.MODERATE, [Category.CLINICAL], detail={}),
          Finding("cg999", "ewas", "methylation", Tier.ROBUST, [Category.AGING], detail={})]
    n = annotate_findings(fs, cache_db=db)
    assert n == 1
    assert fs[0].detail["gnomad_af"] == 0.35
    assert "35.0% of people" in fs[0].description
    assert "gnomad_af" not in fs[1].detail    # CpG marker untouched


def _live_counter(monkeypatch, responses=None):
    """Replace urlopen so a test can count live calls without a network."""
    calls = []
    class _R:
        def __init__(self, body): self._b = body
        def read(self): return self._b
        def __enter__(self): return self
        def __exit__(self, *a): return False
    def fake(req, timeout=None):
        calls.append(req)
        if responses is not None:
            r = responses.pop(0)
            if isinstance(r, Exception):
                raise r
        return _R(b'{"data":{"variant":{"genome":{"af":0.5},"exome":null}}}')
    monkeypatch.setattr("urllib.request.urlopen", fake)
    return calls


def test_cap_stops_live_calls_but_cached_answers_still_serve(tmp_path, monkeypatch):
    """The cap is on what we ask Broad for, not on what we answer — a cached
    variant past the cap must still be annotated, or raising the cap would change
    the report's content rather than just its API usage."""
    db = str(tmp_path / "af.db")
    con = _cache_con(db)
    con.execute("INSERT INTO af VALUES (?,?)", ("9-900-A-G", 0.11)); con.commit(); con.close()
    calls = _live_counter(monkeypatch)
    fs = [Finding(f"1-{i}-A-G", "clinvar", "x", Tier.MODERATE, [Category.CLINICAL], detail={})
          for i in range(5)]
    fs.append(Finding("9-900-A-G", "clinvar", "cached one", Tier.MODERATE, [Category.CLINICAL], detail={}))
    budget = LiveBudget(2)
    n = annotate_findings(fs, cache_db=db, budget=budget)
    assert len(calls) == 2                     # only the allowance was spent
    assert budget.remaining == 0 and budget.spent == 2
    assert fs[5].detail["gnomad_af"] == 0.11   # cached, past the cap, still annotated
    assert n == 3                              # 2 live + 1 cached


def test_over_cap_miss_is_not_cached(tmp_path, monkeypatch):
    """A variant we never asked about must not be recorded as absent — caching that
    would make the cap permanent for that variant across every future report."""
    db = str(tmp_path / "af.db")
    _live_counter(monkeypatch)
    allele_frequency("1-100-A-G", cache_db=db, budget=LiveBudget(0))
    con = _cache_con(db)
    assert con.execute("SELECT count(*) FROM af").fetchone()[0] == 0
    con.close()


def test_429_halts_the_run(tmp_path, monkeypatch):
    """A 429 is the API saying stop. Continuing into it wastes their capacity and
    ours, and looks identical to 'variant unknown' in the output."""
    import urllib.error
    db = str(tmp_path / "af.db")
    err = urllib.error.HTTPError("u", 429, "Too Many Requests", {}, None)
    calls = _live_counter(monkeypatch, responses=[err, err, err, err])
    fs = [Finding(f"1-{i}-A-G", "clinvar", "x", Tier.MODERATE, [Category.CLINICAL], detail={})
          for i in range(4)]
    budget = LiveBudget(10)
    annotate_findings(fs, cache_db=db, budget=budget)
    assert len(calls) == 1 and budget.halted   # stopped after the first refusal
    con = _cache_con(db)
    assert con.execute("SELECT count(*) FROM af").fetchone()[0] == 0   # nothing poisoned
    con.close()
