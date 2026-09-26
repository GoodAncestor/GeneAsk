from geneask.annotators import alphagenome_atlas as atlas


def test_atlas_opt_in_and_configuration(monkeypatch, tmp_path):
    monkeypatch.setenv("ALPHAGENOME_ATLAS_CACHE_DB", str(tmp_path / "cache.db"))
    monkeypatch.delenv("ALPHAGENOME_ATLAS_AVI_FILE", raising=False)
    monkeypatch.delenv("ALPHAGENOME_ATLAS_ENABLED", raising=False)
    assert atlas.query_variant("1-100-A-G")["status"] == "disabled"
    monkeypatch.setenv("ALPHAGENOME_ATLAS_ENABLED", "1")
    monkeypatch.delenv("ALPHA_GENOME_KEY", raising=False)
    assert atlas.query_variant("1-100-A-G")["status"] == "missing_key"
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.setenv("ALPHAGENOME_ATLAS_SCORERS", "")
    assert atlas.query_variant("1-100-A-G")["status"] == "scorers_required"
    assert atlas.query_variant("1-100-AT-G")["status"] == "not_applicable"


def test_atlas_query_contract_and_quota(monkeypatch):
    import sys
    import types
    monkeypatch.setenv("ALPHAGENOME_ATLAS_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    for name in ("alphagenome", "alphagenome.atlas", "alphagenome.data"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    class Matrix:
        def __getitem__(self, key): return .98
    class Metadata:
        @property
        def iloc(self): return self
        def __getitem__(self, key): return {"gene_name": "TAL1"}
    record = types.SimpleNamespace(shape=(1, 1), X=Matrix(), layers={"quantiles": Matrix()},
                                   obs=Metadata(), var=Metadata())
    calls = []
    class Client:
        def query_variant(self, variant, *, requested_scorers):
            calls.append((variant, requested_scorers))
            return {"example": record}
    sys.modules["alphagenome.atlas"].atlas = types.SimpleNamespace(create=lambda key, timeout: Client())
    sys.modules["alphagenome.data"].genome = types.SimpleNamespace(Variant=lambda **kw: kw)
    result = atlas._remote_query("1-100-A-G", ["example"])
    assert result["status"] == "complete" and result["n_tracks"] == 1
    assert result["tracks"][0]["gene_name"] == "TAL1"
    assert calls[0][0]["position"] == 100 and calls[0][1] == ["example"]
    def refusal(self, *a, **kw): raise RuntimeError("RESOURCE_EXHAUSTED")
    monkeypatch.setattr(Client, "query_variant", refusal)
    assert atlas._remote_query("1-100-A-G", ["example"])["status"] == "rate_limited"


def _enable(monkeypatch, tmp_path):
    monkeypatch.setenv("ALPHAGENOME_ATLAS_ENABLED", "1")
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.setenv("ALPHAGENOME_ATLAS_CACHE_DB", str(tmp_path / "cache.db"))
    monkeypatch.delenv("ALPHAGENOME_ATLAS_AVI_FILE", raising=False)
    monkeypatch.delenv("ALPHAGENOME_ATLAS_SCORERS", raising=False)


def _success(vid):
    return {**atlas._base(vid), "status": "complete", "avi_score": .75,
            "tracks": [{"scorer": "AVI_SCORE", "raw_score": .75}],
            "provenance": "alphagenome_atlas_api"}


def test_versioned_cache_offline_and_failure_retry(monkeypatch, tmp_path):
    _enable(monkeypatch, tmp_path)
    calls = []
    def remote(vid, scorers, timeout):
        calls.append((vid, scorers))
        return _success(vid)
    monkeypatch.setattr(atlas, "_run_remote", remote)
    assert atlas.query_variant("chr1-100-A-G")["avi_score"] == .75
    assert calls[0][1] == sorted(atlas.DEFAULT_SCORERS)
    monkeypatch.delenv("ALPHA_GENOME_KEY")
    assert atlas.query_variant("1-100-A-G", offline=True)["cache_hit"]
    assert len(calls) == 1
    assert atlas.query_variant("1-100-A-G", requested_scorers=["AVI_SCORE"], offline=True)["status"] == "offline"
    monkeypatch.setenv("ALPHAGENOME_ATLAS_DATA_VERSION", "future")
    assert atlas.query_variant("1-100-A-G", offline=True)["status"] == "offline"
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.setattr(atlas, "_run_remote", lambda *a: {"status": "timeout"})
    assert atlas.query_variant("1-100-A-G")["status"] == "timeout"
    monkeypatch.setattr(atlas, "_run_remote", remote)
    assert atlas.query_variant("1-100-A-G")["status"] == "complete"
    assert len(calls) == 2


def test_hard_subprocess_timeout_kills_and_reaps(monkeypatch):
    import subprocess
    import sys
    import time
    original_run = subprocess.run
    children = []
    original_popen = subprocess.Popen
    def popen(*args, **kwargs):
        child = original_popen(*args, **kwargs)
        children.append(child)
        return child
    monkeypatch.setattr(subprocess, "Popen", popen)
    def sleeping_worker(command, **kwargs):
        return original_run([sys.executable, "-c", "import time; time.sleep(30)"], **kwargs)
    monkeypatch.setattr(subprocess, "run", sleeping_worker)
    start = time.monotonic()
    result = atlas._run_remote("1-100-A-G", ["AVI_SCORE"], .1)
    assert result["status"] == "timeout"
    assert time.monotonic() - start < 3
    assert children and children[0].poll() is not None


def test_local_scores_survive_remote_failure_and_cap(monkeypatch, tmp_path):
    from geneask.annotators import atlas_avi
    from biocore.providers.base import Finding, Tier, Category
    _enable(monkeypatch, tmp_path)
    records = {f"1-{p}-A-G": {**_success(f"1-{p}-A-G"), "provenance": "alphagenome_atlas_local_avi"}
               for p in (100, 200, 300)}
    def local(variants, *, status):
        status.update(status="ready", available=True)
        return records
    monkeypatch.setattr(atlas_avi, "lookup_many", local)
    calls = []
    def remote(*args):
        calls.append(args)
        return {"status": "timeout"}
    monkeypatch.setattr(atlas, "_run_remote", remote)
    fs = [Finding(vid, "clinvar", "example", Tier.ROBUST, [Category.CLINICAL], detail={}) for vid in records]
    status = {}
    assert atlas.annotate_findings(fs, status=status, max_variants=1) == 3
    assert len(calls) == 1 and status["remote_attempted"] == 1
    assert status["local_hits"] == 3 and status["partial"] == 3
    assert fs[0].detail["alphagenome_atlas"]["remote_status"] == "timeout"
    assert fs[0].tier == Tier.ROBUST


def test_cache_does_not_consume_remote_cap(monkeypatch, tmp_path):
    from biocore.providers.base import Finding, Tier, Category
    _enable(monkeypatch, tmp_path)
    calls = []
    def remote(vid, *args):
        calls.append(vid)
        return _success(vid)
    monkeypatch.setattr(atlas, "_run_remote", remote)
    atlas.query_variant("1-100-A-G")
    fs = [Finding(f"1-{p}-A-G", "variant_lookup", "example", Tier.SPECULATIVE, [Category.CLINICAL], detail={}) for p in (100,200)]
    status = {}
    assert atlas.annotate_findings(fs, status=status, max_variants=1) == 2
    assert calls == ["1-100-A-G", "1-200-A-G"]
    assert status["cache_hits"] == 1 and status["remote_attempted"] == 1


def test_partial_api_results_preserved_but_not_cached(monkeypatch, tmp_path):
    _enable(monkeypatch, tmp_path)
    calls = []
    def remote(vid, *args):
        calls.append(vid)
        return {**_success(vid), "status": "partial", "missing_scorers": ["AVI_SCORE_FEATURE_IMPORTANCE"]}
    monkeypatch.setattr(atlas, "_run_remote", remote)
    assert atlas.query_variant("1-100-A-G")["status"] == "partial"
    assert atlas.query_variant("1-100-A-G")["status"] == "partial"
    assert len(calls) == 2


def test_remote_features_merge_local_avi(monkeypatch, tmp_path):
    _enable(monkeypatch, tmp_path)
    monkeypatch.setattr(atlas, "lookup_local_avi", lambda vid: {**_success(vid), "provenance": "alphagenome_atlas_local_avi"})
    monkeypatch.setattr(atlas, "_run_remote", lambda vid, *args: {
        **atlas._base(vid), "status": "partial", "missing_scorers": ["AVI_SCORE"],
        "tracks": [{"scorer": "AVI_SCORE_FEATURE_IMPORTANCE", "raw_score": .5, "feature_name": "test feature"}], "n_tracks": 1})
    result = atlas.query_variant("1-100-A-G")
    assert result["status"] == "complete" and result["avi_score"] == .75
    assert not result["missing_scorers"] and result["local_avi"]
    assert len(result["tracks"]) == 2


def test_remote_priority_prefers_requested_and_novel_to_resolved(monkeypatch, tmp_path):
    from biocore.providers.base import Finding, Tier, Category
    _enable(monkeypatch, tmp_path)
    calls = []
    def remote(vid, *args):
        calls.append(vid)
        return _success(vid)
    monkeypatch.setattr(atlas, "_run_remote", remote)
    fs = [Finding(f"1-{p}-A-G", "variant_lookup", "example", Tier.SPECULATIVE,
                  [Category.CLINICAL], detail=detail) for p, detail in (
            (100, {"clinical_significance": "Benign"}),
            (200, {"novel_candidate": True}),
            (300, {"research_candidate": True}))]
    status = {}
    atlas.annotate_findings(fs, status=status, max_variants=2)
    assert calls == ["1-300-A-G", "1-200-A-G"]
    assert status["scored"] == 2 and status["skipped"] == 1
