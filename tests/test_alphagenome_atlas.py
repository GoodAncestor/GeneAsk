from geneask.annotators import alphagenome_atlas as atlas


def test_atlas_opt_in_and_configuration(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_ATLAS_ENABLED", raising=False)
    assert atlas.query_variant("1-100-A-G")["status"] == "disabled"
    monkeypatch.setenv("ALPHAGENOME_ATLAS_ENABLED", "1")
    monkeypatch.delenv("ALPHA_GENOME_KEY", raising=False)
    assert atlas.query_variant("1-100-A-G")["status"] == "missing_key"
    monkeypatch.setenv("ALPHA_GENOME_KEY", "unused")
    monkeypatch.delenv("ALPHAGENOME_ATLAS_SCORERS", raising=False)
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
    result = atlas.query_variant("1-100-A-G", requested_scorers=["example"])
    assert result["status"] == "complete" and result["n_tracks"] == 1
    assert result["tracks"][0]["gene_name"] == "TAL1"
    assert calls[0][0]["position"] == 100 and calls[0][1] == ["example"]
    def refusal(self, *a, **kw): raise RuntimeError("RESOURCE_EXHAUSTED")
    monkeypatch.setattr(Client, "query_variant", refusal)
    assert atlas.query_variant("1-100-A-G", requested_scorers=["example"])["status"] == "rate_limited"
