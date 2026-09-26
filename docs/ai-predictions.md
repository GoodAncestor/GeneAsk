# AI prediction integration

AlphaGenome enriches uncertain/conflicting catalogue findings and explicitly
selected `detail.novel_candidate = True` findings. An explicit single-variant
research request may instead set `detail.research_candidate = True`; this grants
scoring eligibility without asserting catalogue absence or biological novelty. Resolved benign/pathogenic
classifications are not overridden by the novel marker. The caller is responsible
for GRCh38 input, genotype quality, and a bounded candidate selection. These
computational results do not change clinical tiers.

Pass a `Pacing` object to `alphagenome_vep.annotate_findings` and export
`pacing.as_dict()`: eligible/scored/failed/skipped/cache_hits/no_result/spent count
unique variant IDs (the function's existing integer return counts enriched
findings). Disabled, missing key/client, partial, and rate-limited runs are
explicit. A successful empty result is not a benign classification.

AlphaGenome details retain the top five tracks, available gene/tissue metadata,
raw and quantile scores, reference assembly, query time and source. Direction is
`unknown`: recommended scorers mix signed differences, absolute differences, and
active-allele scores. The previous universal quantile-sign-to-increase/decrease
mapping was invalid for unsigned scorers. Legacy cached directions are suppressed.
Background quantiles are not patient disease probabilities. Model version is not
invented; `schema_version` versions our summary representation only.

For AlphaMissense, `mirror_status()` checks actual readable schema/rows without a
full scan. Pass `status={}` to `annotate_findings` to collect unique eligible,
scored, not_found, and skipped counts. A missing mirror and a lookup miss are
distinct. A lookup miss does not establish whether a variant is missense.

## AlphaGenome Atlas integration (2026-09-26)

`alphagenome_atlas.query_variant(variant_id, *, requested_scorers=None,
offline=False, cache_db=None, timeout_s=None)` returns a structured status,
tracks, optional `avi_score`, provenance, and cache information. Default explicit
scorers are `AVI_SCORE` and `AVI_SCORE_FEATURE_IMPORTANCE`. Track summaries retain
up to five strongest entries per scorer; raw values are never ranked across
unrelated scorers. Missing requested scorers are explicit. AVI scores describe
research variant impact, not personal disease probability.

`annotate_findings(findings, *, offline=False, status=None, max_variants=None,
cache_db=None)` enriches all valid SNV markers from standard human chromosomes
without changing clinical tiers. The caller must establish GRCh38 coordinates;
this module cannot infer genome assembly from a bare variant string. It returns
enriched finding count, and writes unique-variant coverage to `status`:
`eligible/scored/failed/skipped/cache_hits/local_hits/not_found/partial`, plus
`remote_attempted`, reasons, and local installation status.

API inference is opt-in (`ALPHAGENOME_ATLAS_ENABLED=1`, `ALPHA_GENOME_KEY`). Each
SDK query runs in a subprocess terminated and reaped at its hard timeout
(`ALPHAGENOME_ATLAS_TIMEOUT_S`, default 20 seconds). Report remote calls share a
25-second budget (`ALPHAGENOME_ATLAS_TIME_BUDGET_S`) and maximum 10 remote attempts
(`ALPHAGENOME_ATLAS_MAX_VARIANTS`, maximum 50). Local/cache annotations do not
consume the remote-attempt cap. Remote priority is explicit research requests,
then novel/uncertain/conflicting candidates, then other variants; within each
group higher available local AVI scores come first, with variant-ID tie breaking. Callers with an enclosing subprocess must kill
the whole process group on their own timeout so nested children cannot survive.

Complete API summaries are persisted in `ALPHAGENOME_ATLAS_CACHE_DB` (default
`~/.cache/geneask/alphagenome_atlas.db`). Keys include summary schema version,
`ALPHAGENOME_ATLAS_DATA_VERSION` (default 2026-09), normalized variant, and scorer
set. Cache TTL defaults to 30 days (`ALPHAGENOME_ATLAS_CACHE_TTL_S`). Failures, partial scorer coverage and
missing results are never cached. `offline=True` uses only cache/local scores,
without importing the SDK or contacting an API. Disabled API configuration still
permits existing cache and local records to be read.

### Local indexed AVI scores

Set `ALPHAGENOME_ATLAS_AVI_FILE` to an installed local BGZF/Tabix file or directory
of chromosome shards. `atlas_avi.local_status()` reports missing/invalid indexes
or unsupported schemas. `atlas_avi.lookup_many(variant_ids, path=None,
status=None)` groups indexed seeks into 64kb windows using one open handle per
shard. It returns exact allele matches with `avi_score` and optional `avi_phred`.
A miss is not a benign prediction. Local AVI covers only AVI_SCORE; when requested,
feature importance still requires API/cache coverage. Local scores survive remote
failure, with `missing_scorers` and `remote_status` exposed.

The official 88.5GB archive is described as AVI and Phred-scaled scores in Tabix
format. Its actual column layout has not yet been verified locally because the
download endpoint returned HTTP 500. The adapter intentionally requires named
chromosome, 1-based position, REF, ALT and AVI score columns. It accepts header
aliases or an explicit `ALPHAGENOME_ATLAS_AVI_COLUMNS` comma-separated list supplied
after inspecting the downloaded file. It never guesses column offsets or
silently treats BED starts as 1-based positions. Tests exercise real synthetic
BGZF+TBI files, not an installed full official mirror. Do not report full local
coverage until installation and reference-variant validation succeed.

Sources:
- https://www.alphagenomedocs.com/api/atlas.html
- https://www.alphagenomedocs.com/api/generated/alphagenome.atlas.atlas.AtlasClient.html
- https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/atlas/atlas.py
- https://deepmind.google.com/science/alphagenome/_/download/atlas/avi_scores_snvs_tabix.zip

## Request timeout limitation

Verified against official SDK v0.8.0 source: `DnaClient.score_variant` has no
`timeout` parameter and invokes the streaming RPC without a deadline. We bound
connection readiness with the supported `create(timeout=10.0)` argument. The
report pacing deadline bounds waiting for rate capacity and whether another
request starts; it cannot cancel an in-flight score RPC or the SDK retry loop.
Do not describe this as a hard end-to-end request timeout. A true RPC deadline
requires SDK support or a separately reviewed gRPC channel integration.

Source: https://github.com/google-deepmind/alphagenome/blob/v0.8.0/src/alphagenome/models/dna_client.py

The Atlas `data_version` is a configured cache release label, not a claim that
the API returned a specific server model version. The verified live feature
importance response uses `var.name` for feature labels and has no quantile layer;
raw feature attributions are preserved with their original feature names.

## Commercial and non-commercial output modes

`DNAREPORT_OUTPUT_MODE=noncommercial` preserves the existing research behavior.
`commercial` enables the shared `biocore.licensing` output policy. An unrecognized
value applies commercial restrictions and reports a configuration note; it never
silently enables non-commercial datasets. No earlier output-mode flag was found
in the audited GeneAsk, bio-core or DNA-Report code.

| Data route | Commercial output policy |
| --- | --- |
| Official downloaded AVI SNV / Phred scores | Allowed under the portal's permissive downloadable-artifact category |
| Public AlphaGenome Atlas API, including API-returned AVI and cached API results | Withheld; the app only recognizes the downloaded AVI exception |
| Downloaded AVI feature importance and merged splicing scores | Withheld; portal lists non-commercial use only |
| AlphaGenome API and cached inference | Withheld under non-commercial access/output policy |
| AlphaMissense mirror predictions | Withheld under CC BY-NC-SA 4.0 |

This is the app's source-eligibility policy, not a general assurance about every
possible commercial agreement with the provider. It does not enable separate
Google Cloud commercial contracts. Changing score names or copying API values
into a local cache does not change their provenance or eligibility.

Commercial-mode engine checks happen before restricted API/cache/mirror access.
Atlas uses the actual installed local AVI artifact, with no restricted API
fallback. Model results carry structured license metadata with terms links and
verification date. The shared output filter also removes previously stored
restricted findings, appended prediction text, and derived interpretation text;
the application must filter report-level notes and rebuild derived summaries.
A commercial export must not simply relabel a frozen non-commercial report.

Verified 2026-09-26 against the rendered official download portal and terms:
- https://deepmind.google.com/science/alphagenome/downloads
- https://deepmind.google.com/science/alphagenome/terms
- https://deepmind.google.com/science/alphagenome/output-terms
- https://zenodo.org/records/8208688/files/README.md
