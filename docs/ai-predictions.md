# AI prediction integration

AlphaGenome enriches uncertain/conflicting catalogue findings and explicitly
selected `detail.novel_candidate = True` findings. Resolved benign/pathogenic
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

## AlphaGenome Atlas evaluation (2026-09-25)

Official SDK API is now available:
`alphagenome.atlas.atlas.create(key, timeout=30).query_variant(variant,
requested_scorers=[...])`, returning scorer-name -> AnnData with raw scores,
optional quantile layers and gene/track metadata. `scorer_metadata()` discovers
names and signedness. The optional `alphagenome_atlas.query_variant` adapter uses
this API directly, accepts single-nucleotide GRCh38 substitutions only, and
returns explicit configuration/failure/coverage status. No private endpoint or
invented scorer name is used.

Enable only after installing an SDK version that includes `alphagenome.atlas`,
confirming access terms for the application and validating the desired scorer
names from `scorer_metadata()`. Set `ALPHAGENOME_ATLAS_ENABLED=1`, server-side
`ALPHA_GENOME_KEY`, and `ALPHAGENOME_ATLAS_SCORERS` (comma-separated explicit
names, maximum 20). Call only in an authorized analysis request, not on report
GET. Keep AVI raw scores labeled with their scorer; do not pool them with
molecular effect scores or interpret them as personal disease probability.
This adapter does not silently enable Atlas, replace existing inference, or
claim a production call has succeeded. Add a separately versioned persistent
Atlas cache before using bulk queries; current adapter is for single-variant
evaluation only.

Sources:
- https://www.alphagenomedocs.com/variant_scoring.html
- https://www.alphagenomedocs.com/api/atlas.html
- https://www.alphagenomedocs.com/api/generated/alphagenome.atlas.atlas.AtlasClient.html
- https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/atlas/atlas.py
- https://deepmind.google/science/alphagenome/
