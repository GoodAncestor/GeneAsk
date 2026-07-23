# Provenance — GeneAsk

The ingestion and interpretation code was generalized from the "One Person, Six
Genomes" work (personal `colbyt/genomics` repo): five 23andMe chip tests
(2011–2026) plus a 2014 Complete Genomics whole genome, unified onto GRCh38 with
per-site provenance (4.11M variants).

| GeneAsk module | Generalized from |
|---|---|
| `ingest/from_23andme.py` | 23andme_to_vcf.py — reference-anchored VCF from consumer raw data |
| `ingest/from_complete_genomics.py` | cg_var_to_vcf.py — CG var file → callset |
| `ingest/consensus.py` | build_consensus.py — multi-source provenance-tagged merge |
| `interpret/traits.py` | trait_report.py — genotype + multi-source confidence, now emitting bio-core Findings |
| `interpret/clinvar_screen.py` | ClinVar 157-gene panel screen |
| `scripts/refgt_resolver.py`, `lift_b36_posmap.py` | reference-genotype resolver, build lift (kept as CLI tools) |

Personal data (the unified genome, medical findings) was used as the development
fixture and is NOT committed here; only the generalized code and the public
ClinVar panel are. Evidence tiers map the recovered multi-source confidence:
array+WGS agree or ≥3 sources → ROBUST; single platform → MODERATE; array-vs-WGS
conflict → SPECULATIVE.
