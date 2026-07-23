# GeneAsk

Human genome **variant interpretation** — the sibling to
[MethylAsk](https://github.com/GoodAncestor/MethylAsk) (methylation), both built
on [bio-core](https://github.com/GoodAncestor/bio-core).

MethylAsk answers *"what does my methylation mean"*; GeneAsk answers *"what do my
DNA variants mean"* — a ClinVar clinical screen, trait genotypes with
multi-source confidence, and pharmacogenomics. It ingests 23andMe raw data,
VCFs, and the variant stream from an ONT modBAM (split out by bio-core).

## Scope

GeneAsk holds human variant **knowledge** (ClinVar, PRS, pharmacogenomics, the
trait framework). All organism-agnostic **mechanism** — the provider/report
layer, evidence tiering, coordinate handling — lives in bio-core and is imported.

## Modules

- `geneask.ingest.from_23andme` / `from_complete_genomics` / `consensus` — raw genotype → reference-anchored VCF, multi-source consensus
- `geneask.interpret.traits` — genotype lookup with multi-source confidence → bio-core Findings
- `geneask.interpret.clinvar_screen` — clinical screen against the bundled 157-gene ClinVar panel

## Bundled reference data

- `geneask/data/reference/clinvar_panel_157genes.json.gz` — ACMG SF v3.2 + carrier + pharmacogene panel (157 genes, gene-keyed, release 2026-06-06). Bundled so the clinical screen runs offline.

## The core lesson this encodes

Consumer-array "pathogenic" calls are frequently false positives — in the source
work, 10 of 11 were refuted by whole-genome sequencing. GeneAsk therefore reports
the calling platform and multi-source confidence, and **never presents a
single-array pathogenic hit as robust**. This is the same evidence-tiering
discipline MethylAsk applies to methylation associations.
