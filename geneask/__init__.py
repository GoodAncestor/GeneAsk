# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""GeneAsk: human genome variant interpretation.

The variant-interpretation sibling to MethylAsk (methylation). Answers "what do
my DNA variants mean" — ClinVar clinical screen, trait genotypes with
multi-source confidence, pharmacogenomics. Ingests 23andMe raw data, VCFs, and
the variant stream from an ONT modBAM. Depends on bio-core for shared mechanism
(provider/report layer, evidence tiering); knowledge (ClinVar, PRS, PGx) is here.
"""
__version__ = "0.1.0"
