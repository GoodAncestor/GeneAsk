#!/usr/bin/env python3
"""CLI: annotated-VCF output — run GeneAsk's ClinVar (+AlphaMissense) screen
over a VCF and re-emit it with GA_* INFO fields and provenance in the header,
instead of a human-facing report. New tool (not generalized from the personal
"One Person, Six Genomes" recovery like the other scripts/ here) — see
geneask/interpret/annotated_vcf.py for the annotation logic and field docs.

Usage:
  python annotate_vcf.py --vcf input.vcf[.gz] --out annotated.vcf[.gz] \
      [--alphamissense-db path/to/alphamissense_mirror.db]
"""
import argparse
import json

from geneask.interpret.annotated_vcf import annotate_vcf


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="input VCF or VCF.gz")
    ap.add_argument("--out", required=True, help="output path; gzipped iff it ends .gz")
    ap.add_argument("--alphamissense-db", default=None,
                    help="AlphaMissense mirror SQLite (default: env ALPHAMISSENSE_MIRROR_DB "
                         "or the standard worker mirror path; GA_AM_* stays '.' if absent)")
    return ap.parse_args()


def main():
    a = parse_args()
    stats = annotate_vcf(a.vcf, a.out, alphamissense_db=a.alphamissense_db)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
