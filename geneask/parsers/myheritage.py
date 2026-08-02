# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""MyHeritage DNA raw export. CSV, double-quoted, columns rsid, chrom, pos,
RESULT (concatenated genotype). Detection: 'MyHeritage' in a comment line.
Build 37. Handles the doubled-double-quote field variant via csv.reader.
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult, norm_chrom, norm_allele
from ._helpers import data_rows


class MyHeritageParser(GenotypeParser):
    name = "myheritage"
    default_build = "37"

    def sniff(self, head: list[str]) -> bool:
        return any("myheritage" in l.lower() for l in head[:8])

    def parse(self, path: str) -> ParseResult:
        res = ParseResult(source=self.name, build=self.default_build)
        for f in data_rows(path, ",", quoted=True):
            if not f or f[0].upper().replace('"', '') == "RSID":
                continue
            if len(f) < 4:
                continue
            rsid, chrom, pos, result = f[0], f[1], f[2], f[3].strip()
            try:
                p = int(pos)
            except ValueError:
                continue
            a1 = result[0] if len(result) >= 1 else ""
            a2 = result[1] if len(result) >= 2 else ""
            r = GenotypeRecord(rsid, norm_chrom(chrom), p, norm_allele(a1), norm_allele(a2))
            if r.is_nocall:
                res.n_nocall += 1
            res.records.append(r)
        return res
