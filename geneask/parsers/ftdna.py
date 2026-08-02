# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Family Tree DNA — three export shapes from the same chip data:
  1. CSV, double-quoted, header 'RSID,CHROMOSOME,POSITION,RESULT' (concatenated genotype)
  2. Illumina raw: TAB-delimited, unquoted 'RSID CHROMOSOME POSITION RESULT'
  3. FamFinder: TAB-delimited with separate 'ALLELE1'/'ALLELE2' columns; header
     contains 'famfinder'
All build 37. One parser handles all three by sniffing the header shape.
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult, norm_chrom, norm_allele
from ._helpers import head_lines, data_rows


class FTDNAParser(GenotypeParser):
    name = "ftdna"
    default_build = "37"

    def _shape(self, head: list[str]):
        joined = "\n".join(head).lower()
        hdr = next((l for l in head if l.upper().replace('"', '').startswith("RSID")), "")
        if "famfinder" in joined or "allele1" in hdr.lower():
            return "famfinder"
        if '"' in hdr or ("," in hdr and "\t" not in hdr):
            return "csv"
        if "\t" in hdr and "rsid" in hdr.lower():
            return "illumina"
        return None

    def sniff(self, head: list[str]) -> bool:
        joined = "\n".join(head).lower()
        if "ftdna" in joined or "family tree" in joined or "famfinder" in joined:
            return True
        # header-only signature (some exports have no vendor comment)
        return self._shape(head) is not None and any(
            "rsid,chromosome,position,result" in l.lower().replace('"', '')
            or "rsid\tchromosome\tposition\tresult" in l.lower()
            for l in head)

    def parse(self, path: str) -> ParseResult:
        shape = self._shape(head_lines(path)) or "csv"
        res = ParseResult(source=f"ftdna_{shape}", build=self.default_build)
        if shape == "csv":
            rows = data_rows(path, ",", quoted=True)
        else:
            rows = data_rows(path, "\t")
        for f in rows:
            if not f or f[0].upper().replace('"', '') == "RSID":
                continue
            if shape == "famfinder":
                if len(f) < 5:
                    continue
                rsid, chrom, pos, a1, a2 = f[0], f[1], f[2], f[3], f[4]
            else:
                if len(f) < 4:
                    continue
                rsid, chrom, pos, result = f[0], f[1], f[2], f[3].strip()
                a1 = result[0] if len(result) >= 1 else ""
                a2 = result[1] if len(result) >= 2 else ""
            try:
                p = int(pos)
            except ValueError:
                continue
            r = GenotypeRecord(rsid, norm_chrom(chrom), p, norm_allele(a1), norm_allele(a2))
            if r.is_nocall:
                res.n_nocall += 1
            res.records.append(r)
        return res
