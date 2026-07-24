"""AncestryDNA raw export. Tab-delimited, 5 columns: rsid, chrom, pos, allele1,
allele2. Detection: '#AncestryDNA' in the first comment line. Build 37. No-call
= '0' per allele. Chrom codes 23=X 24=Y 25=PAR->X 26=MT. V1/V2 chips share layout.
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult, norm_chrom, norm_allele
from ._helpers import head_lines, data_rows


class AncestryDNAParser(GenotypeParser):
    name = "ancestrydna"
    default_build = "37"

    def sniff(self, head: list[str]) -> bool:
        return any("ancestrydna" in l.lower() for l in head[:5])

    def parse(self, path: str) -> ParseResult:
        res = ParseResult(source=self.name, build=self.default_build)
        for f in data_rows(path, "\t"):
            if len(f) < 5 or f[0].lower() == "rsid":
                continue
            rsid, chrom, pos, a1, a2 = f[0], f[1], f[2], f[3], f[4]
            try:
                p = int(pos)
            except ValueError:
                continue
            r = GenotypeRecord(rsid, norm_chrom(chrom), p, norm_allele(a1), norm_allele(a2))
            if r.is_nocall:
                res.n_nocall += 1
            res.records.append(r)
        return res
