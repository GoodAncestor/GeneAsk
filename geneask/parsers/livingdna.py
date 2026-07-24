"""Living DNA raw export. TAB-delimited despite a .csv extension. Columns rsid,
chrom, pos, genotype (concatenated). Detection: 'Living DNA' in the first line.
Build 37, forward strand. No-call = '--'. IDs may be AX-/AFFX-prefixed or CHR:POS
positional — passed through as the rsid field.
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult, norm_chrom, norm_allele
from ._helpers import data_rows


class LivingDNAParser(GenotypeParser):
    name = "livingdna"
    default_build = "37"

    def sniff(self, head: list[str]) -> bool:
        return any("living dna" in l.lower() for l in head[:5])

    def parse(self, path: str) -> ParseResult:
        res = ParseResult(source=self.name, build=self.default_build)
        for f in data_rows(path, "\t"):
            if len(f) < 4 or f[0].lower() == "rsid":
                continue
            rsid, chrom, pos, geno = f[0], f[1], f[2], f[3].strip()
            try:
                p = int(pos)
            except ValueError:
                continue
            a1 = geno[0] if len(geno) >= 1 else ""
            a2 = geno[1] if len(geno) >= 2 else ""
            r = GenotypeRecord(rsid, norm_chrom(chrom), p, norm_allele(a1), norm_allele(a2))
            if r.is_nocall:
                res.n_nocall += 1
            res.records.append(r)
        return res
