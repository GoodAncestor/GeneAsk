"""MyHappyGenes (Tempus) genotype export. TAB-delimited, 5 columns: SNP Name,
Chr, Position, 'Allele1 - Forward', 'Allele2 - Forward'. Detection:
'MyHappyGenes' or 'TEMPUS' in a comment. Header claims build 37.1 but real
exports ship GRCh38 coordinates — build is flagged 'unknown' so the caller runs
position-based auto-detection rather than trusting the header (allelix ADR-0021).
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult, norm_chrom, norm_allele
from ._helpers import data_rows


class MyHappyGenesParser(GenotypeParser):
    name = "myhappygenes"
    default_build = "unknown"    # header says 37.1, data is 38 — don't trust it

    def sniff(self, head: list[str]) -> bool:
        j = "\n".join(head[:8]).lower()
        return "myhappygenes" in j or "tempus" in j

    def parse(self, path: str) -> ParseResult:
        res = ParseResult(source=self.name, build=self.default_build)
        res.notes.append("header build claim ignored (known 37.1/38 mismatch); "
                         "run position-based build detection")
        for f in data_rows(path, "\t"):
            if len(f) < 5 or f[0].lower().startswith("snp name") or f[0].lower() == "rsid":
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
