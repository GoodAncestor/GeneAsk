"""Consumer-genotype file parsers — base types.

A parser turns a vendor raw-data export into a stream of normalized
GenotypeRecords, all on the same shape regardless of vendor quirks:
    (rsid, chrom, pos, allele1, allele2)

Vendor formats differ in delimiter, column order, concatenated-vs-split
genotype, chromosome coding, no-call encoding, and build. Each parser owns
those quirks and emits the same record, so downstream (VCF conversion, ClinVar
screen, trait table) is vendor-agnostic.

Format facts (columns/headers/builds/chrom-codes) are just facts — reimplemented
cleanroom from published vendor specs and the `snps` package's documentation.
Nothing here is copied from another implementation.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from abc import ABC, abstractmethod


@dataclass
class GenotypeRecord:
    rsid: str
    chrom: str          # normalized: 1..22, X, Y, MT
    pos: int
    allele1: str        # single base, or "-"/"0" normalized to "" for no-call
    allele2: str        # "" for haploid calls (X/Y/MT in males)

    @property
    def is_nocall(self) -> bool:
        return not self.allele1 and not self.allele2

    @property
    def genotype(self) -> str:
        return (self.allele1 + self.allele2) or "--"


@dataclass
class ParseResult:
    source: str                      # vendor name, e.g. "ancestrydna"
    build: str                       # "36" | "37" | "38" | "unknown" (header claim)
    records: list = field(default_factory=list)
    n_nocall: int = 0
    notes: list = field(default_factory=list)


# canonical chromosome normalization shared across array vendors. Numeric codes
# 23-26 are the Illumina/array convention several vendors use.
_CHR_MAP = {"23": "X", "24": "Y", "25": "X", "26": "MT",  # 25 = PAR -> X
            "XY": "X", "M": "MT", "MT": "MT", "X": "X", "Y": "Y"}


def norm_chrom(c: str) -> str:
    c = c.strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    return _CHR_MAP.get(c, c)


def norm_allele(a: str) -> str:
    """Normalize one allele token; no-call markers ('-','0','.','I','D',...) -> ''.
    Only A/C/G/T survive (indel placeholders I/D can't be coordinate-anchored)."""
    a = a.strip().upper()
    return a if a in ("A", "C", "G", "T") else ""


class GenotypeParser(ABC):
    """One parser per vendor format. `sniff` decides if a file is this format
    (cheap header check); `parse` yields the normalized records."""
    name: str = "base"
    default_build: str = "unknown"

    @abstractmethod
    def sniff(self, head: list[str]) -> bool:
        """Return True if the first ~30 lines look like this vendor's format."""
        ...

    @abstractmethod
    def parse(self, path: str) -> ParseResult:
        ...
