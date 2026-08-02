# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Consumer-genotype parser registry.

detect_parser() sniffs a file against each registered vendor parser (cheap header
check, first match wins) and returns the matching parser, or None. parse_file()
runs it and returns a normalized ParseResult.

Adding a vendor = add one module + append an instance to PARSERS.
"""
from __future__ import annotations
from .base import GenotypeParser, GenotypeRecord, ParseResult
from ._helpers import head_lines
from .twentythreeandme import TwentyThreeAndMeParser
from .ancestrydna import AncestryDNAParser
from .ftdna import FTDNAParser
from .myheritage import MyHeritageParser
from .livingdna import LivingDNAParser
from .myhappygenes import MyHappyGenesParser

# order matters only for ambiguous headers; vendor-comment sniffs are exclusive.
PARSERS: list[GenotypeParser] = [
    TwentyThreeAndMeParser(),
    AncestryDNAParser(),
    MyHeritageParser(),      # MyHeritage before FTDNA: shares the CSV shape but has a vendor comment
    LivingDNAParser(),
    MyHappyGenesParser(),
    FTDNAParser(),           # FTDNA last: its header-only signature is the broadest
]


def detect_parser(path: str) -> GenotypeParser | None:
    head = head_lines(path, 30)
    for p in PARSERS:
        try:
            if p.sniff(head):
                return p
        except Exception:
            continue
    return None


def parse_file(path: str) -> ParseResult | None:
    p = detect_parser(path)
    return p.parse(path) if p else None
