# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Shared line-reading helpers for the vendor parsers."""
from __future__ import annotations
import gzip, csv, io


def open_text(path: str):
    return (gzip.open(path, "rt", encoding="utf-8", errors="replace")
            if path.endswith(".gz")
            else open(path, "r", encoding="utf-8", errors="replace"))


def head_lines(path: str, n: int = 30) -> list[str]:
    out = []
    with open_text(path) as fh:
        for _ in range(n):
            line = fh.readline()
            if not line:
                break
            out.append(line.rstrip("\n\r"))
    return out


def data_rows(path: str, delimiter: str, quoted: bool = False):
    """Yield non-comment, non-blank data rows as lists of fields. Skips lines
    starting with '#' and the column-header line (handled by caller via skip_header)."""
    with open_text(path) as fh:
        if quoted:
            reader = csv.reader(fh, delimiter=delimiter)
            for row in reader:
                if not row or row[0].startswith("#") or not row[0].strip():
                    continue
                yield row
        else:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                yield line.rstrip("\n\r").split(delimiter)
