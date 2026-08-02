"""Annotated-VCF output — VCF in, VCF out, with GeneAsk annotation as declared INFO.

Turns GeneAsk from a report generator into a pipeline step: read any VCF, run it
through GeneAsk's existing per-variant annotation path (the ClinVar clinical
screen — see clinvar_screen.py — layered with AlphaMissense when a mirror is
built), and re-emit the same records with GA_* INFO fields added and matching
##INFO header declarations. Every original record and column survives
unchanged — CHROM/POS/ID/REF/ALT/QUAL/FILTER, the original INFO content,
FORMAT and all sample columns — only fresh INFO keys are appended, and only on
records that actually got a hit. A downstream tool reads the output with any
standard VCF parser (bcftools, pysam, cyvcf2); this is a real VCF, not a
parallel format.

Layers wired by default:
  - ClinVar        always — the full mirror if CLINVAR_MIRROR_DB is built
                    (clinvar_mirror.py), else the bundled 157-gene panel
                    (clinvar_screen.py), same preference order the clinical
                    screen already uses.
  - AlphaMissense   only when a mirror is built (alphamissense.py); silent
                    no-op otherwise, same behaviour as its own annotate_findings().

Deliberately NOT wired here: gnomAD (a per-variant network API call — fine for
the handful of variants in one person's report, wrong for a bulk VCF pipeline
step that could be millions of rows) and GWAS Catalog / CPIC (rsID- or
gene-keyed rather than a single per-ALT-allele scalar, and their free-text
trait/drug-recommendation strings are exactly the kind of value that would
need heavy escaping to survive a comma-delimited INFO list). Both fit the same
(lookup_fn, source_description) shape _clinvar_layer/_alphamissense_layer
already use below, and are natural follow-ons, not a limitation of the design.

Multi-allelic records: every GA_* field is declared Number=A (one value per
ALT allele — VCF's own convention for per-allele annotation, see the spec
sec. 1.6.1) because ClinVar/AlphaMissense key on ('chrom','pos','ref', one
specific alt); collapsing to a single shared value per record would silently
misattribute one allele's annotation to another. An allele with no hit is
represented by the literal '.' at its position in the list, per Number=A
semantics — never a dropped/misaligned field.
"""
from __future__ import annotations
import gzip
import re
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__ as _GENEASK_VERSION
from .clinvar_screen import load_panel, index_by_variant_id

_TOOL = "geneask"

# INFO fields this module can emit, in the (id, number, type, description)
# shape a ##INFO header line needs. Number=A everywhere: see module docstring.
_INFO_FIELDS = [
    ("GA_CLNSIG", "A", "String",
     "GeneAsk/ClinVar clinical significance for this ALT allele "
     "('.' if the allele has no hit in the consulted panel/mirror)"),
    ("GA_CLNGENE", "A", "String",
     "Gene symbol for this ALT allele's ClinVar hit ('.' if none)"),
    ("GA_CLNREVSTAT", "A", "String",
     "ClinVar review status for this ALT allele's hit ('.' if none)"),
    ("GA_CLNSTARS", "A", "Integer",
     "ClinVar gold-star review confidence (0-4) for this ALT allele ('.' if none)"),
    ("GA_AM_CLASS", "A", "String",
     "AlphaMissense pathogenicity class for this ALT allele "
     "('.' if no AlphaMissense mirror was built, or the allele has no hit)"),
    ("GA_AM_SCORE", "A", "Float",
     "AlphaMissense pathogenicity score in [0,1] for this ALT allele ('.' if none)"),
    ("GA_AM_PROT", "A", "String",
     "AlphaMissense protein-variant notation for this ALT allele ('.' if none)"),
]

# VCF String INFO values may not contain the characters that delimit INFO
# fields/list entries ( ; = , ) or whitespace (spec sec. 1.6.1); collapse any
# run of them to a single underscore so e.g. a ClinVar significance of
# "Likely pathogenic, low penetrance" can't be mistaken for two Number=A
# allele entries, or corrupt the surrounding key=value;key=value structure.
# Double-quote isn't a body-line delimiter per spec, but strip it too — an
# upstream source string with a stray '"' has no business surfacing raw in a
# body value either, and it costs nothing to be defensive here.
_UNSAFE_VALUE = re.compile(r'[,;="\s]+')

# a bare (unquoted) token for a VCF meta-line sub-field: none of the chars
# that would otherwise be misread as the next field or corrupt the line.
_BARE_META = re.compile(r'^[^,;="\s]*$')


def _safe_token(s) -> str:
    """Sanitize one annotation value for a Number=A INFO list slot. Falls back
    to '.' (the Number=A 'no value here' marker) for empty/missing input."""
    if s is None:
        return "."
    s = str(s).strip()
    if not s:
        return "."
    return _UNSAFE_VALUE.sub("_", s)


def _escape_meta(v: str) -> str:
    """Backslash-escape a value bound for a double-quoted VCF meta field —
    both ##INFO Description="..." and our own provenance line use this. VCF
    4.2 sec 1.4.2: backslash and double-quote inside a quoted string must be
    backslash-escaped. A header is one line, so embedded newlines are flattened."""
    return v.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").replace("\r", " ")


def _meta_field(value) -> str:
    """Render one <Key=value> sub-value of a VCF meta-information line: bare
    if it needs no escaping, double-quoted + escaped otherwise. This is what
    stops a Description (or a provenance Sources string) containing a comma,
    '=', or quote from corrupting the header line it lives on — requirement
    the annotated-VCF writer must get right, not just the happy path."""
    value = str(value)
    return value if _BARE_META.match(value) else f'"{_escape_meta(value)}"'


def _meta_line(tag: str, **fields) -> str:
    body = ",".join(f"{k}={_meta_field(v)}" for k, v in fields.items())
    return f"##{tag}=<{body}>"


def _info_header_lines(existing_ids: set) -> list[str]:
    return [_meta_line("INFO", ID=i, Number=n, Type=t, Description=d)
            for i, n, t, d in _INFO_FIELDS if i not in existing_ids]


def _provenance_lines(sources: list[str]) -> list[str]:
    """Requirement: a reader must be able to tell what produced this file and
    against what. One line: tool, version, UTC run timestamp, and every
    annotation source actually consulted this run (panel/mirror + release,
    or 'not_available' when a layer's mirror wasn't built)."""
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return [_meta_line("geneaskAnnotateVCF", Tool=_TOOL, Version=_GENEASK_VERSION,
                       RunTimestamp=ts, Sources=";".join(sources))]


def _clinvar_layer(explicit_panel: dict | None):
    """Index + a human-readable description of which ClinVar source was used.
    Mirrors clinvar_screen.index_by_variant_id's own mirror-first preference,
    but (unlike that function) also reports which one it picked, for provenance."""
    if explicit_panel is not None:
        idx = index_by_variant_id(explicit_panel)
        return idx, f"clinvar:explicit_panel ({len(idx)} variants)"
    try:
        from ..annotators.clinvar_mirror import load_panel_from_mirror
        full = load_panel_from_mirror()
    except Exception:
        full = None
    if full:
        return full, f"clinvar:full_mirror ({len(full)} variants, NCBI ClinVar GRCh38)"
    panel = load_panel()
    release = next(iter(panel.values()), {}).get("release", "unknown")
    idx = index_by_variant_id(panel)
    return idx, f"clinvar:bundled_157gene_panel (release={release}, {len(idx)} variants)"


def _alphamissense_layer(db_path: str | None):
    """Returns (lookup_fn_or_None, source_description). lookup_fn_or_None is
    None when no mirror was built — the caller then leaves GA_AM_* at '.'
    without touching sqlite at all, same no-op AlphaMissense already documents
    for annotate_findings()."""
    from ..annotators import alphamissense
    if not alphamissense.mirror_available(db_path):
        return None, "alphamissense:not_available (no mirror built)"
    return (lambda vid: alphamissense.lookup(vid, db_path)), "alphamissense:mirror"


def _chrom_for_lookup(chrom: str) -> str:
    """Annotation sources key variant_id on bare chrom ('1'..'22','X','Y','MT'
    — NCBI/Ensembl style, same as clinvar_mirror.py and to_carried.py). Input
    VCFs may use UCSC 'chr'-prefixed contigs; strip only for the lookup key —
    the record's own CHROM column is never rewritten (round-trip requirement)."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _lookup_alt(chrom: str, pos: str, ref: str, alt: str, clinvar_idx: dict, am_lookup):
    """One ALT allele -> (variant_id, clinvar_record_or_None, alphamissense_record_or_None),
    or None if the allele isn't a plain substitution/indel GeneAsk can key at
    all (no-call '.', spanning-deletion '*', or symbolic/breakend ALT like
    '<DEL>' / 'G]17:1584563]' — none of these have a literal ALT sequence to
    build 'chrom-pos-ref-alt' from)."""
    if not alt or alt in (".", "*") or alt.startswith("<") or "[" in alt or "]" in alt:
        return None
    vid = f"{_chrom_for_lookup(chrom)}-{pos}-{ref}-{alt}"
    cv = clinvar_idx.get(vid)
    am = am_lookup(vid) if am_lookup else None
    return vid, cv, am


def _append_info(info: str, add: dict) -> str:
    """Append GA_* key=value pairs to an existing INFO string. An empty/absent
    INFO ('.') becomes just the new keys; anything already there (including a
    bare FLAG with no '=') is preserved verbatim ahead of ours."""
    parts = [] if info in (".", "") else info.split(";")
    parts.extend(f"{k}={v}" for k, v in add.items())
    return ";".join(parts)


def annotate_vcf(in_path: str, out_path: str, *,
                 clinvar_panel: dict | None = None,
                 alphamissense_db: str | None = None) -> dict:
    """Read in_path (VCF or VCF.gz), annotate with GeneAsk's ClinVar screen
    (+ AlphaMissense when a mirror is available), write out_path in the same
    format (gzipped iff out_path ends '.gz'). Returns run stats:
    {records, alt_alleles, annotated, clinvar_hits, alphamissense_hits, sources}.

    Every input record is re-emitted in its original order with every
    original column untouched. A record GeneAsk can't key at all (symbolic
    ALT, breakend, no-call) or simply has no hit passes through with its INFO
    field byte-identical — no GA_* keys are added unless at least one ALT
    allele actually got a hit (see module docstring on Number=A).

    clinvar_panel overrides the bundled/mirror ClinVar panel (mainly for
    tests); alphamissense_db overrides the AlphaMissense mirror path (default:
    env ALPHAMISSENSE_MIRROR_DB or the standard worker mirror location).
    """
    clinvar_idx, clinvar_src = _clinvar_layer(clinvar_panel)
    am_lookup, am_src = _alphamissense_layer(alphamissense_db)
    sources = [clinvar_src, am_src]

    stats = dict(records=0, alt_alleles=0, annotated=0, clinvar_hits=0,
                 alphamissense_hits=0, sources=sources)

    opener = gzip.open if str(in_path).endswith(".gz") else open
    with opener(in_path, "rt") as fh:
        lines = fh.readlines()

    out_lines: list[str] = []
    existing_info_ids: set = set()
    found_chrom_line = False
    body_start = len(lines)   # no #CHROM line found -> treat as headers-only, no data section
    for i, line in enumerate(lines):
        if line.startswith("##"):
            out_lines.append(line.rstrip("\n"))
            m = re.match(r"##INFO=<ID=([^,]+),", line)
            if m:
                existing_info_ids.add(m.group(1))
        elif line.startswith("#CHROM"):
            found_chrom_line = True
            body_start = i + 1
            break
        else:
            # a body line before any #CHROM header — malformed input; treat
            # everything from here on as data and let the per-record loop
            # below pass it through untouched rather than swallowing it.
            body_start = i
            break

    out_lines.extend(_info_header_lines(existing_info_ids))
    out_lines.extend(_provenance_lines(sources))
    if found_chrom_line:
        out_lines.append(lines[body_start - 1].rstrip("\n"))

    for line in lines[body_start:]:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        stats["records"] += 1
        if len(f) < 8:
            out_lines.append(line)   # too short to be a real record — pass through untouched
            continue
        chrom, pos, vid_col, ref, alt_field, qual, filt, info = f[:8]
        rest = f[8:]
        alts = alt_field.split(",")
        stats["alt_alleles"] += len(alts)

        per_field = {k: [] for k, *_ in _INFO_FIELDS}
        any_hit = False
        for alt in alts:
            hit = _lookup_alt(chrom, pos, ref, alt, clinvar_idx, am_lookup)
            if hit is None:
                for k, *_ in _INFO_FIELDS:
                    per_field[k].append(".")
                continue
            _vid, cv, am = hit
            if cv:
                stats["clinvar_hits"] += 1
                any_hit = True
                per_field["GA_CLNSIG"].append(_safe_token(cv.get("clinical_significance")))
                per_field["GA_CLNGENE"].append(_safe_token(cv.get("gene")))
                per_field["GA_CLNREVSTAT"].append(_safe_token(cv.get("review_status")))
                stars = cv.get("gold_stars")
                per_field["GA_CLNSTARS"].append(str(int(stars)) if stars is not None else ".")
            else:
                per_field["GA_CLNSIG"].append(".")
                per_field["GA_CLNGENE"].append(".")
                per_field["GA_CLNREVSTAT"].append(".")
                per_field["GA_CLNSTARS"].append(".")
            if am:
                stats["alphamissense_hits"] += 1
                any_hit = True
                per_field["GA_AM_CLASS"].append(_safe_token(am.get("am_class")))
                score = am.get("pathogenicity")
                per_field["GA_AM_SCORE"].append(f"{score:.4f}" if score is not None else ".")
                per_field["GA_AM_PROT"].append(_safe_token(am.get("protein_variant")))
            else:
                per_field["GA_AM_CLASS"].append(".")
                per_field["GA_AM_SCORE"].append(".")
                per_field["GA_AM_PROT"].append(".")

        if any_hit:
            stats["annotated"] += 1
            add = {k: ",".join(v) for k, v in per_field.items()}
            info = _append_info(info, add)

        out_lines.append("\t".join([chrom, pos, vid_col, ref, alt_field, qual, filt, info] + rest))

    out_opener = gzip.open if str(out_path).endswith(".gz") else open
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with out_opener(out_path, "wt") as out:
        out.write("\n".join(out_lines) + "\n")
    return stats
