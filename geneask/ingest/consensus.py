#!/usr/bin/env python3
"""Collapse merged 6-sample VCF -> single consensus VCF with provenance INFO.
INFO fields added:
  SRCS   = comma list of sources with a called GT (e.g. A2013,CG)
  NSRC   = number of supporting sources
  SRCMASK= bitmask string over [A2011,A2012,A2013,Av5,Av5ph,CG]
  CONC   = 1 if all calling sources agree on GT (allele-set), 0 otherwise
  PLATF  = ARRAY | WGS | BOTH (which platform types support it)
  CONSGT = consensus genotype (majority; ties -> CG preferred then most-recent array)
The output single sample 'CameronThomson' carries CONSGT as its GT.
"""
import sys
from collections import Counter

SAMPLES = ["A2011", "A2012", "A2013", "Av5", "Av5ph", "CG"]
ARRAY = set(["A2011", "A2012", "A2013", "Av5", "Av5ph"])


def allele_set(rec, s):
    g = rec.samples[s]
    a = g.get('GT')
    if a is None or all(x is None for x in a):
        return None
    alls = rec.alleles
    bases = tuple(sorted(set(alls[x] for x in a if x is not None)))
    return bases, a


def build_consensus(inp: str, outp: str) -> int:
    """Collapse a merged 6-sample VCF into a single provenance-tagged consensus
    VCF (sample 'CameronThomson'). Returns the number of records written.

    pysam is imported inside the function so the module imports without it
    (it is an optional 'ingest' extra); a caller that runs the consensus needs
    pysam installed.
    """
    import pysam
    vin = pysam.VariantFile(inp)
    h = vin.header.copy()
    for tag, num, typ, desc in [
        ("SRCS", ".", "String", "Sources with a called genotype at this site"),
        ("NSRC", "1", "Integer", "Number of sources supporting this site"),
        ("SRCMASK", "1", "String", "Bitmask over A2011,A2012,A2013,Av5,Av5ph,CG"),
        ("CONC", "1", "Integer", "1 if all calling sources agree on genotype allele-set"),
        ("PLATF", "1", "String", "Platform support: ARRAY, WGS, or BOTH"),
        ("CONSGT", "1", "String", "Consensus genotype string")]:
        if tag not in h.info:
            h.info.add(tag, num, typ, desc)
    nh = pysam.VariantHeader()
    for rec in str(h).strip().split('\n'):
        if rec.startswith('##'):
            nh.add_line(rec)
    nh.add_sample("CameronThomson")
    vout = pysam.VariantFile(outp, 'w', header=nh)

    n_written = 0
    for rec in vin:
        if len(rec.alleles) < 2 or rec.alts is None:
            continue
        calls = {}
        for s in SAMPLES:
            r = allele_set(rec, s)
            if r is not None:
                calls[s] = r
        if not calls:
            continue
        ref0 = rec.alleles[0]
        if all(all(b == ref0 for b in calls[s][0]) for s in calls):
            continue
        srcs = [s for s in SAMPLES if s in calls]
        mask = ''.join('1' if s in calls else '0' for s in SAMPLES)
        setcounts = Counter(calls[s][0] for s in srcs)
        top = setcounts.most_common()
        maxn = top[0][1]
        tied = [k for k, v in top if v == maxn]
        if len(tied) == 1:
            consset = tied[0]
        else:
            pref = ['CG', 'Av5', 'A2013', 'A2012', 'A2011', 'Av5ph']
            consset = None
            for p in pref:
                if p in calls and calls[p][0] in tied:
                    consset = calls[p][0]
                    break
            if consset is None:
                consset = tied[0]
        conc = 1 if len(setcounts) == 1 else 0
        has_array = any(s in ARRAY for s in srcs)
        has_wgs = 'CG' in srcs
        platf = 'BOTH' if (has_array and has_wgs) else ('ARRAY' if has_array else 'WGS')
        consgt_tuple = None
        for p in ['CG', 'Av5', 'A2013', 'A2012', 'A2011', 'Av5ph']:
            if p in calls and calls[p][0] == consset:
                consgt_tuple = calls[p][1]
                break
        if consgt_tuple is None:
            consgt_tuple = calls[srcs[0]][1]
        nr = vout.new_record(contig=rec.contig, start=rec.start, stop=rec.stop,
                             alleles=rec.alleles, id=rec.id)
        nr.info['SRCS'] = ','.join(srcs)
        nr.info['NSRC'] = len(srcs)
        nr.info['SRCMASK'] = mask
        nr.info['CONC'] = conc
        nr.info['PLATF'] = platf
        gtstr = '/'.join('.' if x is None else str(x) for x in consgt_tuple)
        nr.info['CONSGT'] = gtstr
        nr.samples['CameronThomson']['GT'] = consgt_tuple
        vout.write(nr)
        n_written += 1
    vin.close()
    vout.close()
    return n_written


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: consensus.py <merged.vcf(.gz)> <out.vcf(.gz)>")
    n = build_consensus(sys.argv[1], sys.argv[2])
    print(f"done: {n} consensus records written")


if __name__ == '__main__':
    main()
