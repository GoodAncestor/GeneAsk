#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (C) 2026 GoodAncestor
"""Build a posmap (rsid -> GRCh37 contig,pos) by lifting a 23andMe b36 file's
positions through an Ensembl NCBI36->GRCh37 chain with CrossMap's Python API."""
import argparse, json
from cmmodule.utils import read_chain_file, map_coordinates

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", required=True)
    a=ap.parse_args()

    maptree, tsize, ssize = read_chain_file(a.chain)
    n_in=n_out=0
    with open(a.infile,'r',errors='replace') as fh, open(a.out,'w') as out:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            f=line.rstrip('\n\r').split('\t')
            if len(f)<4: continue
            rsid,chrom,pos=f[0],f[1],f[2]
            try: p=int(pos)
            except: continue
            n_in+=1
            res=map_coordinates(maptree, chrom, p-1, p, '+')  # 0-based half-open
            if not res or len(res)<2: continue
            tgt=res[1]                        # (chrom, start, end, strand)
            out.write(f"{rsid}\t{tgt[0]}\t{tgt[1]+1}\n")   # 1-based
            n_out+=1
    json.dump(dict(input=n_in, lifted=n_out, dropped=n_in-n_out),
              open(a.stats,'w'), indent=2)
    print(f"lifted {n_out}/{n_in}")

if __name__=='__main__':
    main()
