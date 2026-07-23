#!/usr/bin/env python3
"""Reference-genotype resolver: for trait rsIDs absent from the variant callset,
confirm homozygous-reference status using (1) 23andMe array 0/0 calls and
(2) Complete Genomics reference blocks. Emits the plus-strand hom-ref genotype.

Inputs:
  --rsids           file with one rsid per line (the ABSENT trait SNPs)
  --array-vcfs      comma list of per-source array norm VCFs (contain 0/0 records)
  --cg-tsv          Complete Genomics var file (GRCh37, has ref blocks)
  --fasta37         GRCh37 fasta (to read ref base for CG-confirmed sites)
  --out             output TSV: rsid, chrom37, pos37, ref_base, refgt, confirmed_by
"""
import argparse, subprocess, csv, sys

def run(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rsids', required=True)
    ap.add_argument('--array-vcfs', required=True)
    ap.add_argument('--cg-tsv', required=True)
    ap.add_argument('--fasta37', required=True)
    ap.add_argument('--out', required=True)
    a=ap.parse_args()
    want=[l.strip() for l in open(a.rsids) if l.strip()]
    wantset=set(want)
    # 1) find GRCh37 pos + array 0/0 status from array VCFs
    pos={}   # rsid -> (chrom,pos,ref)
    array_refhom=set()
    for vcf in a.array_vcfs.split(','):
        q=run(f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%GT]\\n' {vcf}")
        for line in q.split('\n'):
            if not line: continue
            f=line.split('\t')
            if len(f)<6: continue
            c,p,rsid,ref,alt,gt=f[:6]
            if rsid in wantset:
                pos[rsid]=(c,p,ref)
                g=gt.replace('|','/').split('/')
                if all(x=='0' for x in g if x not in ('.','')):
                    array_refhom.add(rsid)
    # 2) confirm via CG ref blocks: build interval check per chrom
    # gather needed positions per chrom
    from collections import defaultdict
    bychrom=defaultdict(list)
    for rsid,(c,p,ref) in pos.items():
        bychrom[c].append((int(p),rsid,ref))
    cg_ref=set()
    # single pass over CG file, only ref rows, check overlap
    # build sorted needs
    needs={c:sorted(v) for c,v in bychrom.items()}
    with open(a.cg_tsv) as fh:
        for line in fh:
            if '\tref\t' not in line: continue
            f=line.split('\t')
            if len(f)<7: continue
            chrom=f[3].replace('chr',''); 
            if chrom not in needs: continue
            b=int(f[4]); e=int(f[5])
            for p,rsid,ref in needs[chrom]:
                if b<=p-1 and e>=p:  # 0-based half-open contains 1-based p
                    cg_ref.add(rsid)
    # write
    with open(a.out,'w',newline='') as out:
        w=csv.writer(out,delimiter='\t')
        w.writerow(['rsid','chrom37','pos37','ref_base','refgt','confirmed_by'])
        for rsid in want:
            if rsid in pos:
                c,p,ref=pos[rsid]
                conf=[]
                if rsid in array_refhom: conf.append('array0/0')
                if rsid in cg_ref: conf.append('CG_ref_block')
                refgt=f"{ref}/{ref}" if ref else ''
                w.writerow([rsid,c,p,ref,refgt,'+'.join(conf) if conf else 'unconfirmed'])
            else:
                w.writerow([rsid,'','','','','not_on_array'])
    print("resolver done")

if __name__=='__main__': main()
