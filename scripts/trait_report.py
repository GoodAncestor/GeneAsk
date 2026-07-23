#!/usr/bin/env python3
"""
Trait-analysis framework for the unified multi-source genome.
Given a trait table (rsid, effect_allele, trait, notes) and the unified
callset, report the individual's genotype at each SNP with multi-source
provenance and a confidence flag.

Genotypes are reported on the plus strand as base pairs, matching consumer
raw-data conventions (23andMe/Genetic Genie/Promethease style).

Usage:
  python trait_report.py --vcf unified.GRCh38.vcf.gz --traits traits.csv \
      --out report.csv [--build38-index rsid_pos.tsv]

traits.csv columns (rsid required; others optional):
  rsid, trait, gene, effect_allele, risk_genotype, notes
"""
import argparse, subprocess, sys, csv, json

def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout

def load_callset_by_rsid(vcf, rsids):
    """Return {rsid: dict(chrom,pos,ref,alt,gt_bases,platf,nsrc,conc,srcs,srcmask)}."""
    want=set(rsids)
    out={}
    # single pass over VCF, filter by ID
    q=run(f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t[%GT]\\t%INFO/PLATF\\t%INFO/NSRC\\t%INFO/CONC\\t%INFO/SRCS\\t%INFO/SRCMASK\\n' {vcf}")
    for line in q.split('\n'):
        if not line: continue
        f=line.split('\t')
        if len(f)<11: continue
        chrom,pos,rsid,ref,alt,gt,platf,nsrc,conc,srcs,srcmask=f[:11]
        if rsid not in want: continue
        alts=alt.split(',')
        def base(i):
            if i in ('.',''): return '.'
            if i=='0': return ref
            try: return alts[int(i)-1]
            except: return '.'
        g=gt.replace('|','/').split('/')
        if len(g)==1:
            bases=base(g[0])  # hemizygous
        else:
            bs=[base(x) for x in g]
            if '.' in bs:
                called=[b for b in bs if b!='.']
                bases='/'.join(called+['.']) if called else './.'
            else:
                bases='/'.join(sorted(bs))
        out[rsid]=dict(chrom=chrom,pos=pos,ref=ref,alt=alt,gt_bases=bases,
                       platf=platf,nsrc=int(nsrc),conc=int(conc),srcs=srcs,srcmask=srcmask)
    return out

def confidence(rec):
    if rec is None: return 'ABSENT'
    if rec['platf']=='BOTH' and rec['conc']==1: return 'HIGH (array+WGS agree)'
    if rec['platf']=='BOTH' and rec['conc']==0: return 'CONFLICT (array vs WGS disagree)'
    if rec['nsrc']>=3 and rec['conc']==1: return 'HIGH (multi-source agree)'
    if rec['platf']=='WGS': return 'MODERATE (WGS only)'
    if rec['platf']=='ARRAY': return 'MODERATE (array only)'
    return 'LOW'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--traits', required=True)
    ap.add_argument('--out', required=True)
    a=ap.parse_args()
    traits=list(csv.DictReader(open(a.traits)))
    rsids=[t['rsid'] for t in traits if t.get('rsid')]
    calls=load_callset_by_rsid(a.vcf, rsids)
    rows=[]
    for t in traits:
        rsid=t.get('rsid')
        rec=calls.get(rsid)
        row=dict(t)
        row['genotype']= rec['gt_bases'] if rec else 'not_in_callset'
        row['chrom']= rec['chrom'] if rec else ''
        row['pos_GRCh38']= rec['pos'] if rec else ''
        row['n_sources']= rec['nsrc'] if rec else 0
        row['platform']= rec['platf'] if rec else ''
        row['sources']= rec['srcs'] if rec else ''
        row['confidence']= confidence(rec)
        # carrier of effect allele?
        ea=t.get('effect_allele','').strip()
        if rec and ea:
            row['carries_effect_allele']= 'yes' if ea in rec['gt_bases'].split('/') else 'no'
        else:
            row['carries_effect_allele']=''
        rows.append(row)
    cols=list(rows[0].keys())
    with open(a.out,'w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    n_found=sum(1 for r in rows if r['genotype'] not in ('not_in_callset',))
    print(json.dumps(dict(traits=len(rows), found=n_found, out=a.out)))

if __name__=='__main__':
    main()
