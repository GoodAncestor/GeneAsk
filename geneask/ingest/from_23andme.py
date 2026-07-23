#!/usr/bin/env python3
"""23andMe raw genotype -> VCF, anchoring REF to a reference FASTA.

23andMe format: rsid, chromosome, position, genotype  (plus-strand oriented).
Genotypes: "AA","AG"(SNP), "--"(no-call), single char (haploid Y/MT/male-X),
           "II"/"DD"/"DI"/"I"/"D" (insertion/deletion placeholders — sequence
           not given, so these cannot be turned into coordinate-anchored REF/ALT
           and are exported to a separate unresolved-indel table, not the SNV VCF).
Phased format: rsid, chromosome, position, allele1, allele2 (two columns).
"""
import sys, argparse, gzip

VALID_BASES = set("ACGT")

def open_text(p):
    return gzip.open(p,'rt') if p.endswith('.gz') else open(p,'r',errors='replace')

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-vcf", required=True)            # bgzipped path (we write plain, bgzip outside)
    ap.add_argument("--out-indels", required=True)          # unresolved indel table
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--phased", action="store_true")
    ap.add_argument("--posmap", default=None,
                    help="TSV: rsid<TAB>new_contig<TAB>new_pos to remap coords (e.g. hg18->hg19 lifted). "
                         "Records whose rsid is absent from the map are dropped (unliftable).")
    return ap.parse_args()

def main():
    import pysam  # optional 'ingest' extra; imported lazily
    a = parse_args()
    fa = pysam.FastaFile(a.fasta)
    ref_contigs = set(fa.references)
    # map 23andMe chrom -> fasta contig (Ensembl uses 1..22,X,Y,MT)
    def contig(c):
        if c in ref_contigs: return c
        if c=="MT" and "MT" in ref_contigs: return "MT"
        if c.startswith("chr") and c[3:] in ref_contigs: return c[3:]
        if ("chr"+c) in ref_contigs: return "chr"+c
        return None

    posmap=None
    if a.posmap:
        posmap={}
        with open_text(a.posmap) as ph:
            for ln in ph:
                if not ln.strip(): continue
                q=ln.rstrip('\n').split('\t')
                if len(q)>=3:
                    posmap[q[0]]=(q[1],int(q[2]))

    stats = dict(total=0, written=0, nocall=0, hom_ref=0, het=0, hom_alt=0,
                 haploid=0, multiallelic=0, indel_placeholder=0, ref_mismatch=0,
                 skipped_no_contig=0, bad_allele=0, transitions=0, transversions=0,
                 unliftable=0)
    PURINES={'A','G'}; PYR={'C','T'}
    def is_ts(r,alt):
        return (r in PURINES and alt in PURINES) or (r in PYR and alt in PYR)

    recs=[]  # (contig, pos, rsid, ref, alts_list, gt_string, phased_bool)
    indels=[]
    with open_text(a.infile) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            f=line.rstrip('\n\r').split('\t')
            if a.phased:
                if len(f)<5: continue
                rsid,chrom,pos,al1,al2=f[0],f[1],f[2],f[3].strip(),f[4].strip()
                geno=al1+al2
            else:
                if len(f)<4: continue
                rsid,chrom,pos,geno=f[0],f[1],f[2],f[3].strip()
            stats['total']+=1
            if posmap is not None:
                mp=posmap.get(rsid)
                if mp is None:
                    stats['unliftable']+=1; continue
                ct=contig(mp[0]); p=mp[1]
                if ct is None:
                    stats['skipped_no_contig']+=1; continue
            else:
                ct=contig(chrom)
                if ct is None:
                    stats['skipped_no_contig']+=1; continue
                try: p=int(pos)
                except: continue
            g=geno.upper()
            # no-call
            if g in ('--','','..','00') or set(g)<= {'-','.','0'}:
                stats['nocall']+=1
                recs.append((ct,p,rsid,None,None,'./.',False)); continue
            # indel placeholders
            if set(g) & set('DI'):
                stats['indel_placeholder']+=1
                indels.append((ct,p,rsid,geno)); continue
            # alleles must be valid bases
            alleles=list(g)
            if not all(b in VALID_BASES for b in alleles):
                stats['bad_allele']+=1; continue
            # reference base at this position
            try:
                ref=fa.fetch(ct,p-1,p).upper()
            except Exception:
                stats['skipped_no_contig']+=1; continue
            if ref=='' or ref not in VALID_BASES:
                stats['ref_mismatch']+=1; continue
            haploid = (len(alleles)==1)
            uniq=set(alleles)
            alts=[b for b in sorted(uniq) if b!=ref]
            # build GT
            if not alts:
                gt='0' if haploid else '0/0'
                if haploid: stats['haploid']+=1
                else: stats['hom_ref']+=1
            else:
                if len(alts)==2:
                    stats['multiallelic']+=1
                    code={ref:'0', alts[0]:'1', alts[1]:'2'}
                    gt='/'.join(sorted((code[b] for b in alleles), key=int))
                    for al in alts:
                        if is_ts(ref,al): stats['transitions']+=1
                        else: stats['transversions']+=1
                else:
                    alt=alts[0]
                    code={ref:'0', alt:'1'}
                    if is_ts(ref,alt): stats['transitions']+=1
                    else: stats['transversions']+=1
                    if haploid:
                        gt='1'; stats['haploid']+=1
                    else:
                        gt='/'.join(sorted((code[b] for b in alleles), key=int))
                        if gt=='1/1': stats['hom_alt']+=1
                        else: stats['het']+=1
            recs.append((ct,p,rsid,ref,alts,gt,a.phased))
    fa.close()

    # write VCF (unsorted; bcftools sort will order). Include all contigs present.
    import collections
    # contig lengths
    fa2=pysam.FastaFile(a.fasta)
    clen={c:fa2.get_reference_length(c) for c in fa2.references}
    fa2.close()
    used_contigs=[]
    seen=set()
    for r in recs:
        if r[0] not in seen: seen.add(r[0]); used_contigs.append(r[0])
    with open(a.out_vcf,'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write(f'##source=23andme_to_vcf;sample={a.sample}\n')
        for c in fa2.references if False else used_contigs:
            out.write(f'##contig=<ID={c},length={clen.get(c,0)}>\n')
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##INFO=<ID=SRC,Number=1,Type=String,Description="Source dataset">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}\n')
        for ct,p,rsid,ref,alts,gt,phased in recs:
            if ref is None:  # no-call — emit with symbolic missing (ref unknown); skip from coordinate VCF
                continue
            altf='.' if not alts else ','.join(alts)
            g=gt.replace('/','|') if phased else gt
            out.write(f'{ct}\t{p}\t{rsid}\t{ref}\t{altf}\t.\tPASS\tSRC={a.sample}\tGT\t{g}\n')
            stats['written']+=1
    with open(a.out_indels,'w') as out:
        out.write('chrom\tpos\trsid\tgenotype\n')
        for row in indels: out.write('\t'.join(map(str,row))+'\n')
    import json
    with open(a.out_stats,'w') as out:
        json.dump(stats,out,indent=2)
    print(json.dumps(stats,indent=2))

if __name__=='__main__':
    main()
