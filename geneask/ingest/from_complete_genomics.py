#!/usr/bin/env python3
"""Complete Genomics var-ASM.tsv  ->  GRCh37 VCF (FASTA-anchored, bcftools-normalizable).

CG columns (0-based idx):
 0 locus 1 ploidy 2 allele 3 chrom 4 begin 5 end 6 varType 7 reference
 8 alleleSeq 9 varScoreVAF 10 varScoreEAF 11 varFilter 12 hapLink 13 xRef
 14 alleleFreq 15 alternativeCalls
Coords are 0-based half-open. Strategy: group rows by locus; skip ref/no-call/
no-ref/PAR; reconstruct each haplotype over the union span [B,E] from the FASTA
(seq = fasta[B:begin] + alleleSeq + fasta[end:E]); anchor indels on base B-1;
emit REF/ALT/GT; bcftools norm left-aligns & trims downstream.
"""
import sys, argparse, bz2, gzip, json, re
import pysam

SKIP = {'ref','no-call','no-ref','PAR-called-in-X','no-call-rc','no-call-ri'}
RS = re.compile(r'rs\d+')

def opener(p):
    if p.endswith('.bz2'): return bz2.open(p,'rt')
    if p.endswith('.gz'):  return gzip.open(p,'rt')
    return open(p,'r',errors='replace')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-stats", required=True)
    a=ap.parse_args()

    fa=pysam.FastaFile(a.fasta)
    refset=set(fa.references)
    clen={c:fa.get_reference_length(c) for c in fa.references}
    def contig(c):
        if c in refset: return c
        if c.startswith('chr'):
            b=c[3:]
            b='MT' if b=='M' else b
            if b in refset: return b
        if ('chr'+c) in refset: return 'chr'+c
        return None

    st=dict(loci_total=0, variant_loci=0, written=0, snp=0, ins=0, del_=0, sub=0,
            mnp=0, complex_=0, skipped_type=0, ref_check_fail=0, nocall_allele=0,
            multiallelic=0, hom_alt=0, het=0, hemi=0, with_rsid=0, contig_skip=0,
            span_mismatch=0)

    out=open(a.out_vcf,'w')
    out.write('##fileformat=VCFv4.2\n')
    out.write(f'##source=cg_var_to_vcf;assembly=GS000038910;sample={a.sample}\n')
    out.write('##reference=GRCh37\n')
    for c in fa.references:
        out.write(f'##contig=<ID={c},length={clen[c]}>\n')
    out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    out.write('##INFO=<ID=CGT,Number=1,Type=String,Description="CG varType(s) at locus">\n')
    out.write('##INFO=<ID=CGF,Number=.,Type=String,Description="CG per-allele filter(s)">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}\n')

    def flush(locus_rows):
        # locus_rows: list of column-lists for one locus
        if not locus_rows: return
        st['loci_total']+=1
        vtypes={r[6] for r in locus_rows}
        # skip pure ref/no-call loci
        if vtypes <= SKIP: return
        chrom=locus_rows[0][3]
        ct=contig(chrom)
        if ct is None:
            st['contig_skip']+=1; return
        # collect allele rows keyed by allele id; ignore ref-only sub-rows
        byal={}
        for r in locus_rows:
            al=r[2]
            byal.setdefault(al,[]).append(r)
        st['variant_loci']+=1

        # Determine variant rows (non-skip). 'all' means both haplotypes.
        var_rows=[r for r in locus_rows if r[6] not in SKIP]
        if not var_rows: return
        # union span
        begins=[int(r[4]) for r in var_rows]; ends=[int(r[5]) for r in var_rows]
        B=min(begins); E=max(ends)
        if B<0: return
        # need indel anchor if any allele changes length or is zero-width
        need_anchor=False
        for r in var_rows:
            b,e=int(r[4]),int(r[5]); vt=r[6]
            altseq=r[8]
            reflen=e-b
            if vt in ('ins','del','sub') or len(altseq)!=reflen:
                need_anchor=True
        # collect reference over span from FASTA
        try:
            ref_span=fa.fetch(ct,B,E).upper()
        except Exception:
            st['contig_skip']+=1; return
        anchor=''
        if need_anchor:
            if B-1<0: need_anchor=False
            else:
                anchor=fa.fetch(ct,B-1,B).upper()

        # Build haplotype sequences. Determine which allele ids are variant.
        # Map allele '1','2' -> reconstructed seq over span (or ref if that hap
        # has no variant row => reference). 'all' applies to both.
        hap={}   # allele_id -> seq or None(=nocall)
        present_ids=set()
        has_all=False
        for r in var_rows:
            al=r[2]
            if al=='all': has_all=True
            present_ids.add(al)
        def recon(r):
            b,e=int(r[4]),int(r[5]); altseq=r[8]
            left=fa.fetch(ct,B,b).upper() if b>B else ''
            right=fa.fetch(ct,e,E).upper() if e<E else ''
            return left+altseq.upper()+right
        # per-allele reconstruction; last variant row per allele wins (usually one)
        recon_by_al={}
        for r in var_rows:
            recon_by_al[r[2]]=recon(r)
        # no-call alleles at this locus
        nocall_ids={r[2] for r in locus_rows if r[6] in ('no-call','no-ref','no-call-rc','no-call-ri')}

        # Assemble two haplotypes h1,h2 (or single if 'all'/haploid)
        if has_all:
            s=recon_by_al.get('all')
            haps=[s,s]
        else:
            def hap_seq(idx):
                if idx in recon_by_al: return recon_by_al[idx]
                if idx in nocall_ids: return None    # no-call
                return ref_span                       # reference haplotype
            haps=[hap_seq('1'), hap_seq('2')]
        # if only allele '1' variant exists and no '2' row at all -> hap2 = ref
        # (already handled: hap_seq('2') returns ref_span)

        # Build allele list with anchor
        def full(seq):
            if seq is None: return None
            return (anchor+seq) if need_anchor else seq
        REF=(anchor+ref_span) if need_anchor else ref_span
        if REF=='' :
            return
        alleles=[REF]
        gt=[]
        nocall_flag=False
        for h in haps:
            fs=full(h)
            if fs is None:
                gt.append('.'); nocall_flag=True; continue
            if fs==REF:
                gt.append('0')
            else:
                if fs not in alleles:
                    alleles.append(fs)
                gt.append(str(alleles.index(fs)))
        alts=alleles[1:]
        if not alts:
            return  # no real variant (all ref) — skip
        POS = (B if need_anchor else B+1)  # anchor base 1-based = B; else first span base 1-based = B+1
        # counts
        if len(alts)>1: st['multiallelic']+=1
        # classify
        vt_join=';'.join(sorted(vtypes-SKIP))
        for vt in vtypes:
            if vt=='snp': st['snp']+=1; break
        if 'ins' in vtypes: st['ins']+=1
        if 'del' in vtypes: st['del_']+=1
        if 'sub' in vtypes: st['sub']+=1
        # genotype string
        if has_all:
            gts=gt[0]  # hemizygous/haploid
            st['hemi']+=1
        else:
            # normalize allele order for unphased GT (0/1 not 1/0); keep './x' forms
            if '.' not in gt:
                gt_sorted=sorted(gt, key=int)
            else:
                # put called allele consistently: missing last -> 'x/.'
                gt_sorted=sorted(gt, key=lambda x:(x=='.', x!='.' and int(x) or 0))
            gts='/'.join(gt_sorted)
            if gts in ('1/1','2/2'): st['hom_alt']+=1
            elif '.' in gts: st['nocall_allele']+=1
            else: st['het']+=1
        # rsid
        xref=';'.join(r[13] for r in var_rows if len(r)>13 and r[13])
        m=RS.findall(xref)
        rsid=m[0] if m else '.'
        if rsid!='.': st['with_rsid']+=1
        # join multiple filter values with ',' (';' is the INFO separator);
        # also replace any internal ';' inside a single filter value.
        filt=','.join(sorted({r[11].replace(';',',') for r in var_rows if len(r)>11 and r[11]}))
        vt_join=vt_join.replace(';',',')
        info=f"CGT={vt_join}"
        if filt: info+=f";CGF={filt}"
        out.write(f"{ct}\t{POS}\t{rsid}\t{REF}\t{','.join(alts)}\t.\tPASS\t{info}\tGT\t{gts}\n")
        st['written']+=1

    # stream, grouping consecutive rows by locus id (col 0)
    cur=None; buf=[]
    with opener(a.infile) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('>') or not line.strip():
                continue
            r=line.rstrip('\n\r').split('\t')
            if len(r)<9: continue
            loc=r[0]
            if loc!=cur:
                if buf: flush(buf)
                cur=loc; buf=[r]
            else:
                buf.append(r)
        if buf: flush(buf)
    out.close(); fa.close()
    st['del']=st.pop('del_'); st['sub_']=st.get('sub',0); st['complex']=st.pop('complex_',0)
    json.dump(st, open(a.out_stats,'w'), indent=2)
    print(json.dumps(st,indent=2))

if __name__=='__main__':
    main()
