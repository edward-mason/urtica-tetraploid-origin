
"""gene_divergence.py — common gene-space divergence pipeline (preregistration 2026-10-23).
Usage: python code/gene_divergence.py prepare | blast <TARGET> | summarise
Targets: FH1 FH2 MH1 MH2 DP DA UU
"""
import sys, os, json, re, gzip, time, collections
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from nettle_lib import *

G="data/genomes"; W=os.environ.get("NETTLE_OUTDIR","results/divergence"); os.makedirs(W,exist_ok=True)
TRIM=int(os.environ.get("NETTLE_TRIM","10")); BLASTDIR="results/divergence"
TARGETS={
 "FH1":("FH1/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_genomic.fna.gz","FH1/GCA_048565285.1_U.dioica_female_clone11-4_v1_H1_assembly_report.txt"),
 "FH2":("FH2/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_genomic.fna.gz","FH2/GCA_048565295.1_U.dioica_female_clone11-4_v1_H2_assembly_report.txt"),
 "MH1":("MH1/GCA_054948195.1_U.dioica_male_clone48M_v1_H1_genomic.fna.gz","MH1/GCA_054948195.1_U.dioica_male_clone48M_v1_H1_assembly_report.txt"),
 "MH2":("MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.fna.gz","MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_assembly_report.txt"),
 "DP": ("DP/GCA_964188135.1_drUrtDioi1.1_genomic.fna.gz","DP/GCA_964188135.1_drUrtDioi1.1_assembly_report.txt"),
 "DA": ("DA/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_genomic.fna.gz","DA/GCA_964188125.1_drUrtDioi1.1_alternate_haplotype_assembly_report.txt"),
 "UU": ("UU/GCA_958296335.1_drUrtUren1.1_genomic.fna.gz","UU/GCA_958296335.1_drUrtUren1.1_assembly_report.txt"),
 "BN": ("BN/Urtica_dioica.v1.genome.fasta.gz",None),
}
MATERNAL={1:"MH1",2:"MH1",3:"MH2",4:"MH1",5:"MH2",6:"MH2",7:"MH1",8:"MH1",9:"MH1",10:"MH1",11:"MH1",12:"MH1",13:"MH1"}  # Hirabayashi 2026 Table 2

def prepare():
    """Longest CDS per protein-coding gene from the MH2 RefSeq cds_from_genomic file; length>=300, no internal N."""
    cds=f"{G}/MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_cds_from_genomic.fna.gz"
    rep=assembly_report_map(f"{G}/{TARGETS['MH2'][1]}")
    best={}
    for name,seq in read_fasta(cds):
        if "pseudo=true" in name: continue
        gene=re.search(r"\[gene=([^\]]+)\]",name); prot=re.search(r"\[protein_id=([^\]]+)\]",name)
        loc=re.search(r"\[location=([^\]]+)\]",name); acc=name.split()[0]
        if not (gene and prot and loc): continue
        gid=gene.group(1); L=len(seq)
        if L<300 or L%3!=0 or "N" in seq.upper(): continue
        if gid not in best or L>best[gid][1]:
            sid=acc.split("|")[1] if "|" in acc else acc.split("_cds_")[0].replace("lcl|","")
            sid=re.sub(r"_cds_.*","",sid)
            m=re.search(r"(\d+)\.\.(\d+)",loc.group(1)); start=int(m.group(1)) if m else None
            strand="-" if "complement" in loc.group(1) else "+"
            best[gid]=(seq.upper(),L,prot.group(1),sid,start,strand)
    rows=[]
    with open(f"{W}/query_cds.fa","w") as f:
        for gid,(seq,L,pid,sid,start,strand) in best.items():
            qid=pid
            f.write(f">{qid}\n{seq}\n")
            sn,mol,role=rep.get(sid,(sid,None,None))
            rows.append((qid,gid,sid,chrom_number(sn,mol),start,strand,L))
    df=pd.DataFrame(rows,columns=["qid","gene","mh2_seq","mh2_chr","mh2_start","mh2_strand","cds_len"])
    df.to_csv(f"{W}/query_index.csv",index=False)
    print("queries",len(df),"on chromosomes",df.mh2_chr.notna().sum())

def blast_target(tg):
    fa,rep=TARGETS[tg]
    plain=gunzip_to(f"{G}/{fa}", f"{G}/{tg}.fna")
    db=makeblastdb(plain, f"{G}/blastdb/{tg}", tg)
    t0=time.time()
    out=blastn(f"{W}/query_cds.fa", db, f"{W}/blast_{tg}.tsv", task="dc-megablast", threads=8, max_targets=5, max_hsps=40)
    print(tg,"blast done",f"{time.time()-t0:.0f}s", os.path.getsize(out)//1_000_000,"MB")

def iter_query_groups(path):
    """Stream BLAST tabular output grouped by consecutive qseqid (BLAST writes queries contiguously)."""
    cur=None; buf=[]
    with open(path) as f:
        for line in f:
            q=line[:line.index("\t")]
            if q!=cur:
                if buf: yield cur, buf
                cur=q; buf=[]
            buf.append(line.rstrip("\n").split("\t"))
    if buf: yield cur, buf

def summarise_target(tg, qindex):
    """Per query: best locus, 2nd locus ratio, coverage, anchored sequence."""
    fa,rep=TARGETS[tg]; repm=assembly_report_map(f"{G}/{rep}") if rep else {}
    path=f"{BLASTDIR}/blast_{tg}.tsv"
    rows=[]; anchored={}
    for q, recs in iter_query_groups(path):
        hs=parse_recs(recs); qlen=hs[0]["qlen"]
        loci=cluster_loci(hs)
        if not loci: continue
        b=loci[0]; s2=loci[1]["score"] if len(loci)>1 else 0.0
        sn,mol,role=repm.get(b["sseqid"],(b["sseqid"],None,None))
        c2=loci[1]["cov"] if len(loci)>1 else None; i2=loci[1]["pident"] if len(loci)>1 else None
        rows.append((q,b["sseqid"],chrom_number(sn,mol),b["smin"],b["smax"],b["strand"],b["score"],b["cov"],s2/b["score"] if b["score"] else np.nan,c2,i2,len(loci),b["pident"]))
        anchored[q]=anchored_sequence(b["hsps"],qlen,trim=TRIM)
    res=pd.DataFrame(rows,columns=["qid","sseqid","chr","smin","smax","strand","score","cov","ratio2","cov2","pident2","nloci","pident"])
    res.to_csv(f"{W}/loci_{tg}.csv",index=False)
    with gzip.open(f"{W}/anchored_{tg}.json.gz","wt") as f: json.dump(anchored,f)
    print(tg,"queries with hits",len(res),"single-copy&cov>=0.8:",int(((res["ratio2"]<0.8)&(res["cov"]>=0.8)).sum()),flush=True)

if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="prepare": prepare()
    elif cmd=="blast": blast_target(sys.argv[2])
    elif cmd=="summarise":
        qi=pd.read_csv("results/divergence/query_index.csv")
        for tg in sys.argv[2:]: summarise_target(tg,qi)
