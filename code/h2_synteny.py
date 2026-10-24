
"""h2_synteny.py — reciprocal best BLASTP hits + collinear blocks between U. dioica (UDI, MH2 RefSeq) and outgroups.
Usage: python code/h2_synteny.py blast <SP>   |  python code/h2_synteny.py blocks <SP>
SP in UUR PJU CAN (RAM optional). Preregistered: evalue 1e-10, >=50% mutual coverage, blocks >=5 RBH pairs within 20-gene windows.
"""
import sys, os, subprocess, glob
import numpy as np, pandas as pd
S="results/synteny"; BIN=glob.glob("code/bin/ncbi-blast-*/bin")[0]
def run(cmd):
    r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: raise RuntimeError(r.stderr[-1500:])
def makedb(sp):
    if not os.path.exists(f"{S}/db_{sp}.psq"): run([f"{BIN}/makeblastdb.exe","-in",f"{S}/{sp}.faa","-dbtype","prot","-out",f"{S}/db_{sp}"])
def blastp(q,db,out):
    if os.path.exists(out): return
    run([f"{BIN}/blastp.exe","-query",f"{S}/{q}.faa","-db",f"{S}/db_{db}","-out",out+".part","-evalue","1e-10","-num_threads","8","-max_target_seqs","5","-max_hsps","1",
         "-outfmt","6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"])
    os.replace(out+".part",out)
def rbh(sp):
    cols="qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen".split()
    a=pd.read_csv(f"{S}/UDI_vs_{sp}.tsv",sep="\t",names=cols); b=pd.read_csv(f"{S}/{sp}_vs_UDI.tsv",sep="\t",names=cols)
    for d in (a,b):
        d["qcov"]=(d.qend-d.qstart+1)/d.qlen; d["scov"]=(d.send-d.sstart+1)/d.slen
    a=a[(a.qcov>=0.5)&(a.scov>=0.5)].sort_values("bitscore",ascending=False).drop_duplicates("qseqid")
    b=b[(b.qcov>=0.5)&(b.scov>=0.5)].sort_values("bitscore",ascending=False).drop_duplicates("qseqid")
    m=a.merge(b,left_on=["qseqid","sseqid"],right_on=["sseqid","qseqid"],suffixes=("","_r"))
    pairs=m[["qseqid","sseqid","pident","bitscore"]].rename(columns={"qseqid":"udi","sseqid":"out"})
    pairs.to_csv(f"{S}/rbh_UDI_{sp}.csv",index=False); return pairs
def blocks(sp, min_pairs=5, window=20):
    pairs=pd.read_csv(f"{S}/rbh_UDI_{sp}.csv")
    g1=pd.read_csv(f"{S}/UDI_genes.csv"); g2=pd.read_csv(f"{S}/{sp}_genes.csv")
    # gene order index per chromosome/scaffold
    g1=g1.sort_values(["seqid","start"]); g1["ord"]=g1.groupby("seqid").cumcount()
    g2=g2.sort_values(["seqid","start"]); g2["ord"]=g2.groupby("seqid").cumcount()
    p=pairs.merge(g1[["pid","seqid","chr","ord","start"]],left_on="udi",right_on="pid").merge(g2[["pid","seqid","chr","ord","start"]],left_on="out",right_on="pid",suffixes=("_u","_o"))
    p=p.sort_values(["seqid_u","ord_u"]).reset_index(drop=True)
    # chain: walk pairs in UDI gene order; a pair joins the open block of the same (UDI scaffold, outgroup scaffold)
    # if it lies within `window` genes of that block's last member on BOTH genomes. Linear-time with per-key open blocks.
    su=p.seqid_u.values; so=p.seqid_o.values; ou=p.ord_u.values; oo=p.ord_o.values
    blk=np.full(len(p),-1); bid=0; open_blocks={}   # key (su,so) -> (block id, last ord_u, last ord_o)
    for i in range(len(p)):
        key=(su[i],so[i]); ob=open_blocks.get(key)
        if ob is not None and ou[i]-ob[1]<=window and abs(oo[i]-ob[2])<=window:
            blk[i]=ob[0]; open_blocks[key]=(ob[0],ou[i],oo[i])
        else:
            blk[i]=bid; open_blocks[key]=(bid,ou[i],oo[i]); bid+=1
    p["block"]=blk
    sizes=p.groupby("block").size(); p["block_size"]=p.block.map(sizes)
    p.loc[p.block_size<min_pairs,"block"]=-1
    p.to_csv(f"{S}/blocks_UDI_{sp}.csv",index=False)
    kept=p[p.block>=0]
    print(sp,"RBH pairs",len(pairs),"in blocks",len(kept),"blocks",kept.block.nunique())
    return p
if __name__=="__main__":
    cmd,sp=sys.argv[1],sys.argv[2]
    if cmd=="blast":
        makedb("UDI"); makedb(sp)
        blastp("UDI",sp,f"{S}/UDI_vs_{sp}.tsv"); blastp(sp,"UDI",f"{S}/{sp}_vs_UDI.tsv"); print(sp,"rbh",len(rbh(sp)))
    elif cmd=="blocks": blocks(sp)
