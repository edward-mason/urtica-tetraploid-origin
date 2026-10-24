"""prep_proteomes.py — longest protein per gene and gene-position tables for the H2 synteny analysis.
Inputs: NCBI GFF/protein FASTA (U. dioica MH2 RefSeq, Cannabis cs10, Morus notabilis [scaffold-level; not used downstream]) and Ensembl Rapid Release GFF3/pep (U. urens, Parietaria judaica), all under data/genomes/.
Outputs: results/synteny/<TAG>.faa and <TAG>_genes.csv for TAG in UDI, CAN, MOR, UUR, PJU.
Reconstructed from the session notebook cells of 2026-10-23 (the step was first run inline)."""

import gzip, re, os, collections, pandas as pd
G="data/genomes"
os.makedirs("results/synteny",exist_ok=True)
def read_fa(p):
    name=None; seq=[]
    with gzip.open(p,"rt") as f:
        for l in f:
            if l.startswith(">"):
                if name: yield name,"".join(seq)
                name=l[1:].strip(); seq=[]
            else: seq.append(l.strip())
    if name: yield name,"".join(seq)
def rep_map(path):
    m={}
    for line in open(path,encoding="utf-8",errors="ignore"):
        if line.startswith("#"): continue
        p=line.rstrip("\n").split("\t")
        if len(p)>=7:
            for acc in (p[4],p[6]):
                if acc!="na": m[acc]=(p[0],p[2])
    return m
def chrnum(sn,mol):
    for s in (mol,sn):
        mm=re.fullmatch(r"(?:chr|Chr|SUPER_|LG)?0*(\d{1,2})",str(s))
        if mm: return int(mm.group(1))
    return None

def ncbi_species(tag, gff, faa, rep):
    """NCBI GFF: mRNA/CDS with protein_id; gene position from 'gene' features."""
    rm=rep_map(rep)
    genes={}; prot2gene={}
    with gzip.open(gff,"rt") as f:
        for l in f:
            if l.startswith("#"): continue
            p=l.rstrip("\n").split("\t")
            if len(p)<9: continue
            if p[2]=="gene":
                gid=re.search(r"ID=([^;]+)",p[8]).group(1)
                sn,mol=rm.get(p[0],(p[0],None))
                genes[gid]=(p[0],chrnum(sn,mol),int(p[3]),int(p[4]),p[6])
            elif p[2]=="CDS":
                pid=re.search(r"protein_id=([^;]+)",p[8]); par=re.search(r"Parent=([^;]+)",p[8])
                if pid: prot2gene.setdefault(pid.group(1), None)
                # map through mRNA parent later
    # second pass: mRNA -> gene, CDS -> mRNA
    mrna2gene={}; cds2mrna={}
    with gzip.open(gff,"rt") as f:
        for l in f:
            if l.startswith("#"): continue
            p=l.rstrip("\n").split("\t")
            if len(p)<9: continue
            if p[2] in ("mRNA","transcript"):
                mrna2gene[re.search(r"ID=([^;]+)",p[8]).group(1)]=re.search(r"Parent=([^;]+)",p[8]).group(1)
            elif p[2]=="CDS":
                pid=re.search(r"protein_id=([^;]+)",p[8]); par=re.search(r"Parent=([^;]+)",p[8])
                if pid and par:
                    g=mrna2gene.get(par.group(1), par.group(1))
                    prot2gene[pid.group(1)]=g
    best={}
    for name,seq in read_fa(faa):
        pid=name.split()[0]; g=prot2gene.get(pid)
        if g is None or g not in genes: continue
        if g not in best or len(seq)>len(best[g][1]): best[g]=(pid,seq)
    rows=[]
    with open(f"results/synteny/{tag}.faa","w") as out:
        for g,(pid,seq) in best.items():
            sid,chrn,s,e,strand=genes[g]
            out.write(f">{pid}\n{seq.rstrip('*')}\n"); rows.append((pid,g,sid,chrn,s,e,strand))
    df=pd.DataFrame(rows,columns=["pid","gene","seqid","chr","start","end","strand"]); df.to_csv(f"results/synteny/{tag}_genes.csv",index=False)
    print(tag,len(df),"on chr:",df.chr.notna().sum())
ncbi_species("UDI",f"{G}/MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_genomic.gff.gz",f"{G}/MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_protein.faa.gz",f"{G}/MH2/GCF_054948205.1_U.dioica_male_clone48M_v1_H2_assembly_report.txt")
ncbi_species("MOR",f"{G}/MOR/GCF_000414095.1_ASM41409v2_genomic.gff.gz",f"{G}/MOR/GCF_000414095.1_ASM41409v2_protein.faa.gz",f"{G}/MOR/GCF_000414095.1_ASM41409v2_assembly_report.txt")
ncbi_species("CAN",f"{G}/CAN/GCF_900626175.1_cs10_genomic.gff.gz",f"{G}/CAN/GCF_900626175.1_cs10_protein.faa.gz",f"{G}/CAN/GCF_900626175.1_cs10_assembly_report.txt")
def ensembl_species(tag, gff, pep, rep):
    rm=rep_map(rep); genes={}
    with gzip.open(gff,"rt") as f:
        for l in f:
            if l.startswith("#"): continue
            p=l.rstrip("\n").split("\t")
            if len(p)<9 or p[2]!="gene": continue
            gid=re.search(r"ID=gene:([^;]+)",p[8]).group(1).split(".")[0]
            sn,mol=rm.get(p[0],(p[0],None))
            genes[gid]=(p[0],chrnum(sn,mol) or chrnum(p[0],None),int(p[3]),int(p[4]),p[6])
    best={}
    for name,seq in read_fa(pep):
        pid=name.split()[0]; gm=re.search(r"gene:(\S+)",name)
        g=gm.group(1).split(".")[0] if gm else None
        if g is None or g not in genes: continue
        if g not in best or len(seq)>len(best[g][1]): best[g]=(pid,seq)
    rows=[]
    with open(f"results/synteny/{tag}.faa","w") as out:
        for g,(pid,seq) in best.items():
            sid,chrn,s,e,strand=genes[g]; out.write(f">{pid}\n{seq.rstrip('*')}\n"); rows.append((pid,g,sid,chrn,s,e,strand))
    df=pd.DataFrame(rows,columns=["pid","gene","seqid","chr","start","end","strand"]); df.to_csv(f"results/synteny/{tag}_genes.csv",index=False)
    print(tag,len(df),"on chr:",df.chr.notna().sum(), df.chr.value_counts().sort_index().to_dict())
ensembl_species("UUR",f"{G}/UU_ens/Urtica_urens-GCA_958296335.1-2024_02-genes.gff3.gz",f"{G}/UU_ens/Urtica_urens-GCA_958296335.1-2024_02-pep.fa.gz",f"{G}/UU/GCA_958296335.1_drUrtUren1.1_assembly_report.txt")
ensembl_species("PJU",f"{G}/PJ_ens/Parietaria_judaica-GCA_954871525.1-2024_02-genes.gff3.gz",f"{G}/PJ_ens/Parietaria_judaica-GCA_954871525.1-2024_02-pep.fa.gz",f"{G}/PJ/GCA_954871525.1_drParJuda1.1_assembly_report.txt")
