
"""nettle_lib.py — shared helpers for the U. dioica comparative pipeline (Windows-native, BLAST+ 2.17.0)."""
import gzip, os, re, subprocess, sys, json, collections
import numpy as np, pandas as pd

BLASTBIN = None
def blast_bin():
    global BLASTBIN
    if BLASTBIN is None:
        import glob
        BLASTBIN = glob.glob(os.path.join("code","bin","ncbi-blast-*","bin"))[0]
    return BLASTBIN

def open_any(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def read_fasta(p):
    name=None; seq=[]
    with open_any(p) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name=line[1:].strip(); seq=[]
            else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

def assembly_report_map(path):
    """GenBank/RefSeq accession -> (sequence_name, assigned_molecule, role)."""
    m={}
    with open_any(path) as f:
        for line in f:
            if line.startswith("#"): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<7: continue
            seqname, role, mol, loctype, gb, rel, rs = p[0],p[1],p[2],p[3],p[4],p[5],p[6]
            for acc in (gb, rs):
                if acc and acc!="na": m[acc]=(seqname, mol, role)
    return m

def chrom_number(seqname, mol):
    """Return integer chromosome number if this is a chromosome, else None."""
    for s in (mol, seqname):
        if s is None: continue
        mm=re.fullmatch(r"(?:chr|Chr|CHR|SUPER_|chromosome_?)?0*(\d{1,2})", str(s))
        if mm: return int(mm.group(1))
    return None

def gunzip_to(src, dst):
    if os.path.exists(dst): return dst
    with gzip.open(src,"rb") as fi, open(dst+".part","wb") as fo:
        while True:
            b=fi.read(1<<24)
            if not b: break
            fo.write(b)
    os.replace(dst+".part", dst); return dst

def makeblastdb(fasta, dbprefix, title):
    if os.path.exists(dbprefix+".nsq") or os.path.exists(dbprefix+".00.nsq"): return dbprefix
    cmd=[os.path.join(blast_bin(),"makeblastdb.exe"),"-in",fasta,"-dbtype","nucl","-out",dbprefix,"-title",title,"-parse_seqids","-max_file_sz","4GB"]
    r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: raise RuntimeError(r.stderr[-2000:])
    return dbprefix

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qseq sseq"
def blastn(query, db, out, task="dc-megablast", threads=8, evalue="1e-20", max_targets=10, max_hsps=60, extra=None):
    if os.path.exists(out) and os.path.getsize(out)>0: return out
    cmd=[os.path.join(blast_bin(),"blastn.exe"),"-task",task,"-query",query,"-db",db,"-out",out+".part","-outfmt",OUTFMT,
         "-evalue",evalue,"-num_threads",str(threads),"-max_target_seqs",str(max_targets),"-max_hsps",str(max_hsps),"-dust","no","-soft_masking","false"]
    if extra: cmd+=extra
    r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: raise RuntimeError(r.stderr[-2000:])
    os.replace(out+".part",out); return out

COLS=OUTFMT.split()[1:]
_NUM={"pident":float,"length":int,"mismatch":int,"gapopen":int,"qstart":int,"qend":int,"sstart":int,"send":int,"evalue":float,"bitscore":float,"qlen":int,"slen":int}
_IDX={c:i for i,c in enumerate(COLS)}
def parse_recs(recs):
    """List of split BLAST rows -> list of dicts with numeric fields converted."""
    out=[]
    for r in recs:
        d={c:r[i] for c,i in _IDX.items()}
        if "|" in d["sseqid"]: d["sseqid"]=[x for x in d["sseqid"].split("|") if x][-1]
        for c,t in _NUM.items(): d[c]=t(d[c])
        d["smin"]=min(d["sstart"],d["send"]); d["smax"]=max(d["sstart"],d["send"])
        out.append(d)
    return out

def cluster_loci(hsps, maxgap=50000):
    """hsps: list of dicts (one query). Returns loci sorted by score desc: {sseqid,strand,hsps(list of dicts),score,cov,smin,smax}."""
    groups=collections.defaultdict(list)
    for h in hsps: groups[(h["sseqid"],h["sstrand"])].append(h)
    loci=[]
    for (sid,strand),g in groups.items():
        g.sort(key=lambda h:h["smin"]); cur=[]; last=None
        for h in g:
            if last is not None and h["smin"]-last>maxgap:
                loci.append((sid,strand,cur)); cur=[]; last=None
            cur.append(h); last=h["smax"] if last is None else max(last,h["smax"])
        if cur: loci.append((sid,strand,cur))
    out=[]
    for sid,strand,g in loci:
        qlen=g[0]["qlen"]; covered=np.zeros(qlen+1,dtype=bool); keep=[]
        for h in sorted(g,key=lambda h:-h["bitscore"]):
            a,b=h["qstart"],h["qend"]
            if covered[a:b+1].mean()<=0.2:
                keep.append(h); covered[a:b+1]=True
        out.append({"sseqid":sid,"strand":strand,"hsps":keep,"score":sum(h["bitscore"] for h in keep),"cov":float(covered[1:].mean()),
                    "smin":min(h["smin"] for h in keep),"smax":max(h["smax"] for h in keep),"pident":float(np.mean([h["pident"] for h in keep]))})
    out.sort(key=lambda x:-x["score"])
    return out

def anchored_sequence(hsps, qlen, trim=10):
    """Query-anchored target string: base aligned to each query position; '-' deletion; 'N' uncovered.
    The first/last `trim` aligned query positions of every HSP are discarded unless the HSP reaches the CDS end
    (calibrated on the MH2 self-search: exon-boundary over-extension artefacts)."""
    arr=["N"]*qlen
    for h in hsps:
        q=h["qseq"]; s=h["sseq"].upper(); p=h["qstart"]-1; cols=[]
        for qc,sc in zip(q,s):
            if qc=="-": continue
            cols.append((p,sc)); p+=1
        lo = 0 if h["qstart"]==1 else trim
        hi = len(cols) if h["qend"]==qlen else len(cols)-trim
        for p_,sc in cols[lo:hi]:
            if arr[p_]=="N": arr[p_]=sc
    return "".join(arr)

def pdist(a,b):
    """p-distance over columns where both are ACGT. a,b: numpy <U1 arrays. Returns (mismatches, columns, indel_columns)."""
    A=np.isin(a,list("ACGT")); B=np.isin(b,list("ACGT"))
    both=A&B; n=int(both.sum())
    mm=int((a[both]!=b[both]).sum())
    ind=int(((a=="-")&B).sum()+((b=="-")&A).sum())
    return mm,n,ind
