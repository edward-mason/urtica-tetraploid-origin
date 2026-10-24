
"""Serial downloader with MD5 verification. Usage: python code/download.py <plan.json> <outdir>"""
import sys, os, json, hashlib, requests, time, re
plan=json.load(open(sys.argv[1])); out=sys.argv[2]; os.makedirs(out,exist_ok=True)
log=open(os.path.join(out,"download_log.tsv"),"a")
def dl(url,dest):
    if os.path.exists(dest) and os.path.getsize(dest)>0: return "cached"
    tmp=dest+".part"
    with requests.get(url,stream=True,timeout=600) as r:
        r.raise_for_status()
        with open(tmp,"wb") as f:
            for ch in r.iter_content(1<<22): f.write(ch)
    os.replace(tmp,dest); return "downloaded"
def md5(p):
    h=hashlib.md5()
    with open(p,"rb") as f:
        for ch in iter(lambda: f.read(1<<22), b""): h.update(ch)
    return h.hexdigest()
for key,ent in plan.items():
    d=os.path.join(out,key); os.makedirs(d,exist_ok=True)
    sums={}
    if "md5checksums.txt" in ent["files"]:
        dl(ent["dir"]+"md5checksums.txt", os.path.join(d,"md5checksums.txt"))
        for line in open(os.path.join(d,"md5checksums.txt")):
            parts=line.split()
            if len(parts)==2: sums[os.path.basename(parts[1])]=parts[0]
    for f in ent["files"]:
        if f=="md5checksums.txt": continue
        t0=time.time(); dest=os.path.join(d,f)
        status=dl(ent["dir"]+f,dest)
        m=md5(dest); ok = (sums.get(f)==m) if f in sums else "no_md5"
        log.write("\t".join([time.strftime("%Y-%m-%dT%H:%M:%S"),key,f,ent["dir"]+f,str(os.path.getsize(dest)),m,str(ok),status,f"{time.time()-t0:.0f}s"])+"\n"); log.flush()
        print(key,f,os.path.getsize(dest)//1_000_000,"MB md5_ok=",ok,flush=True)
