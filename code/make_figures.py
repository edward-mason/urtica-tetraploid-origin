
"""make_figures.py — regenerate Fig1–Fig5 from result files. Requires matplotlib; uses plain rcParams if the figure-style helpers are absent."""
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, os
try:
    apply_figure_style(); PL=panel_letter
except NameError:
    mpl.rcParams.update({"font.size":7,"axes.titlesize":7,"axes.labelsize":7,"xtick.labelsize":6,"ytick.labelsize":6,"axes.spines.top":False,"axes.spines.right":False})
    def PL(ax,l): ax.text(-0.15,1.08,l,transform=ax.transAxes,fontweight="bold",fontsize=9)
os.makedirs("paper/figures",exist_ok=True)
df=pd.read_csv("results/divergence/gene_distances.csv")
dDP=pd.read_csv("results/divergence/h1v2_genes_DP.csv"); dDA=pd.read_csv("results/divergence/h1v2_genes_DA.csv")
dDP=dDP[~dDP.paralog_flag.astype(bool)]; dDA=dDA[~dDA.paralog_flag.astype(bool)]
vDP=json.load(open("results/divergence/h1v2_result_DP.json")); vDA=json.load(open("results/divergence/h1v2_result_DA.json"))
rDP=dict(vDP["paralog_filtered"],R_CI=vDP["CI_chromosomes"]["R_CI"],f99_CI=vDP["CI_chromosomes"]["f99_CI"]); rDA=dict(vDA["paralog_filtered"],R_CI=vDA["CI_chromosomes"]["R_CI"],f99_CI=vDA["CI_chromosomes"]["f99_CI"])
bn=json.load(open("results/divergence/h1_BN_bestlocus.json")); cs=json.load(open("results/divergence/copy_structure.json")); cs2=json.load(open("results/tables/copy_structure_v2.json"))
for k,v in cs2.items():
    if isinstance(v,dict) and "frac_1" in v: cs[k]=v
okb=np.ones(len(df),bool)
for t in ["FH1","FH2","MH1","MH2"]: okb&=(df[f"{t}_cov"]>=0.8)&(df[f"{t}_ratio2"]<0.8)
base=df[okb&df.mh2_chr.notna()].copy(); base["b"]=np.fmin(base.d_PAT_FH1,base.d_PAT_FH2); base["W"]=base.d_FH1_FH2
uu=df.loc[okb,"d_FH1_UU"].dropna(); bnok=okb&(df.BN_cov>=0.8); tBN=np.fmin(df.loc[bnok,"d_BN_FH1"],df.loc[bnok,"d_BN_FH2"]).dropna()
bins=np.logspace(-3.3,-0.5,60); cols=["#c0392b","#e67e22","#8e44ad"]
series=[("Within diploid individual (FH1 vs FH2)",base.W,"#9e9e9e"),("Between diploid populations (Bulgaria vs Romania)",base.b,"#4d4d4d"),
        ("Tetraploid Kew, primary (DToL) vs diploid",dDP.t,cols[0]),("Tetraploid Kew, alternate vs diploid",dDA.t,cols[1]),("Bonn assembly (ploidy unresolved) vs diploid",tBN,cols[2]),("U. urens vs U. dioica",uu,"#2471a3")]
# ---- Fig 1
fig,axes=plt.subplots(1,3,figsize=(7.6,3.0),gridspec_kw=dict(width_ratios=[1.6,1,1],wspace=0.6)); ax=axes[0]
for lab,x,c in series:
    x=x[x>0]; ax.hist(np.log10(x),bins=np.log10(bins),histtype="step",density=True,color=c,lw=1.4 if "Kew, primary" in lab else 1.0,label=lab)
ax.set_xlabel("CDS p-distance per gene (log scale; zeros omitted)"); ax.set_ylabel("density of genes"); ax.set_xticks([-3,-2.5,-2,-1.5,-1]); ax.set_xticklabels(["0.1%","0.3%","1%","3%","10%"]); ax.set_xlim(-3.4,-0.4)
ax.set_title("Tetraploid genes diverge from the diploid\nno more than diploid populations do",loc="left"); ax.legend(frameon=False,fontsize=5.5,loc="upper left",handlelength=1.2,labelspacing=0.35); ax.set_ylim(0,ax.get_ylim()[1]*1.45)
labs2=["Kew\nprimary","Kew\nalternate\n(held-out)","Bonn"]
R=[rDP["R"],rDA["R"],bn["R_BN_bestlocus"]]; lo=[rDP["R_CI"][0],rDA["R_CI"][0],np.nan]; hi=[rDP["R_CI"][1],rDA["R_CI"][1],np.nan]
ax=axes[1]
for i,(r_,l,h,c) in enumerate(zip(R,lo,hi,cols)): ax.errorbar(i,r_,yerr=[[r_-l] if l==l else [0],[h-r_] if h==h else [0]],fmt="o",color=c,capsize=3,ms=5)
ax.axhspan(0.67,1.5,color="#e0e0e0",alpha=0.6,lw=0); ax.text(2.4,0.70,"same-lineage range\n(preregistered)",ha="right",va="bottom",fontsize=5.5,color="#555555")
ax.axhline(2,color="#555555",ls=":",lw=0.8); ax.text(2.4,2.04,"allotetraploid criterion (R ≥ 2)",ha="right",va="bottom",fontsize=5.5,color="#555555")
ax.set_xticks(range(3)); ax.set_xticklabels(labs2,fontsize=6); ax.set_ylabel("R (tetraploid / diploid-baseline\nmedian divergence)"); ax.set_ylim(0.5,2.3); ax.set_xlim(-0.5,2.5); ax.set_title("R ≈ 1.1 in both independent\ntetraploid haplotypes",loc="left")
F=[rDP["f99"],rDA["f99"],bn["f99_BN"]]; flo=[rDP["f99_CI"][0],rDA["f99_CI"][0],np.nan]; fhi=[rDP["f99_CI"][1],rDA["f99_CI"][1],np.nan]
ax=axes[2]
for i,(f_,l,h,c) in enumerate(zip(F,flo,fhi,cols)): ax.errorbar(i,100*f_,yerr=[[100*(f_-l)] if l==l else [0],[100*(h-f_)] if h==h else [0]],fmt="o",color=c,capsize=3,ms=5)
ax.axhline(1,color="#9e9e9e",ls="--",lw=0.8); ax.text(2.4,0.93,"1% if identical to baseline",ha="right",va="top",fontsize=5.5,color="#555555")
ax.axhline(20,color="#555555",ls=":",lw=0.8); ax.text(2.4,18.5,"allotetraploid criterion (≥ 20%)",ha="right",va="top",fontsize=5.5,color="#555555")
ax.set_yscale("log"); ax.set_ylim(0.5,30); ax.set_yticks([1,3,10,30]); ax.set_yticklabels(["1%","3%","10%","30%"]); ax.set_xticks(range(3)); ax.set_xticklabels(labs2,fontsize=6)
ax.set_ylabel("genes above 99th percentile\nof diploid baseline"); ax.set_xlim(-0.5,2.5); ax.set_title("1.6–1.8% of genes exceed the diploid\nbaseline's 99th percentile",loc="left")
for ax_,l in zip(axes,"abc"): PL(ax_,l)
fig.subplots_adjust(left=0.07,right=0.99,top=0.82,bottom=0.24); fig.savefig("paper/figures/Fig1_tetraploid_origin.png",dpi=300); fig.savefig("paper/figures/Fig1_tetraploid_origin.pdf"); plt.close(fig)
# ---- Fig 2
order=["FH1","FH2","MH1","MH2","FH1+FH2 (unphased diploid, female)","MH1+MH2 (unphased diploid, male)","DP","DA","DP+DA (tetraploid Kew, both pseudo-haplotypes)","BN"]
names={"FH1":"diploid\n♀ H1","FH2":"diploid\n♀ H2","MH1":"diploid\n♂ H1","MH2":"diploid\n♂ H2","FH1+FH2 (unphased diploid, female)":"diploid ♀\nH1+H2\npooled","MH1+MH2 (unphased diploid, male)":"diploid ♂\nH1+H2\npooled","DP":"tetraploid\nKew prim.","DA":"tetraploid\nKew alt.","DP+DA (tetraploid Kew, both pseudo-haplotypes)":"tetraploid\nKew\npooled","BN":"Bonn\n1.1 Gb"}
cats=[("frac_1","1 locus","#bdbdbd"),("frac_2","2 loci","#7f7f7f"),("frac_3","3 loci","#c0392b"),("frac_4plus","≥4 loci","#7b241c")]
fig,axes=plt.subplots(1,2,figsize=(7.6,3.1),gridspec_kw=dict(width_ratios=[1.9,1],wspace=0.35)); ax=axes[0]; x=np.arange(len(order)); bottom=np.zeros(len(order))
for key,lab,c in cats:
    v=np.array([cs[t][key] for t in order])*100; ax.bar(x,v,bottom=bottom,color=c,label=lab,width=0.7,lw=0); bottom+=v
for i,t in enumerate(order): ax.text(i,bottom[i]+1.5,f"{(cs[t]['frac_3']+cs[t]['frac_4plus'])*100:.0f}%",ha="center",va="bottom",fontsize=6,color="#7b241c")
ax.set_xticks(x); ax.set_xticklabels([names[t] for t in order],fontsize=5); ax.set_ylabel("genes (%)"); ax.set_ylim(0,135); ax.set_yticks([0,20,40,60,80,100]); ax.set_xlim(-0.6,len(order)-0.4)
ax.set_title("Bonn holds ≥3 near-identical copies of 24% of genes (red label):\ntwice the pooled-diploid value (12–13%), three-quarters\nof the pooled-tetraploid value (32%)",loc="left")
ax.legend(frameon=False,fontsize=5.5,loc="upper left",ncol=4,handlelength=1.0,columnspacing=0.8,title="loci per gene at ≥97% identity, ≥80% coverage",title_fontsize=5.5)
ax=axes[1]; okdd=okb&(df.DP_cov>=0.8)&(df.DP_ratio2<0.8)&(df.DA_cov>=0.8)&(df.DA_ratio2<0.8); x1=df.loc[okdd,"d_DA_DP"].dropna()
for lab,xx,c in [("Within tetraploid Kew (primary vs alternate)",x1,"#c0392b"),("Within diploid individual (FH1 vs FH2)",base.W,"#9e9e9e"),("Between diploid populations",base.b,"#4d4d4d")]:
    xx=xx[xx>0]; ax.hist(np.log10(xx),bins=np.log10(bins),histtype="step",density=True,color=c,lw=1.4 if "tetraploid" in lab else 1.0,label=lab)
ax.set_xticks([-3,-2.5,-2,-1.5,-1]); ax.set_xticklabels(["0.1%","0.3%","1%","3%","10%"]); ax.set_xlim(-3.4,-0.4); ax.set_xlabel("CDS p-distance per gene (log scale; zeros omitted)"); ax.set_ylabel("density of genes")
ax.set_title("Tetraploid haplotype pairs diverge\nlike diploid alleles — no\nhomoeologous mode",loc="left"); ax.legend(frameon=False,fontsize=5.5,loc="upper right",handlelength=1.2); ax.set_ylim(0,ax.get_ylim()[1]*1.3)
for ax_,l in zip(axes,"ab"): PL(ax_,l)
fig.subplots_adjust(left=0.08,right=0.99,top=0.82,bottom=0.2); fig.savefig("paper/figures/Fig2_copy_structure.png",dpi=300); fig.savefig("paper/figures/Fig2_copy_structure.pdf"); plt.close(fig)
# ---- Fig 3
cp=pd.read_csv("data/annot/centromere_points.csv"); cp=cp[cp.hap=="H1"].set_index("chr"); lens=json.load(open("data/annot/chrom_lengths.json"))["FH1"]
iv=pd.read_csv("data/annot/tableS8_centromere_intervals.csv"); iv=iv[(iv.hap=="H1")&(iv.len>=100000)]
d3=df[okb].copy(); d3["pi"]=d3[["d_FH1_FH2","d_PAT_FH1","d_PAT_FH2"]].mean(axis=1); d3=d3[d3.pi.notna()&d3.FH1_chr.notna()&d3.FH1_start.notna()]; d3["chr"]=d3.FH1_chr.astype(int)
h3=pd.read_csv("results/divergence/h3_perchr_female.csv").set_index("chr"); cls_col={"polycentric":"#2471a3","acrocentric":"#c0392b","metacentric":"#7d6608"}
fig,axes=plt.subplots(3,5,figsize=(7.4,4.6),sharey=True); axes=axes.ravel(); orderc=[1,2,3,4,5,6,7,10,8,9,11,12,13]
for ax,c in zip(axes,orderc):
    g=d3[d3.chr==c].sort_values("FH1_start"); L=lens[str(c)]; w=(g.FH1_start//2e6).astype(int); med=g.groupby(w).pi.median(); n_=g.groupby(w).size(); xs=(med.index+0.5)*2; kind=cp.loc[c,"class"]; col=cls_col[kind]
    ax.scatter(g.FH1_start/1e6,g.pi*100,s=1.5,color="#bdbdbd",lw=0,rasterized=True); ax.plot(xs[n_>=10],med[n_>=10]*100,color=col,lw=1.4); ivc=iv[iv.chr==c]
    if kind=="acrocentric":
        r_=ivc.loc[ivc.len.idxmax()]; ax.axvspan(r_.start/1e6,r_.end/1e6,color=col,alpha=0.18,lw=0)
        for _,r2 in ivc.drop(r_.name).iterrows(): ax.axvspan(r2.start/1e6,r2.end/1e6,color=col,alpha=0.07,lw=0)
    else:
        for _,r2 in ivc.iterrows(): ax.axvspan(r2.start/1e6,r2.end/1e6,color=col,alpha=0.18,lw=0)
        ax.axvspan(cp.loc[c,"span_start"]/1e6,cp.loc[c,"span_end"]/1e6,color=col,alpha=0.06,lw=0)
    ax.axvline(cp.loc[c,"cen_point"]/1e6,color=col,lw=0.8,ls="--"); ax.set_xlim(0,L/1e6*1.02); ax.set_ylim(0,3.2)
    ax.text(0.03,0.97,f"chr{c} · {kind}\nT = {h3.loc[c,'T']:.2f}",transform=ax.transAxes,va="top",ha="left",fontsize=5.5,color=col)
    if c not in (1,6,11): ax.tick_params(labelleft=False)
    ax.tick_params(labelsize=5.5)
for ax in axes[13:]: ax.axis("off")
axes[14].text(0.0,0.95,"grey dots: genes, π̂ = mean of the three\ndiploid pairwise CDS p-distances\nline: 2-Mb window median (≥10 genes)\ndark shading: RepeatOBserver intervals\n≥100 kb (Hirabayashi 2025 Table S8, H1);\nlight: their span (polycentric) or minor\ncandidates (acrocentric); dashed:\noperational centromere point\nT = median π̂ of the 20 % of genes nearest\nthe point / median of the rest",transform=axes[14].transAxes,fontsize=5.5,va="top")
fig.supxlabel("position on female haplotype 1 (Mb)",fontsize=7); fig.supylabel("CDS nucleotide diversity π̂ (%)",fontsize=7)
fig.suptitle("Diversity is depressed around the predicted centromeric region on every chromosome — most deeply on the predicted-polycentric ones (blue)",fontsize=7.5,x=0.01,ha="left")
fig.subplots_adjust(left=0.07,right=0.99,top=0.9,bottom=0.1,hspace=0.3,wspace=0.12); fig.savefig("paper/figures/Fig3_diversity_landscape.png",dpi=300); fig.savefig("paper/figures/Fig3_diversity_landscape.pdf"); plt.close(fig)
# ---- Fig 4
S="results/synteny"; uu_=pd.read_csv(f"{S}/blocks_UDI_UUR.csv"); uu_=uu_[(uu_.block>=0)&(uu_.chr_o==2)]
fig,axes=plt.subplots(1,3,figsize=(7.4,2.8),gridspec_kw=dict(width_ratios=[1.5,1,1],wspace=0.6)); ax=axes[0]; c8="#c0392b"; c13="#2471a3"
for c,col,lab in [(8,c8,"U. dioica chr8"),(13,c13,"U. dioica chr13")]:
    g=uu_[uu_.chr_u==c]; ax.scatter(g.start_o/1e6,g.start_u/1e6,s=2,color=col,lw=0,label=lab,rasterized=True)
ax.set_xlabel("U. urens chr2 (Mb)"); ax.set_ylabel("U. dioica chr8 / chr13 (Mb)"); ax.set_title("U. urens chr2 is an eight-segment mosaic of\nU. dioica chr8 and chr13",loc="left"); ax.legend(frameon=False,fontsize=6,loc="upper right",markerscale=4); ax.set_xlim(0,33); ax.set_ylim(0,52)
for sp,ax,title in [("PJU",axes[1],"Parietaria"),("CAN",axes[2],"Cannabis")]:
    b=pd.read_csv(f"{S}/blocks_UDI_{sp}.csv"); b=b[b.block>=0]; S8=set(uu_[uu_.chr_u==8].udi); S13=set(uu_[uu_.chr_u==13].udi)
    d8=b[b.udi.isin(S8)].groupby("chr_o").size().sort_values(ascending=False); d13=b[b.udi.isin(S13)].groupby("chr_o").size().sort_values(ascending=False)
    allc=sorted(set(d8.index[:4])|set(d13.index[:4])); y=np.arange(len(allc)); w=0.38
    ax.barh(y+w/2,[d8.get(c,0) for c in allc],height=w,color=c8,label="chr8-derived genes"); ax.barh(y-w/2,[d13.get(c,0) for c in allc],height=w,color=c13,label="chr13-derived genes")
    ax.set_yticks(y); ax.set_yticklabels([f"chr{int(c)}" for c in allc],fontsize=6); ax.set_xlabel("collinear block genes"); ax.set_title(f"{title}: chr8- and chr13-derived\ngenes on disjoint chromosomes",loc="left")
    if sp=="PJU": ax.legend(frameon=False,fontsize=6,loc="lower right")
    ax.invert_yaxis()
for ax_,l in zip(axes,"abc"): PL(ax_,l)
fig.subplots_adjust(left=0.08,right=0.99,top=0.82,bottom=0.2); fig.savefig("paper/figures/Fig4_synteny.png",dpi=300); fig.savefig("paper/figures/Fig4_synteny.pdf"); plt.close(fig)
print("figures written")

# ---- Fig 5 (windowed tracks)
lens=json.load(open("data/annot/chrom_lengths.json"))["MH2"]; Lmax=max(lens.values())/1e6
d5=pd.read_csv("results/divergence/h1v2_genes_DP.csv"); d5a=pd.read_csv("results/divergence/h1v2_genes_DA.csv").set_index("qid")
d5=d5.join(d5a[["t"]].rename(columns={"t":"tA"}),on="qid"); blocks=vDP["paired_scan"]["blocks"]
cols5={"t":"#c0392b","tA":"#e67e22","b":"#4d4d4d","W":"#9e9e9e","d_DA_DP":"#2471a3"}
fig,axes=plt.subplots(13,1,figsize=(7.4,9.5),sharex=True)
for ax,c in zip(axes,range(1,14)):
    g=d5[d5.mh2_chr==c].sort_values("mh2_start"); x=g.mh2_start.values/1e6
    for key,lab in [("W","within diploid (FH1 vs FH2)"),("b","Bulgarian diploid vs female (b)"),("t","tetraploid primary vs female (t)"),("tA","tetraploid alternate vs female"),("d_DA_DP","tetraploid primary vs alternate")]:
        y=pd.Series(g[key].values).rolling(50,center=True,min_periods=25).median().values*100
        ax.plot(x,y,color=cols5[key],lw=1.3 if key=="t" else 0.9,label=lab)
    for bk in blocks:
        if bk["chr"]==c: ax.axvspan(bk["start"]/1e6,bk["end"]/1e6,color="#f5b7b1",alpha=0.5,lw=0)
    ax.axvline(lens[str(c)]/1e6,color="#bbbbbb",lw=0.6); ax.set_ylim(0,1.6); ax.set_yticks([0,0.5,1.0,1.5]); ax.set_xlim(0,Lmax*1.01)
    ax.text(0.995,0.9,f"chr{c}",transform=ax.transAxes,ha="right",va="top",fontsize=6); ax.tick_params(labelsize=5.5)
axes[0].legend(frameon=False,fontsize=5.5,ncol=5,loc="upper left",bbox_to_anchor=(0,1.9),handlelength=1.5,columnspacing=1.0)
axes[-1].set_xlabel("position on male haplotype 2 (Mb); grey tick = chromosome end"); fig.supylabel("50-gene rolling median CDS p-distance (%)",fontsize=7)
fig.suptitle("Along every chromosome the tetraploid tracks the diploid baseline; the two shaded blocks are regions where b (not t) departs",fontsize=7.5,x=0.01,ha="left",y=0.995)
fig.subplots_adjust(left=0.08,right=0.99,top=0.93,bottom=0.05,hspace=0.25)
fig.savefig("paper/figures/Fig5_windowed_tracks.png",dpi=300); fig.savefig("paper/figures/Fig5_windowed_tracks.pdf"); plt.close(fig)
print("figures written")
