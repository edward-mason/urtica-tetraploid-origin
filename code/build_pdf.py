"""build_pdf.py — render paper/manuscript.md (+ figures, + AI disclosure) into a single submission PDF with xhtml2pdf.
Usage (nettle env): python code/build_pdf.py  → paper/Mason_2026_nettle_tetraploid_preprint.pdf (and .html)."""
import re, markdown, os, io, matplotlib
FD=os.path.abspath("paper/fonts").replace("\\","/")  # DejaVu Sans TTFs copied from matplotlib (Bitstream Vera licence)
from xhtml2pdf import pisa
src=open("paper/manuscript.md",encoding="utf-8").read()
# front matter
fm=re.match(r"^---\n(.*?)\n---\n",src,flags=re.S); meta=dict(re.findall(r'^(\w+): "(.*)"$',fm.group(1),flags=re.M)); body=src[fm.end():]
# move AI disclosure in as a section before References
ai=open("paper/AI_DISCLOSURE.md",encoding="utf-8").read()
ai=re.sub(r"^# .*\n","",ai).strip()
i=body.find("# References"); body=body[:i]+"# AI disclosure statement\n\n"+ai+"\n\n# Supplementary material\n\nSupplementary files accompany the preprint as a single archive (see Data availability): S1 preregistration (v1 and v2 with SHA-256 hashes); S2 laboratory notebook; S3 independent reviewer reports and point-by-point responses; S4 claim ledger; S5 validation plan; Table S1 per-gene distances (gene_distances.csv); Table S2 per-gene H1 statistics (h1v2_genes_DP.csv, h1v2_genes_DA.csv); Table S3 per-chromosome H3 statistics; Table S4 synteny blocks; Table S5 copy-structure and sensitivity tables (JSON); code (Snakemake workflow, Python scripts, environment.yml).\n\n"+body[i:]
# figures: replace "# Figure legends" section content with image + legend pairs
figs={"Figure 1":"Fig1_tetraploid_origin.png","Figure 2":"Fig2_copy_structure.png","Figure 3":"Fig3_diversity_landscape.png","Figure 4":"Fig4_synteny.png","Figure 5":"Fig5_windowed_tracks.png"}
def fig_sub(m):
    label=m.group(1); fn=figs[label]; path=os.path.abspath(f"paper/figures/{fn}").replace("\\","/")
    return f'\n\n<div class="fig"><img src="{path}" /></div>\n\n**{label}.'
body=re.sub(r"\n\*\*(Figure \d)\.",fig_sub,body)
# superscripts ^x^ and italics handled by markdown; convert ^...^ to <sup>
body=re.sub(r"\^([^\^\n]{1,40})\^",r"<sup>\1</sup>",body)
html_body=markdown.markdown(body,extensions=["tables","sane_lists"])
title=meta["title"].replace("*Urtica dioica*","<i>Urtica dioica</i>")
css="""
@font-face { font-family: DejaVu; src: url("FDPATH/DejaVuSans.ttf"); }
@font-face { font-family: DejaVu; src: url("FDPATH/DejaVuSans-Bold.ttf"); font-weight: bold; }
@font-face { font-family: DejaVu; src: url("FDPATH/DejaVuSans-Oblique.ttf"); font-style: italic; }
@font-face { font-family: DejaVu; src: url("FDPATH/DejaVuSans-BoldOblique.ttf"); font-weight: bold; font-style: italic; }
@page { size: A4; margin: 2.2cm 2.0cm 2.2cm 2.0cm; @frame footer { -pdf-frame-content: footer; bottom: 0.8cm; height: 0.8cm; } }
body { font-family: DejaVu; font-size: 9.2pt; line-height: 1.35; }
h1.title { font-size: 14.5pt; line-height: 1.25; margin-bottom: 6pt; }
h1 { font-size: 12pt; margin-top: 16pt; margin-bottom: 4pt; border-bottom: 0.5pt solid #999; }
h2 { font-size: 10.2pt; margin-top: 11pt; margin-bottom: 3pt; }
p { margin: 0 0 6pt 0; text-align: justify; }
table { font-size: 7pt; margin: 6pt 0 10pt 0; }
th, td { border: 0.4pt solid #888; padding: 2pt 4pt; vertical-align: top; }
th { background: #eeeeee; }
img { width: 17cm; } div.fig { text-align: center; margin: 10pt 0 4pt 0; }
#footer { font-size: 7.5pt; color: #666; text-align: center; }
""".replace("FDPATH",FD)
front=f'<h1 class="title">{title}</h1>'
html=f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>{front}{html_body}<div id='footer'>Mason — preprint, {meta['date']}. CC BY 4.0.</div></body></html>"
# the markdown already contains the author block; remove the duplicated first markdown H1-less lines? (kept)
open("paper/Mason_2026_nettle_tetraploid_preprint.html","w",encoding="utf-8").write(html)
out="paper/Mason_2026_nettle_tetraploid_preprint.pdf"
with open(out,"wb") as f:
    r=pisa.CreatePDF(io.BytesIO(html.encode("utf-8")),dest=f,encoding="utf-8")
print("errors:",r.err,"->",out,os.path.getsize(out))
