#!/usr/bin/env python3
"""Four-way-core pathway enrichment — clean horizontal lollipop (not dot-grid)."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save
import numpy as np
from pathlib import Path

CDIR = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/four_way_pathways/50_Single")
ONTO = [("kegg", "KEGG", "#E64B35"), ("reactome", "Reactome", "#4DBBD5"),
        ("wikipathways", "WikiPathways", "#00A087"),
        ("biological_process", "GO BP", "#3C5488"),
        ("molecular_function", "GO MF", "#8491B4")]
TOP = {"kegg": 3, "reactome": 3, "wikipathways": 2, "biological_process": 4, "molecular_function": 2}
# drop the most generic GO umbrella terms
STOP = {"biological regulation", "regulation of biological process", "cellular process",
        "positive regulation of biological process", "negative regulation of biological process",
        "positive regulation of cellular process", "negative regulation of cellular process",
        "binding", "protein binding", "regulation of cellular process",
        "single-organism process", "metabolic process"}

def parse(path):
    rows = {}
    for i, ln in enumerate(open(path, encoding="utf-8", errors="replace")):
        if i == 0: continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 11: continue
        t = f[1]
        try: p = float(f[2]); k = int(f[5])
        except ValueError: continue
        if t not in rows or p < rows[t][0]: rows[t] = (p, k)
    items = sorted(rows.items(), key=lambda kv: kv[1][0]); n = len(items)
    out = [(t, k, min(1.0, p*n/(r+1))) for r, (t, (p, k)) in enumerate(items)]
    for i in range(len(out)-2, -1, -1):
        if out[i][2] > out[i+1][2]: out[i] = (out[i][0], out[i][1], out[i+1][2])
    return out

terms = []
for fn, lab, col in ONTO:
    got = [(t, k, fdr, lab, col) for (t, k, fdr) in parse(CDIR / f"{fn}.txt")
           if t.lower() not in STOP and fdr < 0.6]
    terms += got[:TOP[fn]]
# sort all by significance (ascending FDR -> plot smallest FDR at top)
terms.sort(key=lambda x: x[2])
terms = terms[::-1]  # so most significant ends at top after barh

labels = [t[0] if len(t[0]) <= 46 else t[0][:45] + "…" for t in terms]
vals = [-np.log10(max(t[2], 1e-12)) for t in terms]
ks = [t[1] for t in terms]
cols = [t[4] for t in terms]

fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(terms) + 1.0))
y = np.arange(len(terms))
ax.hlines(y, 0, vals, color=cols, lw=2.2, alpha=0.55, zorder=1)
ax.scatter(vals, y, s=[max(30, k*16) for k in ks], c=cols, edgecolor="#222",
           linewidth=0.5, zorder=3)
for yi, k in zip(y, ks):
    ax.text(vals[yi] + 0.12, yi, f"{k}", va="center", fontsize=6.2, color="#333")
ax.axvline(-np.log10(0.05), color="#A0AEC0", ls="--", lw=0.9, zorder=0)
ax.text(-np.log10(0.05), len(terms)-0.3, " FDR 0.05", fontsize=6.3, color="#718096", va="top")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.2)
ax.set_xlabel("−log$_{10}$ FDR   (dot size ∝ genes in term, k)")
ax.set_xlim(0, max(vals) * 1.18)
ax.set_title("Pathway / GO enrichment of the four-way memory core (64 genes)", fontsize=10, pad=10)

from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=7,
                  markeredgecolor="#222", markeredgewidth=0.4, label=l) for _, l, c in ONTO]
ax.legend(handles=handles, title="Ontology", loc="lower right", fontsize=6.8, title_fontsize=7)
save(fig, "Fig5_pathway_fourway")
