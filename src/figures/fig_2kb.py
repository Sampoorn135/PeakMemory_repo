#!/usr/bin/env python3
"""2 kb promoter-proximal: gene-overlap Venn + per-study pathway enrichment."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, DCOLOR, DATASETS, FIGDIR
import numpy as np
from pathlib import Path
from venn import venn
from matplotlib.colors import to_rgba

G = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/genes_2kb")
P = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/pathways_2kb")

# ---- Venn of 2kb promoter gene sets ----
gsets = {d: set(open(G / f"{d}_genes_2kb.txt").read().split()) for d in DATASETS}
fig, ax = plt.subplots(figsize=(6.4, 6.0))
venn({d: gsets[d] for d in DATASETS}, cmap=[to_rgba(DCOLOR[d], 0.45) for d in DATASETS],
     fontsize=8.5, legend_loc="upper left", ax=ax)
core = len(set.intersection(*gsets.values()))
ax.set_title("Promoter-proximal (<2 kb to TSS) memory-gene overlap", fontsize=10.5, pad=12)
ax.text(0.5, -0.02, f"four-way core = {core} genes", transform=ax.transAxes, ha="center",
        fontsize=9, fontweight="bold")
leg = ax.get_legend()
if leg:
    for t, d in zip(leg.get_texts(), DATASETS): t.set_text(f"{d} ({len(gsets[d])})"); t.set_fontsize(8.2)
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig_2kb_gene_venn.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig); print("wrote Fig_2kb_gene_venn.png")

# ---- per-study enrichment lollipops ----
PATHDB = [("kegg", "KEGG", "#E64B35"), ("reactome", "Reactome", "#4DBBD5"),
          ("wikipathways", "WikiPathways", "#00A087")]
STOP = {"biological regulation", "cellular process", "metabolic process",
        "regulation of biological process"}
def parse(path):
    rows = {}
    try:
        for i, ln in enumerate(open(path, encoding="utf-8", errors="replace")):
            if i == 0: continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 11: continue
            t = f[1]
            try: pv = float(f[2]); k = int(f[5])
            except ValueError: continue
            if t not in rows or pv < rows[t][0]: rows[t] = (pv, k)
    except FileNotFoundError:
        return []
    it = sorted(rows.items(), key=lambda kv: kv[1][0]); n = len(it)
    out = [(t, k, min(1, pv*n/(r+1))) for r, (t, (pv, k)) in enumerate(it)]
    for i in range(len(out)-2, -1, -1):
        if out[i][2] > out[i+1][2]: out[i] = (out[i][0], out[i][1], out[i+1][2])
    return out
def top(d, per=2):
    o = []
    for fn, lab, col in PATHDB:
        t = [(x[0], x[1], x[2], col) for x in parse(P / d / f"{fn}.txt")
             if x[0].lower() not in STOP and "podnet" not in x[0].lower()]
        t.sort(key=lambda z: z[2]); o += t[:per]
    o.sort(key=lambda z: z[2]); return o[:6]
def short(t, n=34): return t if len(t) <= n else t[:n-1]+"…"

fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.0))
for ax, d in zip(axes.flat, DATASETS):
    terms = top(d)[::-1]
    if terms:
        y = np.arange(len(terms)); vals = [-np.log10(max(t[2], 1e-12)) for t in terms]
        ax.hlines(y, 0, vals, color=[t[3] for t in terms], lw=2.0, alpha=0.55)
        ax.scatter(vals, y, s=[max(22, t[1]*12) for t in terms], c=[t[3] for t in terms],
                   edgecolor="#222", linewidth=0.4, zorder=3)
        ax.axvline(-np.log10(0.05), color="#A0AEC0", ls="--", lw=0.8)
        ax.set_yticks(y); ax.set_yticklabels([short(t[0]) for t in terms], fontsize=6.4)
        ax.set_xlim(0, max(vals)*1.2 + 0.2)
    ax.set_title(f"{d}  ({len(gsets[d])} promoter genes)", fontsize=8.4)
    ax.tick_params(axis="x", labelsize=6)
from matplotlib.lines import Line2D
h = [Line2D([0],[0],marker="o",color="w",markerfacecolor=c,markersize=6,markeredgecolor="#222",markeredgewidth=0.4,label=l) for _,l,c in PATHDB]
fig.legend(handles=h, title="Ontology", loc="lower center", ncol=3, bbox_to_anchor=(0.5,-0.04), fontsize=7, title_fontsize=7.5)
fig.suptitle("Pathway enrichment of promoter-proximal (<2 kb) memory genes, per study  (−log$_{10}$ FDR; dashed = 0.05)", fontsize=9.8, y=1.0)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig_2kb_enrichment.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig); print("wrote Fig_2kb_enrichment.png")
