#!/usr/bin/env python3
"""Figure 4 — pathway/GO enrichment dot plots for the shared-gene cores."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save
import numpy as np
from pathlib import Path

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison")
GOV = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/gene_overlap")
CORES = [
    ("IMQ ∩ HFSC\n(356 genes)", GOV / "pathways/50_Single"),
    ("Three-way core\n(147 genes)", BASE / "three_way_pathways/50_Single"),
    ("Four-way core\n(64 genes)", BASE / "four_way_pathways/50_Single"),
]
ONTO = [("kegg", "KEGG", "#E64B35"), ("reactome", "Reactome", "#4DBBD5"),
        ("wikipathways", "WikiPathways", "#00A087"), ("biological_process", "GO BP", "#3C5488")]
TOP_PER_ONTO = 3

def parse(path):
    rows = {}
    try:
        for i, ln in enumerate(open(path, encoding="utf-8", errors="replace")):
            if i == 0: continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 11: continue
            t = f[1]
            try: p = float(f[2]); k = int(f[5])
            except ValueError: continue
            if t not in rows or p < rows[t][0]: rows[t] = (p, k)
    except FileNotFoundError:
        return []
    items = sorted(rows.items(), key=lambda kv: kv[1][0])
    n = len(items)
    out = []
    for rank, (t, (p, k)) in enumerate(items, 1):
        fdr = min(1.0, p * n / rank)
        out.append((t, k, fdr))
    # enforce monotone FDR
    for i in range(len(out) - 2, -1, -1):
        if out[i][2] > out[i+1][2]:
            out[i] = (out[i][0], out[i][1], out[i+1][2])
    return out

def short(t, n=42):
    return t if len(t) <= n else t[:n-1] + "…"

fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.6), sharex=False)
for ax, (title, cdir) in zip(axes, CORES):
    rowlabels, xs, sizes, colors = [], [], [], []
    y = 0
    yticks, yticklabels = [], []
    for fn, lab, col in ONTO:
        terms = parse(cdir / f"{fn}.txt")
        terms = [t for t in terms if t[2] < 0.5][:TOP_PER_ONTO]
        for t, k, fdr in terms:
            xs.append(-np.log10(max(fdr, 1e-12)))
            sizes.append(k)
            colors.append(col)
            yticks.append(y); yticklabels.append(short(t))
            y += 1
    yarr = np.array(yticks)
    sc = ax.scatter(xs, yarr, s=[max(22, min(s, 30)*13) for s in sizes], c=colors,
                    edgecolor="#222", linewidth=0.4, alpha=0.92, zorder=3)
    ax.axvline(-np.log10(0.05), color="#A0AEC0", ls="--", lw=0.8, zorder=1)
    ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlabel("−log$_{10}$ FDR")
    ax.set_title(title, fontsize=8.8)
    ax.set_xlim(0, max(xs) * 1.15 if xs else 1)
    ax.grid(axis="x", color="#EDF2F7", lw=0.6, zorder=0)

# legends
from matplotlib.lines import Line2D
onto_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=7,
                       markeredgecolor="#222", markeredgewidth=0.4, label=l)
                for _, l, c in ONTO]
size_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
                       markersize=np.sqrt(max(22, min(k, 30)*13))/2.2, markeredgecolor="#222",
                       markeredgewidth=0.4, label=f"{k}") for k in (3, 10, 30)]
leg1 = fig.legend(handles=onto_handles, title="Ontology", loc="lower center",
                  bbox_to_anchor=(0.32, -0.06), ncol=4, columnspacing=1.0)
fig.legend(handles=size_handles, title="Genes in term (k)", loc="lower center",
           bbox_to_anchor=(0.72, -0.07), ncol=3, columnspacing=1.0)
fig.add_artist(leg1)
fig.text(0.07, 0.97, "a", fontsize=14, fontweight="bold")
fig.suptitle("Pathway / GO enrichment of the shared memory-gene cores (dashed line = FDR 0.05)",
             fontsize=10.5, y=1.0)
save(fig, "Fig4_pathway_enrichment")
