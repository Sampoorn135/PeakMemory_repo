#!/usr/bin/env python3
"""Genomic-feature distribution + enrichment of memory peaks (4 studies)."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, DATASETS, FIGDIR
import json, numpy as np
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap

D = json.load(open("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/genomic_features/genomic_features.json"))
res, gfrac = D["results"], D["genome_frac"]

ORDER = ["Promoter-TSS", "5'UTR", "Exon", "3'UTR", "TTS", "ncRNA/other", "Intron", "Intergenic"]
COL = {"Promoter-TSS": "#E64B35", "5'UTR": "#F0A04B", "Exon": "#F7D154",
       "3'UTR": "#E59866", "TTS": "#B07CC6", "ncRNA/other": "#AEB6BF",
       "Intron": "#4DBBD5", "Intergenic": "#3C5488"}

# ---------- Fig 9: stacked distribution (peaks per dataset + genome reference) ----------
rows = DATASETS + ["genome"]
fig, ax = plt.subplots(figsize=(8.6, 3.8))
y = np.arange(len(rows))[::-1]
left = np.zeros(len(rows))
for cat in ORDER:
    vals = [res[d]["pct"][cat] for d in DATASETS] + [100*gfrac[cat]]
    ax.barh(y, vals, left=left, color=COL[cat], edgecolor="white", linewidth=0.5, label=cat)
    for yi, v, l in zip(y, vals, left):
        if v >= 6:
            ax.text(l + v/2, yi, f"{v:.0f}", ha="center", va="center", fontsize=6.4,
                    color="white", fontweight="bold")
    left += np.array(vals)
ax.set_yticks(y); ax.set_yticklabels(rows, fontsize=8.5)
ax.set_xlim(0, 100); ax.set_xlabel("% of peaks")
ax.set_title("Genomic-feature distribution of memory peaks (StrictMemory 50% Single)", fontsize=10, pad=8)
ax.legend(ncol=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=6.6,
          columnspacing=0.9, handlelength=1.0, handletextpad=0.4)
# annotate distal fraction
for d, yi in zip(DATASETS, y[:len(DATASETS)]):
    ax.text(101, yi, f"distal {res[d]['distal_pct']:.0f}%", va="center", fontsize=6.5, color="#2D3748")
ax.text(-0.085, 1.05, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")
save(fig, "Fig9_genomic_distribution")

# ---------- Fig 10: log2(obs/exp) enrichment heatmap ----------
M = np.array([[res[d]["log2enr"][cat] for d in DATASETS] for cat in ORDER])
vmax = np.nanmax(np.abs(M))
RB = LinearSegmentedColormap.from_list("rb", ["#3C5488", "#9EC9E2", "#FFFFFF", "#F4A582", "#C0392B"])
fig, ax = plt.subplots(figsize=(5.2, 4.4))
im = ax.imshow(M, cmap=RB, vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(DATASETS))); ax.set_xticklabels(DATASETS, rotation=30, ha="right", fontsize=8)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER, fontsize=8)
for i in range(len(ORDER)):
    for j in range(len(DATASETS)):
        v = M[i, j]
        ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=6.8,
                color="white" if abs(v) > vmax*0.6 else "#1A202C")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("log$_2$(obs / exp)", fontsize=7); cb.ax.tick_params(labelsize=6)
ax.set_title("Feature enrichment vs genome\n(red = enriched, blue = depleted)", fontsize=9.5)
ax.text(-0.28, 1.06, "b", transform=ax.transAxes, fontsize=14, fontweight="bold")
save(fig, "Fig10_genomic_enrichment")
print("done")
