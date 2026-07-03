#!/usr/bin/env python3
"""Figure 2 — final memory-peak counts across the four datasets and stringencies."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, DCOLOR, DATASETS, DLABEL
import numpy as np
from pathlib import Path

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
PATHS = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory",
}
CONDS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]
CLABEL = ["50%\nSingle", "50%\nReciprocal", "75%\nSingle", "75%\nReciprocal"]

def n(p):
    return sum(1 for _ in open(p))

counts = {d: [n(PATHS[d] / f"StrictMemory{c}.bed") for c in CONDS] for d in DATASETS}

fig, ax = plt.subplots(figsize=(7.4, 3.9))
x = np.arange(len(CONDS)); w = 0.2
for i, d in enumerate(DATASETS):
    bars = ax.bar(x + (i - 1.5) * w, counts[d], w, color=DCOLOR[d],
                  label=d, edgecolor="white", linewidth=0.5)
    for b, v in zip(bars, counts[d]):
        ax.text(b.get_x() + b.get_width()/2, v + 40, f"{v:,}", ha="center",
                va="bottom", fontsize=5.6, color="#333", rotation=90)

ax.set_xticks(x); ax.set_xticklabels(CLABEL)
ax.set_ylabel("Final memory peaks (n)")
ax.set_xlabel("Stringency variant")
ax.set_ylim(0, max(max(v) for v in counts.values()) * 1.18)
ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.10),
          columnspacing=1.4, handlelength=1.1)
ax.set_title("Acquired memory peaks per dataset (strict, control-subtracted)",
             pad=42, fontsize=10.5)
ax.text(-0.09, 1.04, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")
fig.subplots_adjust(top=0.80)
save(fig, "Fig2_memory_peak_counts")
