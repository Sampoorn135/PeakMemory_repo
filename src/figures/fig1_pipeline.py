#!/usr/bin/env python3
"""Figure 1 — schematic of the BED-based chromatin memory-peak pipeline."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, ACCENT, GREY
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(9.4, 4.8))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc, ec=None, fs=8, tc="white"):
    ec = ec or fc
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=2.5",
                                linewidth=1.2, facecolor=fc, edgecolor=ec, zorder=2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold", zorder=3, linespacing=1.25)

def arrow(x1, y1, x2, y2, color=GREY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=15,
                                 linewidth=1.6, color=color, zorder=1,
                                 shrinkA=2, shrinkB=2))

# main horizontal flow (y ~ 50)
box(1.5, 48, 16, 14, "perturbation\n1…N", "#9F7AEA", fs=8)
box(25, 48, 17, 14, "Step 2\nPerturbation\nconsensus", "#805AD5", fs=7.8)
box(49, 48, 17, 14, "Step 3\nMemory\nmatch", "#2B6CB0", fs=7.8)
box(73, 48, 18, 14, "Step 4\nStrict control\nsubtraction", ACCENT, fs=7.8)
# inputs dropping in from above
box(49, 78, 17, 11, "memory\n(chase endpoint)", "#3182CE", fs=7.5)
box(73, 78, 18, 11, "control\n(baseline)", GREY, fs=7.5)
# output below
box(70, 18, 21, 13, "FINAL\nmemory peaks", "#38A169", fs=9)

# arrows
arrow(17.5, 55, 25, 55)
arrow(42, 55, 49, 55)
arrow(66, 55, 73, 55)
arrow(57.5, 78, 57.5, 62, color="#3182CE")     # memory -> step3
arrow(82, 78, 82, 62, color=GREY)              # control -> step4
arrow(82, 48, 82, 31, color=ACCENT)            # step4 -> final

# captions under steps
ax.text(33.5, 45, "peaks present in\nALL perturbations", ha="center", va="top", fontsize=6.4, color="#553C9A")
ax.text(57.5, 45, "consensus ∩ memory\nunion regions", ha="center", va="top", fontsize=6.4, color="#2C5282")
ax.text(82, 45, "drop if ANY ≥1 bp\noverlap with control", ha="center", va="top", fontsize=6.4, color="#9C4221")
ax.text(80.5, 16, "50% / 75% × Single / Reciprocal", ha="center", va="top", fontsize=6.4,
        color="#22543D", style="italic")

# Step 1 banner
ax.add_patch(FancyBboxPatch((1.5, 68), 40, 7, boxstyle="round,pad=0.4,rounding_size=2",
                            linewidth=0, facecolor="#EDF2F7", zorder=1))
ax.text(21.5, 71.5, "Step 1 — every input BED is sorted + internally merged first",
        ha="center", va="center", fontsize=7, color="#2D3748", fontweight="bold")

# bottom notes
ax.text(34, 27, "Coordinate-only\nno BAM · no counts · no fold-change · no p-values\noverlap geometry computed explicitly in Python",
        ha="center", va="center", fontsize=7, color=GREY, style="italic", linespacing=1.4)
ax.text(50, 4.5, "Biological logic:   reproducibly opened across perturbations   →   still open in memory   →   absent from baseline   =   acquired memory peak",
        ha="center", fontsize=7.3, color="#1A202C")

ax.text(0, 99, "a", fontsize=15, fontweight="bold", va="top")
ax.set_title("BED-based chromatin memory-peak discovery pipeline", fontsize=11, pad=12)
save(fig, "Fig1_pipeline_schematic")
