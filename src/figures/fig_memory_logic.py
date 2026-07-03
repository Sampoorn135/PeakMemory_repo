#!/usr/bin/env python3
"""Simple slide infographic: how a memory peak is defined (+ union operation)."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, FIGDIR
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

GREY = "#A0AEC0"; BASE = "#718096"; PURP = "#9F7AEA"; BLUE = "#2B6CB0"
GREEN = "#38A169"; RED = "#E53E3E"; TRACK = "#E2E8F0"

fig, ax = plt.subplots(figsize=(9.6, 5.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ROWS = [("Day 0  (baseline / control)", 80),
        ("Injury day 1", 65),
        ("Injury day 7", 53),
        ("Memory  (chase)", 40)]
for lab, y in ROWS:
    ax.text(20, y, lab, ha="right", va="center", fontsize=8.2, color="#2D3748", fontweight="bold")

LOCI = [36, 61, 85]          # x centres of the 3 example loci
def track(cx, y):
    ax.plot([cx-8, cx+8], [y, y], color=TRACK, lw=3, solid_capstyle="round", zorder=1)
def peak(cx, y, color, w=13):
    ax.add_patch(FancyBboxPatch((cx-w/2, y-3), w, 6, boxstyle="round,pad=0.2,rounding_size=1.5",
                                facecolor=color, edgecolor="white", lw=0.8, zorder=3))

# draw all tracks
for cx in LOCI:
    for _, y in ROWS:
        track(cx, y)

# --- Locus 1: MEMORY PEAK (open across injury + memory, closed at baseline) ---
cx = LOCI[0]
# baseline closed (no peak). injury D1/D7 + memory open, slightly offset -> union
peak(cx-3, 65, PURP, 11)
peak(cx,   53, PURP, 11)
peak(cx+3, 40, BLUE, 11)
# union region bar
ax.add_patch(FancyBboxPatch((cx-3-5.5-1, 26), (cx+3+5.5)-(cx-3-5.5)+2, 4.5,
            boxstyle="round,pad=0.2,rounding_size=1.5", facecolor=GREEN, edgecolor="white", lw=0.8, zorder=3))
ax.add_patch(FancyArrowPatch((cx, 37), (cx, 31), arrowstyle="-|>", mutation_scale=11, color="#4A5568", lw=1.3))
ax.text(cx, 21, "union →\nmemory peak", ha="center", va="top", fontsize=6.6, color="#22543D", fontweight="bold")
ax.text(cx, 10, "✓ SELECTED", ha="center", fontsize=8.5, color=GREEN, fontweight="bold")

# --- Locus 2: NOT RETAINED (open in injury, CLOSED in memory) ---
cx = LOCI[1]
peak(cx, 65, PURP); peak(cx, 53, PURP)
ax.text(cx, 40, "closed", ha="center", va="center", fontsize=6.2, color=GREY, style="italic")
ax.text(cx, 21, "not open in memory", ha="center", va="top", fontsize=6.6, color="#4A5568")
ax.text(cx, 10, "✗ not retained", ha="center", fontsize=8.5, color="#718096", fontweight="bold")

# --- Locus 3: PRE-EXISTING (open at baseline) ---
cx = LOCI[2]
peak(cx, 80, BASE); peak(cx, 65, PURP); peak(cx, 53, PURP); peak(cx, 40, BLUE)
ax.text(cx, 21, "open at baseline\n→ removed", ha="center", va="top", fontsize=6.6, color="#9B2C2C")
ax.text(cx, 10, "✗ removed", ha="center", fontsize=8.5, color=RED, fontweight="bold")

# legend
ax.add_patch(FancyBboxPatch((2, 92), 4, 3.2, boxstyle="round,pad=0.1,rounding_size=1", facecolor=BLUE, edgecolor="white"))
ax.text(7.5, 93.6, "open region (ATAC peak)", fontsize=7, va="center", color="#2D3748")
ax.plot([55, 60], [93.6, 93.6], color=TRACK, lw=3, solid_capstyle="round")
ax.text(61, 93.6, "closed (no peak)", fontsize=7, va="center", color="#2D3748")

ax.set_title("How a memory peak is defined", fontsize=12.5, pad=10)
ax.text(50, 2.5, "Rule:  open across injury days (consensus)   +   still open in memory (chase)   +   absent at baseline   =   memory peak",
        ha="center", fontsize=7.6, color="#1A202C")
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIGDIR / f"Fig_memory_logic.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig_memory_logic.png")
plt.close(fig)
