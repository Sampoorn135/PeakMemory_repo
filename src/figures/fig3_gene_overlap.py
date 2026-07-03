#!/usr/bin/env python3
"""Figure 3 — cross-dataset gene overlap (pairwise heatmap, peak-vs-gene, UpSet)."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, save, DCOLOR, DATASETS, FIGDIR
import numpy as np
from pathlib import Path
from itertools import combinations
import bisect
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, annotate_set

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
PATHS = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory",
}
COND = "50_Single"

# ---- gene-locus sets + peak sets ----
by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")
genes, peaks = {}, {}
for d in DATASETS:
    bed = load_bed(PATHS[d] / f"StrictMemory{COND}.bed")
    _, loci = annotate_set(bed, by_chrom, t2l, locus_info)
    genes[d] = loci
    pk = {}
    for c, s, e in bed:
        pk.setdefault(c, []).append((s, e))
    for c in pk: pk[c].sort()
    peaks[d] = pk

def pk_overlap(a, b):
    n = 0
    for c, ivs in a.items():
        if c not in b: continue
        bs = [x[0] for x in b[c]]; arr = b[c]
        for s, e in ivs:
            i = bisect.bisect_left(bs, e) - 1; hit = False
            while i >= 0:
                cs, ce = arr[i]
                if ce > s and cs < e: hit = True; break
                if bs[i] < s - 5_000_000: break
                i -= 1
            if hit: n += 1
    return n

# ---- Figure 3: heatmap + peak-vs-gene ----
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9), gridspec_kw={"width_ratios": [1, 1.15]})

# panel A: pairwise shared-gene heatmap
M = np.zeros((4, 4))
for i, x in enumerate(DATASETS):
    for j, y in enumerate(DATASETS):
        M[i, j] = len(genes[x]) if i == j else len(genes[x] & genes[y])
ax = axes[0]
im = ax.imshow(M, cmap="YlGnBu", aspect="equal")
ax.set_xticks(range(4)); ax.set_yticks(range(4))
ax.set_xticklabels(DATASETS, rotation=35, ha="right"); ax.set_yticklabels(DATASETS)
for i in range(4):
    for j in range(4):
        val = int(M[i, j])
        ax.text(j, i, f"{val:,}", ha="center", va="center", fontsize=7.2,
                color="white" if M[i, j] > M.max()*0.55 else "#1A202C",
                fontweight="bold" if i == j else "normal")
ax.set_title("Shared memory genes (pairwise)\ndiagonal = set size", fontsize=8.5)
ax.text(-0.18, 1.06, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.ax.tick_params(labelsize=6)

# panel B: peak vs gene overlap (convergence at gene level)
ax = axes[1]
pairs = list(combinations(DATASETS, 2))
plabels = [f"{x[:3]}–{y[:3]}" for x, y in pairs]
gene_sh = [len(genes[x] & genes[y]) for x, y in pairs]
peak_sh = [pk_overlap(peaks[x], peaks[y]) for x, y in pairs]
xx = np.arange(len(pairs))
ax.bar(xx - 0.2, gene_sh, 0.4, color="#2B6CB0", label="genes shared", edgecolor="white", linewidth=0.4)
ax.bar(xx + 0.2, peak_sh, 0.4, color="#DD6B20", label="peaks shared (≥1 bp)", edgecolor="white", linewidth=0.4)
for i, (g, p) in enumerate(zip(gene_sh, peak_sh)):
    ax.text(i - 0.2, g + 8, str(g), ha="center", va="bottom", fontsize=6, color="#2B6CB0")
    ax.text(i + 0.2, p + 8, str(p), ha="center", va="bottom", fontsize=6, color="#DD6B20")
ax.set_xticks(xx); ax.set_xticklabels(plabels, rotation=30, ha="right")
ax.set_ylabel("Shared (n)")
ax.set_title("Convergence is at the gene level,\nnot at shared coordinates", fontsize=8.5)
ax.legend(loc="upper left")
ax.text(-0.14, 1.06, "b", transform=ax.transAxes, fontsize=14, fontweight="bold")
fig.suptitle(f"Cross-dataset memory-gene overlap ({COND})", fontsize=10.5, y=1.02)
save(fig, "Fig3_gene_overlap")

# ---- Figure 3c: UpSet of the four gene sets ----
from upsetplot import from_contents, UpSet
data = from_contents({d: genes[d] for d in DATASETS})
fig2 = plt.figure(figsize=(8.2, 4.4))
up = UpSet(data, sort_by="cardinality", show_counts=True, element_size=None,
           facecolor="#2D3748", shading_color="#EDF2F7")
up.plot(fig=fig2)
fig2.suptitle(f"Intersections of memory-gene sets across four datasets ({COND})\nrightmost-shared bar = four-way core",
              fontsize=9.5, y=1.0)
for ext in ("png", "pdf"):
    fig2.savefig(FIGDIR / f"Fig3c_gene_upset.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig3c_gene_upset.png")
plt.close(fig2)
print("four-way core genes:", len(set.intersection(*[genes[d] for d in DATASETS])))
