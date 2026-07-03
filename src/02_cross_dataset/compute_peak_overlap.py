#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peak-coordinate overlap across the four datasets for ALL combinations.
Builds a union of all StrictMemory50_Single peaks (merged across datasets); each
merged region's membership = the datasets whose peaks cover it. Reporting the
count of merged regions shared by each pair / triple / the four-way, and a
four-set Venn of peak regions (directly comparable to the gene/TF Venns).
"""
import bisect, json
from pathlib import Path
from itertools import combinations
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak/figures_src")
from fig_style import plt, DCOLOR, DATASETS, FIGDIR
from venn import venn
from matplotlib.colors import to_rgba

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
PATHS = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
OUT = BASE / "cross_dataset_comparison"

def load(p):
    d = {}
    for ln in open(p):
        f = ln.split("\t"); d.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for c in d: d[c].sort()
    return d
peaks = {d: load(PATHS[d]) for d in DATASETS}

# ---- union (merge) of all peaks across datasets -> merged regions ----
merged = {}  # chrom -> list of (s,e)
for d in DATASETS:
    for c, ivs in peaks[d].items():
        merged.setdefault(c, []).extend(ivs)
regions = []  # global list of (chrom, s, e)
for c in merged:
    ivs = sorted(merged[c]); cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e: cur_e = max(cur_e, e)
        else: regions.append((c, cur_s, cur_e)); cur_s, cur_e = s, e
    regions.append((c, cur_s, cur_e))

# ---- per-dataset interval index ----
def index(pk): return {c: ([x[0] for x in v], v) for c, v in pk.items()}
idx = {d: index(peaks[d]) for d in DATASETS}
def hits(ix, c, s, e):
    if c not in ix: return False
    st, arr = ix[c]; i = bisect.bisect_left(st, e) - 1
    while i >= 0:
        cs, ce = arr[i]
        if ce > s and cs < e: return True
        if st[i] < s - 5_000_000: break
        i -= 1
    return False

# ---- membership sets: S[d] = merged-region indices where d has a peak ----
S = {d: set() for d in DATASETS}
for r, (c, s, e) in enumerate(regions):
    for d in DATASETS:
        if hits(idx[d], c, s, e): S[d].add(r)

# ---- all-cases overlap counts ----
rows = []
for combo in list(combinations(DATASETS, 2)) + list(combinations(DATASETS, 3)) + [tuple(DATASETS)]:
    shared = len(set.intersection(*[S[d] for d in combo]))
    rows.append((" ∩ ".join(combo), len(combo), shared))

L = ["# Peak-coordinate overlap across datasets (StrictMemory 50% Single)\n",
     f"Union of all four datasets' memory peaks = **{len(regions):,} merged regions**. "
     "A region is 'shared' by a combination if every member dataset has a peak overlapping it (≥1 bp).\n",
     "| Combination | order | shared peak regions |", "| --- | ---: | ---: |"]
for name, k, n in rows:
    L.append(f"| {name} | {k} | {n} |")
(OUT / "peak_overlap.md").write_text("\n".join(L))

print("union merged regions:", len(regions))
for name, k, n in rows:
    print(f"  {name:42s} {n}")

# ---- four-set Venn of peak regions ----
fig, ax = plt.subplots(figsize=(6.4, 6.0))
venn({d: S[d] for d in DATASETS}, cmap=[to_rgba(DCOLOR[d], 0.45) for d in DATASETS],
     fontsize=8.5, legend_loc="upper left", ax=ax)
core = len(set.intersection(*[S[d] for d in DATASETS]))
ax.set_title("Peak-coordinate overlap across four datasets (50_Single)", fontsize=11, pad=12)
ax.text(0.5, -0.02, f"four-way core = {core} peak regions", transform=ax.transAxes,
        ha="center", fontsize=9, fontweight="bold", color="#1A202C")
leg = ax.get_legend()
if leg:
    for t, d in zip(leg.get_texts(), DATASETS): t.set_text(d); t.set_fontsize(8.5)
for ext in ("png", "pdf"):
    fig.savefig(FIGDIR / f"Fig11_peak_venn.{ext}", dpi=300, bbox_inches="tight")
print("wrote", FIGDIR / "Fig11_peak_venn.png")
plt.close(fig)
