#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pairwise + three-way overlap of the final StrictMemory sets across the three
datasets: IMQ (inflammation), HFSC (hair-follicle wound), distil (GSE221408).
Overlap is computed at TWO levels:
  * peak level   — >=1 bp genomic coordinate overlap (all mm10).
  * gene level   — nearest-TSS gene loci (collapsed isoforms), the biologically
                   meaningful layer (peaks can hit the same gene via different
                   non-overlapping elements).
Outputs counts (pairwise + triple) per stringent condition, plus the three-way
shared gene loci lists.
"""
from __future__ import annotations
import bisect
from pathlib import Path
from itertools import combinations

from annotate_and_overlap import (
    load_tss, build_gene_loci, nearest_tss, load_bed, annotate_set,
)

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
MM10 = BASE / "mm10"
DATASETS = {
    "IMQ":    BASE / "memory_analysis/results/04_final_memory",
    "HFSC":   BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory",
}
CONDS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]
OUT = BASE / "cross_dataset_comparison"
(OUT / "three_way_shared_genes").mkdir(parents=True, exist_ok=True)


# ---- peak-level helpers ---------------------------------------------------
def load_peaks(path):
    d = {}
    for ln in open(path):
        c, s, e = ln.split("\t")[:3]
        d.setdefault(c, []).append((int(s), int(e)))
    for c in d:
        d[c].sort()
    return d

def index(d):
    return {c: ([x[0] for x in v], v) for c, v in d.items()}

def overlaps_any(idx, c, s, e):
    if c not in idx:
        return False
    starts, arr = idx[c]
    i = bisect.bisect_left(starts, e) - 1
    while i >= 0:
        cs, ce = arr[i]
        if ce > s and cs < e:
            return True
        if starts[i] < s - 5_000_000:
            break
        i -= 1
    return False

def n_peaks_shared_by_all(target, others):
    """count peaks in `target` overlapping >=1 peak in EVERY other dataset."""
    idxs = [index(o) for o in others]
    n = 0
    for c, ivs in target.items():
        for s, e in ivs:
            if all(overlaps_any(ix, c, s, e) for ix in idxs):
                n += 1
    return n

def n_peaks_overlapping(a, b):
    ib = index(b); n = 0
    for c, ivs in a.items():
        for s, e in ivs:
            if overlaps_any(ib, c, s, e):
                n += 1
    return n


def jac(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def main():
    by_chrom = load_tss(MM10 / "mm10.tss")
    t2l, locus_info = build_gene_loci(MM10 / "mm10.rna")

    L = ["# Three-Way Overlap of Memory Peaks: IMQ vs HFSC vs distil (mm10)\n",
         "Final `StrictMemory` peaks compared across all three datasets. "
         "Peak level = ≥1 bp coordinate overlap; gene level = nearest-TSS gene "
         "loci (collapsed isoforms).\n"]
    a = L.append

    for cond in CONDS:
        peaks = {k: load_peaks(DATASETS[k] / f"StrictMemory{cond}.bed") for k in DATASETS}
        genes = {}
        for k in DATASETS:
            bed = load_bed(DATASETS[k] / f"StrictMemory{cond}.bed")
            _, loci = annotate_set(bed, by_chrom, t2l, locus_info)
            genes[k] = loci

        a(f"## {cond}\n")
        a("Peaks per set: " + ", ".join(
            f"{k} {sum(len(v) for v in peaks[k].values()):,}" for k in DATASETS) + "  ·  "
          "Gene loci: " + ", ".join(f"{k} {len(genes[k]):,}" for k in DATASETS) + "\n")

        # pairwise
        a("**Pairwise overlap**\n")
        a("| Pair | Peaks shared (A→B / B→A) | Genes shared | Gene Jaccard |")
        a("| --- | --- | ---: | ---: |")
        for x, y in combinations(DATASETS, 2):
            pxy = n_peaks_overlapping(peaks[x], peaks[y])
            pyx = n_peaks_overlapping(peaks[y], peaks[x])
            gshared = len(genes[x] & genes[y])
            a(f"| {x} ∩ {y} | {pxy} / {pyx} | {gshared:,} | {jac(genes[x], genes[y]):.3f} |")
        a("")

        # three-way
        tri_genes = genes["IMQ"] & genes["HFSC"] & genes["distil"]
        tri_peaks = {k: n_peaks_shared_by_all(peaks[k], [peaks[o] for o in DATASETS if o != k])
                     for k in DATASETS}
        a("**All three together (IMQ ∩ HFSC ∩ distil)**\n")
        a(f"- Peaks shared by all three (per dataset): " +
          ", ".join(f"{k} {tri_peaks[k]}" for k in DATASETS) + "\n")
        a(f"- **Gene loci shared by all three: {len(tri_genes):,}**\n")

        reps = sorted({locus_info[l][4] for l in tri_genes if l in locus_info})
        (OUT / "three_way_shared_genes" / f"Shared3way_{cond}.txt").write_text(
            "\n".join(reps) + "\n")
        a("")

    (OUT / "ThreeWay_Overlap_Report.md").write_text("\n".join(L))
    print("wrote", OUT / "ThreeWay_Overlap_Report.md")
    # console summary
    print("Three-way gene loci shared, per condition:")
    for cond in CONDS:
        n = len(open(OUT / "three_way_shared_genes" / f"Shared3way_{cond}.txt").read().split())
        print(f"  {cond}: {n}")


if __name__ == "__main__":
    main()
