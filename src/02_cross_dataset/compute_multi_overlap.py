#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pairwise + N-way overlap of the final StrictMemory sets across FOUR datasets:
IMQ, HFSC, distil, pancreas. Peak level (>=1 bp) and gene level (nearest-TSS
gene loci). Writes the four-way shared gene loci (RefSeq reps) per condition.
"""
from __future__ import annotations
import bisect
from pathlib import Path
from itertools import combinations
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, annotate_set

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
MM10 = BASE / "mm10"
DATASETS = {
    "IMQ":      BASE / "memory_analysis/results/04_final_memory",
    "HFSC":     BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil":   BASE / "distil/distil_analysis/results/04_final_memory",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory",
}
CONDS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]
OUT = BASE / "cross_dataset_comparison"
(OUT / "four_way_shared_genes").mkdir(parents=True, exist_ok=True)


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

def ov(idx, c, s, e):
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

def n_overlap(a, b):
    ib = index(b); n = 0
    for c, ivs in a.items():
        for s, e in ivs:
            if ov(ib, c, s, e):
                n += 1
    return n

def n_shared_all(target, others):
    idxs = [index(o) for o in others]
    n = 0
    for c, ivs in target.items():
        for s, e in ivs:
            if all(ov(ix, c, s, e) for ix in idxs):
                n += 1
    return n

def jac(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def main():
    by_chrom = load_tss(MM10 / "mm10.tss")
    t2l, locus_info = build_gene_loci(MM10 / "mm10.rna")
    names = list(DATASETS)

    L = ["# Four-Way Overlap of Memory Peaks: IMQ · HFSC · distil · pancreas (mm10)\n",
         "Final `StrictMemory` peaks across all four datasets. Peak level = ≥1 bp "
         "coordinate overlap; gene level = nearest-TSS gene loci.\n"]
    a = L.append
    for cond in CONDS:
        peaks = {k: load_peaks(DATASETS[k] / f"StrictMemory{cond}.bed") for k in names}
        genes = {}
        for k in names:
            _, loci = annotate_set(load_bed(DATASETS[k] / f"StrictMemory{cond}.bed"),
                                   by_chrom, t2l, locus_info)
            genes[k] = loci
        a(f"## {cond}\n")
        a("Peaks: " + ", ".join(f"{k} {sum(len(v) for v in peaks[k].values()):,}" for k in names) +
          "  ·  Gene loci: " + ", ".join(f"{k} {len(genes[k]):,}" for k in names) + "\n")
        a("**Pairwise gene overlap**\n")
        a("| Pair | Genes shared | Jaccard | Peaks shared (A→B) |")
        a("| --- | ---: | ---: | ---: |")
        for x, y in combinations(names, 2):
            a(f"| {x} ∩ {y} | {len(genes[x] & genes[y]):,} | {jac(genes[x], genes[y]):.3f} | {n_overlap(peaks[x], peaks[y])} |")
        a("")
        # four-way
        four_genes = set.intersection(*[genes[k] for k in names])
        four_peaks = {k: n_shared_all(peaks[k], [peaks[o] for o in names if o != k]) for k in names}
        a("**All four together**\n")
        a(f"- Peaks shared by all four: " + ", ".join(f"{k} {four_peaks[k]}" for k in names) + "\n")
        a(f"- **Gene loci shared by all four: {len(four_genes):,}**\n")
        reps = sorted({locus_info[l][4] for l in four_genes if l in locus_info})
        (OUT / "four_way_shared_genes" / f"Shared4way_{cond}.txt").write_text("\n".join(reps) + "\n")
        a("")
    (OUT / "FourWay_Overlap_Report.md").write_text("\n".join(L))
    print("wrote", OUT / "FourWay_Overlap_Report.md")
    for cond in CONDS:
        n = len(open(OUT / "four_way_shared_genes" / f"Shared4way_{cond}.txt").read().split())
        print(f"  four-way genes {cond}: {n}")


if __name__ == "__main__":
    main()
