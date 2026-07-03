#!/usr/bin/env python3
"""Sweep the TSS-distance cap: effect on annotation, gene sets, four-way overlap + null."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
import numpy as np
from pathlib import Path
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, nearest_tss

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
GP = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "untitled folder/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
rng = np.random.default_rng(7)
by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")

# precompute per-peak (locus, |dist|) once
peakann = {d: [] for d in DATASETS}
npeaks = {}
for d in DATASETS:
    bed = load_bed(GP[d]); npeaks[d] = len(bed)
    for chrom, s, e in bed:
        hit = nearest_tss(by_chrom, chrom, (s+e)//2)
        if hit:
            peakann[d].append((t2l.get(hit[0], hit[0]), abs(hit[1])))

CAPS = [1_000, 5_000, 10_000, 50_000, 100_000, 500_000, 1_000_000]
print(f"{'cap':>9} | {'assigned%':>9} | " + " ".join(f"{d[:4]:>5}" for d in DATASETS) +
      f" | {'4way':>4} | {'exp':>5} | {'fold':>4} | {'p':>8}")
print("-"*78)
for CAP in CAPS:
    gsets = {d: {l for l, dd in peakann[d] if dd < CAP} for d in DATASETS}
    sizes = {d: len(gsets[d]) for d in DATASETS}
    asg = [sum(1 for _, dd in peakann[d] if dd < CAP) for d in DATASETS]
    apct = 100*sum(asg)/sum(npeaks[d] for d in DATASETS)
    obs4 = len(set.intersection(*[gsets[d] for d in DATASETS]))
    bg = sorted(set().union(*gsets.values())); B = len(bg)
    N = 3000
    nn = np.empty(N, int)
    for i in range(N):
        nn[i] = len(set.intersection(*[set(rng.choice(B, size=sizes[d], replace=False)) for d in DATASETS]))
    exp = nn.mean(); p = (np.sum(nn >= obs4)+1)/(N+1); fold = obs4/exp if exp > 0 else float('inf')
    cap_lab = f"{CAP//1000}kb" if CAP >= 1000 else f"{CAP}bp"
    print(f"{cap_lab:>9} | {apct:8.1f}% | " + " ".join(f"{sizes[d]:>5}" for d in DATASETS) +
          f" | {obs4:>4} | {exp:5.1f} | {fold:4.2f} | {p:8.2e}")
