#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Permutation null tests for the four-way memory convergence.
  GENE null A : resample matched-size gene sets from the accessible-gene
                background (union of genes hit by any dataset) -> controls for
                the fact that all datasets sample the same accessible-gene pool.
  GENE null B : shuffle each dataset's peaks within their chromosome, re-annotate
                to nearest TSS gene -> controls for genomic/gene-density placement.
  TF   null   : resample matched-size enriched-TF sets from the tested-motif
                universe.
Reports observed vs null four-way (and pairwise) overlap, fold, empirical p, z.
"""
import sys, random, math
sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
import numpy as np
from pathlib import Path
from itertools import combinations
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, annotate_set

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
GP = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "untitled folder/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
MOT = BASE / "cross_dataset_comparison/motifs"
rng = np.random.default_rng(7)

def pstats(obs, null):
    null = np.array(null, float)
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    mu, sd = null.mean(), null.std()
    fold = obs / mu if mu > 0 else float("inf")
    z = (obs - mu) / sd if sd > 0 else float("inf")
    return mu, fold, p, z

print("Loading mm10 reference + annotating peaks ...")
by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")

# observed gene sets
gsets = {}
peaks = {}
for d in DATASETS:
    bed = load_bed(GP[d])
    _, loci = annotate_set(bed, by_chrom, t2l, locus_info)
    gsets[d] = loci
    peaks[d] = bed   # list of Interval
obs4 = len(set.intersection(*[gsets[d] for d in DATASETS]))
obs_pair = {f"{x}_{y}": len(gsets[x] & gsets[y]) for x, y in combinations(DATASETS, 2)}

# ---------------- GENE NULL A: resample from accessible-gene background ----------------
background = sorted(set().union(*[gsets[d] for d in DATASETS]))
B = len(background)
sizes = {d: len(gsets[d]) for d in DATASETS}
NA = 5000
null4_A = np.empty(NA, int)
bg = np.array(background, dtype=object)
for i in range(NA):
    sets = [set(rng.choice(B, size=sizes[d], replace=False)) for d in DATASETS]
    null4_A[i] = len(set.intersection(*sets))
# analytic hypergeometric expectation: B * prod(G_d/B)
exp_analytic = B * math.prod(sizes[d] / B for d in DATASETS)

# ---------------- GENE NULL B: genomic peak-shuffle + re-annotate ----------------
# chrom sizes (main chroms)
chromlen = {}
for ln in open(BASE / "mm10/chrom.sizes"):
    p = ln.split()
    if len(p) >= 2 and "_" not in p[0] and p[0] != "chrM":
        chromlen[p[0]] = int(p[1])
# per-chrom TSS positions + locus-int arrays
locus_to_int = {}
def lint(l):
    return locus_to_int.setdefault(l, len(locus_to_int))
chrom_pos, chrom_loc = {}, {}
for c, lst in by_chrom.items():
    chrom_pos[c] = np.array([t for t, _ in lst])
    chrom_loc[c] = np.array([lint(t2l.get(rid, rid)) for _, rid in lst])
# per-dataset peaks grouped by chrom with lengths
pk_by_chrom = {d: {} for d in DATASETS}
for d in DATASETS:
    for iv in peaks[d]:
        pk_by_chrom[d].setdefault(iv[0], []).append(iv[2] - iv[1])
    for c in pk_by_chrom[d]:
        pk_by_chrom[d][c] = np.array(pk_by_chrom[d][c])

def shuffled_gene_set(d):
    loci = set()
    for c, lens in pk_by_chrom[d].items():
        if c not in chrom_pos or c not in chromlen:
            continue
        pos = chrom_pos[c]; locarr = chrom_loc[c]
        hi = np.maximum(chromlen[c] - lens, 1)
        starts = (rng.random(len(lens)) * hi).astype(np.int64)
        centers = starts + lens // 2
        idx = np.searchsorted(pos, centers)
        idx_lo = np.clip(idx - 1, 0, len(pos) - 1)
        idx_hi = np.clip(idx, 0, len(pos) - 1)
        d_lo = np.abs(centers - pos[idx_lo]); d_hi = np.abs(centers - pos[idx_hi])
        nearest = np.where(d_lo <= d_hi, idx_lo, idx_hi)
        loci.update(locarr[nearest].tolist())
    return loci

NB = 1000
null4_B = np.empty(NB, int)
for i in range(NB):
    sets = [shuffled_gene_set(d) for d in DATASETS]
    null4_B[i] = len(set.intersection(*sets))

# ---------------- TF NULL: resample from tested-motif universe ----------------
def enriched(path, q=0.05):
    s = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try: qq = float(f[4])
            except ValueError: continue
            s.add(f[0].split("(")[0].split("/")[0].strip())  # tf name
            if qq >= q:  # still record name in universe but mark non-enriched below
                pass
    return s
def enriched_only(path, q=0.05):
    s = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try: qq = float(f[4])
            except ValueError: continue
            if qq < q: s.add(f[0].split("(")[0].split("/")[0].strip())
    return s
tf_universe = sorted(set().union(*[enriched(MOT/d/"knownResults.txt", q=1.0) for d in DATASETS]))
U = len(tf_universe)
tf_sets = {d: enriched_only(MOT/d/"knownResults.txt") for d in DATASETS}
tf_obs4 = len(set.intersection(*[tf_sets[d] for d in DATASETS]))
tf_sizes = {d: len(tf_sets[d]) for d in DATASETS}
NT = 20000
null4_T = np.empty(NT, int)
for i in range(NT):
    sets = [set(rng.choice(U, size=tf_sizes[d], replace=False)) for d in DATASETS]
    null4_T[i] = len(set.intersection(*sets))
tf_exp_analytic = U * math.prod(tf_sizes[d]/U for d in DATASETS)

# ---------------- report ----------------
print("\n================ GENE four-way overlap ================")
print(f"observed four-way shared genes: {obs4}")
print(f"gene-locus universe (collapsed): {len(locus_to_int)} ; accessible background B = {B}")
print(f"set sizes: " + ", ".join(f"{d}={sizes[d]}" for d in DATASETS))
mu, fold, p, z = pstats(obs4, null4_A)
print(f"  NULL A (resample from accessible background, N={NA}):")
print(f"    expected≈{mu:.2f} (analytic {exp_analytic:.2f}) | fold={fold:.1f}x | z={z:.1f} | p={p:.2e}")
mu, fold, p, z = pstats(obs4, null4_B)
print(f"  NULL B (genomic peak-shuffle, N={NB}):")
print(f"    expected≈{mu:.2f} | fold={fold:.1f}x | z={z:.1f} | p={p:.2e}")

print("\n  pairwise gene overlap vs accessible-background null (analytic expectation):")
for x, y in combinations(DATASETS, 2):
    exp = B * (sizes[x]/B) * (sizes[y]/B)
    print(f"    {x}∩{y}: obs {obs_pair[f'{x}_{y}']:4d}  exp≈{exp:6.1f}  fold {obs_pair[f'{x}_{y}']/exp:4.1f}x")

print("\n================ TF four-way overlap ================")
print(f"observed four-way shared TFs: {tf_obs4}")
print(f"TF universe (unique names tested): {U} ; set sizes: " + ", ".join(f"{d}={tf_sizes[d]}" for d in DATASETS))
mu, fold, p, z = pstats(tf_obs4, null4_T)
print(f"  NULL (resample from tested-motif universe, N={NT}):")
print(f"    expected≈{mu:.2f} (analytic {tf_exp_analytic:.2f}) | fold={fold:.1f}x | z={z:.1f} | p={p:.2e}")
print("\n[caveat] TF motifs are correlated (motif families), so the resampling null is")
print("anti-conservative; treat the TF p-value as an upper-bound-significance estimate.")
