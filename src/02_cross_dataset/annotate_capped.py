#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-annotate memory peaks to genes with a TSS-distance cap, then redo the
gene-overlap + permutation null on the distance-capped gene sets.
A peak is assigned to a gene only if |peak-centre - nearest TSS| < CAP (1000 kb).
Peaks beyond the cap are dropped; only capped assignments are used downstream.
"""
import sys, math
sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
import numpy as np
from pathlib import Path
from itertools import combinations
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, nearest_tss

CAP = 1_000_000           # 1000 kb
BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
GP = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "untitled folder/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
OUT = BASE / "cross_dataset_comparison/genes_1Mb"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")

# ---- capped annotation per dataset ----
gsets, peaks, stats = {}, {}, {}
for d in DATASETS:
    bed = load_bed(GP[d]); peaks[d] = bed
    loci = set(); assigned = 0; rows = []
    for chrom, s, e in bed:
        c = (s + e) // 2
        hit = nearest_tss(by_chrom, chrom, c)
        if hit is None:
            rows.append((chrom, s, e, "NA", "")); continue
        rid, dist = hit
        if abs(dist) < CAP:
            lid = t2l.get(rid, rid); loci.add(lid); assigned += 1
            rows.append((chrom, s, e, rid, dist, lid))
        else:
            rows.append((chrom, s, e, "NA(>cap)", dist))
    gsets[d] = loci
    stats[d] = (len(bed), assigned, len(loci))
    reps = sorted({locus_info[l][4] for l in loci if l in locus_info})
    (OUT / f"{d}_genes_1Mb.txt").write_text("\n".join(reps) + "\n")
    with open(OUT / f"{d}_peak_gene_1Mb.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tnearest_refseq\tdist_to_tss\tgene_locus_rep\n")
        for r in rows:
            rep = locus_info[r[5]][4] if len(r) > 5 and r[5] in locus_info else ""
            fh.write("\t".join(str(x) for x in r[:5]) + f"\t{rep}\n")

print(f"=== distance-capped annotation (cap = {CAP:,} bp = {CAP//1000} kb) ===")
for d in DATASETS:
    npk, asg, ng = stats[d]
    print(f"  {d:9s} peaks={npk:5d}  assigned<cap={asg:5d} ({100*asg/npk:.1f}%)  dropped={npk-asg:4d}  gene loci={ng}")

# ---- overlaps ----
obs4 = len(set.intersection(*[gsets[d] for d in DATASETS]))
print("\n=== gene overlap (distance-capped) ===")
for x, y in combinations(DATASETS, 2):
    print(f"  {x}∩{y}: {len(gsets[x] & gsets[y])}")
print(f"  FOUR-WAY: {obs4}")

# ---- null A: resample from accessible background ----
bg = sorted(set().union(*[gsets[d] for d in DATASETS]))
B = len(bg); sizes = {d: len(gsets[d]) for d in DATASETS}
NA = 5000
nA = np.empty(NA, int)
for i in range(NA):
    nA[i] = len(set.intersection(*[set(rng.choice(B, size=sizes[d], replace=False)) for d in DATASETS]))
expA = nA.mean()

# ---- null B: peak-shuffle + re-annotate (with cap) ----
chromlen = {}
for ln in open(BASE / "mm10/chrom.sizes"):
    p = ln.split()
    if len(p) >= 2 and "_" not in p[0] and p[0] != "chrM":
        chromlen[p[0]] = int(p[1])
locus_to_int = {}
def lint(l): return locus_to_int.setdefault(l, len(locus_to_int))
cpos, cloc = {}, {}
for c, lst in by_chrom.items():
    cpos[c] = np.array([t for t, _ in lst]); cloc[c] = np.array([lint(t2l.get(rid, rid)) for _, rid in lst])
pkc = {d: {} for d in DATASETS}
for d in DATASETS:
    for chrom, s, e in peaks[d]:
        pkc[d].setdefault(chrom, []).append(e - s)
    for c in pkc[d]:
        pkc[d][c] = np.array(pkc[d][c])
def shuf(d):
    loci = set()
    for c, lens in pkc[d].items():
        if c not in cpos or c not in chromlen: continue
        pos = cpos[c]; la = cloc[c]
        starts = (rng.random(len(lens)) * np.maximum(chromlen[c]-lens, 1)).astype(np.int64)
        ctr = starts + lens // 2
        idx = np.searchsorted(pos, ctr)
        lo = np.clip(idx-1, 0, len(pos)-1); hi = np.clip(idx, 0, len(pos)-1)
        dlo = np.abs(ctr-pos[lo]); dhi = np.abs(ctr-pos[hi])
        nn = np.where(dlo <= dhi, lo, hi); dd = np.minimum(dlo, dhi)
        m = dd < CAP
        loci.update(la[nn][m].tolist())
    return loci
NB = 1000
nB = np.array([len(set.intersection(*[shuf(d) for d in DATASETS])) for _ in range(NB)])
expB = nB.mean()

def emp(obs, null):
    null = np.asarray(null, float)
    p = (np.sum(null >= obs)+1)/(len(null)+1)
    z = (obs-null.mean())/null.std() if null.std() > 0 else float("inf")
    return p, z, obs/null.mean() if null.mean() > 0 else float("inf")

pA, zA, fA = emp(obs4, nA); pB, zB, fB = emp(obs4, nB)
print("\n=== NULL on distance-capped four-way gene overlap ===")
print(f"  observed = {obs4}  | accessible background B = {B}")
print(f"  NULL A (resample accessible bg): exp≈{expA:.1f} fold={fA:.2f}x z={zA:.1f} p={pA:.2e}")
print(f"  NULL B (peak-shuffle, capped):   exp≈{expB:.1f} fold={fB:.2f}x z={zB:.1f} p={pB:.2e}")

# summary md
L = [f"# Distance-capped (<{CAP//1000} kb to TSS) gene annotation + overlap + null\n",
     "A peak is assigned to its nearest-TSS gene only if within 1000 kb; only these are used.\n",
     "## Annotation\n| Dataset | peaks | assigned (<cap) | dropped | gene loci |",
     "| --- | ---: | ---: | ---: | ---: |"]
for d in DATASETS:
    npk, asg, ng = stats[d]
    L.append(f"| {d} | {npk:,} | {asg:,} ({100*asg/npk:.1f}%) | {npk-asg} | {ng:,} |")
L += ["\n## Gene overlap\n| Combination | shared genes |", "| --- | ---: |"]
for x, y in combinations(DATASETS, 2):
    L.append(f"| {x} ∩ {y} | {len(gsets[x]&gsets[y])} |")
L.append(f"| **four-way** | **{obs4}** |")
L += ["\n## Four-way null", f"- observed = **{obs4}**; accessible background B = {B}",
      f"- resample-background null: expected ≈ {expA:.1f}, fold {fA:.2f}×, z {zA:.1f}, p {pA:.2e}",
      f"- peak-shuffle null: expected ≈ {expB:.1f}, fold {fB:.2f}×, z {zB:.1f}, p {pB:.2e}"]
(OUT / "summary_1Mb.md").write_text("\n".join(L))
print("\nwrote gene lists + summary to", OUT)
