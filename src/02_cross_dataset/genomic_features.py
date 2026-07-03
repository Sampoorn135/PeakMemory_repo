#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genomic-feature distribution of the memory peaks (all four studies).
Each StrictMemory50_Single peak is assigned to a genomic feature by its centre,
using HOMER's mm10.basic.annotation (a non-overlapping genome tiling), then we
report the % distribution and the log2(observed/expected) enrichment vs the
genome-wide feature coverage. Distal (intron + intergenic) accessible peaks are
the candidate-enhancer fraction.
"""
import bisect, math, json
from pathlib import Path

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
ANN = BASE / "mm10/mm10.basic.annotation"
PATHS = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
DATASETS = list(PATHS)
OUT = BASE / "cross_dataset_comparison/genomic_features"
OUT.mkdir(parents=True, exist_ok=True)

CODE2CAT = {"P": "Promoter-TSS", "5UTR": "5'UTR", "E": "Exon", "I": "Intron",
            "3UTR": "3'UTR", "TTS": "TTS", "N": "Intergenic",
            "ncRNA": "ncRNA/other", "pseudo": "ncRNA/other", "miRNA": "ncRNA/other",
            "snoRNA": "ncRNA/other", "rRNA": "ncRNA/other", "scRNA": "ncRNA/other",
            "snRNA": "ncRNA/other"}
CATS = ["Promoter-TSS", "5'UTR", "Exon", "Intron", "3'UTR", "TTS", "ncRNA/other", "Intergenic"]

# ---- load annotation tiling + genome bp per category ----
seg_start, seg_end, seg_cat = {}, {}, {}
genome_bp = {c: 0 for c in CATS}
with open(ANN, encoding="utf-8", errors="replace") as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        chrom, s, e, code = f[1], f[2], f[3], f[5]
        try:
            s = int(s); e = int(e)
        except ValueError:
            continue
        cat = CODE2CAT.get(code)
        if cat is None:
            continue
        seg_start.setdefault(chrom, []).append(s)
        seg_end.setdefault(chrom, []).append(e)
        seg_cat.setdefault(chrom, []).append(cat)
        genome_bp[cat] += (e - s)
# sort segments per chrom by start
for c in seg_start:
    order = sorted(range(len(seg_start[c])), key=lambda i: seg_start[c][i])
    seg_start[c] = [seg_start[c][i] for i in order]
    seg_end[c] = [seg_end[c][i] for i in order]
    seg_cat[c] = [seg_cat[c][i] for i in order]
# HOMER's intergenic ("N") segments carry placeholder coordinates, so the
# genome baseline is derived from real chromosome sizes: intergenic bp = genome
# minus all annotated (intron/exon/promoter/UTR/TTS/ncRNA) bp.
gsize = 0
for ln in open(BASE / "mm10/chrom.sizes"):
    p = ln.split()
    if len(p) < 2 or "_" in p[0] or p[0] == "chrM":
        continue
    gsize += int(p[1])
nonintergenic = sum(genome_bp[c] for c in CATS if c != "Intergenic")
genome_bp["Intergenic"] = max(1, gsize - nonintergenic)
TOTAL_GENOME = gsize
genome_frac = {c: genome_bp[c] / TOTAL_GENOME for c in CATS}

def feature_of(chrom, center):
    starts = seg_start.get(chrom)
    if not starts:
        return None
    i = bisect.bisect_right(starts, center) - 1
    if i < 0:
        return None
    if seg_end[chrom][i] > center:
        return seg_cat[chrom][i]
    return None

# ---- assign peaks ----
results = {}
for d in DATASETS:
    counts = {c: 0 for c in CATS}
    n = unassigned = 0
    for ln in open(PATHS[d]):
        f = ln.split("\t")
        chrom, s, e = f[0], int(f[1]), int(f[2])
        cat = feature_of(chrom, (s + e) // 2)
        if cat is None:
            unassigned += 1; continue
        counts[cat] += 1; n += 1
    pct = {c: 100 * counts[c] / n for c in CATS}
    enr = {c: math.log2((counts[c]/n) / genome_frac[c]) if counts[c] > 0 else float("nan") for c in CATS}
    results[d] = {"counts": counts, "n": n, "pct": pct, "log2enr": enr,
                  "distal_pct": pct["Intron"] + pct["Intergenic"],
                  "promoter_pct": pct["Promoter-TSS"]}

json.dump({"results": results, "genome_frac": genome_frac}, open(OUT / "genomic_features.json", "w"), indent=2)

# ---- markdown table ----
L = ["# Genomic-feature distribution of memory peaks (StrictMemory 50% Single)\n",
     "Each peak assigned to a feature by its centre (HOMER mm10.basic.annotation). "
     "Distal = Intron + Intergenic ≈ candidate enhancers.\n",
     "## Percent of peaks per feature\n",
     "| Feature | " + " | ".join(DATASETS) + " | genome % |",
     "| --- | " + " | ".join(["---:"]*len(DATASETS)) + " | ---: |"]
for c in CATS:
    L.append(f"| {c} | " + " | ".join(f"{results[d]['pct'][c]:.1f}%" for d in DATASETS) +
             f" | {100*genome_frac[c]:.1f}% |")
L.append(f"| **Distal (enhancer-like)** | " +
         " | ".join(f"**{results[d]['distal_pct']:.1f}%**" for d in DATASETS) + " | — |")
L.append("\n## log2(observed / expected) enrichment vs genome\n")
L.append("| Feature | " + " | ".join(DATASETS) + " |")
L.append("| --- | " + " | ".join(["---:"]*len(DATASETS)) + " |")
for c in CATS:
    L.append(f"| {c} | " + " | ".join(f"{results[d]['log2enr'][c]:+.2f}" for d in DATASETS) + " |")
(OUT / "genomic_features.md").write_text("\n".join(L))

print("=== % distribution ===")
for c in CATS:
    print(f"  {c:14s} " + "  ".join(f"{d}:{results[d]['pct'][c]:5.1f}%" for d in DATASETS))
print("=== distal (intron+intergenic, enhancer-like) ===")
for d in DATASETS:
    print(f"  {d}: {results[d]['distal_pct']:.1f}%  (promoter {results[d]['promoter_pct']:.1f}%)")
print("wrote", OUT / "genomic_features.md")
