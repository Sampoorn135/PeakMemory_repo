#!/usr/bin/env python3
"""Write RefSeq gene lists for all pairwise + three-way memory-gene overlaps."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
from pathlib import Path
from itertools import combinations
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, annotate_set

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
PATHS = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory",
    "pancreas": BASE / "epithelial_pancreas/pancreas_analysis/results/04_final_memory",
}
DATASETS = list(PATHS)
COND = "50_Single"
OUT = BASE / "cross_dataset_comparison/overlap_genes"
OUT.mkdir(parents=True, exist_ok=True)

by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")
gsets = {}
for d in DATASETS:
    _, loci = annotate_set(load_bed(PATHS[d] / f"StrictMemory{COND}.bed"), by_chrom, t2l, locus_info)
    gsets[d] = loci

def reps(loci):
    return sorted({locus_info[l][4] for l in loci if l in locus_info})

summary = []
for combo in list(combinations(DATASETS, 2)) + list(combinations(DATASETS, 3)):
    inter = set.intersection(*[gsets[d] for d in combo])
    name = "_".join(combo)
    (OUT / f"{name}.txt").write_text("\n".join(reps(inter)) + "\n")
    summary.append((name, len(inter)))

for name, n in summary:
    print(f"  {name:28s} {n}")
print("wrote", len(summary), "lists to", OUT)
