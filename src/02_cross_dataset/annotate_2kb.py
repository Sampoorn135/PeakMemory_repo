#!/usr/bin/env python3
"""2 kb TSS-distance gene annotation: per-dataset gene lists + all overlaps."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
from pathlib import Path
from itertools import combinations
from annotate_and_overlap import load_tss, build_gene_loci, load_bed, nearest_tss

CAP = 2_000
BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
GP = {
    "IMQ": BASE / "memory_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "HFSC": BASE / "hfsc/hfsc_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "distil": BASE / "distil/distil_analysis/results/04_final_memory/StrictMemory50_Single.bed",
    "pancreas": BASE / "untitled folder/pancreas_analysis/results/04_final_memory/StrictMemory50_Single.bed",
}
OUT = BASE / "cross_dataset_comparison/genes_2kb"
(OUT / "overlap_lists").mkdir(parents=True, exist_ok=True)

by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")

def reps(loci):
    return sorted({locus_info[l][4] for l in loci if l in locus_info})

gsets = {}
print(f"=== promoter-proximal annotation (|dist to TSS| < {CAP} bp) ===")
for d in DATASETS:
    bed = load_bed(GP[d]); loci = set(); asg = 0
    for chrom, s, e in bed:
        hit = nearest_tss(by_chrom, chrom, (s+e)//2)
        if hit and abs(hit[1]) < CAP:
            loci.add(t2l.get(hit[0], hit[0])); asg += 1
    gsets[d] = loci
    (OUT / f"{d}_genes_2kb.txt").write_text("\n".join(reps(loci)) + "\n")
    print(f"  {d:9s} peaks={len(bed):5d}  promoter-assigned={asg:4d} ({100*asg/len(bed):.1f}%)  gene loci={len(loci)}")

print("\n=== gene overlaps (2 kb) ===")
rows = []
for combo in list(combinations(DATASETS, 2)) + list(combinations(DATASETS, 3)) + [tuple(DATASETS)]:
    inter = set.intersection(*[gsets[d] for d in combo])
    name = "_".join(combo)
    (OUT / "overlap_lists" / f"{name}.txt").write_text("\n".join(reps(inter)) + "\n")
    rows.append((" ∩ ".join(combo), len(combo), len(inter)))
    print(f"  {' ∩ '.join(combo):42s} {len(inter)}")

L = [f"# Promoter-proximal (<{CAP} bp to TSS) gene sets + overlaps\n",
     "Peaks assigned to a gene only if their centre is within 2 kb of a TSS.\n",
     "| Dataset | gene loci (promoter) |", "| --- | ---: |"]
for d in DATASETS:
    L.append(f"| {d} | {len(gsets[d])} |")
L += ["\n## Overlaps\n| Combination | order | shared genes |", "| --- | ---: | ---: |"]
for name, k, n in rows:
    L.append(f"| {name} | {k} | {n} |")
(OUT / "summary_2kb.md").write_text("\n".join(L))
print("\nwrote gene lists + overlaps to", OUT)
