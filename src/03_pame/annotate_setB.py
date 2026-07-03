#!/usr/bin/env python3
"""Nearest-TSS gene annotation of Set B (AP-1-primed, ChIP-bound) regions."""
import sys; sys.path.insert(0, "/sessions/happy-blissful-fermi/mnt/Memory_peak")
from pathlib import Path
from annotate_and_overlap import load_tss, build_gene_loci, nearest_tss

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak")
SETB = BASE / "cross_dataset_comparison/ap1_primed/SetB_memory_AP1motif_bound.bed"
OUT = BASE / "cross_dataset_comparison/ap1_primed"

by_chrom = load_tss(BASE / "mm10/mm10.tss")
t2l, locus_info = build_gene_loci(BASE / "mm10/mm10.rna")

rows = []; loci = set(); loci_100k = set()
bins = {"promoter <2kb": 0, "2-50kb": 0, "50-100kb": 0, ">100kb": 0}
n = 0
for ln in open(SETB):
    f = ln.split("\t"); c, s, e = f[0], int(f[1]), int(f[2]); n += 1
    hit = nearest_tss(by_chrom, c, (s + e) // 2)
    if hit is None:
        continue
    rid, dist = hit; ad = abs(dist); lid = t2l.get(rid, rid)
    rep = locus_info[lid][4] if lid in locus_info else rid
    loci.add(lid)
    if ad <= 100_000:
        loci_100k.add(lid)
    bins["promoter <2kb" if ad < 2000 else "2-50kb" if ad < 50000 else "50-100kb" if ad <= 100000 else ">100kb"] += 1
    rows.append((c, s, e, rep, dist))

def reps(ls):
    return sorted({locus_info[l][4] for l in ls if l in locus_info})

(OUT / "SetB_genes_nearestTSS.txt").write_text("\n".join(reps(loci)) + "\n")
(OUT / "SetB_genes_within100kb.txt").write_text("\n".join(reps(loci_100k)) + "\n")
with open(OUT / "SetB_peak_gene.tsv", "w") as fh:
    fh.write("chrom\tstart\tend\tnearest_refseq\tdist_to_tss\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

print(f"Set B regions: {n}")
print(f"unique nearest genes (loci): {len(loci)}   | within 100 kb: {len(loci_100k)}")
print("distance to nearest TSS:")
for k, v in bins.items():
    print(f"   {k:14s} {v:4d}  ({100*v/n:.0f}%)")
print("wrote SetB_genes_nearestTSS.txt, SetB_genes_within100kb.txt, SetB_peak_gene.tsv")
