#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overlap of HOMER known-motif (TF) enrichment across the three memory datasets.
Enriched TF = known motif with Benjamini q-value < 0.05 in that dataset's
50_Single memory regions. TF identity = name before the first '(' (motifs
collapsed to the TF). Computes pairwise + three-way overlap of TF sets.
"""
from pathlib import Path
from itertools import combinations

BASE = Path("/sessions/happy-blissful-fermi/mnt/Memory_peak/cross_dataset_comparison/motifs")
OUT = BASE
DATASETS = ["IMQ", "HFSC", "distil", "pancreas"]
QCUT = 0.05


def enriched_tfs(path):
    """Return {tf_name: best_q} for motifs with q < QCUT."""
    best = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        next(fh, None)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                q = float(f[4])
            except ValueError:
                continue
            if q >= QCUT:
                continue
            tf = f[0].split("(")[0].split("/")[0].strip()
            if not tf:
                continue
            if tf not in best or q < best[tf]:
                best[tf] = q
    return best


def main():
    tfs = {d: enriched_tfs(BASE / d / "knownResults.txt") for d in DATASETS}
    sets = {d: set(tfs[d]) for d in DATASETS}

    L = ["# Transcription-Factor (Motif) Overlap Across Memory Datasets\n",
         "HOMER known-motif enrichment on each dataset's `StrictMemory50_Single` "
         "regions (`findMotifsGenome.pl -size 200`, vertebrate motifs, GC-matched "
         f"background). **Enriched TF = Benjamini q < {QCUT}**; motifs collapsed to "
         "the TF name.\n"]
    a = L.append

    a("## Enriched TFs per dataset\n")
    a("| Dataset | # enriched TFs (q<0.05) |")
    a("| --- | ---: |")
    for d in DATASETS:
        a(f"| {d} | {len(sets[d])} |")
    a("")

    a("## Pairwise overlap\n")
    a("| Pair | Shared TFs | Jaccard | Shared list |")
    a("| --- | ---: | ---: | --- |")
    for x, y in combinations(DATASETS, 2):
        sh = sorted(sets[x] & sets[y])
        u = len(sets[x] | sets[y])
        jac = len(sets[x] & sets[y]) / u if u else 0
        a(f"| {x} ∩ {y} | {len(sh)} | {jac:.2f} | {', '.join(sh) if sh else '—'} |")
    a("")

    tri = sorted(set.intersection(*[sets[d] for d in DATASETS]))
    a("## All four together (IMQ ∩ HFSC ∩ distil ∩ pancreas)\n")
    a(f"**{len(tri)} TFs enriched in all four memory datasets:**\n")
    a(", ".join(tri) if tri else "_none_")
    a("")
    # union and per-dataset unique
    a("## Dataset-specific (enriched in only one)\n")
    for d in DATASETS:
        others = set().union(*[sets[o] for o in DATASETS if o != d])
        uniq = sorted(sets[d] - others)
        a(f"- **{d}-only ({len(uniq)})**: {', '.join(uniq) if uniq else '—'}")
    a("")

    (OUT / "TF_Overlap_Report_4way.md").write_text("\n".join(L))
    (OUT / "shared_TFs_all4.txt").write_text("\n".join(tri) + "\n")
    for x, y in combinations(DATASETS, 2):
        sh = sorted(sets[x] & sets[y])
        (OUT / f"shared_TFs_{x}_{y}.txt").write_text("\n".join(sh) + "\n")

    print("Enriched TFs:", {d: len(sets[d]) for d in DATASETS})
    print("Pairwise:", {f"{x}∩{y}": len(sets[x] & sets[y]) for x, y in combinations(DATASETS, 2)})
    print(f"Three-way ({len(tri)}):", ", ".join(tri))


if __name__ == "__main__":
    main()
