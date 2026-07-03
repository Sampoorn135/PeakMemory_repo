#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 Gene-level overlap of memory peaks between the two datasets (IMQ vs HFSC)
===============================================================================

Runs AFTER HOMER annotatePeaks.pl has produced annotated tables in
gene_annotation/homer_out/ (one per dataset x stringent condition).

For each stringent condition it:
  1. parses the HOMER table and takes the assigned gene per peak (HOMER's
     "Gene Name" = nearest-TSS gene; optionally restricted to a TSS distance),
  2. builds the set of memory-associated genes for IMQ and for HFSC,
  3. computes the gene-level overlap: shared genes, dataset-unique genes,
     counts, and the gene Jaccard index,
  4. writes the shared / unique gene lists and a markdown report into a NEW
     folder (gene_overlap/).

It also produces a POOLED comparison (union of genes across all four conditions
per dataset). Gene-level overlap is more biologically meaningful than raw peak
coordinate overlap: two datasets can regulate the same gene through different
(non-overlapping) regulatory elements.
===============================================================================
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

VARIANTS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]
DATASETS = ["IMQ", "HFSC"]


def find_col(header: List[str], *needles: str) -> Optional[int]:
    """Return index of the first header cell containing any needle (case-insensitive)."""
    low = [h.lower() for h in header]
    for i, h in enumerate(low):
        if all(n.lower() in h for n in needles):
            return i
    return None


def parse_homer_genes(path: Path, max_dist: Optional[int]) -> Tuple[Set[str], int, int]:
    """
    Parse a HOMER annotated table; return (gene_set, n_peaks, n_assigned).

    max_dist: if set, keep a peak's gene only when |Distance to TSS| <= max_dist.
              if None, keep the nearest gene regardless of distance (HOMER default).
    """
    genes: Set[str] = set()
    n_peaks = 0
    n_assigned = 0
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        gi = find_col(header, "gene", "name")          # "Gene Name"
        di = find_col(header, "distance", "tss")        # "Distance to TSS"
        if gi is None:
            raise ValueError(f"No 'Gene Name' column in {path.name}; header={header[:20]}")
        for row in reader:
            if len(row) <= gi:
                continue
            n_peaks += 1
            gene = row[gi].strip()
            if not gene or gene.upper() in {"NA", "N/A", ""}:
                continue
            if max_dist is not None and di is not None and di < len(row):
                try:
                    if abs(int(float(row[di]))) > max_dist:
                        continue
                except ValueError:
                    pass
            genes.add(gene)
            n_assigned += 1
    return genes, n_peaks, n_assigned


def jaccard(a: Set[str], b: Set[str]) -> float:
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(description="Gene-level overlap of memory peaks across datasets.")
    ap.add_argument("--homer-dir", required=True, type=Path,
                    help="gene_annotation/homer_out containing *.annotated.txt")
    ap.add_argument("--out-dir", required=True, type=Path, help="output folder (created)")
    ap.add_argument("--max-dist", type=int, default=None,
                    help="optional |distance to TSS| cap in bp (default: none = nearest gene)")
    args = ap.parse_args()

    out = args.out_dir
    (out / "gene_lists").mkdir(parents=True, exist_ok=True)
    (out / "shared_genes").mkdir(parents=True, exist_ok=True)

    # ---- load gene sets -------------------------------------------------
    gene_sets: Dict[str, Dict[str, Set[str]]] = {ds: {} for ds in DATASETS}
    peak_counts: Dict[Tuple[str, str], Tuple[int, int]] = {}
    missing = []
    for ds in DATASETS:
        for v in VARIANTS:
            p = args.homer_dir / f"{ds}_StrictMemory{v}.annotated.txt"
            if not p.exists():
                missing.append(p.name)
                continue
            genes, npk, nas = parse_homer_genes(p, args.max_dist)
            gene_sets[ds][v] = genes
            peak_counts[(ds, v)] = (npk, nas)
            (out / "gene_lists" / f"{ds}_{v}_genes.txt").write_text(
                "\n".join(sorted(genes)) + "\n")
    if missing:
        raise SystemExit("Missing HOMER outputs (run the script first):\n  " +
                         "\n  ".join(missing))

    # ---- per-condition overlap -----------------------------------------
    rows = []
    for v in VARIANTS:
        A, B = gene_sets["IMQ"][v], gene_sets["HFSC"][v]
        shared = sorted(A & B)
        (out / "shared_genes" / f"Shared_{v}.txt").write_text("\n".join(shared) + "\n")
        (out / "shared_genes" / f"IMQ_only_{v}.txt").write_text("\n".join(sorted(A - B)) + "\n")
        (out / "shared_genes" / f"HFSC_only_{v}.txt").write_text("\n".join(sorted(B - A)) + "\n")
        rows.append({
            "variant": v, "imq_genes": len(A), "hfsc_genes": len(B),
            "shared": len(shared), "imq_only": len(A - B), "hfsc_only": len(B - A),
            "jaccard": jaccard(A, B),
            "imq_pct": 100.0 * len(A & B) / len(A) if A else 0.0,
            "hfsc_pct": 100.0 * len(A & B) / len(B) if B else 0.0,
        })

    # ---- pooled overlap (union across the 4 conditions) -----------------
    A_pool = set().union(*gene_sets["IMQ"].values())
    B_pool = set().union(*gene_sets["HFSC"].values())
    pooled = {
        "imq_genes": len(A_pool), "hfsc_genes": len(B_pool),
        "shared": len(A_pool & B_pool),
        "imq_only": len(A_pool - B_pool), "hfsc_only": len(B_pool - A_pool),
        "jaccard": jaccard(A_pool, B_pool),
    }
    (out / "shared_genes" / "Shared_POOLED.txt").write_text(
        "\n".join(sorted(A_pool & B_pool)) + "\n")

    # ---- report ---------------------------------------------------------
    L: List[str] = []
    a = L.append
    a("# Gene-Level Overlap of Memory Peaks: IMQ vs HFSC\n")
    a("Final `StrictMemory` peaks were annotated to mm10 genes with HOMER "
      "`annotatePeaks.pl` (nearest-TSS gene)" +
      (f", keeping peaks within ±{args.max_dist:,} bp of a TSS." if args.max_dist else
       " (no distance cut-off).") + "\n")
    a("Overlap here is **gene-level**: two datasets share a gene if a memory peak "
      "in each is assigned to that gene — even if the underlying peaks do not overlap.\n")

    a("## Peaks annotated (per set)\n")
    a("| Dataset | Condition | Peaks | Peaks with a gene | Unique genes |")
    a("| --- | --- | ---: | ---: | ---: |")
    for ds in DATASETS:
        for v in VARIANTS:
            npk, nas = peak_counts[(ds, v)]
            a(f"| {ds} | {v} | {npk:,} | {nas:,} | {len(gene_sets[ds][v]):,} |")
    a("")

    a("## Gene overlap by stringent condition\n")
    a("| Condition | IMQ genes | HFSC genes | Shared | IMQ-only | HFSC-only | "
      "Shared / IMQ | Shared / HFSC | Jaccard |")
    a("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in rows:
        a(f"| {r['variant']} | {r['imq_genes']:,} | {r['hfsc_genes']:,} | "
          f"**{r['shared']:,}** | {r['imq_only']:,} | {r['hfsc_only']:,} | "
          f"{r['imq_pct']:.1f}% | {r['hfsc_pct']:.1f}% | {r['jaccard']:.3f} |")
    a("")

    a("## Pooled overlap (union of genes across all 4 conditions)\n")
    a("| | IMQ | HFSC | Shared | IMQ-only | HFSC-only | Jaccard |")
    a("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    a(f"| genes | {pooled['imq_genes']:,} | {pooled['hfsc_genes']:,} | "
      f"**{pooled['shared']:,}** | {pooled['imq_only']:,} | {pooled['hfsc_only']:,} | "
      f"{pooled['jaccard']:.3f} |")
    a("")
    a("## Files\n")
    a("- `gene_lists/` — full gene list per dataset × condition.")
    a("- `shared_genes/Shared_<condition>.txt` — shared genes per condition; "
      "`Shared_POOLED.txt` — shared across pooled sets.")
    a("- `shared_genes/IMQ_only_*` / `HFSC_only_*` — dataset-specific genes.")
    a("")
    (out / "Gene_Overlap_Report.md").write_text("\n".join(L))
    print(f"Wrote {out/'Gene_Overlap_Report.md'}")
    print("Per-condition shared gene counts:",
          {r["variant"]: r["shared"] for r in rows}, " pooled:", pooled["shared"])


if __name__ == "__main__":
    main()
