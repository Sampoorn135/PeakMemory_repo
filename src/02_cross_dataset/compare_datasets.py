#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 Cross-dataset comparison of final StrictMemory peaks
===============================================================================

Compares the FINAL acquired-memory peaks (results/04_final_memory/StrictMemory*)
between two independent analyses, for each of the four stringent conditions:
    50_Single, 50_Reciprocal, 75_Single, 75_Reciprocal

Both datasets are on the SAME genome assembly (mm10, verified separately), so
peak coordinates are directly comparable.

For each condition we report, using a simple >=1 bp coordinate-overlap rule:
  * peaks in each dataset
  * dataset-A peaks that overlap >=1 dataset-B peak (and %)
  * dataset-B peaks that overlap >=1 dataset-A peak (and %)
  * number of overlapping peak pairs
  * Jaccard index on base pairs: shared_bp / union_bp
Because StrictMemory sets are internally non-overlapping (merged), the bp Jaccard
is exact. The overlapping A-peaks are also written out as a BED per condition.
===============================================================================
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

from memory_peak_pipeline import (
    BedSet, IntervalIndex, Interval, count_intervals, load_bed, logger,
    overlap_bp, setup_logging, write_bed,
)

VARIANTS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]


def total_bp(bedset: BedSet) -> int:
    return sum(iv.length for ivs in bedset.values() for iv in ivs)


def compare(a: BedSet, b: BedSet) -> Tuple[dict, BedSet]:
    """Compare two StrictMemory sets; return stats dict and A-peaks overlapping B."""
    b_index = IntervalIndex(b)
    a_index = IntervalIndex(a)

    a_overlapping: BedSet = {}
    n_a_hit = 0
    n_pairs = 0
    shared_bp = 0
    for chrom, ivs in a.items():
        for A in ivs:
            hits = b_index.query(A)
            if hits:
                n_a_hit += 1
                n_pairs += len(hits)
                shared_bp += sum(overlap_bp(A, h) for h in hits)
                a_overlapping.setdefault(chrom, []).append(Interval(chrom, A.start, A.end))

    # B peaks overlapping any A (other direction)
    n_b_hit = 0
    for chrom, ivs in b.items():
        for B in ivs:
            if a_index.query(B):
                n_b_hit += 1

    na, nb = count_intervals(a), count_intervals(b)
    union_bp = total_bp(a) + total_bp(b) - shared_bp
    stats = {
        "a_total": na,
        "b_total": nb,
        "a_overlapping_b": n_a_hit,
        "b_overlapping_a": n_b_hit,
        "a_pct": (100.0 * n_a_hit / na) if na else 0.0,
        "b_pct": (100.0 * n_b_hit / nb) if nb else 0.0,
        "pairs": n_pairs,
        "shared_bp": shared_bp,
        "jaccard_bp": (shared_bp / union_bp) if union_bp else 0.0,
    }
    return stats, a_overlapping


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare StrictMemory peaks across two datasets.")
    ap.add_argument("--dir-a", required=True, type=Path, help="results/04_final_memory of dataset A")
    ap.add_argument("--dir-b", required=True, type=Path, help="results/04_final_memory of dataset B")
    ap.add_argument("--name-a", default="A")
    ap.add_argument("--name-b", default="B")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    setup_logging(True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    isect_dir = args.out_dir / "overlapping_peaks"
    isect_dir.mkdir(exist_ok=True)

    rows = []
    for v in VARIANTS:
        a = load_bed(args.dir_a / f"StrictMemory{v}.bed")
        b = load_bed(args.dir_b / f"StrictMemory{v}.bed")
        logger.info("Comparing %s ...", v)
        stats, a_ov = compare(a, b)
        write_bed(a_ov, isect_dir / f"Overlap_{v}_{args.name_a}_in_{args.name_b}.bed")
        rows.append((v, stats))
        logger.info("  %s: %s=%d (%d overlap, %.1f%%) | %s=%d (%d overlap, %.1f%%) | Jaccard=%.4f",
                    v, args.name_a, stats["a_total"], stats["a_overlapping_b"], stats["a_pct"],
                    args.name_b, stats["b_total"], stats["b_overlapping_a"], stats["b_pct"],
                    stats["jaccard_bp"])

    # markdown report
    L = []
    a_ = L.append
    a_(f"# Cross-Dataset Memory-Peak Overlap: {args.name_a} vs {args.name_b}\n")
    a_("Final `StrictMemory` peaks compared by genomic coordinate (both mm10). ")
    a_("Overlap rule: **≥ 1 bp** shared.\n")
    a_(f"- Dataset A = **{args.name_a}**")
    a_(f"- Dataset B = **{args.name_b}**\n")
    a_("## Overlap by stringent condition\n")
    a_(f"| Condition | {args.name_a} peaks | {args.name_b} peaks | "
       f"{args.name_a}∩{args.name_b} (A-side) | {args.name_a}∩{args.name_b} (B-side) | "
       f"Overlapping pairs | Jaccard (bp) |")
    a_("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for v, s in rows:
        a_(f"| {v} | {s['a_total']:,} | {s['b_total']:,} | "
           f"{s['a_overlapping_b']:,} ({s['a_pct']:.1f}%) | "
           f"{s['b_overlapping_a']:,} ({s['b_pct']:.1f}%) | "
           f"{s['pairs']:,} | {s['jaccard_bp']:.4f} |")
    a_("")
    a_("## How to read this\n")
    a_(f"- **A-side overlap** = number of {args.name_a} peaks that hit ≥1 {args.name_b} peak.")
    a_(f"- **B-side overlap** = number of {args.name_b} peaks that hit ≥1 {args.name_a} peak.")
    a_("- The two sides differ because one peak in one set can overlap several in the other,")
    a_("  and peak widths differ between the datasets.")
    a_("- **Jaccard (bp)** = shared base pairs / union base pairs (0 = none, 1 = identical).")
    a_(f"- The actual overlapping {args.name_a} peaks per condition are in `overlapping_peaks/`.")
    a_("")
    (args.out_dir / "CrossDataset_Overlap_Report.md").write_text("\n".join(L))
    logger.info("Wrote report -> %s", args.out_dir / "CrossDataset_Overlap_Report.md")


if __name__ == "__main__":
    main()
