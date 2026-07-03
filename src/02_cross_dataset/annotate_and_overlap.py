#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 Local (no-HOMER) peak->gene annotation + cross-dataset gene overlap
===============================================================================

Uses the HOMER mm10 GENOME data the user supplied (mm10/mm10.tss + mm10/mm10.rna)
purely as a coordinate reference. NO HOMER binary is run; everything is Python.

Why this is equivalent to HOMER's nearest-gene call:
  * mm10.tss gives every RefSeq transcript's TSS (interval centre) + strand.
  * Each peak is assigned to the gene whose TSS is nearest the peak centre —
    the same "Gene Name = nearest TSS" logic HOMER's annotatePeaks.pl uses by
    default.

Gene vs transcript:
  RefSeq lists many transcripts per gene. We collapse transcripts into GENE LOCI
  by merging same-strand overlapping gene bodies (mm10.rna), so the overlap is
  computed at gene level, not isoform level. Each locus is labelled by a
  representative RefSeq accession + its merged coordinates. (Gene SYMBOLS need
  HOMER's separate data/accession map, which was not provided; symbols can be
  added later without changing any counts.)

Outputs (new folder, default gene_overlap/):
  annotated/<DS>_<variant>.annotated.tsv   per-peak nearest gene + distance
  gene_lists/<DS>_<variant>_genes.txt      gene loci hit by that set
  shared_genes/Shared_<variant>.txt        loci shared by IMQ & HFSC (+ unique)
  Gene_Overlap_Report.md                   summary tables
===============================================================================
"""
from __future__ import annotations

import argparse
import bisect
from pathlib import Path
from typing import Dict, List, Optional, Tuple

VARIANTS = ["50_Single", "50_Reciprocal", "75_Single", "75_Reciprocal"]
DATASETS = ["IMQ", "HFSC"]


# ---------------------------------------------------------------------------
# Load TSS reference (mm10.tss): RefSeqID  chrom  start  end  strand(0=+,1=-)
# TSS = centre of the [start,end] promoter interval.
# ---------------------------------------------------------------------------
def load_tss(path: Path) -> Dict[str, List[Tuple[int, str]]]:
    by_chrom: Dict[str, List[Tuple[int, str]]] = {}
    tss_of: Dict[str, Tuple[str, int, int]] = {}
    n = 0
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        rid, chrom, s, e, strand = f[0], f[1], int(f[2]), int(f[3]), f[4]
        tss = (s + e) // 2
        by_chrom.setdefault(chrom, []).append((tss, rid))
        tss_of[rid] = (chrom, tss, 0 if strand == "0" else 1)
        n += 1
    for c in by_chrom:
        by_chrom[c].sort()
    print(f"[INFO] TSS reference: {n} transcripts across {len(by_chrom)} chromosomes")
    load_tss.tss_of = tss_of  # type: ignore[attr-defined]
    return by_chrom


# ---------------------------------------------------------------------------
# Collapse transcripts -> gene loci using gene bodies (mm10.rna):
#   RefSeqID  chrom  start  end  strand(+/-)  exon_structure...
# Merge overlapping bodies on the SAME strand; each merged span = one gene locus.
# ---------------------------------------------------------------------------
def build_gene_loci(rna_path: Path) -> Tuple[Dict[str, str], Dict[str, Tuple[str, int, int, str]]]:
    bodies: List[Tuple[str, str, int, int, str]] = []  # (rid,chrom,start,end,strand)
    for ln in open(rna_path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        rid, chrom, s, e, strand = f[0], f[1], int(f[2]), int(f[3]), f[4]
        bodies.append((rid, chrom, s, e, strand))

    # group by (chrom, strand), sort, merge overlapping intervals into loci
    from collections import defaultdict
    groups: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
    for rid, chrom, s, e, strand in bodies:
        groups[(chrom, strand)].append((s, e, rid))

    transcript_to_locus: Dict[str, str] = {}
    locus_info: Dict[str, Tuple[str, int, int, str]] = {}  # locus_id -> (chrom,start,end,strand)
    locus_n = 0
    for (chrom, strand), ivs in groups.items():
        ivs.sort()
        cur_s, cur_e, members = None, None, []
        def flush():
            nonlocal locus_n
            if not members:
                return
            locus_n += 1
            lid = f"L{locus_n}"
            rep = sorted(members, key=lambda r: r)[0]  # deterministic representative
            locus_info[lid] = (chrom, cur_s, cur_e, strand, rep)  # type: ignore
            for m in members:
                transcript_to_locus[m] = lid
        for s, e, rid in ivs:
            if cur_s is None:
                cur_s, cur_e, members = s, e, [rid]
            elif s <= cur_e:  # overlap -> same locus
                cur_e = max(cur_e, e)
                members.append(rid)
            else:
                flush()
                cur_s, cur_e, members = s, e, [rid]
        flush()
    print(f"[INFO] Collapsed {len(bodies)} transcripts -> {len(locus_info)} gene loci")
    return transcript_to_locus, locus_info


# ---------------------------------------------------------------------------
# Nearest-TSS assignment
# ---------------------------------------------------------------------------
def nearest_tss(by_chrom, chrom: str, center: int) -> Optional[Tuple[str, int]]:
    """Return (RefSeqID, signed_distance) of nearest TSS on chrom, or None."""
    arr = by_chrom.get(chrom)
    if not arr:
        return None
    positions = [p for p, _ in arr]
    i = bisect.bisect_left(positions, center)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            tss, rid = arr[j]
            d = center - tss
            if best is None or abs(d) < abs(best[1]):
                best = (rid, d)
    return best


def load_bed(path: Path) -> List[Tuple[str, int, int]]:
    out = []
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def annotate_set(bed: List[Tuple[str, int, int]], by_chrom, t2l, locus_info):
    """Return (rows, gene_locus_set). rows = per-peak annotation tuples."""
    rows = []
    loci = set()
    for chrom, s, e in bed:
        c = (s + e) // 2
        hit = nearest_tss(by_chrom, chrom, c)
        if hit is None:
            rows.append((chrom, s, e, "NA", "", "NA"))
            continue
        rid, dist = hit
        lid = t2l.get(rid, f"solo:{rid}")
        if lid not in locus_info and lid.startswith("solo:"):
            # transcript had no body record; treat as its own locus
            locus_info[lid] = (chrom, c, c, ".", rid)
        rep = locus_info[lid][4] if lid in locus_info else rid
        rows.append((chrom, s, e, rid, dist, lid))
        loci.add(lid)
    return rows, loci


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mm10-dir", required=True, type=Path)
    ap.add_argument("--imq-dir", required=True, type=Path)
    ap.add_argument("--hfsc-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--max-dist", type=int, default=None,
                    help="optional |distance to TSS| cap (bp). Default: none (nearest).")
    args = ap.parse_args()

    by_chrom = load_tss(args.mm10_dir / "mm10.tss")
    t2l, locus_info = build_gene_loci(args.mm10_dir / "mm10.rna")

    out = args.out_dir
    for sub in ("annotated", "gene_lists", "shared_genes"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    src = {"IMQ": args.imq_dir, "HFSC": args.hfsc_dir}
    gene_sets: Dict[str, Dict[str, set]] = {ds: {} for ds in DATASETS}
    peakn: Dict[Tuple[str, str], Tuple[int, int]] = {}

    for ds in DATASETS:
        for v in VARIANTS:
            bed = load_bed(src[ds] / f"StrictMemory{v}.bed")
            rows, loci = annotate_set(bed, by_chrom, t2l, locus_info)
            if args.max_dist is not None:
                keep = set()
                frows = []
                for r in rows:
                    if r[3] != "NA" and r[4] != "" and abs(int(r[4])) <= args.max_dist:
                        keep.add(r[5]); frows.append(r)
                    else:
                        frows.append((r[0], r[1], r[2], "NA", "", "NA"))
                rows, loci = frows, keep
            # write per-peak annotation
            with open(out / "annotated" / f"{ds}_{v}.annotated.tsv", "w") as fh:
                fh.write("chrom\tstart\tend\tnearest_refseq\tdist_to_tss\tgene_locus\tlocus_rep\tlocus_coord\n")
                for chrom, s, e, rid, dist, lid in rows:
                    if lid in locus_info:
                        lc, ls, le, lstr, rep = locus_info[lid]
                        coord = f"{lc}:{ls}-{le}({lstr})"
                    else:
                        rep, coord = rid, ""
                    fh.write(f"{chrom}\t{s}\t{e}\t{rid}\t{dist}\t{lid}\t{rep}\t{coord}\n")
            # gene list = representative RefSeq per locus
            reps = sorted({locus_info[l][4] for l in loci if l in locus_info})
            (out / "gene_lists" / f"{ds}_{v}_genes.txt").write_text("\n".join(reps) + "\n")
            gene_sets[ds][v] = loci
            n_assigned = sum(1 for r in rows if r[5] != "NA")
            peakn[(ds, v)] = (len(bed), n_assigned)

    # ---- per-condition overlap (by gene locus) -------------------------
    rows_rep = []
    for v in VARIANTS:
        A, B = gene_sets["IMQ"][v], gene_sets["HFSC"][v]
        shared = A & B
        def reps_of(s):
            return sorted({locus_info[l][4] for l in s if l in locus_info})
        (out / "shared_genes" / f"Shared_{v}.txt").write_text("\n".join(reps_of(shared)) + "\n")
        (out / "shared_genes" / f"IMQ_only_{v}.txt").write_text("\n".join(reps_of(A - B)) + "\n")
        (out / "shared_genes" / f"HFSC_only_{v}.txt").write_text("\n".join(reps_of(B - A)) + "\n")
        rows_rep.append({
            "v": v, "a": len(A), "b": len(B), "shared": len(shared),
            "ao": len(A - B), "bo": len(B - A), "jac": jaccard(A, B),
            "ap": 100*len(shared)/len(A) if A else 0, "bp": 100*len(shared)/len(B) if B else 0,
        })

    A_pool = set().union(*gene_sets["IMQ"].values())
    B_pool = set().union(*gene_sets["HFSC"].values())
    pooled = {"a": len(A_pool), "b": len(B_pool), "shared": len(A_pool & B_pool),
              "ao": len(A_pool - B_pool), "bo": len(B_pool - A_pool),
              "jac": jaccard(A_pool, B_pool)}
    (out / "shared_genes" / "Shared_POOLED.txt").write_text(
        "\n".join(sorted({locus_info[l][4] for l in (A_pool & B_pool) if l in locus_info})) + "\n")

    # ---- report ---------------------------------------------------------
    L: List[str] = []
    a = L.append
    a("# Gene-Level Overlap of Memory Peaks: IMQ vs HFSC (mm10)\n")
    a("Final `StrictMemory` peaks were annotated to the **nearest TSS** using the "
      "supplied HOMER mm10 reference (`mm10.tss` / `mm10.rna`), entirely in Python "
      "(no HOMER binary). Transcript isoforms were collapsed into **gene loci** "
      "(same-strand overlapping gene bodies), so overlap is gene-level." +
      (f" Peaks were kept only within ±{args.max_dist:,} bp of a TSS." if args.max_dist else
       " No distance cut-off (nearest gene always assigned).") + "\n")
    a("> Gene identifiers are representative **RefSeq accessions** (this HOMER genome "
      "folder has no symbol map). Counts are final; symbols can be relabelled later.\n")

    a("## Peaks annotated\n")
    a("| Dataset | Condition | Peaks | Assigned a gene | Gene loci hit |")
    a("| --- | --- | ---: | ---: | ---: |")
    for ds in DATASETS:
        for v in VARIANTS:
            npk, nas = peakn[(ds, v)]
            a(f"| {ds} | {v} | {npk:,} | {nas:,} | {len(gene_sets[ds][v]):,} |")
    a("")

    a("## Gene overlap by stringent condition\n")
    a("| Condition | IMQ genes | HFSC genes | Shared | IMQ-only | HFSC-only | "
      "Shared/IMQ | Shared/HFSC | Jaccard |")
    a("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in rows_rep:
        a(f"| {r['v']} | {r['a']:,} | {r['b']:,} | **{r['shared']:,}** | "
          f"{r['ao']:,} | {r['bo']:,} | {r['ap']:.1f}% | {r['bp']:.1f}% | {r['jac']:.3f} |")
    a("")
    a("## Pooled overlap (union of genes across all 4 conditions)\n")
    a("| | IMQ | HFSC | Shared | IMQ-only | HFSC-only | Jaccard |")
    a("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    a(f"| gene loci | {pooled['a']:,} | {pooled['b']:,} | **{pooled['shared']:,}** | "
      f"{pooled['ao']:,} | {pooled['bo']:,} | {pooled['jac']:.3f} |")
    a("")
    a("## Files\n")
    a("- `annotated/` — per-peak nearest gene + signed distance to TSS.")
    a("- `gene_lists/` — gene loci (representative RefSeq) hit by each set.")
    a("- `shared_genes/Shared_*` — shared loci per condition and pooled; `*_only_*` = dataset-specific.")
    a("")
    (out / "Gene_Overlap_Report.md").write_text("\n".join(L))
    print("[INFO] per-condition shared loci:", {r["v"]: r["shared"] for r in rows_rep},
          "pooled:", pooled["shared"])
    print("[INFO] wrote", out / "Gene_Overlap_Report.md")


if __name__ == "__main__":
    main()
