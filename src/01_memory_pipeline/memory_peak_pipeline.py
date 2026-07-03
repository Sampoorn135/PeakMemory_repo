#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 BED-Based Chromatin Memory Peak Discovery Pipeline
===============================================================================

PURPOSE
-------
Identify chromatin "memory" peaks using ONLY genomic interval coordinates from
BED files. The pipeline is deliberately *coordinate-only*: it never reads BAM
files, read counts, signal tracks, fold-changes, p-values, or any differential
accessibility statistic. Every decision is made from interval overlap geometry,
implemented explicitly in Python so each step is transparent and auditable.

BIOLOGICAL LOGIC (the question this answers)
--------------------------------------------
A "memory" peak is a regulatory element (an open-chromatin / ATAC peak) that:

  1. is OPENED REPRODUCIBLY across every perturbation state
     (i.e. it is a real, shared consequence of the perturbation, not noise);
  2. is STILL OPEN in the final "memory"/chase state
     (the element retained its open state after the perturbation was removed);
  3. is NOT already open in the baseline/"control" state
     (so it represents an *acquired* change, not a pre-existing open region).

Steps 2-4 below encode exactly these three biological criteria.

INPUT (generic pipeline)
------------------------
A working directory containing:
    control.bed     -> baseline / control state
    memory.bed      -> final chase / memory state
    sample1.bed ...  -> perturbation states (variable count, >=1)
Files are auto-detected by name; nothing is hard-coded.

OUTPUTS
-------
results/
  01_preprocessed/      sorted + internally-merged copy of every input BED
  02_consensus/         the 4 perturbation-consensus sets
  03_memory_candidates/ memory regions matched to each consensus set
  04_final_memory/      memory regions after strict control subtraction (FINAL)
  reports/MemoryPeak_Report.md   full QC report

The 4 variants come from two orthogonal choices:
  - threshold:    0.50 or 0.75   (how much overlap counts as "the same peak")
  - directionality: Single        (overlap must cover >= threshold of the
                                    *reference* peak only)
                    Reciprocal     (overlap must cover >= threshold of *both*
                                    peaks -- a stricter, symmetric criterion)

COORDINATE MODEL
----------------
BED is 0-based, half-open: an interval [start, end) has length (end - start)
and the bp overlap of [s1,e1) and [s2,e2) is max(0, min(e1,e2) - max(s1,s2)).
All arithmetic below uses this convention.

Author: pipeline generated for chromatin memory analysis.
===============================================================================
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import logging
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np  # noqa: F401  (kept available for downstream numeric work)


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logger = logging.getLogger("memory_peak")


def setup_logging(verbose: bool = True) -> None:
    """Configure terminal logging so every major action is reported."""
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="[%(levelname)s] %(message)s",
    )


# ---------------------------------------------------------------------------
# Core data structure
# ---------------------------------------------------------------------------
@dataclass
class Interval:
    """A single genomic interval (BED 0-based, half-open [start, end))."""
    chrom: str
    start: int
    end: int
    extra: Tuple[str, ...] = field(default_factory=tuple)  # preserved columns 4+

    @property
    def length(self) -> int:
        return self.end - self.start


# A "BedSet" is simply a chromosome-keyed dict of sorted Interval lists.
BedSet = Dict[str, List[Interval]]


# ---------------------------------------------------------------------------
# BED input / output
# ---------------------------------------------------------------------------
def _open_maybe_gzip(path: Path):
    """Open plain or gzipped text transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_bed(path: Path) -> BedSet:
    """
    Load a BED file into a chromosome-keyed dict of Intervals.

    - Skips browser/track/comment header lines.
    - Keeps columns 4+ as 'extra' so they can be preserved where meaningful.
    - Intervals are NOT sorted/merged here; that is the job of preprocess().
    """
    bedset: BedSet = {}
    n = 0
    with _open_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                # Non-numeric coordinate -> not a data row, skip.
                continue
            if end <= start:
                # Zero/negative-length interval is meaningless; skip defensively.
                continue
            extra = tuple(cols[3:]) if len(cols) > 3 else tuple()
            bedset.setdefault(chrom, []).append(Interval(chrom, start, end, extra))
            n += 1
    logger.info("Loaded %s  (%d intervals, %d chromosomes)", path.name, n, len(bedset))
    return bedset


def write_bed(bedset: BedSet, path: Path, keep_extra: bool = False) -> int:
    """
    Write a BedSet to disk, sorted by (chrom, start, end).

    keep_extra=False writes a clean 3-column BED. We default to 3 columns for
    any interval whose coordinates were *constructed* by the pipeline (consensus
    regions, union memory regions, merged intervals), because per-peak metadata
    from the source rows no longer applies to the new coordinates.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(path, "wt") as fh:
        for chrom in _sorted_chroms(bedset.keys()):
            for iv in sorted(bedset[chrom], key=lambda x: (x.start, x.end)):
                if keep_extra and iv.extra:
                    fh.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\t" + "\t".join(iv.extra) + "\n")
                else:
                    fh.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\n")
                count += 1
    logger.info("Wrote %s  (%d intervals)", path.name, count)
    return count


def _sorted_chroms(chroms) -> List[str]:
    """Natural-ish chromosome ordering: chr1, chr2, ... chr10, ... chrX, chrY, chrM."""
    def keyfn(c: str):
        body = c[3:] if c.lower().startswith("chr") else c
        if body.isdigit():
            return (0, int(body), "")
        return (1, 0, body)
    return sorted(chroms, key=keyfn)


def count_intervals(bedset: BedSet) -> int:
    return sum(len(v) for v in bedset.values())


# ---------------------------------------------------------------------------
# Overlap primitives  (explicit, transparent)
# ---------------------------------------------------------------------------
def overlap_bp(a: Interval, b: Interval) -> int:
    """Base pairs shared by two intervals on the SAME chromosome (0 if none)."""
    if a.chrom != b.chrom:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def passes_single(ref: Interval, other: Interval, threshold: float) -> bool:
    """
    'Single' criterion: overlap covers >= threshold of the REFERENCE peak.
        overlap / len(ref) >= threshold
    Directional: only the reference peak must be sufficiently covered.
    """
    ov = overlap_bp(ref, other)
    if ov == 0:
        return False
    return ov / ref.length >= threshold


def passes_reciprocal(ref: Interval, other: Interval, threshold: float) -> bool:
    """
    'Reciprocal' criterion: overlap covers >= threshold of BOTH peaks.
        overlap / len(ref)   >= threshold  AND
        overlap / len(other) >= threshold
    Symmetric and stricter than 'single'.
    """
    ov = overlap_bp(ref, other)
    if ov == 0:
        return False
    return (ov / ref.length >= threshold) and (ov / other.length >= threshold)


def interval_passes(ref: Interval, other: Interval, threshold: float, reciprocal: bool) -> bool:
    """Dispatch to the chosen overlap criterion."""
    return (passes_reciprocal if reciprocal else passes_single)(ref, other, threshold)


# ---------------------------------------------------------------------------
# Fast overlap candidate lookup
# ---------------------------------------------------------------------------
class IntervalIndex:
    """
    Per-chromosome sorted index for fast overlap queries.

    For each chromosome we store intervals sorted by start, plus a running
    'max end so far' is NOT used; instead we bound the search using the maximum
    interval length on that chromosome so a simple bisect window is exhaustive.
    This keeps the logic explicit and dependency-free (no intervaltree needed).
    """

    def __init__(self, bedset: BedSet):
        self._by_chrom: Dict[str, List[Interval]] = {}
        self._starts: Dict[str, List[int]] = {}
        self._max_len: Dict[str, int] = {}
        for chrom, ivs in bedset.items():
            s = sorted(ivs, key=lambda x: x.start)
            self._by_chrom[chrom] = s
            self._starts[chrom] = [iv.start for iv in s]
            self._max_len[chrom] = max((iv.length for iv in s), default=0)

    def query(self, q: Interval) -> List[Interval]:
        """Return all indexed intervals that overlap q by >= 1 bp."""
        chrom = q.chrom
        if chrom not in self._by_chrom:
            return []
        starts = self._starts[chrom]
        ivs = self._by_chrom[chrom]
        max_len = self._max_len[chrom]
        # Any overlapping interval must start before q.end and end after q.start.
        # Its start is therefore in [q.start - max_len, q.end).
        lo = bisect.bisect_left(starts, q.start - max_len)
        hi = bisect.bisect_left(starts, q.end)
        hits: List[Interval] = []
        for i in range(lo, hi):
            iv = ivs[i]
            if iv.end > q.start:  # confirm true overlap
                hits.append(iv)
        return hits


# ---------------------------------------------------------------------------
# Step 1: Preprocessing  (sort + merge internally overlapping intervals)
# ---------------------------------------------------------------------------
def sort_and_merge(bedset: BedSet) -> BedSet:
    """
    Sort intervals and merge those that overlap or touch WITHIN a single file.

    Biological purpose: a single regulatory element can be represented by
    several overlapping/adjacent peak calls. Collapsing them to one interval
    removes duplicate representations and makes all downstream overlap logic
    deterministic. Merged coordinates span the union of the merged peaks.

    Note: merged intervals are emitted as bare coordinates (extra columns from
    the originals no longer describe the new merged span).
    """
    out: BedSet = {}
    for chrom, ivs in bedset.items():
        s = sorted(ivs, key=lambda x: (x.start, x.end))
        merged: List[Interval] = []
        for iv in s:
            if merged and iv.start <= merged[-1].end:  # overlap or touch
                last = merged[-1]
                if iv.end > last.end:
                    merged[-1] = Interval(chrom, last.start, iv.end)
            else:
                merged.append(Interval(chrom, iv.start, iv.end))
        out[chrom] = merged
    return out


# ---------------------------------------------------------------------------
# Step 2: Perturbation consensus
# ---------------------------------------------------------------------------
def build_consensus(
    perturbations: Sequence[BedSet],
    pert_names: Sequence[str],
    threshold: float,
    reciprocal: bool,
) -> BedSet:
    """
    Build a consensus set of peaks PRESENT IN ALL perturbation files.

    Algorithm (deterministic, anchored on a reference file):
      * The FIRST perturbation file is the reference.
      * For each reference peak R:
          - For every OTHER perturbation file F, search for peaks overlapping R
            that satisfy the chosen criterion relative to R.
          - If even one file F contributes NO qualifying peak, R is dropped:
            the element is not reproducibly present across all perturbations.
          - Otherwise the cluster = R plus every qualifying peak from every
            other file.
      * Consensus coordinates for a retained cluster:
            start = min(start) over all matched peaks
            end   = max(end)   over all matched peaks

    With exactly two perturbation files (reference + one other) the "anchor on
    reference" rule and a pairwise "consecutive matched peaks" rule are
    identical, so this faithfully implements the spec while generalising to any
    number of perturbation files.

    Biological meaning: retains only chromatin elements that respond the same
    way in every perturbation replicate/condition -> a robust, shared signal.
    """
    if not perturbations:
        return {}

    ref_set = perturbations[0]
    other_indices = [IntervalIndex(b) for b in perturbations[1:]]
    other_names = list(pert_names[1:])

    consensus: BedSet = {}
    n_ref = count_intervals(ref_set)
    retained = 0

    for chrom, ref_ivs in ref_set.items():
        for R in ref_ivs:
            matched: List[Interval] = [R]
            present_in_all = True
            for idx in other_indices:
                # candidate peaks in this other file overlapping R
                candidates = idx.query(R)
                qualifying = [c for c in candidates if interval_passes(R, c, threshold, reciprocal)]
                if not qualifying:
                    present_in_all = False
                    break
                matched.extend(qualifying)
            if not present_in_all:
                continue
            start = min(iv.start for iv in matched)
            end = max(iv.end for iv in matched)
            consensus.setdefault(chrom, []).append(Interval(chrom, start, end))
            retained += 1

    # Collapse any consensus intervals that now overlap each other (different
    # reference peaks can expand into overlapping spans). This keeps the
    # consensus set a clean, non-redundant interval set.
    consensus = sort_and_merge(consensus)

    crit = "Reciprocal" if reciprocal else "Single"
    logger.info(
        "  consensus(threshold=%.2f, %s): %d reference peaks -> %d retained -> %d after merge "
        "(other files: %s)",
        threshold, crit, n_ref, retained, count_intervals(consensus),
        ", ".join(other_names) if other_names else "(none)",
    )
    return consensus


# ---------------------------------------------------------------------------
# Step 3: Compare consensus peaks against memory
# ---------------------------------------------------------------------------
def match_memory(
    consensus: BedSet,
    memory: BedSet,
    threshold: float,
    reciprocal: bool,
) -> BedSet:
    """
    For each consensus peak, find ALL memory peaks satisfying the overlap
    criterion (not just the best/first), then build ONE union region:

        start = min(consensus.start, all matching memory starts)
        end   = max(consensus.end,   all matching memory ends)

    A consensus peak with no qualifying memory peak yields no output.

    Biological meaning: keep the perturbation-shared element only if it is also
    open in the final memory state; the union span records the full extent of
    the retained-open region.
    """
    mem_index = IntervalIndex(memory)
    out: BedSet = {}
    n_consensus = count_intervals(consensus)
    n_with_match = 0

    for chrom, cons_ivs in consensus.items():
        for C in cons_ivs:
            candidates = mem_index.query(C)
            # For 'single', the criterion is overlap >= threshold of the
            # CONSENSUS peak (the reference here). For 'reciprocal', both the
            # consensus peak and the memory peak must be >= threshold covered.
            matches = [m for m in candidates if interval_passes(C, m, threshold, reciprocal)]
            if not matches:
                continue
            start = min([C.start] + [m.start for m in matches])
            end = max([C.end] + [m.end for m in matches])
            out.setdefault(chrom, []).append(Interval(chrom, start, end))
            n_with_match += 1

    out = sort_and_merge(out)  # union regions may overlap; collapse to clean set
    crit = "Reciprocal" if reciprocal else "Single"
    logger.info(
        "  memory match(threshold=%.2f, %s): %d consensus peaks -> %d matched -> %d regions",
        threshold, crit, n_consensus, n_with_match, count_intervals(out),
    )
    return out


# ---------------------------------------------------------------------------
# Step 4: Strict control subtraction
# ---------------------------------------------------------------------------
def filter_control(memory_regions: BedSet, control: BedSet) -> Tuple[BedSet, int]:
    """
    STRICT rule: if a memory region overlaps ANY control peak by >= 1 bp,
    discard the entire region. No thresholds, no partial subtraction.

    Biological meaning: any overlap with the baseline state means the element
    was already (partly) open before the perturbation, so it cannot be an
    *acquired* memory peak.

    Returns (surviving_regions, n_removed).
    """
    control_index = IntervalIndex(control)
    out: BedSet = {}
    removed = 0
    for chrom, ivs in memory_regions.items():
        for iv in ivs:
            hits = control_index.query(iv)  # query() already requires >=1 bp overlap
            if hits:
                removed += 1
                continue
            out.setdefault(chrom, []).append(iv)
    logger.info(
        "  control filter: removed %d / %d regions overlapping control -> %d final",
        removed, count_intervals(memory_regions), count_intervals(out),
    )
    return out, removed


# ---------------------------------------------------------------------------
# Length statistics
# ---------------------------------------------------------------------------
def length_stats(bedset: BedSet) -> Dict[str, float]:
    """count, mean, median, min, max of interval lengths."""
    lengths = [iv.length for ivs in bedset.values() for iv in ivs]
    if not lengths:
        return {"count": 0, "mean": 0.0, "median": 0.0, "min": 0, "max": 0}
    return {
        "count": len(lengths),
        "mean": statistics.mean(lengths),
        "median": statistics.median(lengths),
        "min": min(lengths),
        "max": max(lengths),
    }


# ---------------------------------------------------------------------------
# The 4 variants
# ---------------------------------------------------------------------------
VARIANTS: List[Tuple[str, float, bool]] = [
    ("50_Single", 0.50, False),
    ("50_Reciprocal", 0.50, True),
    ("75_Single", 0.75, False),
    ("75_Reciprocal", 0.75, True),
]


# ---------------------------------------------------------------------------
# File auto-detection (generic pipeline; nothing hard-coded)
# ---------------------------------------------------------------------------
def detect_inputs(input_dir: Path) -> Tuple[Path, Path, List[Path]]:
    """
    Auto-detect control.bed, memory.bed, and all perturbation sample BEDs.

    - control.bed and memory.bed are matched by exact stem (case-insensitive).
    - Every other *.bed that is not control/memory is treated as a perturbation
      sample. Files matching 'sample*' sort first and in natural numeric order.
    """
    beds = sorted(input_dir.glob("*.bed"))
    control = memory = None
    samples: List[Path] = []
    for p in beds:
        stem = p.stem.lower()
        if stem == "control":
            control = p
        elif stem == "memory":
            memory = p
        else:
            samples.append(p)

    def sample_key(p: Path):
        m = re.search(r"(\d+)", p.stem)
        return (0, int(m.group(1))) if m else (1, p.stem)

    samples.sort(key=sample_key)

    if control is None:
        raise FileNotFoundError(f"control.bed not found in {input_dir}")
    if memory is None:
        raise FileNotFoundError(f"memory.bed not found in {input_dir}")
    if not samples:
        raise FileNotFoundError(f"no perturbation sample*.bed found in {input_dir}")
    return control, memory, samples


# ---------------------------------------------------------------------------
# Main generic pipeline
# ---------------------------------------------------------------------------
def run_pipeline(input_dir: Path, results_dir: Path) -> dict:
    """Execute the full generic pipeline on an input directory."""
    logger.info("=" * 70)
    logger.info("BED-BASED CHROMATIN MEMORY PEAK PIPELINE")
    logger.info("=" * 70)

    control_path, memory_path, sample_paths = detect_inputs(input_dir)
    logger.info("Detected control : %s", control_path.name)
    logger.info("Detected memory  : %s", memory_path.name)
    logger.info("Detected %d perturbation files: %s",
                len(sample_paths), ", ".join(p.name for p in sample_paths))

    d_pre = results_dir / "01_preprocessed"
    d_cons = results_dir / "02_consensus"
    d_mem = results_dir / "03_memory_candidates"
    d_final = results_dir / "04_final_memory"
    d_rep = results_dir / "reports"
    for d in (d_pre, d_cons, d_mem, d_final, d_rep):
        d.mkdir(parents=True, exist_ok=True)

    # ---- Load + Step 1 preprocess every input -----------------------------
    logger.info("-" * 70)
    logger.info("STEP 1: Preprocessing (sort + internal merge)")
    raw_inputs: Dict[str, BedSet] = {}
    pre_inputs: Dict[str, BedSet] = {}

    def _ingest(label: str, path: Path) -> BedSet:
        raw = load_bed(path)
        pre = sort_and_merge(raw)
        write_bed(pre, d_pre / f"{label}.bed")
        raw_inputs[label] = raw
        pre_inputs[label] = pre
        logger.info("  %-10s %d raw -> %d after sort+merge", label,
                    count_intervals(raw), count_intervals(pre))
        return pre

    control = _ingest("control", control_path)
    memory = _ingest("memory", memory_path)
    pert_sets: List[BedSet] = []
    pert_labels: List[str] = []
    for i, sp in enumerate(sample_paths, start=1):
        label = f"sample{i}"
        pert_sets.append(_ingest(label, sp))
        pert_labels.append(label)

    # ---- Step 2 consensus, Step 3 memory, Step 4 control (per variant) -----
    summary = {
        "input_dir": str(input_dir),
        "results_dir": str(results_dir),
        "control_name": control_path.name,
        "memory_name": memory_path.name,
        "sample_names": [p.name for p in sample_paths],
        "n_perturbations": len(sample_paths),
        "peak_counts": {label: count_intervals(pre_inputs[label]) for label in pre_inputs},
        "raw_peak_counts": {label: count_intervals(raw_inputs[label]) for label in raw_inputs},
        "consensus": {},
        "memory": {},
        "removed": {},
        "final": {},
        "final_stats": {},
    }

    for name, thr, recip in VARIANTS:
        logger.info("-" * 70)
        logger.info("VARIANT %s  (threshold=%.2f, %s)", name, thr,
                    "Reciprocal" if recip else "Single")

        logger.info("STEP 2: Building Consensus%s", name)
        cons = build_consensus(pert_sets, pert_labels, thr, recip)
        write_bed(cons, d_cons / f"Consensus{name}.bed")
        summary["consensus"][name] = count_intervals(cons)

        logger.info("STEP 3: Matching Consensus%s against memory", name)
        mem = match_memory(cons, memory, thr, recip)
        write_bed(mem, d_mem / f"Memory{name}.bed")
        summary["memory"][name] = count_intervals(mem)

        logger.info("STEP 4: Strict control subtraction for %s", name)
        final, removed = filter_control(mem, control)
        write_bed(final, d_final / f"StrictMemory{name}.bed")
        summary["removed"][name] = removed
        summary["final"][name] = count_intervals(final)
        summary["final_stats"][name] = length_stats(final)

    # ---- QC report --------------------------------------------------------
    logger.info("-" * 70)
    logger.info("Generating QC report")
    write_report(summary, d_rep / "MemoryPeak_Report.md")
    logger.info("DONE. Results in %s", results_dir)
    return summary


# ---------------------------------------------------------------------------
# QC report
# ---------------------------------------------------------------------------
def write_report(summary: dict, path: Path) -> None:
    """Render the full markdown QC report required by the spec."""
    L: List[str] = []
    a = L.append

    a("# Chromatin Memory Peak Discovery — QC Report\n")
    a("This report is generated entirely from BED interval coordinates. ")
    a("No BAM files, read counts, signal tracks, fold-changes, or p-values were used.\n")

    # Input summary
    a("## Input Summary\n")
    a(f"- Input directory: `{summary['input_dir']}`")
    a(f"- Control file: `{summary['control_name']}`")
    a(f"- Memory file: `{summary['memory_name']}`")
    a(f"- Perturbation files detected: **{summary['n_perturbations']}**")
    for s in summary["sample_names"]:
        a(f"    - `{s}`")
    a("")

    # Peak counts per input
    a("## Peak Counts (per input file)\n")
    a("Counts are after Step 1 sort + internal merge; the raw column is the count before merging.\n")
    a("| Role | File | Peaks (merged) | Peaks (raw) |")
    a("| --- | --- | ---: | ---: |")
    role_to_name = {"control": summary["control_name"], "memory": summary["memory_name"]}
    for i, s in enumerate(summary["sample_names"], start=1):
        role_to_name[f"sample{i}"] = s
    for label in summary["peak_counts"]:
        a(f"| {label} | `{role_to_name.get(label, '')}` | "
          f"{summary['peak_counts'][label]:,} | {summary['raw_peak_counts'][label]:,} |")
    a("")

    # Consensus statistics
    a("## Consensus Statistics\n")
    a("Peaks reproducibly present across **all** perturbation files.\n")
    a("| Consensus set | Peaks |")
    a("| --- | ---: |")
    for name, _, _ in VARIANTS:
        a(f"| Consensus{name} | {summary['consensus'][name]:,} |")
    a("")

    # Memory candidate statistics
    a("## Memory Candidate Statistics\n")
    a("Consensus peaks that are also open in the memory state (union regions).\n")
    a("| Memory candidate set | Regions |")
    a("| --- | ---: |")
    for name, _, _ in VARIANTS:
        a(f"| Memory{name} | {summary['memory'][name]:,} |")
    a("")

    # Control filtering statistics
    a("## Control Filtering Statistics\n")
    a("Strict rule: a memory region overlapping **any** control peak by ≥1 bp is removed entirely.\n")
    a("| Variant | Candidates | Removed by control overlap | Surviving (final) |")
    a("| --- | ---: | ---: | ---: |")
    for name, _, _ in VARIANTS:
        cand = summary["memory"][name]
        rem = summary["removed"][name]
        fin = summary["final"][name]
        a(f"| {name} | {cand:,} | {rem:,} | {fin:,} |")
    a("")

    # Final peak counts
    a("## Final Peak Counts\n")
    a("| Final set | Peaks |")
    a("| --- | ---: |")
    for name, _, _ in VARIANTS:
        a(f"| StrictMemory{name} | {summary['final'][name]:,} |")
    a("")

    # Length statistics
    a("## Peak Length Statistics (final sets)\n")
    a("| Final set | Count | Mean (bp) | Median (bp) | Min (bp) | Max (bp) |")
    a("| --- | ---: | ---: | ---: | ---: | ---: |")
    for name, _, _ in VARIANTS:
        st = summary["final_stats"][name]
        a(f"| StrictMemory{name} | {st['count']:,} | {st['mean']:.1f} | "
          f"{st['median']:.1f} | {st['min']:,} | {st['max']:,} |")
    a("")

    a("## Method Notes\n")
    a("- **Single** criterion: overlap ≥ threshold of the *reference* peak only.")
    a("- **Reciprocal** criterion: overlap ≥ threshold of *both* peaks.")
    a("- Consensus is anchored on the first perturbation file; a reference peak ")
    a("  is retained only if every other perturbation file contributes a qualifying overlap.")
    a("- Memory regions are the union span of each consensus peak and *all* its ")
    a("  qualifying memory peaks (not the single best/first overlap).")
    a("- BED coordinates are 0-based, half-open; overlap = max(0, min(ends) − max(starts)).")
    a("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L))
    logger.info("Wrote QC report -> %s", path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="BED-based chromatin memory peak discovery (coordinate-only)."
    )
    ap.add_argument("--input-dir", required=True, type=Path,
                    help="Directory containing control.bed, memory.bed, sample*.bed")
    ap.add_argument("--results-dir", required=True, type=Path,
                    help="Directory to write results/ subfolders")
    ap.add_argument("--quiet", action="store_true", help="Reduce logging")
    args = ap.parse_args()

    setup_logging(verbose=not args.quiet)
    run_pipeline(args.input_dir, args.results_dir)


if __name__ == "__main__":
    main()
