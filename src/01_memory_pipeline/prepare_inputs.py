#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 Replicate -> Condition consensus, and role assignment for the memory pipeline
===============================================================================

This dataset ships 8 raw ATAC-seq peak BED files = 4 biological conditions,
each with 2 replicates:

    GSM5229061_D6_Ctrl_..._Rep1_peaks.bed.gz   \
    GSM5229062_D6_Ctrl_..._Rep2_peaks.bed.gz   /  -> D6_Ctrl
    GSM5229063_D6_IMQ_..._Rep1_peaks.bed.gz    \
    GSM5229064_D6_IMQ_..._Rep2_peaks.bed.gz    /  -> D6_IMQ
    GSM5229065_D30_Ctrl_..._Rep1_peaks.bed.gz  \
    GSM5229066_D30_Ctrl_..._Rep2_peaks.bed.gz  /  -> D30_Ctrl
    GSM5229067_D30_PI_..._Rep1_peaks.bed.gz    \
    GSM5229068_D30_PI_..._Rep2_peaks.bed.gz    /  -> D30_PI

STEP A — collapse replicates to one BED per condition
-----------------------------------------------------
We collapse each condition's TWO replicates into ONE reproducible peak set using
*exactly the same consensus logic* applied to perturbations (build_consensus),
treating the two replicates as the "perturbation files". The chosen rule is
**50% Single** (overlap >= 50% of the reference/Rep1 peak). A peak is kept only
if it reproduces across both replicates; the kept coordinate is the union span
of the matched replicate peaks. This is the standard "reproducible peaks across
replicates" idea, expressed in the same coordinate-only language as the rest of
the pipeline.

STEP B — assign biological roles (user-specified mapping)
---------------------------------------------------------
    D6_Ctrl  (baseline arm, perturbation 1) -> sample1.bed
    D6_IMQ   (inflammation,  perturbation 2) -> sample2.bed
    D30_Ctrl (time control)                  -> control.bed
    D30_PI   (memory state)                  -> memory.bed

The collapsed role BEDs are written to <inputs_dir>/ where the generic pipeline
auto-detects them. Nothing about the threshold or mapping is hidden: both are
declared here in plain sight.
===============================================================================
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

from memory_peak_pipeline import (
    BedSet,
    build_consensus,
    count_intervals,
    load_bed,
    logger,
    setup_logging,
    sort_and_merge,
    write_bed,
)

# Replicate-merge rule (chosen by the user): 50% Single
REP_MERGE_THRESHOLD = 0.50
REP_MERGE_RECIPROCAL = False

# Condition -> role mapping (user-specified biological design)
CONDITION_TO_ROLE: Dict[str, str] = {
    "D6_Ctrl": "sample1",    # perturbation 1 (baseline arm)
    "D6_IMQ": "sample2",     # perturbation 2 (inflammation)
    "D30_Ctrl": "control",   # baseline / time control
    "D30_PI": "memory",      # final memory state
}

ROLE_DESCRIPTION = {
    "sample1": "D6_Ctrl  (perturbation 1)",
    "sample2": "D6_IMQ   (perturbation 2)",
    "control": "D30_Ctrl (baseline / control)",
    "memory": "D30_PI   (memory state)",
}

# Parse condition + replicate from filenames like
# GSM5229061_D6_Ctrl_ATACseq_WT_Rep1_peaks.bed.gz
FNAME_RE = re.compile(
    r"GSM\d+_(?P<cond>D\d+_[A-Za-z0-9]+)_ATACseq_WT_(?P<rep>Rep\d+)_peaks\.bed(?:\.gz)?$"
)


def discover_conditions(raw_dir: Path) -> Dict[str, List[Path]]:
    """Group raw replicate BEDs by biological condition."""
    groups: Dict[str, List[Path]] = {}
    for p in sorted(raw_dir.iterdir()):
        m = FNAME_RE.search(p.name)
        if not m:
            continue
        cond = m.group("cond")
        groups.setdefault(cond, []).append(p)
    for cond in groups:
        groups[cond].sort(key=lambda x: x.name)  # Rep1 before Rep2
    return groups


def collapse_replicates(rep_paths: List[Path], cond: str) -> Tuple[BedSet, dict]:
    """
    Collapse a condition's replicates into one reproducible BedSet using the
    same build_consensus logic used for perturbations (50% Single).
    """
    rep_sets: List[BedSet] = []
    rep_labels: List[str] = []
    info = {"condition": cond, "replicates": []}
    for rp in rep_paths:
        raw = load_bed(rp)
        pre = sort_and_merge(raw)  # sort + internal merge each replicate first
        rep_sets.append(pre)
        rep_labels.append(rp.name)
        info["replicates"].append({"file": rp.name,
                                    "raw": count_intervals(raw),
                                    "merged": count_intervals(pre)})

    if len(rep_sets) == 1:
        logger.warning("  %s has a single replicate; using it as-is.", cond)
        merged = rep_sets[0]
    else:
        merged = build_consensus(
            rep_sets, rep_labels,
            threshold=REP_MERGE_THRESHOLD, reciprocal=REP_MERGE_RECIPROCAL,
        )
    info["reproducible_peaks"] = count_intervals(merged)
    return merged, info


def prepare(raw_dir: Path, inputs_dir: Path) -> dict:
    """Build the 4 role BEDs (control/memory/sample1/sample2) from 8 raw files."""
    inputs_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=" * 70)
    logger.info("STEP A/B: Replicate consensus (50%% Single) + role assignment")
    logger.info("=" * 70)

    groups = discover_conditions(raw_dir)
    if not groups:
        raise FileNotFoundError(f"No recognizable GSM replicate BEDs in {raw_dir}")

    logger.info("Detected %d conditions: %s", len(groups), ", ".join(sorted(groups)))
    report = {"conditions": [], "mapping": {}}

    for cond, reps in sorted(groups.items()):
        if cond not in CONDITION_TO_ROLE:
            logger.warning("Condition %s has no role mapping; skipping.", cond)
            continue
        role = CONDITION_TO_ROLE[cond]
        logger.info("-" * 70)
        logger.info("Condition %-9s (%d replicates) -> role '%s' [%s]",
                    cond, len(reps), role, ROLE_DESCRIPTION[role])
        merged, info = collapse_replicates(reps, cond)
        info["role"] = role
        out_path = inputs_dir / f"{role}.bed"
        write_bed(merged, out_path)
        logger.info("  %s -> %s (%d reproducible peaks)",
                    cond, out_path.name, count_intervals(merged))
        report["conditions"].append(info)
        report["mapping"][role] = {"condition": cond,
                                   "file": f"{role}.bed",
                                   "reproducible_peaks": count_intervals(merged)}

    return report


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Collapse ATAC replicates to condition consensus and assign roles."
    )
    ap.add_argument("--raw-dir", required=True, type=Path,
                    help="Directory with the 8 raw GSM replicate BED(.gz) files")
    ap.add_argument("--inputs-dir", required=True, type=Path,
                    help="Output directory for control/memory/sample*.bed")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    setup_logging(verbose=not args.quiet)
    prepare(args.raw_dir, args.inputs_dir)


if __name__ == "__main__":
    main()
