#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 GSE221408 (distal wound-memory) role assignment for the memory pipeline
===============================================================================

Dataset: epidermal ATAC-seq, IDR conservative peaks (mm10), from
"Tissue memory relies on stem cell priming in distal undamaged areas"
(Nat Cell Biol 2023). Timeline (from the paper):
  0w     = homeostasis / baseline, before the first wound
  1w pw1 = 1 week post FIRST wound  (first wound healing)
  8w pw1 = 8 weeks post first wound = reset homeostasis = MEMORY / chase state
  1w pw2 = 1 week post SECOND wound (second wound healing / re-challenge)

User-specified role mapping (use only the required set):
  0w     -> control.bed   (baseline; strict subtraction)
  8wpw1  -> memory.bed    (memory / chase endpoint)
  1wpw1  -> sample1.bed   (perturbation 1 = first wound healing)
  1wpw2  -> sample2.bed   (perturbation 2 = second wound healing)
  0w_AC, 8wpw1_MEM        -> NOT USED

IDR conservative peaks are already replicate-reproducible, so there is no
replicate merge — each role file is just Step-1 sort + internal merge.
Then the standard pipeline builds the perturbation consensus (peaks shared by
BOTH wound healings), keeps those still open in memory (8wpw1), and strictly
removes anything overlapping baseline (0w).
===============================================================================
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

from memory_peak_pipeline import (
    BedSet, build_consensus, count_intervals, load_bed, logger,
    setup_logging, sort_and_merge, write_bed,
)

MERGE_THRESHOLD = 0.50          # only used if a condition ever has >1 file
MERGE_RECIPROCAL = False

CONDITION_TO_ROLE: Dict[str, str] = {
    "0w": "control",     # baseline
    "8wpw1": "memory",   # memory / chase endpoint
    "1wpw1": "sample1",  # perturbation 1 (first wound healing)
    "1wpw2": "sample2",  # perturbation 2 (second wound healing)
}
ROLE_DESCRIPTION = {
    "control": "0w     (baseline)",
    "memory":  "8wpw1  (memory / chase endpoint)",
    "sample1": "1wpw1  (first wound healing)",
    "sample2": "1wpw2  (second wound healing)",
}
UNUSED_CONDITIONS = {
    "0w_AC": "week-0 away-control region, not used",
    "8wpw1_MEM": "8w distal-memory-region file, not used (8wpw1 chosen as memory)",
}

# GSE221408_ATAC_<cond>-idr.conservative_peak.narrowPeak
FNAME_RE = re.compile(r"GSE\d+_ATAC_(?P<cond>.+?)-idr\.conservative_peak\.narrowPeak$")


def discover(raw_dir: Path) -> Dict[str, List[Path]]:
    groups: Dict[str, List[Path]] = {}
    for p in sorted(raw_dir.iterdir()):
        if not p.is_file():
            continue
        m = FNAME_RE.search(p.name)
        if m:
            groups.setdefault(m.group("cond"), []).append(p)
    return groups


def build_role(files: List[Path]) -> BedSet:
    sets = [sort_and_merge(load_bed(f)) for f in files]
    if len(sets) == 1:
        return sets[0]
    return build_consensus(sets, [f.name for f in files],
                           threshold=MERGE_THRESHOLD, reciprocal=MERGE_RECIPROCAL)


def prepare(raw_dir: Path, inputs_dir: Path) -> dict:
    inputs_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=" * 70)
    logger.info("GSE221408 distal wound-memory: role assignment")
    logger.info("=" * 70)
    groups = discover(raw_dir)
    if not groups:
        raise FileNotFoundError(f"No GSE221408 narrowPeak files in {raw_dir}")
    logger.info("Detected conditions: %s", ", ".join(sorted(groups)))

    for cond, reason in UNUSED_CONDITIONS.items():
        if cond in groups:
            logger.info("Skipping %-10s -> %s", cond, reason)

    missing = [c for c in CONDITION_TO_ROLE if c not in groups]
    if missing:
        logger.warning("MISSING required conditions: %s", ", ".join(missing))

    report = {"mapping": {}, "missing": missing}
    for cond, files in sorted(groups.items()):
        if cond not in CONDITION_TO_ROLE:
            continue
        role = CONDITION_TO_ROLE[cond]
        logger.info("-" * 70)
        logger.info("Condition %-10s -> role '%s' [%s]", cond, role, ROLE_DESCRIPTION[role])
        merged = build_role(files)
        out_path = inputs_dir / f"{role}.bed"
        write_bed(merged, out_path)
        logger.info("  %s -> %s (%d peaks)", cond, out_path.name, count_intervals(merged))
        report["mapping"][role] = {"condition": cond, "peaks": count_intervals(merged)}
    return report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--inputs-dir", required=True, type=Path)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    setup_logging(verbose=not args.quiet)
    prepare(args.raw_dir, args.inputs_dir)


if __name__ == "__main__":
    main()
