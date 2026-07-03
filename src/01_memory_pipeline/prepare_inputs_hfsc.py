#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 HFSC (GSE165312) role assignment for the memory pipeline
===============================================================================

Dataset: hair-follicle stem cell wound-memory ATAC peaks (mm10), one BED per
condition (NO replicates this time):

    GSE165312_Ctrl_D0_peaks_...bed.gz      unwounded control, day 0   (UNUSED)
    GSE165312_Ctrl_D80_peaks_...bed.gz     unwounded control, day 80  -> control
    GSE165312_Dremel_D3_peaks_...bed.gz    wound, day 3               -> sample1
    GSE165312_Dremel_D7_peaks_...bed.gz    wound, day 7               -> sample2
    GSE165312_Dremel_D80_peaks_...bed.gz   wound, day 80 (memory)     -> memory

Biological logic (identical framework to the first dataset):
  * Perturbation states = the acute wound timepoints Dremel_D3 and Dremel_D7.
    Their consensus (peaks present in BOTH, built INSIDE the main pipeline using
    the same overlap logic) defines reproducibly wound-opened elements.
  * Memory state = Dremel_D80 (still open long after the wound).
  * Control = Ctrl_D80 ONLY (user choice): the time-matched unwounded control.
    Any memory region overlapping a Ctrl_D80 peak (>=1 bp) is removed, so only
    *acquired* (wound-specific, non-baseline) memory peaks survive.

There are no replicates to collapse here, so "merge the inputs with the same
logic" reduces to the Step-1 sort + internal merge applied to each single file.
If a condition ever had multiple files, this script would collapse them with the
same build_consensus(50% Single) rule used previously (kept for generality).
===============================================================================
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

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

# Same replicate/merge rule as before, used only if a condition has >1 file.
MERGE_THRESHOLD = 0.50
MERGE_RECIPROCAL = False

# Condition -> role mapping (user-specified). Ctrl_D0 intentionally omitted.
CONDITION_TO_ROLE: Dict[str, str] = {
    "Ctrl_D80": "control",    # time-matched unwounded control (strict subtraction)
    "Dremel_D3": "sample1",   # perturbation 1 (acute wound, day 3)
    "Dremel_D7": "sample2",   # perturbation 2 (acute wound, day 7)
    "Dremel_D80": "memory",   # memory state (wound, day 80)
}
ROLE_DESCRIPTION = {
    "control": "Ctrl_D80   (unwounded day-80 control)",
    "sample1": "Dremel_D3  (perturbation 1)",
    "sample2": "Dremel_D7  (perturbation 2)",
    "memory":  "Dremel_D80 (memory state)",
}
# Conditions present but deliberately NOT used, with reason (for the log/audit).
UNUSED_CONDITIONS = {"Ctrl_D0": "day-0 baseline not used; Ctrl_D80 chosen as control"}

# Parse condition from filenames like GSE165312_Dremel_D3_peaks_SPM9_Blktclean.bed.gz
FNAME_RE = re.compile(
    r"GSE\d+_(?P<cond>Ctrl_D\d+|Dremel_D\d+)_peaks_.*\.bed(?:\.gz)?$"
)


def discover_conditions(raw_dir: Path) -> Dict[str, List[Path]]:
    groups: Dict[str, List[Path]] = {}
    for p in sorted(raw_dir.iterdir()):
        if not p.is_file():
            continue
        m = FNAME_RE.search(p.name)
        if not m:
            continue
        groups.setdefault(m.group("cond"), []).append(p)
    for cond in groups:
        groups[cond].sort(key=lambda x: x.name)
    return groups


def build_role(files: List[Path]) -> BedSet:
    """Single file -> sort+merge; multiple files -> consensus (same logic)."""
    sets = [sort_and_merge(load_bed(f)) for f in files]
    if len(sets) == 1:
        return sets[0]
    return build_consensus(sets, [f.name for f in files],
                           threshold=MERGE_THRESHOLD, reciprocal=MERGE_RECIPROCAL)


def prepare(raw_dir: Path, inputs_dir: Path) -> dict:
    inputs_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=" * 70)
    logger.info("HFSC role assignment (control = Ctrl_D80 only)")
    logger.info("=" * 70)

    groups = discover_conditions(raw_dir)
    if not groups:
        raise FileNotFoundError(f"No recognizable GSE165312 BEDs in {raw_dir}")
    logger.info("Detected %d conditions: %s", len(groups), ", ".join(sorted(groups)))

    for cond, reason in UNUSED_CONDITIONS.items():
        if cond in groups:
            logger.info("Skipping condition %-9s -> %s", cond, reason)

    report = {"mapping": {}}
    for cond, files in sorted(groups.items()):
        if cond not in CONDITION_TO_ROLE:
            continue
        role = CONDITION_TO_ROLE[cond]
        logger.info("-" * 70)
        logger.info("Condition %-10s (%d file) -> role '%s' [%s]",
                    cond, len(files), role, ROLE_DESCRIPTION[role])
        merged = build_role(files)
        out_path = inputs_dir / f"{role}.bed"
        write_bed(merged, out_path)
        logger.info("  %s -> %s (%d peaks)", cond, out_path.name, count_intervals(merged))
        report["mapping"][role] = {"condition": cond, "file": f"{role}.bed",
                                   "peaks": count_intervals(merged)}
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description="Assign HFSC GSE165312 files to pipeline roles.")
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--inputs-dir", required=True, type=Path)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    setup_logging(verbose=not args.quiet)
    prepare(args.raw_dir, args.inputs_dir)


if __name__ == "__main__":
    main()
