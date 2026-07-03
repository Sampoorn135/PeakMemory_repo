#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Role assignment for the epithelial_pancreas injury time-course (mm10, MACS2
narrowPeak, one file per condition; no replicates -> sort+merge only).

Mapping (user-confirmed):
    CTRL  -> control.bed   (baseline)
    Day1  -> sample1.bed   (perturbation 1, acute injury)
    Day7  -> sample2.bed   (perturbation 2)
    Day28 -> memory.bed    (recovery / memory endpoint)

A memory peak = reproducibly opened in both acute timepoints (Day1 & Day7),
still open at Day28, and absent from the untreated CTRL baseline.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
from memory_peak_pipeline import (
    load_bed, sort_and_merge, write_bed, count_intervals, logger, setup_logging,
)

CONDITION_TO_ROLE = {
    "CTRL":  "control",
    "Day1":  "sample1",
    "Day7":  "sample2",
    "Day28": "memory",
}
ROLE_DESC = {
    "control": "CTRL  (baseline)",
    "sample1": "Day1  (perturbation 1, acute)",
    "sample2": "Day7  (perturbation 2)",
    "memory":  "Day28 (memory endpoint)",
}
FNAME_RE = re.compile(r"GSM\d+_ATAC\.(?P<cond>CTRL|Day\d+)\.peaks\.narrowPeak$")


def prepare(raw_dir: Path, inputs_dir: Path):
    inputs_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=" * 60)
    logger.info("epithelial_pancreas role assignment")
    logger.info("=" * 60)
    found = {}
    for p in sorted(raw_dir.iterdir()):
        m = FNAME_RE.search(p.name)
        if m:
            found[m.group("cond")] = p
    logger.info("Detected conditions: %s", ", ".join(sorted(found)))
    for cond, role in CONDITION_TO_ROLE.items():
        if cond not in found:
            logger.warning("  MISSING %s (role %s)", cond, role)
            continue
        pre = sort_and_merge(load_bed(found[cond]))
        out = inputs_dir / f"{role}.bed"
        write_bed(pre, out)
        logger.info("  %-6s -> %s (%d peaks) [%s]", cond, out.name,
                    count_intervals(pre), ROLE_DESC[role])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--inputs-dir", required=True, type=Path)
    args = ap.parse_args()
    setup_logging(True)
    prepare(args.raw_dir, args.inputs_dir)


if __name__ == "__main__":
    main()
