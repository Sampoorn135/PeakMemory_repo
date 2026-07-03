# Chromatin Memory Peak Discovery — QC Report

This report is generated entirely from BED interval coordinates. 
No BAM files, read counts, signal tracks, fold-changes, or p-values were used.

## Input Summary

- Input directory: `distil/distil_analysis/inputs`
- Control file: `control.bed`
- Memory file: `memory.bed`
- Perturbation files detected: **2**
    - `sample1.bed`
    - `sample2.bed`

## Peak Counts (per input file)

Counts are after Step 1 sort + internal merge; the raw column is the count before merging.

| Role | File | Peaks (merged) | Peaks (raw) |
| --- | --- | ---: | ---: |
| control | `control.bed` | 35,195 | 35,195 |
| memory | `memory.bed` | 34,644 | 34,644 |
| sample1 | `sample1.bed` | 49,078 | 49,078 |
| sample2 | `sample2.bed` | 42,200 | 42,200 |

## Consensus Statistics

Peaks reproducibly present across **all** perturbation files.

| Consensus set | Peaks |
| --- | ---: |
| Consensus50_Single | 36,370 |
| Consensus50_Reciprocal | 36,138 |
| Consensus75_Single | 30,899 |
| Consensus75_Reciprocal | 29,138 |

## Memory Candidate Statistics

Consensus peaks that are also open in the memory state (union regions).

| Memory candidate set | Regions |
| --- | ---: |
| Memory50_Single | 27,686 |
| Memory50_Reciprocal | 27,527 |
| Memory75_Single | 20,021 |
| Memory75_Reciprocal | 18,981 |

## Control Filtering Statistics

Strict rule: a memory region overlapping **any** control peak by ≥1 bp is removed entirely.

| Variant | Candidates | Removed by control overlap | Surviving (final) |
| --- | ---: | ---: | ---: |
| 50_Single | 27,686 | 24,590 | 3,096 |
| 50_Reciprocal | 27,527 | 24,461 | 3,066 |
| 75_Single | 20,021 | 18,004 | 2,017 |
| 75_Reciprocal | 18,981 | 17,129 | 1,852 |

## Final Peak Counts

| Final set | Peaks |
| --- | ---: |
| StrictMemory50_Single | 3,096 |
| StrictMemory50_Reciprocal | 3,066 |
| StrictMemory75_Single | 2,017 |
| StrictMemory75_Reciprocal | 1,852 |

## Peak Length Statistics (final sets)

| Final set | Count | Mean (bp) | Median (bp) | Min (bp) | Max (bp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| StrictMemory50_Single | 3,096 | 671.8 | 641.0 | 179 | 1,915 |
| StrictMemory50_Reciprocal | 3,066 | 669.8 | 638.0 | 179 | 1,915 |
| StrictMemory75_Single | 2,017 | 649.4 | 622.0 | 179 | 1,915 |
| StrictMemory75_Reciprocal | 1,852 | 642.4 | 615.5 | 179 | 1,915 |

## Method Notes

- **Single** criterion: overlap ≥ threshold of the *reference* peak only.
- **Reciprocal** criterion: overlap ≥ threshold of *both* peaks.
- Consensus is anchored on the first perturbation file; a reference peak 
  is retained only if every other perturbation file contributes a qualifying overlap.
- Memory regions are the union span of each consensus peak and *all* its 
  qualifying memory peaks (not the single best/first overlap).
- BED coordinates are 0-based, half-open; overlap = max(0, min(ends) − max(starts)).
