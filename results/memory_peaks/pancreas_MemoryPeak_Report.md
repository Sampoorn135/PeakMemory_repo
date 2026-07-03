# Chromatin Memory Peak Discovery — QC Report

This report is generated entirely from BED interval coordinates. 
No BAM files, read counts, signal tracks, fold-changes, or p-values were used.

## Input Summary

- Input directory: `epithelial_pancreas/pancreas_analysis/inputs`
- Control file: `control.bed`
- Memory file: `memory.bed`
- Perturbation files detected: **2**
    - `sample1.bed`
    - `sample2.bed`

## Peak Counts (per input file)

Counts are after Step 1 sort + internal merge; the raw column is the count before merging.

| Role | File | Peaks (merged) | Peaks (raw) |
| --- | --- | ---: | ---: |
| control | `control.bed` | 140,295 | 140,295 |
| memory | `memory.bed` | 130,155 | 130,155 |
| sample1 | `sample1.bed` | 164,695 | 164,695 |
| sample2 | `sample2.bed` | 145,344 | 145,344 |

## Consensus Statistics

Peaks reproducibly present across **all** perturbation files.

| Consensus set | Peaks |
| --- | ---: |
| Consensus50_Single | 79,991 |
| Consensus50_Reciprocal | 75,264 |
| Consensus75_Single | 55,514 |
| Consensus75_Reciprocal | 42,429 |

## Memory Candidate Statistics

Consensus peaks that are also open in the memory state (union regions).

| Memory candidate set | Regions |
| --- | ---: |
| Memory50_Single | 46,787 |
| Memory50_Reciprocal | 42,366 |
| Memory75_Single | 28,090 |
| Memory75_Reciprocal | 18,670 |

## Control Filtering Statistics

Strict rule: a memory region overlapping **any** control peak by ≥1 bp is removed entirely.

| Variant | Candidates | Removed by control overlap | Surviving (final) |
| --- | ---: | ---: | ---: |
| 50_Single | 46,787 | 43,526 | 3,261 |
| 50_Reciprocal | 42,366 | 39,283 | 3,083 |
| 75_Single | 28,090 | 27,237 | 853 |
| 75_Reciprocal | 18,670 | 18,124 | 546 |

## Final Peak Counts

| Final set | Peaks |
| --- | ---: |
| StrictMemory50_Single | 3,261 |
| StrictMemory50_Reciprocal | 3,083 |
| StrictMemory75_Single | 853 |
| StrictMemory75_Reciprocal | 546 |

## Peak Length Statistics (final sets)

| Final set | Count | Mean (bp) | Median (bp) | Min (bp) | Max (bp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| StrictMemory50_Single | 3,261 | 400.1 | 362.0 | 151 | 2,441 |
| StrictMemory50_Reciprocal | 3,083 | 395.0 | 357.0 | 151 | 2,441 |
| StrictMemory75_Single | 853 | 364.2 | 335.0 | 151 | 1,352 |
| StrictMemory75_Reciprocal | 546 | 349.2 | 324.0 | 151 | 1,352 |

## Method Notes

- **Single** criterion: overlap ≥ threshold of the *reference* peak only.
- **Reciprocal** criterion: overlap ≥ threshold of *both* peaks.
- Consensus is anchored on the first perturbation file; a reference peak 
  is retained only if every other perturbation file contributes a qualifying overlap.
- Memory regions are the union span of each consensus peak and *all* its 
  qualifying memory peaks (not the single best/first overlap).
- BED coordinates are 0-based, half-open; overlap = max(0, min(ends) − max(starts)).
