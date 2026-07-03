# Chromatin Memory Peak Discovery — QC Report

This report is generated entirely from BED interval coordinates. 
No BAM files, read counts, signal tracks, fold-changes, or p-values were used.

## Input Summary

- Input directory: `hfsc/hfsc_analysis/inputs`
- Control file: `control.bed`
- Memory file: `memory.bed`
- Perturbation files detected: **2**
    - `sample1.bed`
    - `sample2.bed`

## Peak Counts (per input file)

Counts are after Step 1 sort + internal merge; the raw column is the count before merging.

| Role | File | Peaks (merged) | Peaks (raw) |
| --- | --- | ---: | ---: |
| control | `control.bed` | 32,236 | 32,236 |
| memory | `memory.bed` | 32,976 | 32,976 |
| sample1 | `sample1.bed` | 36,222 | 36,222 |
| sample2 | `sample2.bed` | 35,430 | 35,430 |

## Consensus Statistics

Peaks reproducibly present across **all** perturbation files.

| Consensus set | Peaks |
| --- | ---: |
| Consensus50_Single | 29,034 |
| Consensus50_Reciprocal | 28,957 |
| Consensus75_Single | 24,779 |
| Consensus75_Reciprocal | 24,243 |

## Memory Candidate Statistics

Consensus peaks that are also open in the memory state (union regions).

| Memory candidate set | Regions |
| --- | ---: |
| Memory50_Single | 20,674 |
| Memory50_Reciprocal | 19,759 |
| Memory75_Single | 14,515 |
| Memory75_Reciprocal | 12,333 |

## Control Filtering Statistics

Strict rule: a memory region overlapping **any** control peak by ≥1 bp is removed entirely.

| Variant | Candidates | Removed by control overlap | Surviving (final) |
| --- | ---: | ---: | ---: |
| 50_Single | 20,674 | 18,467 | 2,207 |
| 50_Reciprocal | 19,759 | 17,645 | 2,114 |
| 75_Single | 14,515 | 13,525 | 990 |
| 75_Reciprocal | 12,333 | 11,508 | 825 |

## Final Peak Counts

| Final set | Peaks |
| --- | ---: |
| StrictMemory50_Single | 2,207 |
| StrictMemory50_Reciprocal | 2,114 |
| StrictMemory75_Single | 990 |
| StrictMemory75_Reciprocal | 825 |

## Peak Length Statistics (final sets)

| Final set | Count | Mean (bp) | Median (bp) | Min (bp) | Max (bp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| StrictMemory50_Single | 2,207 | 1019.6 | 919.0 | 307 | 4,668 |
| StrictMemory50_Reciprocal | 2,114 | 990.5 | 905.0 | 307 | 4,668 |
| StrictMemory75_Single | 990 | 1072.0 | 986.0 | 307 | 4,668 |
| StrictMemory75_Reciprocal | 825 | 1017.5 | 947.0 | 307 | 3,351 |

## Method Notes

- **Single** criterion: overlap ≥ threshold of the *reference* peak only.
- **Reciprocal** criterion: overlap ≥ threshold of *both* peaks.
- Consensus is anchored on the first perturbation file; a reference peak 
  is retained only if every other perturbation file contributes a qualifying overlap.
- Memory regions are the union span of each consensus peak and *all* its 
  qualifying memory peaks (not the single best/first overlap).
- BED coordinates are 0-based, half-open; overlap = max(0, min(ends) − max(starts)).
