# Chromatin Memory Peak Discovery — QC Report

This report is generated entirely from BED interval coordinates. 
No BAM files, read counts, signal tracks, fold-changes, or p-values were used.

## Input Summary

- Input directory: `memory_analysis/inputs`
- Control file: `control.bed`
- Memory file: `memory.bed`
- Perturbation files detected: **2**
    - `sample1.bed`
    - `sample2.bed`

## Peak Counts (per input file)

Counts are after Step 1 sort + internal merge; the raw column is the count before merging.

| Role | File | Peaks (merged) | Peaks (raw) |
| --- | --- | ---: | ---: |
| control | `control.bed` | 105,195 | 105,195 |
| memory | `memory.bed` | 107,340 | 107,340 |
| sample1 | `sample1.bed` | 113,838 | 113,838 |
| sample2 | `sample2.bed` | 156,734 | 156,734 |

## Consensus Statistics

Peaks reproducibly present across **all** perturbation files.

| Consensus set | Peaks |
| --- | ---: |
| Consensus50_Single | 87,690 |
| Consensus50_Reciprocal | 81,702 |
| Consensus75_Single | 75,380 |
| Consensus75_Reciprocal | 57,859 |

## Memory Candidate Statistics

Consensus peaks that are also open in the memory state (union regions).

| Memory candidate set | Regions |
| --- | ---: |
| Memory50_Single | 76,544 |
| Memory50_Reciprocal | 72,785 |
| Memory75_Single | 56,579 |
| Memory75_Reciprocal | 47,119 |

## Control Filtering Statistics

Strict rule: a memory region overlapping **any** control peak by ≥1 bp is removed entirely.

| Variant | Candidates | Removed by control overlap | Surviving (final) |
| --- | ---: | ---: | ---: |
| 50_Single | 76,544 | 74,219 | 2,325 |
| 50_Reciprocal | 72,785 | 70,971 | 1,814 |
| 75_Single | 56,579 | 55,438 | 1,141 |
| 75_Reciprocal | 47,119 | 46,504 | 615 |

## Final Peak Counts

| Final set | Peaks |
| --- | ---: |
| StrictMemory50_Single | 2,325 |
| StrictMemory50_Reciprocal | 1,814 |
| StrictMemory75_Single | 1,141 |
| StrictMemory75_Reciprocal | 615 |

## Peak Length Statistics (final sets)

| Final set | Count | Mean (bp) | Median (bp) | Min (bp) | Max (bp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| StrictMemory50_Single | 2,325 | 329.8 | 292.0 | 86 | 2,169 |
| StrictMemory50_Reciprocal | 1,814 | 291.5 | 250.0 | 86 | 1,790 |
| StrictMemory75_Single | 1,141 | 304.4 | 259.0 | 86 | 2,169 |
| StrictMemory75_Reciprocal | 615 | 236.0 | 183.0 | 86 | 1,790 |

## Method Notes

- **Single** criterion: overlap ≥ threshold of the *reference* peak only.
- **Reciprocal** criterion: overlap ≥ threshold of *both* peaks.
- Consensus is anchored on the first perturbation file; a reference peak 
  is retained only if every other perturbation file contributes a qualifying overlap.
- Memory regions are the union span of each consensus peak and *all* its 
  qualifying memory peaks (not the single best/first overlap).
- BED coordinates are 0-based, half-open; overlap = max(0, min(ends) − max(starts)).
