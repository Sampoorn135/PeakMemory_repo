# Cross-Dataset Memory-Peak Overlap: IMQ vs HFSC

Final `StrictMemory` peaks compared by genomic coordinate (both mm10). 
Overlap rule: **≥ 1 bp** shared.

- Dataset A = **IMQ**
- Dataset B = **HFSC**

## Overlap by stringent condition

| Condition | IMQ peaks | HFSC peaks | IMQ∩HFSC (A-side) | IMQ∩HFSC (B-side) | Overlapping pairs | Jaccard (bp) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 50_Single | 2,325 | 2,207 | 11 (0.5%) | 11 (0.5%) | 11 | 0.0020 |
| 50_Reciprocal | 1,814 | 2,114 | 5 (0.3%) | 5 (0.2%) | 5 | 0.0008 |
| 75_Single | 1,141 | 990 | 2 (0.2%) | 2 (0.2%) | 2 | 0.0003 |
| 75_Reciprocal | 615 | 825 | 1 (0.2%) | 1 (0.1%) | 1 | 0.0001 |

## How to read this

- **A-side overlap** = number of IMQ peaks that hit ≥1 HFSC peak.
- **B-side overlap** = number of HFSC peaks that hit ≥1 IMQ peak.
- The two sides differ because one peak in one set can overlap several in the other,
  and peak widths differ between the datasets.
- **Jaccard (bp)** = shared base pairs / union base pairs (0 = none, 1 = identical).
- The actual overlapping IMQ peaks per condition are in `overlapping_peaks/`.
