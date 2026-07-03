# Gene-Level Overlap of Memory Peaks: IMQ vs HFSC (mm10)

Final `StrictMemory` peaks were annotated to the **nearest TSS** using the supplied HOMER mm10 reference (`mm10.tss` / `mm10.rna`), entirely in Python (no HOMER binary). Transcript isoforms were collapsed into **gene loci** (same-strand overlapping gene bodies), so overlap is gene-level. No distance cut-off (nearest gene always assigned).

> Gene identifiers are representative **RefSeq accessions** (this HOMER genome folder has no symbol map). Counts are final; symbols can be relabelled later.

## Peaks annotated

| Dataset | Condition | Peaks | Assigned a gene | Gene loci hit |
| --- | --- | ---: | ---: | ---: |
| IMQ | 50_Single | 2,325 | 2,325 | 2,055 |
| IMQ | 50_Reciprocal | 1,814 | 1,814 | 1,661 |
| IMQ | 75_Single | 1,141 | 1,141 | 1,064 |
| IMQ | 75_Reciprocal | 615 | 615 | 593 |
| HFSC | 50_Single | 2,207 | 2,207 | 1,983 |
| HFSC | 50_Reciprocal | 2,114 | 2,114 | 1,906 |
| HFSC | 75_Single | 990 | 990 | 931 |
| HFSC | 75_Reciprocal | 825 | 825 | 781 |

## Gene overlap by stringent condition

| Condition | IMQ genes | HFSC genes | Shared | IMQ-only | HFSC-only | Shared/IMQ | Shared/HFSC | Jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50_Single | 2,055 | 1,983 | **356** | 1,699 | 1,627 | 17.3% | 18.0% | 0.097 |
| 50_Reciprocal | 1,661 | 1,906 | **265** | 1,396 | 1,641 | 16.0% | 13.9% | 0.080 |
| 75_Single | 1,064 | 931 | **94** | 970 | 837 | 8.8% | 10.1% | 0.049 |
| 75_Reciprocal | 593 | 781 | **41** | 552 | 740 | 6.9% | 5.2% | 0.031 |

## Pooled overlap (union of genes across all 4 conditions)

| | IMQ | HFSC | Shared | IMQ-only | HFSC-only | Jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| gene loci | 2,055 | 1,983 | **356** | 1,699 | 1,627 | 0.097 |

## Files

- `annotated/` — per-peak nearest gene + signed distance to TSS.
- `gene_lists/` — gene loci (representative RefSeq) hit by each set.
- `shared_genes/Shared_*` — shared loci per condition and pooled; `*_only_*` = dataset-specific.
