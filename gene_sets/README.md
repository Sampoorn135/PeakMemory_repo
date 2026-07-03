# Gene sets

Curated gene lists, signatures, and region (BED) files produced by the Peak Memory pipeline,
gathered here in one place for reuse. Code that generates them is in [`../src/`](../src);
the narrative is in [`../docs/Methods_and_Results.md`](../docs/Methods_and_Results.md).

All mouse files are mm10 and use gene symbols (or RefSeq where noted); human files are hg38 / HGNC symbols.

## 01_memory_genes_per_dataset
Genes near each dataset's StrictMemory peaks (nearest-TSS / 2 kb window), per dataset and per
stringency variant (`50/75` × `Single/Reciprocal`). Datasets: IMQ, HFSC, distil, pancreas.

## 02_memory_genes_shared
Cross-dataset comparisons of the memory-gene lists: `Shared_*` (common to datasets),
`*_only_*` (dataset-specific), and `Shared_POOLED`. `*_SYMBOLS` files are symbol-mapped versions.
(Key result: this overlap is **not** enriched above an accessibility-matched background.)

## 03_pame_core
The PAME funnel core. Set B = memory ∩ AP-1 motif ∩ AP-1 binding (467 regions);
Set A / **PAME** = Set B ∩ H3K27ac-active (209 regions / 202 genes).
- `SetA_AP1primed_H3K27ac.bed`, `SetB_memory_AP1motif_bound.bed`, `IMQ_memory_AP1motif.bed` — regions
- `SetA_gene_symbols.txt` (202 mouse genes), `SetA_genes_nearestTSS.txt`, `SetA_genes_within100kb.txt`
- `SetA_human_signature.txt` (183 human orthologs), `SetA_mouse2human_orthologs.tsv`
- `SetA_peak_gene.tsv`, `SetB_peak_gene.tsv` — per-peak gene assignment
- `STRING_gene_list.txt` — input to the STRING network

## 04_pame_pathways
Ontology / pathway enrichment of the PAME (Set A) genes: `GO_BP/CC/MF.csv`, `KEGG.csv`,
`Reactome.csv`, `WikiPathways.csv`, and the curated `SetA_top_pathways.tsv`
(adhesion / junction / cytoskeleton signal).

## 05_aging_dynamics
Behaviour of PAME enhancers in the epidermal aging atlas (GSE288730).
- `Age_gained_genes.txt` / `Age_lost_genes.txt` / `Stable_genes.txt` and matching `*_peaks.bed`
- `SetA_age_dynamic_peak_groups.csv`, `SetA_all_atlas_overlap_genes.txt`, `SetA_aging_peaks.csv`
- `persample_setA.csv` — per-animal (sample-level) scores used for the rigorous test
- `dar_overlap/` — overlap of PAME/memory sets with age differentially-accessible regions
- `motifs/` — HOMER AP-1 motif enrichment for AgedVsAdult and AgedVsYoung DARs

## 06_cancer_signature
Human / cancer arm.
- `SetA_AP1primed_hg38.bed` — 116 PAME regions lifted mm10→hg38; `liftover_drops.txt` — dropped regions
- `atac_signal_long.csv` — PAME accessibility across normal / AK / cSCC (GSE277274)
- `PAME_human_signature_183genes.txt` — full human PAME signature
- `TCGA_14gene_AP1_signature.txt` — the 14-gene prognostic signature (TCGA-HNSC; HR ≈ 1.50, Cox p = 0.014)
