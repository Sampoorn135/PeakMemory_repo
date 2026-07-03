# Data manifest

This project is a **meta-analysis of published data**. None of the primary sequencing
data is generated here; everything starts from processed peak calls, count matrices, or
signal tracks downloaded from public repositories. Large files are **not** committed to
the repo (see `.gitignore`); this manifest records exactly what to download and where it
sits in the analysis so the pipeline can be reproduced.

## Genome build
All mouse analyses use **mm10**. The human cancer arm uses **hg38** (PAME peaks are
lifted mm10 → hg38). Mouse → human gene mapping uses orthology (mousipy / HGNC).

## Mouse ATAC-seq datasets (memory-peak discovery)

| Short name | Accession | Tissue / model | Role in pipeline |
|---|---|---|---|
| **IMQ** | GSE (GSM5229061–GSM5229068) | Imiquimod psoriasiform skin inflammation; epidermal stem cells | Memory-peak source + PAME source. Control = D6/D30 Ctrl; perturbation = D6 IMQ; memory/chase = D30 post-inflammation (PI). 8 replicate peak BEDs (`GSM5229061…68_*_peaks.bed.gz`). |
| **HFSC** | GSE165312 | Hair-follicle stem cells, Dremel skin wound | Cross-dataset comparison. Ctrl_D80 → control; Dremel_D3/D7 → perturbations; Dremel_D80 → memory. |
| **distil** | GSE221408 | Distal (unwounded) wound-memory epidermis | Cross-dataset comparison (IDR conservative narrowPeak per condition). |
| **pancreas** | GSE (epithelial pancreas injury time-course) | Pancreatic epithelium injury/recovery | Cross-dataset comparison. CTRL → control; Day1/Day7 → perturbations; Day28 → memory. |

> Replicate → pipeline-role assignment is done by the `prepare_inputs*.py` scripts
> (`src/01_memory_pipeline/`); each script documents the exact filename → role mapping.

## Supporting tracks / references (PAME construction)

| Resource | Accession / source | Use |
|---|---|---|
| Mouse keratinocyte **H3K27ac** ChIP-seq (bigWig) | **GSE86900** | "Active enhancer" filter — defines which AP-1-primed memory peaks are active (Set A / PAME). |
| **UniBind** AP-1 robust TFBS (mm10) | UniBind (Robust collection) | "AP-1 bound" filter — keeps memory∩motif peaks that overlap experimental AP-1 binding. |
| mm10 genome FASTA + `mm10.tss`, `chrom.sizes` | UCSC / HOMER | Motif scanning, nearest-TSS annotation. (`mm10_kept/` holds the small TSS/chrom files used.) |
| HOMER known vertebrate motif PWMs | HOMER install | `findMotifsGenome.pl` motif enrichment; AP-1 PWM for direct scanning. |

## Aging dataset (time / aging-epigenetics arm)

| Resource | Accession | Notes |
|---|---|---|
| Mouse aging single-nucleus ATAC atlas | **GSE288730** (Lu, Zhou & Cao 2025; "Organism-wide cellular dynamics and epigenomic remodeling in mammalian aging") | We subset the **skin** matrix to *interfollicular epidermal basal cells* (`subset_epdsc.py`) → `EpdSC_peak_count.h5ad` = 11,326 nuclei × 1,341,077 peaks, mm10, ages 1/5/21 mo. The `.h5ad` is ~GB-scale and **not committed**. |

## Cancer / pre-cancer arm

| Resource | Accession | Notes |
|---|---|---|
| Human skin normal → actinic keratosis (AK) → cSCC ATAC | **GSE277274** | Pre-cancer priming test; PAME peaks lifted to hg38 and quantified over these bigWigs (`quantify_atac.py` via remote range-read). |
| **TCGA-HNSC** RNA-seq (STAR counts) + clinical | GDC Data Portal (project TCGA-HNSC) | Survival validation. One folder per case under `data/raw/` (STAR `augmented_star_gene_counts.tsv`), plus `clinical.tsv`, `follow_up.tsv`, `gdc_sample_sheet*.tsv`, `metadata.cart*.json`. Built into `raw_counts_matrix.csv` / `normalized_logCPM.csv` → 521-patient cohort. **Not committed** (≈0.4 GB). |

## What *is* committed (small, in `results/`)
Final peak BEDs for the four memory sets, the 209-region PAME BED (`SetA_AP1primed_H3K27ac.bed`),
its mouse (202) and human (183) gene lists, the lifted hg38 PAME BED, the per-stage QC
and overlap/enrichment markdown reports, and the aging peak table. These let a reader
inspect every result without re-downloading raw data.
