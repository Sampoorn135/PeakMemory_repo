# Peak Memory — Methods and Results

*A cross-tissue meta-analysis of epigenetic injury/inflammation memory, its conserved AP-1 regulatory program, and its fate in aging and cancer.*

Sampoorn Pandey

---

## Summary

Inflammation and injury leave a **chromatin memory**: ATAC-seq peaks that open during a
perturbation and remain accessible after it resolves. We re-derived memory peaks identically
in four independent mouse ATAC-seq datasets and asked whether the memory converges across
tissues, what its core regulatory unit is, and what becomes of it in aging and cancer.

The exact memory **coordinates do not overlap** across tissues (0 of ~10,600 union regions
shared by all four), and the modest **gene-level overlap is fully explained by shared
chromatin accessibility** (≤1× over an accessibility-matched null). What *is* convergent is the
**transcription-factor program**: an **AP-1 / KLF / TEAD / CTCF** motif signature shared ~10×
more than chance (z ≈ 25, p < 5×10⁻⁵). Distilling memory peaks in epidermal stem cells to
**PAME — Primed AP-1 Memory Enhancers** (memory ∩ AP-1 motif ∩ AP-1 binding ∩ H3K27ac;
209 regions / 202 genes) yields elements that sit on the **cell-adhesion / junction /
cytoskeleton (EMT/invasion)** machinery. In aging these enhancers are **relatively protected**
from the genome-wide accessibility decline (directionally specific but a weak effect at the
animal level), and in cancer a **14-gene AP-1 signature is prognostic in TCGA-HNSC
(HR ≈ 1.50, Cox p = 0.014)**.

---

## 1. Methods

### 1.1 Datasets

All mouse analyses use genome build **mm10**; the human cancer arm uses **hg38** (PAME peaks
lifted mm10 → hg38). Mouse → human gene mapping uses orthology (`mousipy`, HGNC/HCOP/BioMart).
This is a meta-analysis of published, processed data; no primary sequencing was generated. The
full accession manifest is in [`DATA.md`](DATA.md).

**Memory-peak discovery (mouse ATAC-seq):**

| Short name | Accession | Tissue / model | Control → perturbation → memory |
|---|---|---|---|
| **IMQ** | GSM5229061–GSM5229068 | Imiquimod psoriasiform skin inflammation; epidermal stem cells | D6/D30 Ctrl → D6 IMQ → D30 post-inflammation |
| **HFSC** | GSE165312 | Hair-follicle stem cells, Dremel skin wound | Ctrl_D80 → Dremel_D3/D7 → Dremel_D80 |
| **distil** | GSE221408 | Distal (unwounded) wound-memory epidermis | per-condition IDR conservative narrowPeak |
| **pancreas** | (epithelial pancreas injury) | Pancreatic epithelium injury/recovery | CTRL → Day1/Day7 → Day28 |

**Supporting tracks for PAME:** mouse keratinocyte **H3K27ac** ChIP-seq (GSE86900, MK1/MK2,
build confirmed mm10); **UniBind** robust AP-1 TFBS (mm10); HOMER vertebrate motif PWMs;
UCSC/HOMER mm10 TSS and chrom-size references.

**Aging:** mouse aging single-nucleus ATAC atlas **GSE288730** (Lu, Zhou & Cao 2025), skin
matrix subset to interfollicular epidermal basal cells → 11,326 nuclei × 1,341,077 peaks,
mm10, three ages.

**Cancer / pre-cancer:** human normal → actinic keratosis (AK) → cSCC ATAC-seq **GSE277274**;
**TCGA-HNSC** RNA-seq (STAR counts) + clinical from the GDC.

### 1.2 Memory-peak discovery pipeline

The pipeline (`src/01_memory_pipeline/`) operates on **BED interval coordinates only** — no BAM
files, read counts, signal, fold-changes, or p-values — so a "memory peak" is defined purely by
reproducible presence/absence across conditions. For each dataset:

1. **Sort + internal merge** of every input peak file (control, perturbation replicates, memory).
2. **Consensus** = peaks reproducibly present across **all** perturbation replicates. Four
   variants are computed by overlap rule: 50% vs 75% reciprocal-fraction, and Single vs
   Reciprocal direction. The primary analysis set is **50_Single**.
3. **Memory candidates** = consensus perturbation peaks that are also open in the memory state.
4. **Strict control filtering** = any memory candidate overlapping *any* control peak by ≥1 bp is
   removed entirely (i.e. we keep only regions that were closed at baseline, opened on
   perturbation, and stayed open). The survivors are the **StrictMemory** peaks.

### 1.3 Cross-dataset convergence and permutation nulls

Convergence (`src/02_cross_dataset/`) was tested at three levels on the StrictMemory50_Single
sets: (i) **coordinates** — peaks merged into a union; a region is "shared" if every member
dataset has a peak overlapping it (≥1 bp); (ii) **genes** — peaks annotated to the nearest TSS
using the HOMER mm10 reference, with transcript isoforms collapsed to gene loci (no distance
cut-off); (iii) **TF motifs** — HOMER known-motif enrichment (`findMotifsGenome.pl -size 200`,
GC-matched background), enriched TF = Benjamini **q < 0.05**, motifs collapsed to TF name.

Significance of the four-way overlaps was assessed by **permutation nulls** (`null_test.py`):
for genes, (A) resampling from the accessible-gene background and (B) genomic peak-shuffle with
re-annotation; for TF motifs, resampling from the tested-motif universe (N = 20,000
permutations). The whole-genome null is reported only to show why it is *inappropriate* here.

### 1.4 PAME construction (epidermal stem cells, IMQ)

Starting from the 2,325 IMQ StrictMemory50_Single peaks (`src/03_pame/`):

- **+ AP-1 motif** — HOMER AP-1(bZIP) PWM present in the peak → 656 regions.
- **+ AP-1 binding (Set B)** — overlaps a UniBind robust AP-1 ChIP TFBS → **467 regions**.
- **+ active enhancer (Set A / PAME)** — H3K27ac signal above background, where "active" = mean
  H3K27ac over the region ≥ the 95th percentile of 2,000 random genomic regions (threshold
  1.059), read directly from the GSE86900 bigWig (pure-Python reader) → **209 regions**.

The 258 Set B regions not in Set A are "poised/primed-only" (memory + AP-1 motif + AP-1 bound
but H3K27ac-low). Set A regions were annotated to nearest TSS and GO/pathway enrichment run with
HOMER `annotatePeaks.pl` / `findGO.pl` (mouse).

### 1.5 Aging analysis

Set A/B and the AP-1-bound regions were intersected with the GSE288730 epidermal-basal-cell
atlas (`src/04_aging/`). Because single-cell ATAC has many non-independent cells per animal,
the primary test was **pseudobulk / sample-level** (animal as the experimental unit) with an
accessibility-matched permutation null; per-cell tests are reported alongside but treated as
anti-conservative (pseudoreplication). Three complementary analyses were run: pseudobulk
log2(Aged/Young), a Mann-Whitney rank-skew of per-peak log2FC (set vs rest), and per-cell
Set-A enrichment by age, plus the animal-level Young-vs-Aged comparison.

### 1.6 Cancer / pre-cancer

PAME peaks were lifted **mm10 → hg38** with `pyliftover` (UCSC chain); uniquely-mapped regions
were retained (`src/05_cancer_precancer/`). The 202 mouse PAME genes were mapped to human
orthologs (`mousipy`) → a **183-gene human PAME signature**. For TCGA-HNSC
(`src/06_tcga_survival/`), STAR gene counts were filtered (zero-count and low-expression genes
removed) and converted to log-CPM; per-gene Cox proportional-hazards and Kaplan-Meier tests
were run, and signature **module scores** (mean z-scored expression) were tested with a Cox
model on overall survival. The pre-cancer GSE277274 quantification (lifted PAME regions over
normal/AK/cSCC ATAC bigWigs) is set up in `quantify_atac.py`.

---

## 2. Results

### 2.1 Memory peaks are recovered in all four tissues

The strict pipeline yields a few thousand high-confidence memory peaks per dataset
(StrictMemory50_Single):

| Dataset | Candidates | Removed by control overlap | **StrictMemory peaks** |
|---|---:|---:|---:|
| IMQ | 76,544 | 74,219 | **2,325** |
| HFSC | 20,674 | 18,467 | **2,207** |
| distil | 27,686 | 24,590 | **3,096** |
| pancreas | 46,787 | 43,526 | **3,261** |

The strict control filter is aggressive (it removes >90% of candidates in skin), leaving only
regions that were genuinely closed at baseline and remained open after the perturbation
resolved.

### 2.2 Memory coordinates do not overlap across tissues

Merging all four StrictMemory sets gives a union of **10,604 regions**. **No region is shared by
all four datasets**, and three-way sharing is essentially nil (best three-way = 3 regions). The
strongest pairwise overlap is HFSC ∩ distil = 149 regions — both skin, as expected.

| Combination | Order | Shared regions |
|---|---:|---:|
| HFSC ∩ distil | 2 | 149 |
| IMQ ∩ distil | 2 | 45 |
| distil ∩ pancreas | 2 | 36 |
| IMQ ∩ pancreas | 2 | 25 |
| HFSC ∩ pancreas | 2 | 23 |
| IMQ ∩ HFSC | 2 | 11 |
| HFSC ∩ distil ∩ pancreas | 3 | 3 |
| IMQ ∩ distil ∩ pancreas | 3 | 1 |
| **all four** | 4 | **0** |

### 2.3 Gene overlap is real but explained by shared accessibility

Annotating to nearest TSS, the four sets hit 2,055 (IMQ), 1,983 (HFSC), 2,522 (distil) and 2,615
(pancreas) gene loci, with **64 genes shared by all four**. Against a naïve whole-genome null
that looks ~37× enriched — but that null is wrong. When the correct **accessible-gene background**
is used, the observed 64 is *below* expectation:

| Null model | Expected four-way genes | Fold | z | p |
|---|---:|---:|---:|---:|
| A — accessible-gene background | 84.5 | 0.8× | −2.4 | 0.995 |
| B — peak-shuffle + re-annotate | 149.5 | 0.4× | −9.1 | ~1.0 |

All six pairwise gene overlaps are likewise at or below their accessible-background expectation
(≈0.5–0.7×). The apparent gene convergence is an artifact of all four datasets sampling memory
genes from the same accessible-gene pool — **not** a shared memory-gene program.

### 2.4 The transcription-factor program is convergent (~10× over chance)

HOMER motif enrichment gives 149 (IMQ), 141 (HFSC), 167 (distil) and 174 (pancreas) enriched TFs
(q < 0.05). **68 TFs are enriched in all four datasets** — and here the convergence is genuine:

| Null model | Expected four-way TFs | Fold | z | p |
|---|---:|---:|---:|---:|
| Resample from tested-motif universe (450 TFs; N = 20,000) | 6.7 | **10.2×** | **25.3** | **< 5×10⁻⁵** |

The shared program is dominated by **AP-1 (Fos/Jun/Fosl/Atf/BATF/Fra), KLF (KLF1/3/4/5/6/…),
TEAD (TEAD1–4) and CTCF**, plus p53/p63/p73 and Sp1/2. The four-way shared set was never
approached in 20,000 permutations. (Caveat: motifs within a family are correlated, so this null
is anti-conservative — treat p as an upper bound; the z ≈ 25 margin keeps the conclusion robust.)
The three-way skin overlap (IMQ ∩ HFSC ∩ distil) is even tighter at 86 shared TFs.

**Central claim:** injury/inflammation memory converges at the **trans-regulatory (TF) level** —
a conserved AP-1/KLF/TEAD/CTCF enhancer-memory program — implemented through largely
dataset-specific *cis*-elements and target genes.

### 2.5 PAME: distilling memory to active AP-1 enhancers

In epidermal stem cells the memory peaks funnel cleanly to an active AP-1 compartment:

| Step | Filter | Regions |
|---|---|---:|
| Memory (StrictMemory50_Single) | retained-accessible after perturbation | 2,325 |
| + AP-1 motif | HOMER AP-1(bZIP) PWM | 656 |
| **Set B** + AP-1 binding | overlaps UniBind AP-1 TFBS | **467** |
| **Set A / PAME** + H3K27ac active | active enhancer | **209** |

PAME regions are strongly enriched for active-enhancer status: **45% of AP-1-primed regions are
H3K27ac-active vs 5% of random background (8.9× enrichment; hypergeometric p ≈ 10⁻⁹⁴)**, with
median H3K27ac ~3× background in Set B. They are predominantly **distal enhancers** (9% promoter,
66% 2–50 kb, 14% 50–100 kb, 11% >100 kb), annotating to **202 unique genes** (180 within 100 kb).

### 2.6 PAME genes converge on adhesion / junction / cytoskeleton (EMT/invasion)

The 202 PAME genes are enriched for **cell adhesion, cell junctions and the actin cytoskeleton** —
exactly the machinery dismantled during EMT and epithelial-cancer invasion, and AP-1 (FOS/JUN) is
the canonical driver of that program:

| Database | Top specific term | p | Example genes |
|---|---|---:|---|
| KEGG | Adherens junction | 1.9×10⁻⁶ | Ptpn6, Pard3, Wasf2, Nectin1, Sorbs1, Actn1, Lmo7 |
| Reactome | Cell–cell communication | 1.4×10⁻⁵ | Col17a1, Lims1, Pard3, Nectin1, Cadm1, Ptpn6, Actn1 |
| Reactome | Cell-junction organization | 1.5×10⁻⁵ | Cadm1, Nectin1, Pard3, Actn1, Col17a1, Lims1 |
| GO-MF | actin / cytoskeletal-protein binding | 1.1×10⁻⁵ | Myo1d, Mrtfa, Phactr1, Arhgef7, Lmo7 (12 total) |
| KEGG | Axon guidance / Sema3A–Plexin | 5.9×10⁻⁴ | Nrp1, Plxna1, Slit3, Pard3, Gnai3 |

Notable drivers include **Col17a1** (hemidesmosome collagen linking to EpdSC aging and skin SCC),
the adhesion/polarity tumor suppressors **Cadm1, Nectin1, Pard3**, and actin/migration effectors
**Wasf2 (WAVE2), Lims1 (PINCH), Actn1, Mrtfa, Mical2, Arhgef7**. (Broad developmental GO terms
also rank highly, as expected for any ~200-gene set near tissue-identity genes; the specific
adhesion/cytoskeleton signal is the robust, interpretable finding.)

### 2.7 Aging: PAME enhancers resist the genome-wide accessibility decline

Across three complementary tests, PAME/AP-1-primed enhancers are **skewed toward higher
accessibility in older epidermal stem cells** while background regions decline:

| Analysis | Test | Set A result |
|---|---|---|
| Pseudobulk log2(Aged/Young) | accessibility-matched permutation | **+0.091** vs null −0.060, **p = 5×10⁻⁴** |
| Rank skew toward age-gain | Mann-Whitney, set vs rest | **AUC 0.598, p = 1.8×10⁻⁷** |
| Per-cell Set-A enrichment | by age | median 2.44 → 2.91 → 2.97; KW **p = 1.5×10⁻¹¹** |

The direction is specific (background *loses* accessibility with age, ≈ −0.06; the primed set
*gains*, +0.06 to +0.09) and reflects a broad population shift rather than a rare expanding clone.
**Important rigor caveat:** when the *animal* is the experimental unit (the correct way to claim an
aging effect), the increase is **not significant** (Young vs Aged p = 0.15; Cohen's d = 0.07 per
cell) and is non-monotonic (peaks at Adult). The very small per-cell p-values are inflated by
pseudoreplication. So this is a directionally consistent, region-class-specific trend, but a weak
and statistically unconfirmed age effect — convincing for its *consistency and direction*, not its
magnitude (~1.08× pseudobulk).

### 2.8 Cancer: a 14-gene AP-1 signature is prognostic in TCGA-HNSC

Of the 209 PAME regions, **116 (55%) lifted uniquely to hg38**. The 202 mouse genes mapped to a
**183-gene human signature**; 114 were measurable in TCGA-HNSC. Per-gene Cox analysis found 16
genes nominally associated with overall survival (P < 0.05; e.g. MSANTD3, FAM53B, PRR7, CA12,
GADD45B, UBE2L3, ACTN1, MYOF), though none survive FDR correction. Tested as **module scores** in
the TCGA-HNSC cohort (477 patients, 220 deaths):

| Signature | Genes | Hazard ratio (95% CI) | Cox p | Log-rank p |
|---|---:|---|---:|---:|
| 106-gene AP-1 signature | 106 | 1.31 (0.88–1.93) | 0.18 | 0.15 |
| **14-gene AP-1 signature** | 14 | **1.50 (1.09–2.08)** | **0.014** | 0.058 |

Higher activity of the focused **14-gene AP-1-primed signature predicts worse overall survival**
(HR ≈ 1.50, Cox p = 0.014; concordance 0.54). This provides a direct, prognostic link from the
mouse injury-memory enhancers to human squamous-cancer outcome. The pre-cancer
(normal → AK → cSCC) accessibility-trend test over GSE277274 is set up but exploratory (small n
per state; bulk ATAC mixes tumor/stroma/immune).

---

## 3. Limitations

- **Coordinate vs program.** The convergent signal is at the TF-motif level; shared coordinates
  are ~0 and shared genes are explained by accessibility. Claims should stay at the regulatory-
  program level, not "shared memory loci".
- **Motif-null correlation.** The 10× TF enrichment uses a null that treats motifs as independent;
  family correlation makes it anti-conservative (upper bound on significance).
- **Aging effect is weak.** Significant at the (pseudoreplicated) cell level but not at the animal
  level; only three age bins; modest effect size.
- **Cross-species loss.** Only 116/209 PAME regions survive liftOver; the human signature drops
  non-conserved enhancers. The TCGA result is gene-expression-based, not chromatin-based.
- **No FDR-significant single genes** in TCGA; the signal is at the aggregate-signature level.
- **Pre-cancer arm is exploratory** (small n, bulk composition).

---

## 4. Data and code availability

All accessions are listed in [`DATA.md`](DATA.md). Analysis code is in `src/` (organized by stage
`01_memory_pipeline` → `06_tcga_survival`); small final outputs and per-stage reports are in
`results/`. Large inputs (h5ad atlas, TCGA matrices, bigWigs) are not committed and are
reproducible from the manifest.
