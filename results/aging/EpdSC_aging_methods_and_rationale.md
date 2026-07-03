# Aging epidermal-stem-cell analysis — full workflow, verification & rationale

## 1. Data source and upstream processing (done by the atlas authors, not us)
- **Atlas:** GSE288730 — *"Organism-wide cellular dynamics and epigenomic remodeling in mammalian aging"* (Lu, Zhou & Cao, Rockefeller University; 2025; PMID 40463164).
- **Scale:** ~7 million single nuclei, 21 tissues, both sexes, 3 ages → 536 main cell types / 1,828 subtypes, on a shared set of **~1.3 M cis-regulatory elements (peaks), mm10**.
- **Assay:** **EasySci-ATAC-seq** — a combinatorial-indexing single-**nucleus** ATAC method (PMID 38036784), sequenced on Illumina NovaSeq 6000; C57BL/6 wild-type mice.
- **What the authors did:** isolated nuclei → **DAPI-gated/sorted intact single nuclei** (the "DAPI" in their library names; this is nuclear-stain gating, **not** antibody/surface-marker FACS) → EasySci combinatorial indexing → aligned reads, called a common ~1.3 M-peak set, built the **cell × peak count matrix** → per-nucleus QC (reads/nucleus, promoter ratio) → unsupervised clustering of accessibility → **annotated cell types**. The skin organ is "back skin" (sample GSM8774016).
- **Important:** the label *"Interfollicular epidermal basal cells"* is the atlas's **computational** annotation (chromatin-based clustering), not a sorted/marker-purified population. We took these labels as given.
- **We did not** do alignment, peak-calling, clustering, or annotation. We started from the authors' processed cell × peak matrix.

## 2. From the atlas to our working dataset (what we did)
- Downloaded the full **skin** matrix, then subset to the epidermal-stem-cell compartment with `subset_epdsc.py` (regex match `interfollicular epidermal basal`, low-RAM backed read, keeping **all** such nuclei across every age and sex and **all** peaks).
- Output: `EpdSC_peak_count.h5ad` = **11,326 nuclei × 1,341,077 peaks**.

## 3. Cell verification & timepoints (what / how / why)
| Check | How | Result | Why |
|---|---|---|---|
| Cell-type purity | `Main_cell_type` unique values | **1** value: "Interfollicular epidermal basal cells" | confirm clean subset, no contaminating types |
| Identity match | compare to IMQ study | interfollicular **basal** layer = the epidermal stem/progenitor compartment we used in IMQ | ensures cross-dataset comparison is like-for-like |
| Timepoints | `Age` field ↔ atlas design | **Young = 1 mo (4,968), Adult = 5 mo (3,474), Aged = 21 mo (2,884)**; balanced M/F; 21 skin mice (8 / 7 / 6) | defines the developmental vs aging windows; checks power per group |
| Genome build | max peak coordinate vs mm10/mm39 | **mm10** (chr1 max-end 195.31 M > mm39 195.15 M, < mm10 195.47 M) | coordinate-compatible with Set A/B → no liftOver |
| Count sanity | sample of X matrix | uint32 counts (0–10), median ~3,600 counts & ~2,540 peaks per nucleus; QC fields (`total_reads`, `promoter_ratio`) present | confirms real counts (not pre-binarized) and adequate depth |

Verdict: suitable for downstream analysis.

## 4. Foundation for the analyses (what / how)
- **Map our regions onto the atlas peak space** (both mm10) by coordinate overlap (bounded-window interval search): Set A 209 → **224** atlas peaks; Set B 467 → **497**; UniBind AP-1 TFBS → **420**.
- **Pseudobulk by age:** summed the count matrix per peak within each age group (memory-safe CSR block reads), then **CPM-normalized** per age (removes differences in cell number and depth).
- Why: lets us ask, for the *same* genomic loci defined in the IMQ epidermis, how accessibility behaves across the atlas's age axis.

## 5. The four analyses — what, how, which, why
1. **Pseudobulk accessibility vs age (A1) + AP-1 binding (A2).** *How:* per-peak log2(fold-change) between ages, compared to an **accessibility-matched permutation null** (random peaks matched on baseline CPM, 2,000 permutations). *Which sets:* Set A, Set B, AP-1-bound. *Why:* tests whether the primed region **class** changes accessibility beyond what background regions of the same baseline accessibility do — the matched null is essential because these are high-accessibility enhancers.
2. **Rank-based skew (A3).** *How:* Mann-Whitney U of per-peak log2FC, set vs all other peaks (reports AUC). *Why:* uses every peak (no arbitrary fold-change cutoff), so it is properly powered — the hard-threshold DAR count (~500 genome-wide) was too sparse to test overlap.
3. **Single-cell + sample-level (A4).** *How:* a per-nucleus accessibility score over the Set-A peaks (depth-normalized), tested across ages two ways — **cells as the unit** (the "usual" single-cell test) and **mice as the unit** (per-sample means, 21 mice). *Why:* the single-cell view shows heterogeneity, but cells from one mouse are not independent, so the **animal-level** test is the statistically valid one for an age effect (avoids pseudoreplication).
4. **Developmental vs aging decomposition.** *How:* split the age axis into 1→5 mo (development) and 5→21 mo (true aging) and run A1's matched-null test on each. *Why:* "Young = 1 month" is juvenile, so a single Young→Aged contrast conflates maturation with aging.

## 6. Results (refined)
- **Development (1→5 mo):** Set A +0.122 vs background +0.092 (p = 0.06, **n.s.**) — a **global** developmental opening; primed regions open along with the rest of the genome, not specifically.
- **Aging (5→21 mo):** background **declines** (−0.075) but Set A/B decline **less** (−0.03 / −0.02), significantly above background (**p = 0.015 / 0.001**) → **relative maintenance: AP-1-primed enhancers resist the genome-wide accessibility loss of aging.**
- **Single-cell vs animal-level:** cells-as-unit p = 5×10⁻¹⁰ but effect size negligible (Cohen's d = 0.07); **mice-as-unit p = 0.15 (n.s.)** — the cell-level p is pseudoreplication-inflated.

## 7. Honest conclusion
The defensible statement is **not** "AP-1-primed enhancers gain accessibility with age," but: *during postnatal maturation the epidermal genome opens globally (primed regions included, non-specifically); during true aging the genome loses accessibility and the AP-1-primed memory enhancers are relatively protected from that loss (peak-class p ≈ 0.001–0.015).* This is a subtle, region-class-level effect; it is **not** robust as a per-animal increase (n.s.), and the very small single-cell p-values overstate it. Confirming it would need more aged mice and/or an orthogonal aging cohort. Files: `AP1primed_aging_EpdSC.png`, `SetA_persample_age.png`, `EpdSC_aging_results.md`, `SetA_aging_peaks.csv`, `persample_setA.csv`.
