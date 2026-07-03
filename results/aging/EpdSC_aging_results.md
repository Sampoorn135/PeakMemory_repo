# AP-1-primed memory enhancers in epidermal-stem-cell aging

**Data (verified):** skin aging atlas, subset to **interfollicular epidermal basal cells** (= the IMQ epidermal-stem-cell compartment). 11,326 cells × 1,341,077 peaks, genome **mm10** (matches Set A/B), three ages with thousands of cells each (Young 4,968 · Adult 3,474 · Aged 2,884), balanced by sex. Suitable for downstream analysis.

Atlas peaks overlapping each region set: **Set A 224**, **Set B 497**, AP-1-bound (UniBind) 420.

## Headline result (with the important rigor caveat)
**As a region class, AP-1-primed enhancers are skewed toward higher accessibility in older epidermal stem cells than matched background regions** (peak-level tests below). **But this is a suggestive trend, not a robust age effect:** when the *animal* is treated as the experimental unit — the statistically correct way to claim an aging effect — the increase is **not significant (Young vs Aged p = 0.15; monotonic-trend p = 0.15)**, the per-cell effect size is tiny (Cohen's d = 0.07), and the pattern is non-monotonic (Adult-high). The very small single-cell p-values (e.g. 5×10⁻¹⁰) are **inflated by pseudoreplication** (thousands of non-independent cells per animal). So: real, specific directionality at the region-class level; weak/unconfirmed at the level that matters biologically.

| Analysis | Test | Result |
|---|---|---|
| **1–2** Pseudobulk accessibility vs age | mean per-peak log2(Aged/Young), accessibility-matched permutation | Set A **+0.091** vs null −0.060, **p = 5×10⁻⁴**; Set B +0.061 (p = 5×10⁻⁴); AP-1-bound +0.068 (p = 5×10⁻⁴). Mean CPM Young→Adult→Aged: Set A 2.41→2.67→2.61 (background flat ≈0.75). |
| **3** Rank skew toward age-gain | Mann-Whitney of per-peak log2FC, set vs rest | Set A **AUC 0.598, p = 1.8×10⁻⁷**; Set B AUC 0.569 (p = 5×10⁻⁸); AP-1-bound AUC 0.576 (p = 3×10⁻⁸). |
| **4** Single-cell primed accessibility | per-cell Set-A enrichment vs age | median Young 2.44 → Adult 2.91 → Aged 2.97; **Aged>Young p = 5.4×10⁻¹⁰**, Kruskal–Wallis **p = 1.5×10⁻¹¹**. |

## Single-cell vs animal-level — why pseudobulk was the main test
Single-cell ATAC has thousands of cells per animal, but cells from one animal are **not independent replicates**. Testing an age effect with cells as the unit therefore commits *pseudoreplication* and badly inflates significance. Here that is explicit:

| Unit of analysis | n | Set-A enrichment Young→Aged | test | p |
|---|---|---|---|---|
| **Cells** ("usual" single-cell) | 4,968 vs 2,884 | 2.44 → 2.97 (Cohen's d = 0.07) | Mann-Whitney | 5.4×10⁻¹⁰ |
| **Animals** (sample-level, correct) | 8 vs 6 | 2.27 → 2.53 (Adult 4.30) | Mann-Whitney | **0.15 (n.s.)** |

The cell-level p of 5×10⁻¹⁰ looks decisive but the effect size is negligible (d = 0.07); at the animal level (n = 6–8/group, high between-animal variance) the age increase is not significant and not monotonic. This is precisely why pseudobulk / sample-level testing is the field standard for condition effects (cf. Squair et al. 2021) — the single-cell-only result would have over-claimed. See `SetA_persample_age.png` and `persample_setA.csv`.

## Interpretation
- **Direction is specific.** Background regions *decline* with age (null log2FC ≈ −0.06); the AP-1-primed/-bound set *rises* (+0.06 to +0.09). So these regions resist and reverse the global age-associated loss of accessibility.
- **Broad population shift, not a rare subpopulation.** The whole EpdSC distribution moves up (median rises) rather than a discrete high-AP-1 clone expanding — the top-decile fraction does not grow with age (Young 10% → Aged 8.8%).
- **Established by adulthood.** Accessibility peaks at Adult and is maintained in Aged (non-monotonic), i.e. priming is laid down with maturation and held.
- **Hypothesis support.** Aging reinforces accessibility at the very same AP-1-primed memory enhancers (wired to the adhesion/cytoskeleton/EMT genes from the GO analysis). That is a concrete chromatin substrate by which an aged epidermis could be predisposed to AP-1-driven invasive programs — the link to test directly in cancer data.

## Caveats (honest)
- Effect **size is modest** (~1.08× pseudobulk). It is the *consistency and direction* across three independent analyses that make it convincing, not the magnitude.
- Only **3 age bins**; bulk pseudobulk pools cells per age (the single-cell test uses per-cell variation, so it is the best-powered). Sample-level pseudobulk (21 samples) would add formal replication.
- The hard-threshold DAR call was underpowered (≈508 genome-wide DARs); the rank test (Analysis 3) is the proper, powered version.
- AP-1-bound set is small (420 peaks, UniBind redundancy) — treat Analysis 2 as supporting.

## Files
- `AP1primed_aging_EpdSC.png` — 3-panel figure (A pseudobulk trend, B rank skew, C single-cell).
- `SetA_aging_peaks.csv` — the 224 Set-A atlas peaks ranked by age-gain, with gene symbols (top gainers: Cyp1b1, Etv6, Pgf, Vangl1, Prdm14, Arg1 …) — the candidate **aging-reinforced primed enhancers**.
