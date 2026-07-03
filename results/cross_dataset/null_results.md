# Permutation null tests — four-way memory convergence

Tests whether the observed four-way overlaps exceed chance. Peaks/genes:
StrictMemory 50% Single. One-sided p = P(null ≥ observed).

## Peak (coordinate) level
Observed four-way shared peak regions = **0** (best pairwise 149 / ~10,600 union).
Effectively no shared coordinates — far below any chance expectation.

## Gene level (nearest-TSS gene loci)
Observed four-way shared genes = **64**. Set sizes: IMQ 2055, HFSC 1983, distil 2522, pancreas 2615.

| Null model | Expected | Fold | z | p (enrichment) |
| --- | ---: | ---: | ---: | ---: |
| A — resample from accessible-gene background (B = 6,828) | 84.5 | 0.8× | −2.4 | 0.995 |
| B — genomic peak-shuffle + re-annotate | 149.5 | 0.4× | −9.1 | ~1.0 |

Pairwise gene overlaps are likewise **at or below** their accessible-background
expectation (observed/expected ≈ 0.5–0.7× for all six pairs).

**Interpretation:** once you control for the fact that all four datasets draw
memory genes from the *same accessible-gene pool*, the gene-level overlap is **not
enriched — if anything slightly depleted.** The apparent "gene convergence" is
explained by shared chromatin accessibility, not by a shared memory-gene program.
(Against a naïve whole-genome null it would look ~37× enriched — that null is
inappropriate here.)

## TF-motif level
Observed four-way shared enriched TFs (q<0.05) = **68**. Set sizes: IMQ 149,
HFSC 141, distil 167, pancreas 174; tested-motif universe = 450 TFs.

| Null model | Expected | Fold | z | p |
| --- | ---: | ---: | ---: | ---: |
| resample from tested-motif universe (N=20,000) | 6.7 | **10.2×** | **25.3** | **< 5×10⁻⁵** |

**Interpretation:** the shared TF program (**AP-1 / KLF / TEAD / CTCF**) is ~10×
more shared than chance and never approached in 20,000 permutations — highly
significant. Caveat: motifs within a family are correlated, so this null is
anti-conservative; treat as an upper bound on significance, but the margin
(z = 25) is large enough that the conclusion is robust.

## Bottom line
Convergence across these four injury-memory datasets is **specifically at the
trans-regulatory (TF) level**, not at the level of shared coordinates (≈0) or
shared genes (explained by accessibility). The defensible central claim is a
**conserved AP-1/KLF/TEAD enhancer-memory regulatory program**, implemented
through largely dataset-specific cis-elements and target genes.
