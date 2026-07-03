# AP-1-primed regions in epidermal stem cells — Set A (with H3K27ac)

**Cell type:** interfollicular epidermal stem cells (IMQ skin dataset, mm10)

## The funnel
| Step | Filter | Regions |
|---|---|---|
| Memory peaks (StrictMemory50_Single) | retained-accessible after perturbation | 2,325 |
| + AP-1 motif | HOMER AP-1(bZIP) PWM in the peak | 656 |
| **Set B** = + AP-1 binding | overlaps a UniBind AP-1 ChIP TFBS | **467** |
| **Set A** = + H3K27ac active | H3K27ac signal above background | **209** |

The 258 Set B regions **not** in Set A are "**poised/primed-only**" (memory + AP-1 motif + AP-1 bound, but H3K27ac-low) — a latent compartment distinct from the 209 **active** AP-1-primed enhancers in Set A.

## H3K27ac source
GSE86900 — **mouse keratinocyte** (MK1/MK2) H3K27ac ChIP-seq (Sinha lab). Read directly from the
GEO bigWig; genome build **confirmed mm10** (chr1 = 195,471,971), so coordinates match Set B exactly.
"Active" = mean H3K27ac over the region ≥ the 95th percentile of 2,000 random genomic regions (threshold 1.059).

## Result (sanity check passed)
- **45% of AP-1-primed regions are H3K27ac-active** vs **5% of random background → 8.9× enrichment**.
- Median H3K27ac signal 3× higher in Set B than background.
- Hypergeometric enrichment **p ≈ 10⁻⁹⁴**.
This confirms AP-1-primed memory regions are strongly biased toward active epidermal enhancers — not a random subset.

## Set A gene annotation (nearest TSS)
- 209 regions → **202 unique nearest genes** (180 within 100 kb).
- Location: 9% promoter (<2 kb), 66% 2–50 kb, 14% 50–100 kb, 11% >100 kb → predominantly **distal enhancers**.

## Files
- `SetA_AP1primed_H3K27ac.bed` — 209 active AP-1-primed regions (col4 = H3K27ac signal)
- `SetB_H3K27ac_signal.tab` — all 467 Set B regions with H3K27ac signal + active flag (0/1)
- `background_H3K27ac_signal.tab` — 2,000 background regions (threshold reference)
- `SetA_peak_gene.tsv`, `SetA_genes_within100kb.txt`, `SetA_genes_nearestTSS.txt`

## Next step (for the cancer-priming hypothesis)
Convert the Set A gene list (currently RefSeq IDs) → gene symbols → human orthologs to build the
**AP-1-primed gene signature**, then test in TCGA (e.g., HNSC/SKCM/cervical squamous) whether high
signature activity in tumors associates with worse survival.
