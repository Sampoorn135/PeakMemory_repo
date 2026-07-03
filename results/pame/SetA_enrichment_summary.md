# Set A — gene symbols + pathway enrichment

**Input:** 209 AP-1-primed + H3K27ac-active regions → **202 unique gene symbols**
(via HOMER `annotatePeaks.pl mm10`; enrichment via HOMER `findGO.pl mouse`).
Full list: `SetA_gene_symbols.txt`. Full table: `SetA_top_pathways.tsv`.

## Headline: the active AP-1-primed genes converge on CELL ADHESION / JUNCTIONS / CYTOSKELETON

| Database | Top specific term | p | Genes |
|---|---|---|---|
| KEGG | **Adherens junction** | 1.9×10⁻⁶ | Ptpn6, Pard3, Wasf2, Nectin1, Sorbs1, Actn1, Lmo7 |
| Reactome | **Cell–cell communication** | 1.4×10⁻⁵ | Col17a1, Lims1, Pard3, Nectin1, Cadm1, Ptpn6, Actn1 |
| Reactome | **Cell-junction organization** | 1.5×10⁻⁵ | Cadm1, Nectin1, Pard3, Actn1, Col17a1, Lims1 |
| Reactome | Cytoskeletal remodeling & cell spreading (IPP/PINCH) | 6.5×10⁻⁴ | Actn1, Lims1 |
| Reactome | Nectin/Necl trans-heterodimerization | 9.7×10⁻⁴ | Cadm1, Nectin1 |
| Reactome | Cell–extracellular-matrix interactions | 4.1×10⁻³ | Lims1, Actn1 |
| GO-MF | **actin binding** / cytoskeletal-protein binding | 1.1×10⁻⁵ | Myo1d, Mrtfa, Phactr1, Arhgef7, Lmo7 … (12) |
| KEGG | Axon guidance / Sema3A–Plexin (guidance & repulsion) | 5.9×10⁻⁴ | Nrp1, Plxna1, Slit3, Pard3, Gnai3 |

Broader GO terms (developmental process p≈6×10⁻¹¹, anatomical-structure development, protein binding) are also top-ranked, consistent with memory regions sitting near tissue-identity/developmental genes.

## Why this matters for the cancer-priming hypothesis
Adherens junctions, nectin/cadherin adhesion, the PINCH–ILK–parvin (IPP) focal-adhesion complex, actin remodeling and cell-spreading are exactly the machinery that is **dismantled during EMT and epithelial-cancer invasion** — and **AP-1 (FOS/JUN) is the canonical transcriptional driver of that invasive program.** So the regions that "remember" accessibility, carry an AP-1 motif, are AP-1-bound, and are *already active enhancers* in keratinocytes are wired into the adhesion→cytoskeleton→migration axis. That is a concrete, testable basis for "AP-1 priming → loss of adhesion control → worse prognosis."

Driver genes worth flagging:
- **Col17a1** — hemidesmosome collagen; its loss/asymmetric partitioning governs **epidermal stem-cell aging** and marks **skin SCC** → bridges this set to the aging-atlas arm.
- **Cadm1, Nectin1, Pard3** — adhesion/polarity genes that act as **tumor suppressors**; their dysregulation drives carcinoma invasion.
- **Wasf2 (WAVE2), Lims1 (PINCH), Actn1, Mrtfa (MRTF-A), Mical2, Pkn2, Arhgef7** — actin/migration machinery; MRTF-A and MICAL2 are established **EMT/invasion** effectors.
- **Grhl1** (epidermal differentiation TF), **Cxcl1** (AP-1-driven inflammatory chemokine), **Fgfr2, Rarb, Znrf3** (epithelial growth / Wnt control).

## Caveats
- Enrichment used HOMER's default genome-wide background. Broad GO terms (development, "binding") are partly expected for any ~200-gene set; the **specific adhesion/junction/cytoskeleton signal is the robust, interpretable finding**. A matched background (all memory-near genes, or accessible epidermal genes) would tighten the broad terms — I can re-run that if you want.
- `msigdb`/`cosmic` (oncogenic gene sets) were empty in this HOMER mouse install, so direct cancer-gene-set enrichment isn't shown here; the human-ortholog → TCGA test is the way to get that directly.

## Files
`SetA_gene_symbols.txt` (202 symbols) · `SetA_top_pathways.tsv` · `seta_go/SetA.annotated.txt` (per-peak symbol annotation) · `seta_go/GO_refseq/` & `GO_symbols/` (full HOMER output, open `geneOntology.html`).
