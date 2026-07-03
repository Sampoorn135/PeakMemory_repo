# Resume: Precancer → Cancer progression analysis (Set A)

## Goal
Map mouse epidermal-SC memory enhancers (Set A: AP1-primed + H3K27ac, mm10) onto the
human normal → actinic keratosis (AK) → cutaneous SCC (cSCC) progression, and test
whether these memory enhancers/genes open or close along progression.

## Input (already in this folder)
- `SetA_AP1primed_H3K27ac.bed` — 209 mm10 regions (col4 = H3K27ac signal)
- `SetA_gene_symbols.txt` — 202 mouse gene symbols

## DONE
- **Gene ortholog mapping (mouse→human):** complete.
  - `SetA_mouse2human_orthologs.tsv` — per-gene mapping (homology type, confidence, source)
  - `SetA_human_signature.txt` — **183 human genes** (161 one2one + 5 one2many + 10 HCOP-recovered)
  - 25 mouse genes have no human ortholog — all are predicted Rik/Gm loci (18), miRNAs (6),
    1 pseudogene. No protein-coding gene lost. Tool: `mousipy` (bundles HCOP + BioMart tables).

## TODO (needs internet — set sandbox allowlist to "All domains" BEFORE starting the session)
1. **Peak liftOver mm10 → hg38**
   - Chain: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz
   - Tool: `pyliftover` (already pip-installable). Keep uniquely-mapped regions; report drop rate.
2. **Chromatin track — GSE277274 (ATAC-seq, 13 samples: normal/AK/cSCC, human, bigWig)**
   - Tar: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE277274&format=file  (~1 GB, BW)
   - Quantify mean ATAC signal of each lifted hg38 region per sample (pyBigWig).
   - Build region×sample matrix → test normal vs AK vs cSCC trend (Jonckheere-Terpstra /
     Kruskal–Wallis), per-region and as an aggregate memory-enhancer score.
3. **Gene track — test the 183-gene signature across states**
   - RNA option: m6A "input" samples (GSE277275) or external AK/cSCC expression (E-GEOD-45216).
   - ssGSEA/GSVA score per sample → trend across normal/AK/cSCC.

## Sibling accessions (same SuperSeries, PRJNA1161632)
- GSE277274 = ATAC-seq (USE THIS)  | GSE277275 = m6A-seq | GSE277273 = 850K methylation

## Caveats to keep in mind
- Bulk ATAC = tumor+stroma+immune mix; a "closing" peak may be composition change (anchor in scATAC if needed).
- Small n per state (5 AK / 5 cSCC / 10 normal-ish) → treat trend as exploratory.
- Cross-species liftOver loses non-conserved enhancers; report how many of 209 survive.
- Confirm bigWig genome build = hg38 before overlaying lifted coords.
