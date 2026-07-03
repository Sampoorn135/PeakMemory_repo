# Pathway / GO Enrichment of the Shared Memory Genes (IMQ ∩ HFSC)

Genes shared between the two datasets' final memory peaks were tested for pathway / ontology enrichment with HOMER `findGO.pl` (mouse), which uses a genome-wide gene background and the hypergeometric test. **FDR** below is a Benjamini–Hochberg correction computed across all terms in each ontology; raw *p* is HOMER's hypergeometric p-value.

> Main set = **50_Single shared genes (356)** — identical to the pooled shared set. *k* = shared genes in the term.

## Enriched pathways — 50_Single shared set (356 genes)

### KEGG  (0 terms at FDR<0.05 of 229 tested)

| Term | k | raw p | FDR | Example genes |
| --- | ---: | ---: | ---: | --- |
| Axon guidance | 10 | 3.3e-04 | 6.1e-02 | Lrrc4, Bmp7, Bmpr1b, Efna5, Sema4a, Robo2, Ppp3ca, Pard3 |
| Sphingolipid signaling pathway | 8 | 5.3e-04 | 6.1e-02 | S1pr3, Sgms1, Cers6, Pik3r1, Ppp2r5c, Gnai3, Nsmaf, Asah2 |
| cAMP signaling pathway | 9 | 3.1e-03 | 1.6e-01 | Adrb2, Atp2b1, Ppp1cb, Pde4d, Pik3r1, Gnai3, Creb3l1, Fxyd2 |
| cGMP-PKG signaling pathway | 8 | 4.1e-03 | 1.6e-01 | Gnai3, Creb3l1, Fxyd2, Ppp3ca, Adrb2, Ppif, Atp2b1, Ppp1cb |
| Proteoglycans in cancer | 9 | 4.2e-03 | 1.6e-01 | Itgav, Ppp1cb, Pik3r1, Gab1, Iqgap1, Tiam1, Mmp2, Ezr |
| Adherens junction | 5 | 4.8e-03 | 1.6e-01 | Iqgap1, Lmo7, Map3k7, Nectin1, Pard3 |
| Regulation of actin cytoskeleton | 9 | 5.4e-03 | 1.6e-01 | Ezr, Tiam1, Diaph1, Iqgap1, Fgf7, Pik3r1, Ppp1cb, Itga6 |
| Amoebiasis | 6 | 5.6e-03 | 1.6e-01 | C8g, Serpinb3b, Lama3, Pik3r1, Hspb1, C8b |
| Adrenergic signaling in cardiomyocytes | 7 | 7.2e-03 | 1.8e-01 | Adrb2, Ppp1cb, Atp2b1, Fxyd2, Gnai3, Ppp2r5c, Creb3l1 |
| Small cell lung cancer | 5 | 8.7e-03 | 2.0e-01 | Itga6, Rarb, Itgav, Pik3r1, Lama3 |

### Reactome  (4 terms at FDR<0.05 of 505 tested)

| Term | k | raw p | FDR | Example genes |
| --- | ---: | ---: | ---: | --- |
| Cell junction organization | 7 | 1.1e-04 | 3.2e-02 | Lama3, Itga6, Cadm1, Cldn1, Lims1, Pard3, Nectin1 |
| Molecules associated with elastic fibres | 5 | 1.3e-04 | 3.2e-02 | Bmp7, Fbln2, Ltbp1, Ltbp3, Itgav |
| Cell-Cell communication | 8 | 2.3e-04 | 3.9e-02 | Cadm1, Lama3, Itga6, Nectin1, Pard3, Cldn1, Lims1, Iqgap1 |
| Elastic fibre formation | 5 | 3.6e-04 | 4.5e-02 | Bmp7, Fbln2, Ltbp1, Ltbp3, Itgav |
| MET activates PI3K/AKT signaling | 2 | 4.0e-03 | 3.4e-01 | Gab1, Pik3r1 |
| Nectin/Necl  trans heterodimerization | 2 | 4.0e-03 | 3.4e-01 | Cadm1, Nectin1 |
| Laminin interactions | 3 | 5.0e-03 | 3.4e-01 | Col18a1, Itga6, Lama3 |
| Regulation of gene expression by Hypoxia-inducible Factor | 2 | 5.6e-03 | 3.4e-01 | Cited2, Hif1a |
| PI-3K cascade:FGFR2 | 3 | 6.5e-03 | 3.4e-01 | Fgf7, Gab1, Pik3r1 |
| Cell-cell junction organization | 4 | 6.7e-03 | 3.4e-01 | Cadm1, Nectin1, Pard3, Cldn1 |

### WikiPathways  (2 terms at FDR<0.05 of 88 tested)

| Term | k | raw p | FDR | Example genes |
| --- | ---: | ---: | ---: | --- |
| XPodNet - protein-protein interactions in the podocyte expanded by STRING | 35 | 1.4e-04 | 1.3e-02 | Bmp7, Ppp3ca, Stt3b, Htra1, Palld, Map4k1, Diaph1, Cldn1 |
| PodNet: protein-protein interactions in the podocyte | 17 | 6.9e-04 | 3.0e-02 | Ezr, Arhgap24, Hif1a, Lims1, Robo2, Pik3r1, Utrn, Igsf5 |
| Signal Transduction of S1P Receptor | 3 | 1.3e-02 | 3.8e-01 | Asah2, S1pr3, Gnai3 |
| IL-3 Signaling Pathway | 6 | 2.6e-02 | 5.7e-01 | Pik3r1, Mmp2, Gab1, Hspb1, Stat6, Syk |
| Lung fibrosis | 4 | 5.0e-02 | 6.0e-01 | Bmp7, Fgf7, Mmp2, Atp11a |
| Nuclear Receptors | 3 | 5.5e-02 | 6.0e-01 | Rarb, Rora, Nr4a2 |
| FAS pathway and Stress induction of HSP regulation | 3 | 5.5e-02 | 6.0e-01 | Hspb1, Mapkapk3, Map3k7 |
| G13 Signaling Pathway | 3 | 5.8e-02 | 6.0e-01 | Diaph1, Iqgap1, Ppp1cb |
| Splicing factor NOVA regulated synaptic proteins | 3 | 7.0e-02 | 6.0e-01 | Efna5, Atp2b1, Cadm1 |
| IL-6 signaling Pathway | 5 | 7.5e-02 | 6.0e-01 | Map3k7, Ppp2r5c, Hspb1, Pik3r1, Gab1 |

### GO Biological Process  (563 terms at FDR<0.05 of 4604 tested)

| Term | k | raw p | FDR | Example genes |
| --- | ---: | ---: | ---: | --- |
| biological regulation | 232 | 3.8e-18 | 1.7e-14 | Gm14137, Fam3b, C8b, Zfp169, Rbpjl, Nfil3, Zbtb7a, Sptssb |
| regulation of biological process | 224 | 6.4e-17 | 1.5e-13 | Sgms1, Nrxn2, Bcl9, Cldn1, Tead4, Trim71, Mak, Scgb1a1 |
| cellular process | 259 | 1.1e-16 | 1.7e-13 | Tead4, Bcl9, Cldn1, Nrxn2, Sgms1, Mak, Trim71, Emp1 |
| positive regulation of biological process | 146 | 1.7e-16 | 1.9e-13 | Golga2, Lmo7, Cited2, Fgf7, Snx18, Pkp4, Efna5, Mapre2 |
| multicellular organism development | 109 | 4.6e-15 | 3.9e-12 | Sh3pxd2a, Cnot2, Pard3, Robo2, Ltbp1, Cryaa, Ggnbp1, Syk |
| positive regulation of cellular process | 134 | 5.1e-15 | 3.9e-12 | Atxn7, Klf6, Tank, Camta1, Gab1, Lpar3, Bmp7, Tiam1 |
| regulation of cellular process | 211 | 6.9e-15 | 4.5e-12 | Zfp385a, Mtss1, Oas1f, Rbm20, Usp47, Mir21c, Sox8, Pde4d |
| regulation of multicellular organismal process | 89 | 1.2e-14 | 6.7e-12 | Fgf7, Efna5, Golga2, Cited2, Adrb2, Nfkbiz, Klf4, Mapre2 |
| anatomical structure development | 126 | 2.3e-14 | 1.2e-11 | Trim71, Hif1a, Dnajb6, Nectin1, Nrxn2, Bcl9, Rubie, Cldn1 |
| cell differentiation | 99 | 3.8e-14 | 1.6e-11 | Palld, Kidins220, Flrt3, Bmpr1b, Ezr, Pde4d, Arid5b, Tenm4 |

### GO Molecular Function  (42 terms at FDR<0.05 of 831 tested)

| Term | k | raw p | FDR | Example genes |
| --- | ---: | ---: | ---: | --- |
| protein binding | 195 | 1.5e-20 | 9.0e-18 | Tmed10, Snx27, Rarb, Tnks1bp1, Ccnl1, Lims1, Ppp2r5c, Cpeb2 |
| binding | 249 | 2.2e-20 | 9.0e-18 | Col18a1, Ltbp1, Pheta1, Ggnbp1, Kctd3, Klf3, Rbms2, Gpr155 |
| protein domain specific binding | 27 | 2.2e-06 | 4.8e-04 | Trak1, Tmem88, Adgrl2, Igsf5, Cited2, Skp1, Syk, Arhgap29 |
| cation binding | 75 | 2.3e-06 | 4.8e-04 | Rnf144b, Lims1, Kalrn, Map3k20, Nudt3, Trp63, Cdc7, Mmp2 |
| metal ion binding | 73 | 3.0e-06 | 5.0e-04 | P3h2, Fah, Asah2, Gmpr, Ece1, Ocm, Zfp90, Hivep2 |
| molecular function regulator activity | 49 | 3.9e-06 | 5.2e-04 | Il9, Nsmaf, Arhgap29, Il33, Trib2, Efna5, Gdf11, Htra1 |
| enzyme binding | 53 | 4.4e-06 | 5.2e-04 | Tank, Tiam1, Pstpip1, Iqgap1, Stat6, Nkain1, Fmnl2, Map3k7 |
| small molecule binding | 99 | 8.4e-06 | 8.7e-04 | Stt3b, Rarb, Lims1, Trp63, Adrb2, Prickle2, Tut7, Cd93 |
| transcription coregulator binding | 10 | 1.0e-05 | 9.7e-04 | Rora, Klf4, Hif1a, Vgll4, Chd6, Hnf1b, Stat6, Map3k7 |
| PDZ domain binding | 9 | 1.5e-05 | 1.1e-03 | Atp2b1, Igsf5, Adgrl2, Slc22a12, Synj2, Kidins220, Tmem88, Cadm1 |

## How enrichment holds across stringency

Number of significantly enriched terms (FDR < 0.05) in each shared set:

| Ontology | 50_Single (356) | 50_Reciprocal (265) | 75_Single (94) | 75_Reciprocal (41) |
| --- | ---: | ---: | ---: | ---: |
| KEGG | 0 | 0 | 7 | 16 |
| Reactome | 4 | 2 | 0 | 0 |
| WikiPathways | 2 | 2 | 0 | 2 |
| GO Biological Process | 563 | 395 | 145 | 245 |
| GO Molecular Function | 42 | 41 | 2 | 0 |

## Top pathway per condition (most significant KEGG / Reactome)

| Condition | Top KEGG (k, FDR) | Top Reactome (k, FDR) |
| --- | --- | --- |
| 50_Single | Axon guidance (10, FDR 6.1e-02) | Cell junction organization (7, FDR 3.2e-02) |
| 50_Reciprocal | Axon guidance (8, FDR 9.0e-02) | Molecules associated with elastic fibres (5, FDR 1.3e-02) |
| 75_Single | Axon guidance (8, FDR 2.0e-05) | Role of phospholipids in phagocytosis (2, FDR 2.8e-01) |
| 75_Reciprocal | Osteoclast differentiation (4, FDR 5.2e-03) | Role of phospholipids in phagocytosis (2, FDR 8.4e-02) |

## Files

- `pathways/<condition>/` — full HOMER findGO output for every ontology (`biological_process.txt`, `kegg.txt`, `reactome.txt`, `wikipathways.txt`, `molecular_function.txt`, `msigdb.txt`, `geneOntology.html`, …).
