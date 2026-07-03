# Transcription-Factor (Motif) Overlap Across Memory Datasets

HOMER known-motif enrichment on each dataset's `StrictMemory50_Single` regions (`findMotifsGenome.pl -size 200`, vertebrate motifs, GC-matched background). **Enriched TF = Benjamini q < 0.05**; motifs collapsed to the TF name.

## Enriched TFs per dataset

| Dataset | # enriched TFs (q<0.05) |
| --- | ---: |
| IMQ | 149 |
| HFSC | 141 |
| distil | 167 |
| pancreas | 174 |

## Pairwise overlap

| Pair | Shared TFs | Jaccard | Shared list |
| --- | ---: | ---: | --- |
| IMQ ∩ HFSC | 102 | 0.54 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, EHF, EKLF, ELF3, ELF5, ERG, ETS1, ETV1, ETV4, Emx2, En1, FOXK1, FOXK2, Fli1, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GABPA, GRHL2, HEB, HLF, HOXA3, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, Lhx1, Lhx3, Lhx6, MafA, MafB, MafF, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFE2L2, NFIL3, NFkB-p65-Rel, Npas4, Nrf2, Oct11, Oct2, Oct4, Oct6, PBX2, Pdx1, Rfx6, STAT4, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, Zic3, c-Jun-CRE, p53, p63, p73 |
| IMQ ∩ distil | 103 | 0.48 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, E2F3, EBF1, EKLF, ERG, Egr1, Egr2, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, Fox:Ebox, FoxL2, Foxa2, Foxa3, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GLIS3, GRHL2, HIC1, HLF, HOXA3, Hoxc10, Hoxc9, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, Lhx1, Lhx3, MafA, MafB, MafK, Maz, NF-E2, NF1, NF1-halfsite, NFAT, NFAT:AP1, NFE2L2, NFIL3, Nkx6.1, Nrf2, Oct11, Oct6, PBX2, Pdx1, RUNX-AML, Rfx6, SCL, STAT4, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD1, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, c-Jun-CRE, p53, p63, p73 |
| IMQ ∩ pancreas | 90 | 0.39 | AP-1, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Bcl11a, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, EBF1, EKLF, Emx2, En1, FOXA1, FOXK1, FOXK2, Fos, Fosl2, Fox:Ebox, FoxL2, Foxa2, Foxa3, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GATA3, GLIS3, Gata4, Gata6, HLF, HOXA3, Jun-AP1, JunB, JunD, KLF1, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, MITF, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT, NFAT:AP1, NFE2L2, NFIL3, Nkx6.1, Nrf2, OCT:OCT-short, PBX2, Pdx1, RUNX-AML, Rfx6, SCL, STAT4, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD1, TEAD2, TEAD3, TEAD4, THRb, TRPS1, Tcf12, Tgif2, Tlx?, Unknown-ESC-element, ZNF519, Zic2, Zic3, c-Jun-CRE |
| HFSC ∩ distil | 105 | 0.52 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CEBP:AP1, CREB5, CTCF, CTCF-SatelliteElement, Chop, DLX1, DLX2, DLX5, Dlx3, EKLF, ERG, EWS:ERG-fusion, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GATA, GRHL2, Gsx2, HLF, HOXA3, Isl1, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, LHX9, Lhx1, Lhx2, Lhx3, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFATC2, NFE2L2, NFIL3, NFkB-p65, NFkB2-p52, Nrf2, Oct11, Oct6, PBX2, Pdx1, Rbpj1, Rfx6, STAT4, Smad2, Smad3, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, ZEB1, ZNF416, c-Jun-CRE, p53, p63, p73 |
| HFSC ∩ pancreas | 91 | 0.41 | AP-1, Ascl1, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, CEBP, CEBP:AP1, CRE, CREB5, CTCF, CTCF-SatelliteElement, Chop, DLX1, DLX2, DLX5, Dlx3, EKLF, Emx2, En1, FOXK1, FOXK2, FOXP1, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GATA, Gsx2, HLF, HOXA3, Isl1, Jun-AP1, JunB, JunD, KLF1, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, LHX9, Lhx2, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFATC2, NFE2L2, NFIL3, NFkB-p65, Nrf2, PBX2, Pdx1, RFX, RORg, RORgt, Rfx2, Rfx6, STAT1, STAT4, STAT5, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tgif2, Tlx?, USF1, Unknown-ESC-element, Zic3, c-Jun-CRE |
| distil ∩ pancreas | 110 | 0.48 | AP-1, AR-halfsite, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, CEBP, CEBP:AP1, CREB5, CTCF, CTCF-SatelliteElement, Chop, DLX1, DLX2, DLX5, Dlx3, EBF1, EBF2, EKLF, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, Fox:Ebox, FoxL2, Foxa2, Foxa3, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GATA, GLIS3, GRE, Gsx2, HLF, HOXA3, HRE, Hoxb4, IRF2, IRF8, ISRE, Isl1, Jun-AP1, JunB, JunD, KLF1, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, LHX9, Lhx2, MafA, MafB, MafK, Mef2a, Mef2b, Mef2c, Mef2d, NF-E2, NF1, NF1-halfsite, NF1:FOXA1, NFAT, NFAT:AP1, NFATC2, NFE2L2, NFIL3, NFkB-p65, Nkx6.1, Nrf2, PBX2, PGR, Pdx1, RUNX-AML, Rfx6, SCL, STAT4, STAT6, Six1, Six2, Six4, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD1, TEAD2, TEAD3, TEAD4, Tcf12, Tcf21, Tgif2, Tlx?, Unknown-ESC-element, WT1, c-Jun-CRE |

## All four together (IMQ ∩ HFSC ∩ distil ∩ pancreas)

**68 TFs enriched in all four memory datasets:**

AP-1, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, EKLF, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, HLF, HOXA3, Jun-AP1, JunB, JunD, KLF1, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFE2L2, NFIL3, Nrf2, PBX2, Pdx1, Rfx6, STAT4, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tgif2, Tlx?, Unknown-ESC-element, c-Jun-CRE

## Dataset-specific (enriched in only one)

- **IMQ-only (19)**: Bcl6, ETS:RUNX, Gata1, Gata2, HIF-1a, HIF2a, Hoxa10, Hoxd11, Hoxd9, Meis1, PAX5, Ptf1a, REST-NRSF, Snail1, Tbet, Tbr1, Tbx5, Tgif1, ZNF341
- **HFSC-only (10)**: ELF1, ETS, EWS:FLI1-fusion, Elf4, Elk1, Etv2, SPDEF, Tbx20, Usf2, ZEB2
- **distil-only (25)**: CDX4, Cdx2, Elk4, HOXB13, Hoxa11, Hoxa13, Hoxc13, Hoxc6, Hoxd13, IRF1, MYRF, SOX1, Smad4, Sox15, Sox17, Sox2, Sox21, Sox3, Sox4, Sox7, Sox9, Unknown, ZNF189, ZNF322, Zfp809
- **pancreas-only (42)**: ARE, Ap4, Atoh1, Atoh7, BHLHA15, BMYB, Brn2, CEBP:CEBP, COUP-TFII, EAR2, EBF, Erra, Hand2, LXRE, MRE, MYB, Mesp1, Myf5, MyoD, MyoG, NeuroD1, NeuroG2, Olig2, PR, Pitx1:Ebox, RORa, RUNX, RUNX1, RUNX2, Reverb, Rfx1, Srebp1a, TCF4, THRa, Twist, Twist2, X-box, ZBTB18, ZNF91, Zac1, Zic, Znf263
