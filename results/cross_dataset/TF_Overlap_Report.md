# Transcription-Factor (Motif) Overlap Across Memory Datasets

HOMER known-motif enrichment on each dataset's `StrictMemory50_Single` regions (`findMotifsGenome.pl -size 200`, vertebrate motifs, GC-matched background). **Enriched TF = Benjamini q < 0.05**; motifs collapsed to the TF name.

## Enriched TFs per dataset

| Dataset | # enriched TFs (q<0.05) |
| --- | ---: |
| IMQ | 149 |
| HFSC | 141 |
| distil | 167 |

## Pairwise overlap

| Pair | Shared TFs | Jaccard | Shared list |
| --- | ---: | ---: | --- |
| IMQ ∩ HFSC | 102 | 0.54 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, EHF, EKLF, ELF3, ELF5, ERG, ETS1, ETV1, ETV4, Emx2, En1, FOXK1, FOXK2, Fli1, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GABPA, GRHL2, HEB, HLF, HOXA3, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, Lhx1, Lhx3, Lhx6, MafA, MafB, MafF, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFE2L2, NFIL3, NFkB-p65-Rel, Npas4, Nrf2, Oct11, Oct2, Oct4, Oct6, PBX2, Pdx1, Rfx6, STAT4, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, Zic3, c-Jun-CRE, p53, p63, p73 |
| IMQ ∩ distil | 103 | 0.48 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, E2F3, EBF1, EKLF, ERG, Egr1, Egr2, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, Fox:Ebox, FoxL2, Foxa2, Foxa3, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GLIS3, GRHL2, HIC1, HLF, HOXA3, Hoxc10, Hoxc9, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, Lhx1, Lhx3, MafA, MafB, MafK, Maz, NF-E2, NF1, NF1-halfsite, NFAT, NFAT:AP1, NFE2L2, NFIL3, Nkx6.1, Nrf2, Oct11, Oct6, PBX2, Pdx1, RUNX-AML, Rfx6, SCL, STAT4, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD1, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, c-Jun-CRE, p53, p63, p73 |
| HFSC ∩ distil | 105 | 0.52 | AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CEBP:AP1, CREB5, CTCF, CTCF-SatelliteElement, Chop, DLX1, DLX2, DLX5, Dlx3, EKLF, ERG, EWS:ERG-fusion, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GATA, GRHL2, Gsx2, HLF, HOXA3, Isl1, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, LHX9, Lhx1, Lhx2, Lhx3, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFATC2, NFE2L2, NFIL3, NFkB-p65, NFkB2-p52, Nrf2, Oct11, Oct6, PBX2, Pdx1, Rbpj1, Rfx6, STAT4, Smad2, Smad3, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, ZEB1, ZNF416, c-Jun-CRE, p53, p63, p73 |

## All three together (IMQ ∩ HFSC ∩ distil)

**86 TFs enriched in all three memory datasets:**

AP-1, AP-2alpha, AP-2gamma, Atf1, Atf2, Atf3, Atf4, Atf7, BATF, BORIS, Bach1, Bach2, Brn1, CEBP, CREB5, CTCF, CTCF-SatelliteElement, Chop, EKLF, ERG, Emx2, En1, FOXK1, FOXK2, Fos, Fosl2, FoxL2, Foxf1, Foxo1, Foxo3, Fra1, Fra2, GRHL2, HLF, HOXA3, Jun-AP1, JunB, JunD, KLF1, KLF10, KLF14, KLF17, KLF3, KLF5, KLF6, Klf15, Klf4, Klf9, Lhx1, Lhx3, MafA, MafB, MafK, NF-E2, NF1, NF1-halfsite, NFAT:AP1, NFE2L2, NFIL3, Nrf2, Oct11, Oct6, PBX2, Pdx1, Rfx6, STAT4, Sox10, Sox6, Sp1, Sp2, Sp5, Stat3, Stat3+il21, TEAD, TEAD2, TEAD3, TEAD4, Tcf12, Tcfcp2l1, Tgif2, Tlx?, Unknown-ESC-element, c-Jun-CRE, p53, p63, p73

## Dataset-specific (enriched in only one)

- **IMQ-only (30)**: Bcl11a, Bcl6, ETS:RUNX, FOXA1, GATA3, Gata1, Gata2, Gata4, Gata6, HIF-1a, HIF2a, Hoxa10, Hoxd11, Hoxd9, MITF, Meis1, OCT:OCT-short, PAX5, Ptf1a, REST-NRSF, Snail1, THRb, TRPS1, Tbet, Tbr1, Tbx5, Tgif1, ZNF341, ZNF519, Zic2
- **HFSC-only (20)**: Ascl1, CRE, ELF1, ETS, EWS:FLI1-fusion, Elf4, Elk1, Etv2, FOXP1, RFX, RORg, RORgt, Rfx2, SPDEF, STAT1, STAT5, Tbx20, USF1, Usf2, ZEB2
- **distil-only (45)**: AR-halfsite, CDX4, Cdx2, EBF2, Elk4, GRE, HOXB13, HRE, Hoxa11, Hoxa13, Hoxb4, Hoxc13, Hoxc6, Hoxd13, IRF1, IRF2, IRF8, ISRE, MYRF, Mef2a, Mef2b, Mef2c, Mef2d, NF1:FOXA1, PGR, SOX1, STAT6, Six1, Six2, Six4, Smad4, Sox15, Sox17, Sox2, Sox21, Sox3, Sox4, Sox7, Sox9, Tcf21, Unknown, WT1, ZNF189, ZNF322, Zfp809
