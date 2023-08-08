# Pourmorady et al - Code for Data Analysis
![OE](https://github.com/arielpourmorady/Pourmorady_etal/assets/70916908/eb65d636-d6cc-40b5-8402-c6de752cb008)

## Introduction
This repository is divided into several directories and subdirectories that contain scripts for analyzing data associated with the paper "Olfactory receptor mRNAs act as “selfish” non-coding RNAs that enforce transcriptional singularity" (Pourmorady et al., 202X).

All analysis, with the exception of pymol generated models, was performed in RStudio.

Directories are divided by methodology and may contain subdirectories that are specific to certain panels within figures. 

To independently analyze this data, we provide all input and annotation files which are attached to the GEO superseries generated for this paper.

Figure|Panel|Parent Directory|Folder
---|---|---|---
Fig. 1|A|multiome|multiome_partI.input_processing
Fig. 1|B|multiome|multiome_partII.pseudotime_projection
Fig. 1|C|multiome|multiome_partII.pseudotime_projection
Fig. 1|D|multiome|multiome_partIII.pseudotime_plots
Fig. 1|E|multiome|multiome_partIV.accessibility_perGI_tileplot
Fig. 1|F|multiome|multiome_partIII.pseudotime_plots
Fig. 1|G|multiome|multiome_partIII.pseudotime_plots
Fig. 1|H|multiome|multiome_partIII.pseudotime_plots
Fig. 1|I|multiome|multiome_partIII.pseudotime_plots
Fig. 1|J|multiome|multiome_partIII.pseudotime_plots
Fig. 1|K|multiome|multiome_partI.input_processing
Fig. 1|L|multiome|multiome_partV.GI_Determinism
Fig. 1|M|multiome|multiome_partV.GI_Determinism
Fig. 1|N|multiome|multiome_partV.GI_Determinism
Fig. 2|C|dipc|dipc.ptI.enhancer.preferences.revised
Fig. 2|D|dipc|dipc.ptI.enhancer.preferences.revised
Fig. 2|E|dipc|dipc.ptII.GIHDeterminism
Fig. 2|F|dipc|dipc.ptIII.inactive.hub.extraction
Fig. 2|G|dipc|dipc.ptIV.hubs.per.cell
Fig. 2|I|dipc|dipc.ptV.GIH_ORcomp_contact_specificity
Fig. 2|J|dipc|dipc.ptV.GIH_ORcomp_contact_specificity
Fig. 2|K|dipc|dipc.ptV.GIH_ORcomp_contact_specificity
Fig. 2|L|dipc|dipc.ptV.GIH_ORcomp_contact_specificity
Fig. 3|A|liquidhic|liquid.hic.heatmaps
Fig. 3|B|liquidhic|liquid.hic.LOS
Fig. 3|C|hichip|
Fig. 3|D|hichip|
Fig. 3|E|hichip|
Fig. 3|F|hichip|
Fig. 4|B|rna|rna.final
Fig. 4||hic|hic
Fig. 4|C|rna|rna.final
Fig. 4||hic|hic
Fig. 4|H|rna|ompttatetop2.mor28icretdt.rna
Fig. 4|I|hic|ompttatetop2.mor28icretdt.hic
Fig. 4|J|hic|ompttatetop2.mor28icretdt.hic
Fig. 4|K|hic|ompttatetop2.mor28icretdt.hic
Fig. 4|L|hic|ompttatetop2.mor28icretdt.hic
Fig. 4|M|hic|ompttatetop2.mor28icretdt.hic
Fig. 4|N|hic|ompttatetop2.mor28icretdt.hic
Fig. 5|A|rna|rna.final
Fig. 5||hic|hic
Fig. 5|C|rna|rna.final
Fig. 5|D|rna|rna.final
Fig. 5||hic|hic
Supp. Fig. 1|A|rna|rna.final
Supp. Fig. 1|B|rna|rna.final
Supp. Fig. 1|C|rna|rna.final
Supp. Fig. 1|D|multiome|multiome_partI.input_processing
Supp. Fig. 1|E|multiome|multiome_partI.input_processing
Supp. Fig. 3|C|dipc|dipc.ptVI.counting.enhancers.by.distance
Supp. Fig. 3|D|dipc|dipc.ptVI.counting.enhancers.by.distance
Supp. Fig. 3|E|dipc|dipc.ptIII.inactive.hub.extraction
Supp. Fig. 3|F|dipc|dipc.ptIII.inactive.hub.extraction
Supp. Fig. 4|B|dipc|dipc.ptIII.inactive.hub.extraction
Supp. Fig. 4|C|dipc|dipc.ptIII.inactive.hub.extraction
Supp. Fig. 4|D|dipc|dipc.ptV.GIH_ORcomp_contact_specificity
Supp. Fig. 5|B|liquidhic|liquid.hic.LOS
Supp. Fig. 5|C|liquidhic|liquid.hic.LOS
Supp. Fig. 6|A|hic|hic
Supp. Fig. 6|B|hic|hic
Supp. Fig. 6|C|hic|hic
Supp. Fig. 6|D|hic|hic
Supp. Fig. 6|E|dipc|dipc.ptVII.inactive.hub.extraction.tan
Supp. Fig. 6|F|dipc|dipc.ptVII.inactive.hub.extraction.tan/dipc.ptVIII.tanPCAandGIHsize
Supp. Fig. 7||hic|ompttatetop2.mor28icretdt.hic
Supp. Fig. 8|B|rna|rna.final
Supp. Fig. 8|C|rna|rna.final
Supp. Fig. 8|D|hic|hic

## Input Files
[Dip-C]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fdipc%2Etar%2Egz)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fdipc%2Etar%2Egz

[Hi-C]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fhic%2Etar%2Egz)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fhic%2Etar%2Egz

[HiChIP]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fhichip%2Etar%2Egz)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fhichip%2Etar%2Egz 

[Liquid Hi-C]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fliquidhic%2Etar%2Egz)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fliquidhic%2Etar%2Egz

[Raw NGS Data]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file
