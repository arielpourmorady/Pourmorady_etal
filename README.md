# Pourmorady et al - Code for Data Analysis
![OE](https://github.com/arielpourmorady/Pourmorady_etal/assets/70916908/eb65d636-d6cc-40b5-8402-c6de752cb008)

## Introduction
This repository is divided into several directories and subdirectories that contain scripts for analyzing data associated with the paper "RNA-mediated symmetry breaking enables singular olfactory receptor choice" (Pourmorady et al., 2023).

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
Fig. 2|B|dipc|dipc.ptII.enhancer.preferences
Fig. 2|C|dipc|dipc.ptIV.inactive.hub.extraction
Fig. 2|E|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Fig. 2|F|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Fig. 2|G|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Fig. 2|H|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Fig. 3|A|liquidhic|liquid.hic.heatmaps
Fig. 3|B|liquidhic|liquid.hic.heatmaps.replicates
Fig. 3|C|hichip|
Fig. 3|D|hichip|
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
Fig. 5|B|rna|rna.final
Fig. 5||hic|hic
Supp Fig. 1|A|multiome_P2nc|multiome_partI.input_processing
Supp Fig. 1|B|multiome_P2nc|multiome_partII.pseudotime_projection
Supp Fig. 1|C|multiome_P2nc|multiome_partII.pseudotime_projection
Supp Fig. 1|D|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|E|multiome_P2nc|multiome_partIV.accessibility_perGI_tileplot
Supp Fig. 1|F|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|G|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|H|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|I|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|J|multiome_P2nc|multiome_partIII.pseudotime_plots
Supp Fig. 1|K|multiome_P2nc|multiome_partI.input_processing
Supp Fig. 1|L|multiome_P2nc|multiome_partV.GI_Determinism
Supp Fig. 1|M|multiome_P2nc|multiome_partV.GI_Determinism
Supp Fig. 1|N|multiome_P2nc|multiome_partV.GI_Determinism
Supp. Fig. 2|A|rna|rna.final
Supp. Fig. 2|B|rna|rna.final
Supp. Fig. 2|C|rna|rna.final
Supp. Fig. 2|D|multiome|multiome_partI.input_processing
Supp. Fig. 2|E|multiome|multiome_partI.input_processing
Supp. Fig. 4|C|dipc|dipc.ptI.counting.enhancers.by.distance
Supp. Fig. 4|D|dipc|dipc.ptI.counting.enhancers.by.distance
Supp. Fig. 5|B|dipc|dipc.ptII.enhancer.preferences
Supp. Fig. 5|C|dipc|dipc.ptIII.GIHDeterminism
Supp. Fig. 5|D|dipc|dipc.ptV.hubs.per.cell
Supp. Fig. 5|E|dipc|dipc.ptIV.inactive.hub.extraction
Supp. Fig. 6|B|dipc|dipc.ptVI.EqualHubComparison.QC
Supp. Fig. 6|C|dipc|dipc.ptVI.EqualHubComparison.QC
Supp. Fig. 6|D|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Supp. Fig. 6|E|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Supp. Fig. 6|F|dipc|dipc.ptVII.GIH_ORcomp_contact_specificity
Supp. Fig. 7|B|liquidhic|liquid.hic.LOS
Supp. Fig. 7|C|liquidhic|liquid.hic.LOS
Supp. Fig. 7|E|hichip|
Supp. Fig. 7|F|hichip|
Supp. Fig. 7|G|hichip|
Supp. Fig. 7|H|hichip|
Supp. Fig. 8|A|hic|hic
Supp. Fig. 8|B|hic|hic
Supp. Fig. 8|C|hic|hic
Supp. Fig. 8|D|hic|hic
Supp. Fig. 8|E|dipc|dipc.ptVIII.inactive.hub.extraction.tan/dipc.ptIX.tanPCAandGIHsize
Supp. Fig. 8|F|dipc|dipc.ptVIII.inactive.hub.extraction.tan/dipc.ptIX.tanPCAandGIHsize
Supp. Fig. 9|A|hic|ompttatetop2.mor28icretdt.hic
Supp. Fig. 9|B|hic|ompttatetop2.mor28icretdt.hic
Supp. Fig. 10|B|rna|rna.final
Supp. Fig. 10|C|rna|rna.final
Supp. Fig. 10|D|hic|hic
Supp. Fig. 10|G|rna|rna.OMPM71nc_v_OMPtTAtetOGFP
Supp. Fig. 10|H|rna|rna.OMPM71nc_v_OMPGFP
Supp. Fig. 10|I|rna|rna.OMPM71nc_v_OMPGFPOMPtTA


## Input Files
[Dip-C]([url]([https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fdipc%2Etar%2Egz](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230380))): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file&file=GSE230380%5FPourmorady%5Fetal%5Fassociated%5Ffiles%5Fdipc%2Etar%2Egz

[Raw NGS Data]([url](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file)): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE230380&format=file
