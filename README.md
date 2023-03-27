# Pourmorady_etal

1. Figure 1 - multiome

2. Figure 2
 - 2A/B - pymol
 - 2C - dipc.enhancer_preferences
 - 2D - dipc.hubs.per.cell_2D
 - 2E - dipc.inactive.hub.extraction
 - 2F - dipc.GIH_ORcomp_contact_specificity
 - 2G - figure_2G.dipc.OR_to_GIH.contact_specificity
 
3. Figure 3
 - 3A - liquid.hic.heatmaps
 - 3B - liquid.hic.heatmaps.replicates
 - 3C - hichip
 
4. Figure 5
 - RNAseq - rna.knockouts
 - HiC - hic.p2.ko
 
Data Analysis & Calculations

Multiome Experiments
Data processing with Signac
Signac 1.9.0 developed by the Satija Lab was used for all analyses which was performed in R 4.1.3.

Data input and quality control
All data was loaded and processed following the vignette described by the Satija Lab online (https://stuartlab.org/signac/articles/pbmc_multiomic.html). Briefly, outputs from cellranger-arc count, filtered_feature_bc_matrix.h5 (RNA) and atac_fragments.tsv.gz (ATAC), were loaded and combined into a Seurat object, and Chromatin Assay, respectively. The multiome was quality controlled by retaining cells with a TSS enrichment > 1, nucleosome signal < 2, ATAC counts between 1,000 and 10,000, and RNA counts between 1,000 and 25,000. Peaks were called using macs2 on the entire multiome, using the fragment file as the input to the peak calling algorithm, CallPeaks, with the following arguments: -g 2.7e9 -f BED –nomodel –extsize 200 –shift -100, and removing peaks overlapping annotated genomic blacklist regions for the mm10 genome.

Gene expression and DNA accessibility processing, and cell annotation
Gene expression UMI count data was normalized using SCTransform and we performed PCA using the RunPCA function in Seurat. Dimension reduction was performed on the “peaks” Chromatin Assay using LSI, with the following functions FindTopFeatures (arguments: min.cutoff = 5), RunTFIDF, and RunSVD. Joint UMAP was performed using the FindMultiModalNeighbors function to calculate closest neighbors using a weighted combination of RNA and ATAC similarities. Clusters were identified using weighted nearest neighbors data and were annotated by their expression of known marker genes for cells in the olfactory neuronal lineage [17, 18].

mOSN cCRE Identification
mOSN marker genes were identified by using the FindMarkers function on RNA data from the mOSN cluster. Only genes detected in > 25% of cells were extracted. To identify cCREs > 2000 bp but < 500,000 bp away from the TSS of top mOSN marker genes we performed peak-to-gene linkage using the LinkPeaks function with the following arguments: peak.assay = "peaks", expression.assay = "SCT", distance = 5e+05, min.distance = 2000, min.cells = 10, n_sample = 200, pvalue_cutoff = 0.05, score_cutoff = 0.2.

scATAC tile plots & percent accessibility per enhancer
The TilePlot function was used to generate tile plots of the number of fragments within a window surrounding Greek Islands and a cCRE. The FeatureMatrix function was used to extract fragments within a bed file of Greek Islands and mOSN cCREs. The percent accessibility per enhancer was measured by calculating the frequency with which all mOSNs had at least one fragment (Tn5 insertion event) within each enhancer.

Pseudotime projection with Monocle3 & verification with expression of neuronal marker genes
We converted our MOE multiome Seurat object to a CellDataSet using as.cell_data_set( ), and projected it into low-dimensional space via UMAP. We ordered cells and chose the beginning-most node in the GBC cluster as the root node for trajectory generation. We then verified the accuracy of the trajectory by examining the SCT normalized expression of known marker genes: GBC – Ascl1, INP1/2 – Neurod1, INP3 – Tex15, iOSN – Gap43, and mOSN – Carns1.

Enhancer accessibility and max OR transcript levels over pseudotime
The FeatureMatrix function was used to extract fragment counts within bed files of three enhancer groups: Greek Islands, mOSN cCREs, and Lhx2 and EBF1 bound peaks in mOSNs identified in Monahan et al. 2017 [11], and were normalized to the total number of ATAC counts per cell. Total enhancer accessibility per cell was extracted by measuring the total number of normalized fragments per enhancer (Fig. 1F) or for all enhancers (Fig. 1E, 1G-H) within an enhancer group, which was then averaged for all cells sharing the same pseudotime rounded to the nearest integer. To measure the accessibility of only active GIs, we measured the average normalized fragment count of only active GIs per cell, which was then plotted across pseudotime.


 
Dip-C Experiments
Hierarchical clustering of Greek Island spatial relationships
Refer to Fig. S2A for schematic. To identify active and inactive enhancer complexes, the function ‘dip-c pd’ was used to generate symmetric matrices, measured in particle radii (p.r.), of the pairwise distance relationships between all Greek Islands as well as the active OR allele. We scaled distance matrices and performed hierarchical clustering based off the euclidean relationships between all elements in a matrix using the hclust function from the R package “dendextend.” We cut dendrograms at a height of 2.5 p.r. using the cutree function, to extract all branches with an average radius of the same distance, which was verified in Fig. S2H. 

Identification of Active GIH, Inactive GIHs, Inactive GIH most similar to Active GIH, and OR genes within GIHs
The active GIH was identified as the branch containing the active OR allele by hierarchical clustering (see above). Inactive GIHs were identified as branches containing greater than the average size of the active GIH (> 5.59 GIs, Fig. 2D top), or within ± 1SD the size of the average active GIH (2.71 – 8.47 GIs, Fig. 2D bottom). For CSS comparison of active and inactive hubs, only the inactive GIH in every cell that had the most similar concentration of GIs as the active GIH in the same cell was used for CSS analysis. This inactive GIH was identified by counting the number of GIs in a cell’s active GIH, and finding another GIH (dendrogram branch) containing the most similar number of GIs (Fig. S2A5-9). The similar size and density of active GIH and the most similar inactive GIH was verified in Fig. S2B-C. OR genes contained within GIHs were identified by using the dip-c pd function to extract all ORs within 5 p.r. of all GIs participating in a GIH.

Contact specificity score in Dip-C
Haplotype imputed single-cell contacts were extracted as described in [14] and were binned at 50-kb resolution and normalized to contacts per million. To compute single-cell CSS, contacts within a 1-Mb radius surrounding two regions of interest were computed, between all regions of interest in a cell, and summed across cells, generating a 2-Mb-by-2-Mb matrix of contacts. Contact specificity scores were generated by normalizing the contacts within each 50-kb-by-50-kb bin to the total number of contacts in the matrix. The normalized contacts in the 50-kb-by-50-kb bin at the center of the matrix was measured to generate the CSS for the interaction of interest. To fairly compare and perform statistics on structures that don’t contain an identical number of elements, we generated scaled CSS values where the contacts in a 50-kb-by-50-kb bin at the center of the matrix in an individual cell was divided by all contacts in the 2-Mb by 2-Mb matrix for the desired region in all cells.
 
Bulk Hi-C and Liquid Hi-C Experiments

In situ Hi-C, Liquid Hi-C, and HiChIP contact analysis & contact specificity score calculation (CSS)
Contacts were extracted at 50-kb resolution from cooler files and normalized to contacts per billion for comparison between libraries. To compute CSS, contacts within a 1-Mb radius surrounding two regions of interest were computed, between all regions of interest in one condition and summed, generating a 2-Mb-by-2-Mb matrix of contacts. Contact specificity scores were generated by normalizing the contacts within each 50-kb-by-50-kb bin to the total number of contacts in the matrix. The normalized contacts in the 50-kb-by-50-kb bin at the center of the matrix was measured to generate the CSS for the interaction of interest. 

Loss of Structure in Liquid Hi-C
Loss of structure calculation was performed as described in Belaghzal et al. 2021[7] with the following modifications. To examine increasing loss of structure with increased lengths of pre-digestion, normalized contacts made within a 50-kb bin were divided by the total number of contacts made by the same bin to its surrounding 3-Mb radius in cis, generating the loss of structure score (LOS) for that bin. Successful Liquid Hi-C was verified by observing higher levels of LOS for compartment A vs. compartment B chromatin, that increased with longer pre-digestion. Compartment profiles were generated by performing eigendecomposition on 50-kb binned cooler files with cooltools.

