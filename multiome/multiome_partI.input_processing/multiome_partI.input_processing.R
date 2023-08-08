# Integrating ATAC with scRNA-seq data https://satijalab.org/signac/articles/moe_rna_II_vignette.html

library(dplyr)
library(tidyr)
library(reshape2)
library(Signac)
library(Seurat)
library(SeuratDisk)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

set.seed(1234)
options(scipen=999)


# load the RNA and ATAC data
counts <- Read10X_h5("/data/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/data/outs/atac_fragments.tsv.gz"
fragments <- CreateFragmentObject(fragpath)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA data
moe <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
moe[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
moe

## Quality Control

DefaultAssay(moe) <- "ATAC"

moe <- NucleosomeSignal(moe)
moe <- TSSEnrichment(moe)

VlnPlot(
  object = moe,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
moe <- subset(
  x = moe,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
moe

# call peaks using MACS2
peaks <- CallPeaks(moe, 
                   macs2.path = "/usr/local/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(moe),
  features = peaks,
  cells = colnames(moe)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
moe[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# Gene expression data processing

DefaultAssay(moe) <- "RNA"
moe <- SCTransform(moe)
moe <- RunPCA(moe)


DefaultAssay(moe) <- "peaks"
moe <- FindTopFeatures(moe, min.cutoff = 5)
moe <- RunTFIDF(moe)
moe <- RunSVD(moe)

## Joint UMAP Visualiazation
# build a joint neighbor graph using both assays
moe <- FindMultiModalNeighbors(
  object = moe,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

moe <- FindClusters(moe, 
                    graph.name = "wknn",
                    resolution = 0.4)

# build a joint UMAP visualization
moe <- RunUMAP(
  object = moe,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(moe, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

######################################################
######################################################
########## Validation of Cluster Identities ########## 
######################################################
######################################################

DefaultAssay(moe) <- "SCT"
moe.markers <- FindAllMarkers(moe, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
top.markers <- moe.markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

View(moe.markers)
# Cluster 0 - Omp & Cnga2 = mOSN
# Cluster 10 - Ascl1 = GBC
# Cluster 4 - Neurod1 = INP
# Cluter 4 - Gap43 - iOSN

plot1 <- DimPlot(moe, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

DefaultAssay(moe_rna) <- 'SCT'
plot2 <- FeaturePlot(moe, reduction = "umap", features = "Ascl1", order = TRUE)
plot2
plot3 <- FeaturePlot(moe, reduction = "umap", features = "Neurod1", order = TRUE)
plot3
plot4 <- FeaturePlot(moe, reduction = "umap", features = "Gap43", order = TRUE)
plot4
plot5 <- FeaturePlot(moe, reduction = "umap", features = "Omp", order = TRUE)
plot5


######################################################
######################################################
############ Link Peaks to Top OSN Genes ############# 
######################################################
######################################################

DefaultAssay(moe) <- "RNA"

cluster0.markers <- FindMarkers(moe, ident.1 = 0, min.pct = 0.25)

DefaultAssay(moe) <- "peaks"

moe <- RegionStats(moe,  genome = BSgenome.Mmusculus.UCSC.mm10)
moe <- LinkPeaks(
  object = moe,
  peak.assay = "peaks",
  expression.assay = "SCT",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = 2000,
  min.cells = 10,
  genes.use = row.names(cluster0.markers),
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.2,
  gene.id = FALSE,
  verbose = TRUE
)

links <- Links(moe) %>% 
  as.data.frame() %>% 
  dplyr::filter(score > 0)
links_df <- links[order(-links$score),] %>%
  dplyr::select(peak) %>%
  dplyr::rename("loc" = "peak")


# Fraction Accessibility of GIs vs. cCREs

mOSN_moe <- subset(moe, idents = 0)

gi_bed <- as.data.frame(read.table("/data/multiome/Greek_Islands.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

fraction_accessibility <- function(bed){
  df <- FeatureMatrix(
    fragments,
    features = bed[,1],
    cells = Cells(mOSN_moe),
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE) %>% 
    as.data.frame()
  df$loc <- rownames(df)
  df <- melt(df)
  df$value[df$value > 0] <- 1
  mOSN_cellcount <- unique(df$variable) %>% length()
  df <- df %>% group_by(loc) %>%
    summarise(fraction_acc = sum(value)/mOSN_cellcount)
  return(df)
}

gi_fraction_accessibility <- fraction_accessibility(gi_bed)
gi_fraction_accessibility$type <- "GIs"

link_fraction_accessibility <- fraction_accessibility(links_df)
link_fraction_accessibility$type <- "links"

fraction_accessibility_df <- rbind(gi_fraction_accessibility, link_fraction_accessibility)

ggplot(fraction_accessibility_df, aes(x = type, y = fraction_acc)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.2) +
  ylab("Percent accessibility \n per enhancer (%)") + 
  theme_classic()

##### ##### ##### ##### ##### 
##### Length Normalized ##### 
##### ##### ##### ##### ##### 

gi_fraction_accessibility_distnorm <- gi_fraction_accessibility %>%
  separate(loc, sep = "-", into = c("chr", "loc1", "loc2")) %>%
  mutate(length = as.numeric(loc2) - as.numeric(loc1)) %>% 
  mutate(distnorm_fraction_acc = 1000*fraction_acc/length) 

link_fraction_accessibility_distnorm <- link_fraction_accessibility %>%
  separate(loc, sep = "[.]", into = c("chr", "loc1", "loc2")) %>%
  mutate(length = as.numeric(loc2) - as.numeric(loc1)) %>% 
  mutate(distnorm_fraction_acc = 1000*fraction_acc/length) 

fraction_accessibility_df_distnorm <- rbind(gi_fraction_accessibility_distnorm, link_fraction_accessibility_distnorm)

ggplot(df3, aes(x = type, y = distnorm_fraction_acc)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.2) +
  ylab("Percent accessibility \n per enhancer (%)") + 
  theme_classic()

gi_fraction_accessibility_distnorm %>%
  group_by(type) %>%
  summarise(mean = 100*mean(distnorm_fraction_acc),
            sd = 100*sd(distnorm_fraction_acc),
            n = n())

link_fraction_accessibility_distnorm %>%
  group_by(type) %>%
  summarise(mean = 100*mean(distnorm_fraction_acc),
            sd = 100*sd(distnorm_fraction_acc),
            n = n())

t.test(gi_fraction_accessibility_distnorm$distnorm_fraction_acc ,
       link_fraction_accessibility_distnorm$distnorm_fraction_acc)


#write.table(links_df, file = "/data/finalpaper_August2023/multiome_partI.input_processing/mOSN_cCREs.txt", sep = "\t")
#saveRDS(moe, file = "/data/finalpaper_August2023/multiome_partI.input_processing/moe.rds")
#saveRDS(mOSN_moe, file = "/data/finalpaper_August2023/multiome_partI.input_processing/mOSN_moe.rds")


########################
########################
### NATURE REVISIONS ###
########################
########################

# Clustering MOE based off of RNA data alone

# # Note that all operations below are performed on the RNA assay Set and verify that the
# # default assay is RNA
# DefaultAssay(moe_rna) <- "RNA"
# DefaultAssay(moe_rna)
# # perform visualization and clustering steps
# moe_rna <- NormalizeData(moe_rna)
# moe_rna <- FindVariableFeatures(moe_rna)
# moe_rna <- ScaleData(moe_rna)
# moe_rna <- RunPCA(moe_rna, verbose = FALSE)
# moe_rna <- FindNeighbors(moe_rna, dims = 1:30)
# moe_rna <- FindClusters(moe_rna, resolution = 0.8, verbose = FALSE)
# moe_rna <- RunUMAP(moe_rna, dims = 1:30)
# 
# plot6 <- DimPlot(moe_rna, repel = TRUE, group.by = ident('seurat_clusters') ,label = TRUE) + NoLegend()
# plot7 <- DimPlot(moe_rna, repel = TRUE, group.by = ident('wknn_res.0.4'), label = TRUE) + NoLegend()
# plot8 <- DimPlot(moe, label = TRUE, repel = TRUE, group.by = ident('wknn_res.0.4'), reduction = "umap") + NoLegend() #UMAP projection of SCT/peak weighted nearest neighbor
# plot6 + plot7 + plot8

# TilePlots of all cCREs

GIs_2000bp_df <- fraction_accessibility_df %>% 
  dplyr::filter(type == "GIs") %>%
  separate(loc, sep = "-", into = c("chr", "loc1", "loc2")) %>%
  mutate(loc1 = as.numeric(loc1) - 1000,
         loc2 = as.numeric(loc2) + 1000,
         loc_spaced = paste(chr, loc1, loc2, sep="-")) %>%
  dplyr::select(loc_spaced, type)

links_2000bp_df <- fraction_accessibility_df %>% 
  dplyr::filter(type == "links") %>%
  separate(loc, sep = "[.]", into = c("chr", "loc1", "loc2")) %>%
  mutate(loc1 = as.numeric(loc1) - 1000,
         loc2 = as.numeric(loc2) + 1000,
         loc_spaced = paste(chr, loc1, loc2, sep="-")) %>%
  dplyr::select(loc_spaced, type)

tileplot_GIs_fn <- function(loc){
  tile_plot <- TilePlot(
    object = mOSN_moe,
    region = loc,
    idents = 0,
    tile.cells = 200
  )
  ggsave(paste("/data/finalpaper_revisions/multiome_partI.input_processing/GI_tileplots/", loc, ".pdf", sep = ""), tile_plot)
  return(tile_plot)
}

tileplot_links_fn <- function(loc){
  tile_plot <- TilePlot(
    object = mOSN_moe,
    region = loc,
    idents = 0,
    tile.cells = 200
  )
  ggsave(paste("/data/finalpaper_revisions/multiome_partI.input_processing/cCRE_tileplots/", loc, ".pdf", sep = ""), tile_plot)
  return(tile_plot)
}

for (i in c(1:nrow(GIs_2000bp_df))){
  tileplot_GIs_fn(as.character(GIs_2000bp_df[i,1]))
}

for (i in c(1:nrow(links_2000bp_df))){
  tileplot_links_fn(as.character(links_2000bp_df[i,1]))
}


#### Adding OR promoter to Accessibility Boxpolot in mOSNs

or_promoter_bed <- as.data.frame(read.table("/data/annotations/OR-xscript-Soria+UCSC.mm10.GeneID.300bp-promoter.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

or_promoter_fraction_accessibility <- fraction_accessibility(or_promoter_bed)
or_promoter_fraction_accessibility$type <- "or_promoters"

fraction_accessibility_df <- rbind(gi_fraction_accessibility, link_fraction_accessibility, or_promoter_fraction_accessibility)
fraction_accessibility_df$type <- factor(fraction_accessibility_df$type, levels = c("links", "GIs", "or_promoters"))

ggplot(fraction_accessibility_df, aes(x = type, y = fraction_acc)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.2) +
  ylab("Percent accessibility \n per site (%)") + 
  theme_classic()

fraction_accessibility_df %>% head()
fraction_accessibility_df %>% group_by(type) %>%
  summarise(mean = mean(fraction_acc)*100,
            sd = sd(fraction_acc)*100)
