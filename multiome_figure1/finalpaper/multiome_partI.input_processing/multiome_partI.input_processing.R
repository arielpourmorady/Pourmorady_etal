# Integrating ATAC with scRNA-seq data https://satijalab.org/signac/articles/pbmc_vignette.html

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


# Before doing peak calling, you have to use Seurat
# to cluster your cells and then assign cell types
# to your data! This is important.

### AT THIS TIME RUN SEURAT_CLUSTERING.R

#celltypes <- moe_rna@active.ident
#moe <- AddMetaData(object = moe,
#                   metadata = celltypes,
#                   col.name = "celltypes")
#Idents(moe) <- "celltypes"

# call peaks using MACS2
#peaks <- CallPeaks(moe, 
#                   macs2.path = "/usr/local/bin/macs2",
#                   group.by = "celltypes")

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

mean(gi_fraction_accessibility$fraction_acc)
mean(link_fraction_accessibility$fraction_acc)

t.test(link_fraction_accessibility$fraction_acc , gi_fraction_accessibility$fraction_acc)

######################################################
######################################################
############### ATAC Tiling for cCREs ################ 
######################################################
######################################################

DefaultAssay(mOSN_moe) <- "ATAC"

#chr9-49880001-49885000 (boundaries are adjusted so that enhnacer is in center)
tile_plot1 <- TilePlot(
  object = mOSN_moe,
  region = "chr9-49883700-49886300",
  idents = 0,
  tile.cells = 200
)
tile_plot1

#H
tile_plot2 <- TilePlot(
  object = mOSN_moe,
  region = c("chr14-52546909-52549295"),
  idents = 0,
  tile.cells = 200,
)

#Lipsi
tile_plot4 <- TilePlot(
  object = mOSN_moe,
  region = "chr2-37126864-37129314",
  idents = 0,
  tile.cells = 200,
)

tile_plot1 + tile_plot2 + tile_plot4


#write.table(links_df, file = "/data/finalpaper/multiome_partI.input_processing/mOSN_cCREs.txt", sep = "\t")
#saveRDS(moe, file = "/data/finalpaper/multiome_partI.input_processing/moe.rds")
#saveRDS(mOSN_moe, file = "/data/finalpaper/multiome_partI.input_processing/mOSN_moe.rds")

