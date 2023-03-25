library(dplyr)
library(tidyverse)
library("DESeq2")
library("BiocParallel")
library("GenomicAlignments")
library("GenomicFeatures")
library("GenomicRanges")
library("Rsamtools")
library("AnnotationDbi")
library("rtracklayer")
library("org.Mm.eg.db")
library("genefilter")
library("ggrepel")
library(reshape2)
library("RColorBrewer")
library("pheatmap")
library(gridExtra)
library(fgsea)
library(cowplot)
library(grid)
library(patchwork)

theme_set(theme_cowplot())
filter_genes <- function(x){
  x <- intersect(x, rownames(dds))
  return(x)
}
remove_ORs <- function(x) {
  x <- setdiff(x,Olfr)
  return(x)
}
firstColToRowName <- function(x) {
  rownames(x) <- x[,1]
  x[,1] <- NULL
  return(x)
}
#Load reference
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes)
exonsByGene <- exonsBy(genes_txdb, by="gene")
sampleTable <- read.delim("/media/storageE/ariel/R/finalpaper/rna.development.p2.ko/sampletable.txt",header =T)
bamfiles <- filenames <- file.path(sampleTable$data)
bamfiles <- BamFileList(bamfiles)
rna <- summarizeOverlaps(features=exonsByGene, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=FALSE, preprocess.reads=invertStrand)
colData(rna) <- DataFrame(sampleTable)
colnames(rna)<-sampleTable$library

Atf5_reps <- c('Atf5_rep2', 'Atf5_rep3')
Icam_reps <- c('Icam_rep1', 'Icam_rep2')
Mash1_reps <- c('Mash_rep1', 'Mash_rep2')
Mash1Ngn_reps <- c('MashNgn_rep1', 'MashNgn_rep2')
Ngn_reps <- c('Ngn_rep1', 'Ngn_rep2')
OMP_reps <- c('OMP_rep1', 'OMP_rep2')

P2_reps <- c('OMPtTA_tetOP2_rep1',
             'gg8tTA_tetOP2_WT_rep1',
             'gg8tTA_tetOP2_KO_rep1',
             'OMPtTA_tetOP2_KO_rep1',
             'gg8tTA_tetOP2_KO_rep1',
             'gg8tTA_tetOP2_KO_rep2',
             'OMPtTA_tetOP2_rep2',
             'OMPtTA_tetOP2_KO_rep2',
             'OMPtTA_tetOP2_KO_rep3',
             'gg8tTA_tetOP2_KO_rep3')

neuronal_lineage <- c(Atf5_reps, Icam_reps, Mash1_reps, Mash1Ngn_reps, Ngn_reps, OMP_reps) 
neuronal_P2_lineage <- c(neuronal_lineage, P2_reps)



# set field for comparing samples 
dds <- DESeqDataSet(rna, design = ~ condition)  # adjust field used for design as needed
# Set levels
#dds$tissue <- relevel(dds$tissue, ref="Ngn") # adjust ref as needed
dds$condition <- factor(dds$condition, levels=c("Icam","Mash1","Mash1Ngn","Ngn","Atf5", "OMP", 
                                                'gg8tTA_tetOP2_WT',
                                                'gg8tTA_tetOP2_KO',
                                                'OMPtTA_tetOP2',
                                                'OMPtTA_tetOP2_KO')) 
# remove genes with fewer than 10 mapped reads across all conditions
dds <- dds[ rowSums(counts(dds)) > 0, ]
# differential expressed genes
dds <- DESeq(dds)

# Finding differentially expressed genes



# Gene set analysis
Olfr <- filter_genes(scan(file = "/media/storageA/kevin/Ldb1_Paper/annotation/Olfr-geneID-for-DESeq2.txt", what = character()))
Pcdh <- filter_genes(scan(file = "/media/storageA/kevin/annotation/genes.Pcdh.txt", what = character()))
Synapsis <- filter_genes(scan(file="/media/storageA/kevin/annotation/genes.synapsis.txt",what = character()))
nonOlfr <- setdiff(rownames(dds),Olfr)
FishOR <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-FishOR.txt", what = character()))
nonFishOR <- setdiff(Olfr,FishOR)
Meitoic_phenotype <- filter_genes(scan(file = "/media/storageA/kevin/annotation/genes.meiotic_phenotype.txt", what = character()))
OMP_High <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes_OMP-Enriched.txt", what = character()))
Icam_High <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes_Icam-Enriched.txt", what = character()))
Ngn_High <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes_Ngn-Enriched.txt", what = character() ))
ORs_and_differential_genes <- list(Icam_High,Ngn_High,OMP_High,Olfr)
ORs_and_differential_genes <- Reduce(union,ORs_and_differential_genes)
plain_genes <- setdiff(rownames(dds),ORs_and_differential_genes)
TAAR <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.TAAR.txt", what = character()))
zone5 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone5.txt", what = character()))
zone4.5 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone4.5.txt", what = character()))
zone4 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone4.txt", what = character()))
zone3.5 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone3.5.txt", what = character()))
zone3 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone3.txt", what = character()))
zone2.5 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone2.5.txt", what = character()))
zone2 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone2.txt", what = character()))
zone1.5 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone1.5.txt", what = character()))
zone1 <- filter_genes(scan(file ="/media/storageA/kevin/annotation/genes.ORs-zone1.txt", what = character()))
zone1_markers <- filter_genes(scan(file="/media/storageA/kevin/annotation/genes.Zone1-markers-from-Lisa.padj-05.txt", what = character()))
zone5_markers <- filter_genes(scan(file="/media/storageA/kevin/annotation/genes.Zone5-markers-from-Lisa.padj-05.txt", what = character()))
DNA_binding <- filter_genes(scan(file="/media/storageA/kevin/annotation/genes.DNA-binding_SeqSpecific.txt", what = character()))
nucleosome <- filter_genes(scan(file="/media/storageA/kevin/annotation/genes.GO_CC_MM_NUCLEOSOME.txt",what=character()))

#annotation
ORclusters <- scan(file ="/media/storageA/kevin/annotation/ORclusters.txt", what = character())
ORs.byCluster <- read.delim("/media/storageA/kevin/annotation/ORs-by-cluster.txt", header= FALSE)
colnames(ORs.byCluster) <- c("OR","Cluster")
ORs.byCluster.bed <- read.delim("/media/storageA/kevin/annotation/ORs-by-cluster.bed", header= FALSE)
colnames(ORs.byCluster.bed) <- c("chr","left","right","OR","Cluster")
cluster.info <- data.frame(ORs.byCluster.bed[,c(5,1:3)])
rownames(cluster.info) <- ORs.byCluster.bed$OR
ORs.byZone <- read.delim("/media/storageA/kevin/annotation/ORs-by-zone.txt", header= FALSE)
colnames(ORs.byZone) <- c("OR","zone")
zone.info <- data.frame(zone = ORs.byZone[,2])
rownames(zone.info) <- ORs.byZone$OR

# generate normalized counts for QC
# Normalizing transformation of count data and export normalized counts (normalized counts per gene)
rld <- rlog(dds, blind =T)
rlogMat <- assay(rld)
rlog_rna <- as.data.frame(rlogMat)
#write.table(as.data.frame(rlogMat), file="Developmental_Proteome_RNASeq.12Aug19.rlogMat.txt")

#### QC
##PCA analysis

plotPCA(rld[nonOlfr, neuronal_lineage], intgroup=c("condition"))
plotPCA(rld[nonOlfr, neuronal_P2_lineage], intgroup=c("condition"))
# 
# res <- results(dds, contrast = c("condition", 'OMPtTA_tetOM71_KO_10ug', 'OMP'))
# resSig <- subset(res, padj < 0.05)
# resSig <- resSig[rownames(resSig) %in% nonOlfr,]
# topgenes <- head(resSig[order(-resSig$log2FoldChange), ], 100) %>% rownames()

geneCounts_Krt5 <- plotCounts(dds[,neuronal_lineage], gene ="Krt5", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ascl1 <- plotCounts(dds[,neuronal_lineage], gene ="Ascl1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Neurog1 <- plotCounts(dds[,neuronal_lineage], gene ="Neurog1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Neurod1 <- plotCounts(dds[,neuronal_lineage], gene ="Neurod1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Tex15 <- plotCounts(dds[,neuronal_lineage], gene ="Tex15", intgroup = c("condition"),returnData = TRUE)
geneCounts_Gap43 <- plotCounts(dds[,neuronal_lineage], gene ="Gap43", intgroup = c("condition"),returnData = TRUE)
geneCounts_Omp <- plotCounts(dds[,neuronal_lineage], gene ="Omp", intgroup = c("condition"),returnData = TRUE)
geneCounts_Atf5 <- plotCounts(dds[,neuronal_lineage], gene ="Atf5", intgroup = c("condition"),returnData = TRUE)

a <- ggplot(geneCounts_Krt5, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Krt5 expression")

b <- ggplot(geneCounts_Ascl1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ascl1 expression")

c <- ggplot(geneCounts_Neurog1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Neurog1 expression")

d <- ggplot(geneCounts_Neurod1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Neurod1 expression")

e <- ggplot(geneCounts_Tex15, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Tex15 expression")

f <- ggplot(geneCounts_Gap43, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Gap43 expression")

g <- ggplot(geneCounts_Omp, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Omp expression")

h <- ggplot(geneCounts_Atf5, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Atf5 expression")

grid.arrange(a, b, c, d, 
             e, f, h, g, ncol = 4)


geneCounts_Lhx2 <- plotCounts(dds[,neuronal_lineage], gene ="Lhx2", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ebf1 <- plotCounts(dds[,neuronal_lineage], gene ="Ebf1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ldb1 <- plotCounts(dds[,neuronal_lineage], gene ="Ldb1", intgroup = c("condition"),returnData = TRUE)

i <- ggplot(geneCounts_Lhx2, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Lhx2 expression")

j <- ggplot(geneCounts_Ebf1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ebf1 expression")

k <- ggplot(geneCounts_Ldb1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ldb1 expression")

i + j + k


