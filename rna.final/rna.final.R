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

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_finalized/rna.final") 

#######################################
############## FUNCTIONS ##############
#######################################

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

#######################################
######## EXTRACTING RAW COUNTS ########
#######################################

#Load reference
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes)
exonsByGene <- exonsBy(genes_txdb, by="gene")

# Only do this one time to generate table rds file of raw counts

# # Developmental Files
# sampleTable_dev <- read.delim("rna.final/sampletable_complete.txt",header =T)
# bamfiles_dev <- filenames <- file.path(sampleTable_dev$data)
# bamfiles_dev <- BamFileList(bamfiles_dev)
# rna_dev <- summarizeOverlaps(features=exonsByGene, reads=bamfiles, mode="Union", singleEnd=FALSE, 
#                          ignore.strand=FALSE, preprocess.reads=invertStrand)
# 
# 
# colData(rna_dev) <- DataFrame(sampleTable_dev)
# colnames(rna_dev)<-sampleTable$library
# saveRDS(rna_dev, file = "rna.final/rna_dev.rds")
# 
# # mutantsant Files
# sampleTable_mutants <- read.delim("rna.final/sampletable_complete_ptII.txt",header =T)
# bamfiles_mutants <- filenames <- file.path(sampleTable_mutants$data)
# bamfiles_mutants <- BamFileList(bamfiles_mutants)
# rna_mutants <- summarizeOverlaps(features=exonsByGene, reads=bamfiles, mode="Union", singleEnd=FALSE, 
#                              ignore.strand=FALSE, preprocess.reads=invertStrand)
# 
# 
# colData(rna_mutants) <- DataFrame(sampleTable_mutants)
# colnames(rna_mutants)<-sampleTable$library
# saveRDS(rna_mutants, file = "rna.final/rna_mutants.rds")

#######################################
########## LOADING RAW COUNTS #########
#######################################

rna_dev <- readRDS(file = "rna.final/rna_dev.rds")
rna_mutants <- readRDS(file = "rna.final/rna_mutants.rds")

setequal(rowData(rna_dev)@rownames, rowData(rna_mutants)@rownames) ## TRUE
idx <- match(rowData(rna_dev)@rownames, rowData(rna_mutants)@rownames)
rna <- cbind(rna_dev, rna_mutants[idx,])

P2_reps <- c('OMPtTA_tetOP2_rep1',
             'gg8tTA_tetOP2_WT_rep1',
             'gg8tTA_tetOP2_KO_rep1',
             'OMPtTA_tetOP2_KO_rep1',
             'gg8tTA_tetOP2_KO_rep1',
             'gg8tTA_tetOP2_KO_rep2',
             'OMPtTA_tetOP2_rep2',
             'OMPtTA_tetOP2_KO_rep2',
             'OMPtTA_tetOP2_KO_rep3',
             'gg8tTA_tetOP2_KO_rep3',
             'gg8tTA_tetOP2_WT_rep2',
             'gg8tTA_tetOP2_WT_rep3')

tetOGFP_reps <- c('gg8tTA_tetOGFP_rep1',
                  'gg8tTA_tetOGFP_rep2',
                  'OMPtTA_tetOGFP_rep1',
                  'OMPtTA_tetOGFP_rep2')

M71_reps <- c('OMPtTA_tetOM71_KO_rep1',
              'OMPtTA_tetOM71_KO_rep2',
              'OMPtTA_tetOM71_KO_rep3')

Atf5_reps <- c('Atf5_rep2', 'Atf5_rep3')
Icam_reps <- c('Icam_rep1', 'Icam_rep2')
Mash1_reps <- c('Mash_rep1', 'Mash_rep2')
Mash1Ngn_reps <- c('MashNgn_rep1', 'MashNgn_rep2')
Ngn_reps <- c('Ngn_rep1', 'Ngn_rep2')
OMP_reps <- c('OMP_rep1', 'OMP_rep2')

neuronal_lineage <- c(Icam_reps,Mash1_reps, Mash1Ngn_reps, Ngn_reps, Atf5_reps, OMP_reps)
neuronal_lineage_P2_reps <- c(neuronal_lineage, P2_reps)

#######################################
########### DESEQ GENERATION ##########
#######################################

# set field for comparing samples 
dds <- DESeqDataSet(rna, design = ~ condition)  # adjust field used for design as needed
# Set levels
#dds$tissue <- relevel(dds$tissue, ref="Ngn") # adjust ref as needed
dds$condition <- factor(dds$condition, levels=c('Icam',
                                                'Mash1',
                                                'Mash1Ngn',
                                                'Ngn',
                                                'Atf5',
                                                'OMP',
                                                'gg8tTA_tetOP2_WT',
                                                'gg8tTA_tetOP2_KO',
                                                'OMPtTA_tetOP2',
                                                'OMPtTA_tetOP2_KO',
                                                'gg8tTA_tetOGFP',
                                                'OMPtTA_tetOGFP',
                                                'OMPtTA_tetOM71_KO')) 

# remove genes with fewer than 10 mapped reads across all conditions
dds <- dds[rowSums(counts(dds)) > 0, ]
# differential expressed genes
dds <- DESeq(dds)

#######################################
############ DATA ANALYSIS ############
#######################################

# Gene set analysis
Olfr <- filter_genes(scan(file = "/media/storageA/kevin/Ldb1_Paper/annotation/Olfr-geneID-for-DESeq2.txt", what = character()))
nonOlfr <- setdiff(rownames(dds),Olfr)

# generate normalized counts for QC
# Normalizing transformation of count data and export normalized counts (normalized counts per gene)
rld <- rlog(dds, blind =T)
rlogMat <- assay(rld)
rlog_rna <- as.data.frame(rlogMat)

##PCA analysis

plotPCA(rld[nonOlfr, neuronal_lineage_P2_reps], intgroup=c("condition"))
plotPCA(rld[nonOlfr, ], intgroup=c("condition"))


## OR logFC Comparison against tetO-GFP

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOP2", "OMPtTA_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]

resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_a <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

a <- ggplot(resLFC_a %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOP2 vs. OMPtTA_tetOGFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")
a

resLFC <- lfcShrink(dds, contrast = c("condition", "gg8tTA_tetOP2_WT", "OMPtTA_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_b <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

b <- ggplot(resLFC_b %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("gg8tTA_tetOP2_WT vs. OMPtTA_tetOGFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOP2_KO", "OMPtTA_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_c <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

c <- ggplot(resLFC_c %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOP2_KO vs. OMPtTA_tetOGFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "gg8tTA_tetOP2_KO", "gg8tTA_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_d <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

d <- ggplot(resLFC_d %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("gg8tTA_tetOP2_KO vs. gg8tTA_tetOGFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOM71_KO", "OMPtTA_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_e <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

e <- ggplot(resLFC_e %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOM71_KO vs. OMPtTA_tetOGFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")



## OR logFC Comparison against same developmental stage

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOP2", "OMP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]

resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_f <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

f <- ggplot(resLFC_f %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOP2 vs. OMP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "gg8tTA_tetOP2_WT", "OMP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_g <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

g <- ggplot(resLFC_g %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("gg8tTA_tetOP2_WT vs. OMP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOP2_KO", "OMP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_h <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

h <- ggplot(resLFC_h %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOP2_KO vs. OMP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "gg8tTA_tetOP2_KO", "Atf5"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_i <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

i <- ggplot(resLFC_i %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("gg8tTA_tetOP2_KO vs. Atf5") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOM71_KO", "OMP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_j <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

j <- ggplot(resLFC_j %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  ylim(-10,15) + 
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOM71_KO vs. OMP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")


a + b + c + d + e + f + g + h + i + j + plot_layout(nrow = 2)


# Developmental Gene Counts

geneCounts_Krt5 <- plotCounts(dds[,neuronal_lineage], gene ="Krt5", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ascl1 <- plotCounts(dds[,neuronal_lineage], gene ="Ascl1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Neurog1 <- plotCounts(dds[,neuronal_lineage], gene ="Neurog1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Neurod1 <- plotCounts(dds[,neuronal_lineage], gene ="Neurod1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Tex15 <- plotCounts(dds[,neuronal_lineage], gene ="Tex15", intgroup = c("condition"),returnData = TRUE)
geneCounts_Gap43 <- plotCounts(dds[,neuronal_lineage], gene ="Gap43", intgroup = c("condition"),returnData = TRUE)
geneCounts_Omp <- plotCounts(dds[,neuronal_lineage], gene ="Omp", intgroup = c("condition"),returnData = TRUE)
geneCounts_Atf5 <- plotCounts(dds[,neuronal_lineage], gene ="Atf5", intgroup = c("condition"),returnData = TRUE)

k <- ggplot(geneCounts_Krt5, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Krt5 expression")

l <- ggplot(geneCounts_Ascl1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ascl1 expression")

m <- ggplot(geneCounts_Neurog1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Neurog1 expression")

n <- ggplot(geneCounts_Neurod1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Neurod1 expression")

o <- ggplot(geneCounts_Tex15, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Tex15 expression")

p <- ggplot(geneCounts_Gap43, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Gap43 expression")

q <- ggplot(geneCounts_Omp, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Omp expression")

r <- ggplot(geneCounts_Atf5, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Atf5 expression")

k + l + m + n + o + p + q + r + plot_layout(ncol = 4)


geneCounts_Lhx2 <- plotCounts(dds[,neuronal_lineage], gene ="Lhx2", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ebf1 <- plotCounts(dds[,neuronal_lineage], gene ="Ebf1", intgroup = c("condition"),returnData = TRUE)
geneCounts_Ldb1 <- plotCounts(dds[,neuronal_lineage], gene ="Ldb1", intgroup = c("condition"),returnData = TRUE)

s <- ggplot(geneCounts_Lhx2, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("Lhx2 expression")

t <- ggplot(geneCounts_Ebf1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ebf1 expression")

u <- ggplot(geneCounts_Ldb1, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("Ldb1 expression")

s + t + u

