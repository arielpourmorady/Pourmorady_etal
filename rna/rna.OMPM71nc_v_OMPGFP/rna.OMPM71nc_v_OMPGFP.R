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
library(topGO)
library(KEGGREST)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/rna/rna.OMPM71nc_v_OMPGFP/") 

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

# Developmental Files
sampleTable <- read.delim("sampletable.txt",header =T)
bamfiles <- filenames <- file.path(sampleTable$data)
bamfiles <- BamFileList(bamfiles)
rna <- summarizeOverlaps(features=exonsByGene, reads=bamfiles, mode="Union", singleEnd=FALSE,
                         ignore.strand=FALSE, preprocess.reads=invertStrand)


colData(rna) <- DataFrame(sampleTable)
colnames(rna)<-sampleTable$library
saveRDS(rna, file = "rna.rds")

#######################################
########### DESEQ GENERATION ##########
#######################################

# set field for comparing samples 
dds <- DESeqDataSet(rna, design = ~ condition)  # adjust field used for design as needed
# Set levels
#dds$tissue <- relevel(dds$tissue, ref="Ngn") # adjust ref as needed
dds$condition <- factor(dds$condition, levels=c('OMP_GFP',
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

##########################################################
############ Differential Expression Analysis ############
##########################################################

resLFC <- lfcShrink(dds, contrast = c("condition", "OMPtTA_tetOM71_KO", "OMP_GFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
resLFC_sig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj < 0.05) %>% mutate(sig = "yes")
resLFC_nonsig <- resLFC %>% as.data.frame() %>% dplyr::filter(padj >= 0.05) %>% mutate(sig = "no")
resLFC_nonsig_NA <- resLFC %>% as.data.frame() %>% dplyr::filter(is.na(padj)) %>% mutate(sig = "no")
resLFC_d <- rbind(resLFC_sig, resLFC_nonsig, resLFC_nonsig_NA)

ggplot(resLFC_d %>% arrange(sig), aes(x = baseMean, y = log2FoldChange)) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_hline(yintercept = 0, linetype="dotted", color = "red", size=1.5) +  
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOM71_KO vs. OMP_GFP") + 
  ylab("log2FC mut/ctrl") + 
  xlab("mean of normalized counts")

ggplot(resLFC_d, aes(x = log2FoldChange, y = -log(padj))) + 
  geom_point(size = 0.5, aes(color = sig)) + 
  scale_color_manual(values = c("darkslategrey", "red")) + 
  geom_vline(xintercept = 0, linetype="dotted", color = "black", size = 1) +  
  theme_classic() +
  theme(legend.position = "none") + 
  ggtitle("OMPtTA_tetOM71_KO vs. OMP_GFP") + 
  ylab("-log10(padj)") + 
  xlab("log2FoldChange")


resLFC_d$comparison <- "OMPtTA_tetOM71_KO vs. OMP_GFP"
write.table(resLFC_d, file = "SuppFig10h_MAPlot.txt")
