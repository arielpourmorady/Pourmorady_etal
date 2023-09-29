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

setwd("/media/storageE/ariel/R/finalpaper_August2023/rna/ompttatetop2.mor28icretdt.rna/") 

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
# sampleTable_dev <- read.delim("sampletable_dev.txt",header =T)
# bamfiles_dev <- filenames <- file.path(sampleTable_dev$data)
# bamfiles_dev <- BamFileList(bamfiles_dev)
# rna_dev <- summarizeOverlaps(features=exonsByGene, reads=bamfiles_dev, mode="Union", singleEnd=FALSE,
#                          ignore.strand=FALSE, preprocess.reads=invertStrand)
# 
# 
# colData(rna_dev) <- DataFrame(sampleTable_dev)
# colnames(rna_dev)<-sampleTable_dev$library
# saveRDS(rna_dev, file = "rna_dev.rds")
# 
# # mutantsant Files
# sampleTable_mutants <- read.delim("sampletable_RYG.txt",header =T)
# bamfiles_mutants <- filenames <- file.path(sampleTable_mutants$data)
# bamfiles_mutants <- BamFileList(bamfiles_mutants)
# rna_mutants <- summarizeOverlaps(features=exonsByGene, reads=bamfiles_mutants, mode="Union", singleEnd=FALSE,
#                              ignore.strand=FALSE, preprocess.reads=invertStrand)
# 
# 
# colData(rna_mutants) <- DataFrame(sampleTable_mutants)
# colnames(rna_mutants)<-sampleTable_mutants$library
# saveRDS(rna_mutants, file = "rna_mutants.rds")

#######################################
########## LOADING RAW COUNTS #########
#######################################

rna_dev <- readRDS(file = "rna_dev.rds")
rna_mutants <- readRDS(file = "rna_mutants.rds")

setequal(rowData(rna_dev)@rownames, rowData(rna_mutants)@rownames) ## TRUE
idx <- match(rowData(rna_dev)@rownames, rowData(rna_mutants)@rownames)
rna <- cbind(rna_dev, rna_mutants[idx,])

#############################
########## DESEQ ############
#############################

Atf5_reps <- c('Atf5_rep2', 'Atf5_rep3')
Icam_reps <- c('Icam_rep1', 'Icam_rep2')
Mash1_reps <- c('Mash_rep1', 'Mash_rep2')
Mash1Ngn_reps <- c('MashNgn_rep1', 'MashNgn_rep2')
Ngn_reps <- c('Ngn_rep1', 'Ngn_rep2')
OMP_reps <- c('OMP_rep1', 'OMP_rep2')

mor28icretdTom_ompttatetop2 <- c('mor28iCretdTom_OMPtTAtetOP2_green_rep1',
                                 'mor28iCretdTom_OMPtTAtetOP2_red_rep1',
                                 'mor28iCretdTom_OMPtTAtetOP2_yellow_rep1',
                                 'mor28iCretdTom_OMPtTAtetOP2_green_rep2',
                                 'mor28iCretdTom_OMPtTAtetOP2_red_rep2',
                                 'mor28iCretdTom_OMPtTAtetOP2_yellow_rep2',
                                 'mor28iCretdTom_OMPtTAtetOP2_green_rep3',
                                 'mor28iCretdTom_OMPtTAtetOP2_red_rep3',
                                 'mor28iCretdTom_OMPtTAtetOP2_yellow_rep3')

neuronal_lineage <- c(Atf5_reps, Icam_reps, Mash1_reps, Mash1Ngn_reps, Ngn_reps, OMP_reps) 
neuronal_RYG_lineage <- c(neuronal_lineage, mor28icretdTom_ompttatetop2)



# set field for comparing samples 
dds <- DESeqDataSet(rna, design = ~ condition)  # adjust field used for design as needed
# Set levels
dds$condition <- factor(dds$condition, levels=c("Icam","Mash1","Mash1Ngn","Ngn","Atf5", "OMP",
                                                'mor28iCretdTom_OMPtTAtetOP2_green',
                                                'mor28iCretdTom_OMPtTAtetOP2_red',
                                                'mor28iCretdTom_OMPtTAtetOP2_yellow')) 
# remove genes with fewer than 10 mapped reads across all conditions
dds <- dds[ rowSums(counts(dds)) > 0, ]
# differential expressed genes
dds <- DESeq(dds)


################################
########## ANALYSIS ############
################################

# Gene set analysis
Olfr <- filter_genes(scan(file = "/media/storageA/kevin/Ldb1_Paper/annotation/Olfr-geneID-for-DESeq2.txt", what = character()))
nonOlfr <- setdiff(rownames(dds),Olfr)

# generate normalized counts for QC
# Normalizing transformation of count data and export normalized counts (normalized counts per gene)
rld <- rlog(dds, blind =T)
rlogMat <- assay(rld)
rlog_rna <- as.data.frame(rlogMat)

#write.table(rlog_rna, file = "mor28iCretdTom_OMPtTAtetOP2_normalizedCounts.txt", sep = "\t", row.names = TRUE, col.names = TRUE) 
#This file only contains rld normalized reads for the mor28icreomptTAtetOP2 mice


### Calculate unnormalized FPKM to generate TPM
fpkm_unscaled <- as.data.frame(fpkm(dds, robust = FALSE)) %>% rownames_to_column(var="GeneID")
fpkm_unscaled <- fpkm_unscaled %>% gather(geno,fpkm,-GeneID) %>% separate(geno,sep="_rep",into=c("geno","rep"))

### Calculate TPM
tpm <- fpkm_unscaled %>% group_by(geno,rep) %>% summarize(total_fpkm = sum(fpkm)) %>% 
  full_join(fpkm_unscaled) %>% mutate(tpm = (fpkm / total_fpkm) *1000000) %>% 
  dplyr::select(GeneID,geno,rep, tpm)

neuronal_RYG_lineage_norm_counts <- tpm %>% filter(GeneID %in% Olfr)

neuronal_RYG_lineage_norm_counts_Olfr17 <- neuronal_RYG_lineage_norm_counts %>%
  filter(GeneID == "Olfr17")
neuronal_RYG_lineage_norm_counts_Olfr17$OR <- "P2"

neuronal_RYG_lineage_norm_counts_Olfr1507 <- neuronal_RYG_lineage_norm_counts %>%
  filter(GeneID == "Olfr1507")
neuronal_RYG_lineage_norm_counts_Olfr1507$OR <- "mor28"

neuronal_RYG_lineage_norm_counts_inactive <- neuronal_RYG_lineage_norm_counts %>%
  filter(GeneID != "Olfr1507", GeneID != "Olfr17")
neuronal_RYG_lineage_norm_counts_inactive$OR <- "inactive"





neuronal_RYG_lineage_norm_counts <- rbind(neuronal_RYG_lineage_norm_counts_Olfr17, 
                                          neuronal_RYG_lineage_norm_counts_Olfr1507, 
                                          neuronal_RYG_lineage_norm_counts_inactive) %>%
  group_by(geno, rep, GeneID, OR) %>%
  summarise(tpm = mean(tpm)) %>%
  group_by(geno, rep, OR) %>%
  summarise(sum_counts = sum(tpm),
            mean_counts = median(tpm)) %>%
  mutate(norm_counts = mean_counts/sum_counts)

neuronal_RYG_lineage_norm_counts$geno <- factor(neuronal_RYG_lineage_norm_counts$geno, 
                                                levels = c("Icam","Mash","MashNgn","Ngn","Atf5", "OMP", 
                                                           'mor28iGFP',
                                                           'mor28iCretdTom_OMPtTAtetOP2_green',
                                                           'mor28iCretdTom_OMPtTAtetOP2_red',
                                                           'mor28iCretdTom_OMPtTAtetOP2_yellow'))

neuronal_RYG_lineage_norm_counts <- neuronal_RYG_lineage_norm_counts %>%
  filter(geno == 'mor28iCretdTom_OMPtTAtetOP2_green'|  geno == 'mor28iCretdTom_OMPtTAtetOP2_red'| geno == 'mor28iCretdTom_OMPtTAtetOP2_yellow')



################################################################################################################
a <- ggplot(neuronal_RYG_lineage_norm_counts, aes(x = geno, y = mean_counts, colour = OR)) + 
  geom_jitter(stat = "identity", position = position_dodge(width = 0.4)) + 
  ggtitle("TPM of OR Expression") +
  theme()+ 
  theme_light() + 
  ylab("TPM") + 
  scale_y_log10() + 
  scale_color_manual(values = c("black", "firebrick3", "green4"))
a

neuronal_RYG_lineage_norm_counts_supplement <- neuronal_RYG_lineage_norm_counts %>% dplyr::select(geno, rep, OR, mean_counts)
write.table(neuronal_RYG_lineage_norm_counts_supplement, file = "/media/storageE/ariel/R/finalpaper_August2023/rna/ompttatetop2.mor28icretdt.rna/Fig4h_ORSwitchingTPM.txt", sep = "\t")


geneCounts_Olfr17 <- plotCounts(dds[,mor28icretdTom_ompttatetop2], gene ="Olfr17", intgroup = c("condition"),returnData = TRUE)
geneCounts_Olfr1507 <- plotCounts(dds[,mor28icretdTom_ompttatetop2], gene ="Olfr1507", intgroup = c("condition"),returnData = TRUE)


t <- ggplot(geneCounts_Olfr17, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("P2 expression")

u <- ggplot(geneCounts_Olfr1507, aes(x = condition, y = count)) +
  scale_y_log10() +  geom_point(cex = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ggtitle("mor28 expression")

 t + u
