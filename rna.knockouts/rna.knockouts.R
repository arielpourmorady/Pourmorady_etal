library(GenomicFeatures)
library(DESeq2)
library(Rsamtools) #needed for BamFileList
library(rtracklayer) #needed to import GTF as GRanges object
library(GenomicAlignments) #used for summarizeOverlaps
library( "RColorBrewer" ) #for making heatmaps
library(gridExtra)
library("pheatmap") #used for heatmaps
library(ggrepel) #for genes points on plot
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
library(patchwork)
#edit 
#setwd("/path")

#load in annotation
#which GTF file should I use for the best OR annotation, Soria annotation not as acurate as Tan's
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
#previously used this gtf: /seq/mm10/iGenome/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf
genes_txdb <- makeTxDbFromGRanges(genes) # creates database
exonsByGene <- exonsBy(genes_txdb, by="gene") #GRangesList object, run exonsByGene$Olfr1507 to see the annotated exons
gene_width <- data.frame(sum(width(exonsByGene)))

#read in my data
sampleTable <- read.delim("/media/storageE/ariel/RNA/sample_table_knockouts_OMP_230129.txt",header =T)
#sampleTable <- sampleTable[1:2,]

bamfiles <- filenames <- file.path(sampleTable$data) #creates character vector with paths to bam files
bamfiles <- BamFileList(bamfiles, yieldSize=1000000)
sampleTable

# not sure if this is accurate
se <- summarizeOverlaps(features=exonsByGene, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE,
                        preprocess.reads=invertStrand)#,
#preprocess.reads=invertStrand) #for truseq data include preprocess.reads=invertStrand


colData(se) <-DataFrame(sampleTable) #is this additional metadata necessary?
colnames(se)<-sampleTable$geno

#one factor design formula is simple: ~ condition
#to do multifactor design (for example to control for a batch), then use design formula ~ batch + condition
#se$date<-factor(se$date)
dds <- DESeqDataSet(se, design = ~ geno)
dds$geno <- relevel(dds$geno, ref="OMPGFP") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 0.0, ] #should I use O reads instead?

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

#############################
# get FPKM AND TPM from DDS #
#############################
### Calculate FPKM
fpkm_scaled <- as.data.frame(fpkm(dds, robust = TRUE)) # robus normalizes on DEseq "scaling factor"
### Calculate unnormalized FPKM to generate TPM
fpkm_unscaled <- as.data.frame(fpkm(dds, robust = FALSE)) %>% rownames_to_column(var="GeneID")
fpkm_unscaled <- fpkm_unscaled %>% gather(geno,fpkm,-GeneID) #%>% separate(condition,sep=".",into=c("geno","rep"))

#fpkm_unscaled_a <- fpkm_unscaled[grepl("^[^.]*[.][^.]*$", fpkm_unscaled$condition),] 
#fpkm_unscaled_a <- within(fpkm_unscaled_a, condition<-data.frame(do.call('rbind', strsplit(as.character(condition), '.', fixed=TRUE)))) %>% as.data.frame()
#fpkm_unscaled_a %>% rename("geno" = "condition.X1")
#fpkm_unscaled_a %>% head()
#separate(data = fpkm_unscaled_a, col = condition, sep='.', into=c("geno", "rep")) 

### Calculate TPM
tpm <- fpkm_unscaled %>% group_by(geno) %>% summarize(total_fpkm = sum(fpkm)) %>% 
  full_join(fpkm_unscaled) %>% mutate(tpm = (fpkm / total_fpkm) *1000000) %>% 
  dplyr::select(GeneID,geno,tpm)

#tpm.mean <- tpm %>% separate(condition,sep="_2",into=c("condition","rep")) %>% select(-rep) %>% group_by(condition,GeneID) %>% summarize(mean_tpm = mean(tpm))
#tpm.table <- tpm %>% group_by(condition,GeneID) %>% dplyr::select(GeneID,condition,tpm)%>% spread(condition,tpm)
#tpm.mean.table <- tpm.mean %>% group_by(condition,GeneID) %>% dplyr::select(GeneID,condition,mean_tpm)%>% spread(condition,mean_tpm)

df <- tpm %>% filter(geno == "OMPGFP") %>%
  filter(GeneID %in% Olfr) %>% as.data.frame()
write.table(df, file = "/media/storageE/ariel/R/finalpaper/rna.knockouts/or_expression.txt", sep = "\t")

#### Ariel

"OMPtTA_tetOP2"
"OMPtTA_tetOP2_KO"
"OMPtTA_tetOM71_KO"
"OMPGFP"

resLFC <- lfcShrink(dds, contrast = c("geno", "OMPtTA_tetOGFP", "OMPGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "OMPtTA_tetOGFP \n vs. OMPGFP")

resLFC <- lfcShrink(dds, contrast = c("geno", "OMPtTA_tetOP2", "OMPGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "OMPtTA_tetOP2 \n vs. OMPGFP")

resLFC <- lfcShrink(dds, contrast = c("geno", "OMPtTA_tetOP2_KO", "OMPGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "OMPtTA_tetOP2_KO \n vs. OMPGFP")

resLFC <- lfcShrink(dds, contrast = c("geno", "OMPtTA_tetOM71_KO", "OMPGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "OMPtTA_tetOM71_KO \n vs. OMPGFP")






resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2_KO", "gg8tta_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2_KO vs. gg8tTA_tetOGFP") 

resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2_KO", "gg8tTA_tetOP2tg"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2_KO vs. gg8tTA_tetOP2tg", ylim = c(-7,7))

resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2tg", "gg8tta_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% Olfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2tg vs. gg8tta_tetOGFP") 

resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2_KO", "gg8tta_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% nonOlfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2_KO vs. gg8tTA_tetOGFP")

resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2_KO", "gg8tTA_tetOP2tg"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% nonOlfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2_KO vs. gg8tTA_tetOP2tg")

resLFC <- lfcShrink(dds, contrast = c("geno", "gg8tTA_tetOP2tg", "gg8tta_tetOGFP"), type="normal")
resLFC <- resLFC[row.names(resLFC) %in% nonOlfr,]
plotMA(resLFC, main = "gg8tTA_tetOP2tg vs. gg8tta_tetOGFP")
