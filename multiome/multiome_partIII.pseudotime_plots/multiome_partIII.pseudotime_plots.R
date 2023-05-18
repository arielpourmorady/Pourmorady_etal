library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyverse)
library(cowplot)
library(reshape2)

set.seed(1234)
options(scipen=999)

fragpath <- "/data/outs/atac_fragments.tsv.gz"
fragments <- CreateFragmentObject(fragpath)

neurons <- readRDS("/data/finalpaper/multiome_partII.pseudotime_projection/neurons.rds")

#################################################################################
#################################################################################
#################################################################################
#################################################################################


##################################################################################
##################################################################################  
############## Adding Maximal OR Expression per Cell to Metadata ################# 
##################################################################################
##################################################################################

# You are writing code to add top expressed OR to metadata

olfr_anno <- read.table("/data/annotations/Olfr-CDS-ucsc.bed") %>%
  as.data.frame() %>%
  dplyr::select(V4, V1) %>%
  dplyr::rename("olfr" = "V4", "chr" = "V1")

rna <- GetAssayData(object = neurons, 
                    assay = "SCT",
                    slot = "data") 

olfrs <- rna %>% 
  row.names() %>% 
  as.data.frame() %>% 
  dplyr::filter(grepl("Olfr", .))  

# Extracting the most highly expressed Olfr per Cell
df <- rna %>% as.data.frame() 
maxOR_exp_percell <- df %>% mutate(genes = row.names(df)) %>% melt() %>% 
  dplyr::filter(genes %in% olfrs$.) %>%
  group_by(variable) %>%
  dplyr::filter(value == max(value)) %>%
  group_by(variable) %>%
  summarise(value = unique(value))
row.names(maxOR_exp_percell) <- maxOR_exp_percell$variable
maxOR_exp_percell <- maxOR_exp_percell %>% dplyr::select(value)

neurons <- AddMetaData(object = neurons,
                   metadata = maxOR_exp_percell$value,
                   col.name = "ORtranscript")


#########################################################################################
######################################################################################### 
###### Adding Greek Island/Lhx2_Ebf1/mOSN_cCRE Accessibility per Cell to Metadata #######
#########################################################################################
#########################################################################################

lhx2_ebf_bed <- as.data.frame(read.table("/data/annotations/OMP-Ebf+Lhx2.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

gi_bed <- as.data.frame(read.table("/data/multiome/Greek_Islands.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

mOSN_cCREs_bed <- as.data.frame(read.table("/data/finalpaper/multiome_partI.input_processing/mOSN_cCREs.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE))

extract_peaks <- function(bed){
  df <- FeatureMatrix(
    fragments,
    features = bed[,1],
    cells = Cells(neurons),
    process_n = 2000,
    sep = c("-", "-"),
    verbose = TRUE) %>% as.data.frame()
  df$loc <- rownames(df)
  df <- melt(df) %>%
    dplyr::filter(value > 0) %>%
    group_by(variable) %>%
    summarise(Sum_Acc_percell = sum(value),
              Count = n()) %>%
    mutate(Acc_perEnh_percell = Sum_Acc_percell/Count) %>%
    dplyr::rename("cell" = "variable")
  neurons_cells <- as.data.frame(Cells(neurons)) %>% dplyr::rename("cell" = "Cells(neurons)")
  df <- left_join(neurons_cells, df)
  df[is.na(df)] <- 0
  row.names(df) <-df$cell
  return(df)
}

gi_peaks <- extract_peaks(gi_bed)
lhx2_ebf_peaks <- extract_peaks(lhx2_ebf_bed)
mOSN_cCREs_peaks <- extract_peaks(mOSN_cCREs_bed)

neurons <- AddMetaData(
  object = neurons,
  metadata = gi_peaks$Sum_Acc_percell,
  col.name = "Sum_gi_peaks"
)

neurons <- AddMetaData(
  object = neurons,
  metadata = gi_peaks$Acc_perEnh_percell,
  col.name = "Acc_perEnh_percell_gi_peaks"
)

neurons <- AddMetaData(
  object = neurons,
  metadata = lhx2_ebf_peaks$Sum_Acc_percell,
  col.name = "Sum_lhx2_ebf_peaks"
)

neurons <- AddMetaData(
  object = neurons,
  metadata = mOSN_cCREs_peaks$Sum_Acc_percell,
  col.name = "Sum_mOSN_cCREs_peaks"
)

neurons@meta.data$Sum_gi_peaks <- neurons@meta.data$Sum_gi_peaks/neurons@meta.data$nCount_ATAC
neurons@meta.data$Sum_lhx2_ebf_peaks <- neurons@meta.data$Sum_lhx2_ebf_peaks/neurons@meta.data$nCount_ATAC
neurons@meta.data$Sum_mOSN_cCREs_peaks <- neurons@meta.data$Sum_mOSN_cCREs_peaks/neurons@meta.data$nCount_ATAC

#########################################
#########################################
###### Pseudotime Plots for Paper #######
#########################################
#########################################


metadata <- data.frame(neurons@meta.data)


time <- metadata %>% dplyr::select(pseudotime, Sum_gi_peaks, Acc_perEnh_percell_gi_peaks, Sum_lhx2_ebf_peaks, Sum_mOSN_cCREs_peaks, ORtranscript)
time$pseudotime <- round(time$pseudotime)
time[is.na(time)] <- 0

df <- time %>% 
  group_by(pseudotime) %>% 
  summarise(
    mean = mean(Sum_gi_peaks),
    se = sqrt(var(Sum_gi_peaks) / length(Sum_gi_peaks)),
    mean_norm = mean(ORtranscript),
    se_norm = sqrt(var(ORtranscript) / length(ORtranscript))
  )

a <- ggplot(df, aes(x = pseudotime)) + 
  geom_line(aes(y = mean), colour = "black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour = "black") +
  geom_line(aes(y = mean_norm*0.0008222), colour = "blue3") +
  geom_errorbar(aes(ymin=(mean_norm*0.0008222)-(se_norm*0.0008222), ymax=(mean_norm*0.0008222)+(se_norm*0.0008222)), colour = "blue3") +
  scale_y_continuous(
    name = "Cummulative Greek Island Accessibility",
    sec.axis = sec_axis(trans=~./0.0008222, name="OR Expression")
  ) + 
  theme_classic() + theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "blue3")
  )
a


df <- time %>% 
  group_by(pseudotime) %>% 
  summarise(
    mean = mean(Sum_gi_peaks),
    se = sqrt(var(Sum_gi_peaks) / length(Sum_gi_peaks)),
    mean_norm = mean(Acc_perEnh_percell_gi_peaks),
    se_norm = sqrt(var(Acc_perEnh_percell_gi_peaks) / length(Acc_perEnh_percell_gi_peaks))
  )

b <- ggplot(df, aes(x = pseudotime)) + 
  geom_line(aes(y = mean), colour = "black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour = "black") +
  geom_line(aes(y = mean_norm*0.001208), colour = "springgreen4") +
  geom_errorbar(aes(ymin=(mean_norm*0.001208)-(se_norm*0.001208), ymax=(mean_norm*0.001208)+(se_norm*0.001208)), colour = "springgreen4") +
  scale_y_continuous(
    name = "Cummulative Greek Island Accessibility",
    sec.axis = sec_axis(trans=~./0.001208, name="Accessibility per Active Greek Island")
  ) + 
  theme_classic() + theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "springgreen4")
  )
b

df <- time %>% 
  group_by(pseudotime) %>% 
  summarise(
    mean = mean(Sum_gi_peaks),
    se = sqrt(var(Sum_gi_peaks) / length(Sum_gi_peaks)),
    mean_norm = mean(Sum_lhx2_ebf_peaks),
    se_norm = sqrt(var(Sum_lhx2_ebf_peaks) / length(Sum_lhx2_ebf_peaks))
  )

c <- ggplot(df, aes(x = pseudotime)) + 
  geom_line(aes(y = mean), colour = "black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour = "black") +
  geom_line(aes(y = mean_norm*0.02108), colour = "firebrick2") +
  geom_errorbar(aes(ymin=(mean_norm*0.02108)-(se_norm*0.02108), ymax=(mean_norm*0.02108)+(se_norm*0.02108)), colour = "firebrick2") +
  scale_y_continuous(
    name = "Cummulative Greek Island Accessibility",
    sec.axis = sec_axis(trans=~./0.1208, name="Cummulative Lhx2/EBF1 Accessibility")
  ) + 
  theme_classic() + theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "firebrick2")
  )
c

df <- time %>% 
  group_by(pseudotime) %>% 
  summarise(
    mean = mean(Sum_gi_peaks),
    se = sqrt(var(Sum_gi_peaks) / length(Sum_gi_peaks)),
    mean_norm = mean(Sum_mOSN_cCREs_peaks),
    se_norm = sqrt(var(Sum_mOSN_cCREs_peaks) / length(Sum_mOSN_cCREs_peaks))
  )

d <- ggplot(df, aes(x = pseudotime)) + 
  geom_line(aes(y = mean), colour = "black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour = "black") +
  geom_line(aes(y = mean_norm*0.3637101), colour = "maroon2") +
  geom_errorbar(aes(ymin=(mean_norm*0.3637101)-(se_norm*0.3637101), ymax=(mean_norm*0.3637101)+(se_norm*0.3637101)), colour = "maroon2") +
  scale_y_continuous(
    name = "Cummulative Greek Island Accessibility",
    sec.axis = sec_axis(trans=~./0.3637101, name="Cummulative top OSN Enhancer Accessibility")
  ) + 
  theme_classic() + theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "maroon2")
  )
d

d | c | a | b

#saveRDS(neurons, file = "/data/finalpaper/multiome_partIII.pseudotime_plots/neurons.rds")



