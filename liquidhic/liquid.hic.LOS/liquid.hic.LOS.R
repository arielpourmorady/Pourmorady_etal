library(data.table)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(ggpubr)
library(vctrs)
library(tidyr)
library(stringr)
library(readr) # read_delim
library(tibble)
library(cowplot)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(R.utils)
library(populationPDXdesign)
library(multiHiCcompare)
library(patchwork)

# set directory to current directory
# this directory should contain all files in the GitHub folder

working_directory <- "/media/storageE/ariel/R/finalpaper_August2023/liquidhic/liquid.hic.LOS/"
setwd(working_directory) 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

## Script to extract eigen score ##
#java -jar juicer_tools.jar eigenvector <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
# execute the above script in bash

## Juicer Merged Files with KR Normalization ##
options(scipen = 999)

dump_dir <- "/media/storageD/distiller/dump/"
bin_size <- 50000

## FUNCTIONS ##
format_trans_path <- function(name){
  path <- str_c(dump_dir,"trans.",name,".",bin_size,".txt.gz")
  return(path)
}
format_cis_path <- function(name){
  path <- str_c(dump_dir,"cis.",name,".",bin_size,".txt.gz")
  return(path)
}

# Loading .txt.gz file from dump
gg8.p2.ctrl.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.merge"
gg8.p2.ctrl.trans_hic_path <-format_trans_path(gg8.p2.ctrl.name)
gg8.p2.ctrl.cis_hic_path <-format_cis_path(gg8.p2.ctrl.name)
gg8.p2.ctrl.total_count <- 362670259

gg8.p2.5m.name <- "gg8tTA.tetOP2.het.control.5minPreDig.LiquidHiC.merge"
gg8.p2.5m.trans_hic_path <-format_trans_path(gg8.p2.5m.name)
gg8.p2.5m.cis_hic_path <-format_cis_path(gg8.p2.5m.name)
gg8.p2.5m.total_count <- 266530739

gg8.p2.30m.name <- "gg8tTA.tetOP2.het.control.30minPreDig.Liquid.HiC.210807.merge"
gg8.p2.30m.trans_hic_path <-format_trans_path(gg8.p2.30m.name)
gg8.p2.30m.cis_hic_path <-format_cis_path(gg8.p2.30m.name)
gg8.p2.30m.total_count <- 170046887

gg8.p2.60m.name <- "gg8tTA.tetOP2.het.control.60minPreDig.Liquid.HiC.210807.merge"
gg8.p2.60m.trans_hic_path <-format_trans_path(gg8.p2.60m.name)
gg8.p2.60m.cis_hic_path <-format_cis_path(gg8.p2.60m.name)
gg8.p2.60m.total_count <- 200315002

gg8tta.tetop2.h27ac.HiChIP.IP.name <- "gg8tTAtetOP2.H3K27ac.HiChIP.AP151.AP152_230921"
gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path <-format_trans_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path <-format_cis_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.total_count <- 298853841

## ESTABLISH BINS ##

bed_to_bait <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(bait = str_c(chrNo,"_",V2)) %>% select(bait)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

bed_to_prey <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(prey = str_c(chrNo,"_",V2)) %>% select(prey)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

All_Bins_bait <- bed_to_bait("mm10_assembled.50kb.bed")
OR_Clusters_50kb_bait <- bed_to_bait("ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% 
  mutate(type="OR_Clusters")
Islands_50kb_bait <- bed_to_bait("Greek_Islands.50kb.mm10.bed") %>% 
  mutate(type="Islands")

All_Bins_prey <- bed_to_prey("mm10_assembled.50kb.bed")
OR_Clusters_50kb_prey <- bed_to_prey("ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% 
  mutate(type="OR_Clusters")
Islands_50kb_prey <- bed_to_prey("Greek_Islands.50kb.mm10.bed") %>% 
  mutate(type="Islands")


# All Bins ... Juicer KR Normalized Merged Files


libraries <- c("gg8.p2.ctrl",
               "gg8.p2.5m",
               "gg8.p2.30m",
               "gg8.p2.60m",
               "gg8.p2.hichip") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(gg8.p2.ctrl.trans_hic_path,
                    gg8.p2.5m.trans_hic_path,
                    gg8.p2.30m.trans_hic_path,
                    gg8.p2.60m.trans_hic_path,
                    gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path, , control.acid.r1.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.30.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(gg8.p2.ctrl.cis_hic_path,
                  gg8.p2.5m.cis_hic_path,
                  gg8.p2.30m.cis_hic_path,
                  gg8.p2.60m.cis_hic_path,
                  gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(gg8.p2.ctrl.total_count,
                 gg8.p2.5m.total_count,
                 gg8.p2.30m.total_count,
                 gg8.p2.60m.total_count,
                 gg8tta.tetop2.h27ac.HiChIP.IP.total_count)


###################################
#### These Hi-C Files are Large and
#### therefore need to be filtered
#### to load the whole thing into R
###################################

hic_contacts = data.frame()
for (i in c(1:5)) {
  cis_all = fread(cis_hic_path[i], col.names = c("bait", "prey", "contact")) 
  cis_short <- filter(cis_all, cis_all$bait == cis_all$prey)
  cis_short$arch = 'cis_short'
  cis_long <- filter(cis_all, cis_all$bait != cis_all$prey)
  cis_long$arch = 'cis_long'
  hic_contacts_loop <- rbind(cis_short, cis_long)
  hic_contacts_loop$geno <- print(libraries[i])
  hic_contacts_loop[is.na(hic_contacts_loop)] <- 0
  totalcount <- hic_contacts_loop %>% group_by(geno) %>% summarise(sum = sum(contact)) # You should activate this when you ARE NOT FILTERING
  hic_contacts_loop$norm <- hic_contacts_loop$contact/totalcount$sum # You should activate this when you ARE NOT FILTERING
  hic_contacts_loop <- select(hic_contacts_loop, -contact)
  hic_contacts <- rbind(hic_contacts, hic_contacts_loop)
  rm(cis_short, cis_long, cis_all, hic_contacts_loop)
}



# Eigen vs LOS gg8-tTA TETOP2 HET

## Relevant Bins

chr2_bins_bait <- dplyr::filter(All_Bins_bait, grepl("2_",bait))
chr2_bins_bait <- dplyr::filter(chr2_bins_bait, !grepl("12_",bait))
chr2_bins_bait <- chr2_bins_bait %>%
  tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) %>%
  mutate(bait_loc = as.numeric(bait_loc)) 
chr2_bins_bait_filter <- chr2_bins_bait %>%
  dplyr::filter(bait_loc > 3e6, bait_loc < (182100000 - 3e6))

chr2_bins_bait_filter_downsample <- sample_n(chr2_bins_bait_filter, 50)

LOS_function <- function(hic_data){
  master_df = data.frame()
  for(i in chr2_bins_bait_filter$bait){
    loc <- strsplit(i, split = "_")
    loc <- as.numeric(loc[[1]][2])
    df <- hic_data
    df <- df %>% filter(bait == i) %>%
      tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) %>%
      tidyr::separate('prey', into = c('prey_chr','prey_loc'), sep = '_', remove = F) %>%
      mutate(bait_loc = as.numeric(bait_loc), prey_loc = as.numeric(prey_loc))
    df_filt <- df %>%
      filter(prey_loc > (loc - 3e6), prey_loc < (loc + 3e6)) %>%
      group_by(bait) %>% summarise(near_norm = sum(norm)) %>%
      dplyr::mutate(fraction_norm = near_norm/sum(df$norm))
    master_df <- rbind(master_df, df_filt)
  }
  return(master_df)
}

###########
control <- hic_contacts %>% filter(geno == "gg8.p2.ctrl", bait %in% chr2_bins_bait$bait)
Five.Min <- hic_contacts %>% filter(geno == "gg8.p2.5m", bait %in% chr2_bins_bait$bait)
Thirty.Min <- hic_contacts %>% filter(geno == "gg8.p2.30m", bait %in% chr2_bins_bait$bait)
Sixty.Min <- hic_contacts %>% filter(geno == "gg8.p2.60m", bait %in% chr2_bins_bait$bait)

control.LOS <- LOS_function(control)
Five.Min.LOS <- LOS_function(Five.Min)
Thirty.Min.LOS <- LOS_function(Thirty.Min)
Sixty.Min.LOS <- LOS_function(Sixty.Min)

LOS.5min <- merge(control.LOS, Five.Min.LOS, by = "bait")
LOS.5min$LOS <- (LOS.5min$fraction_norm.x - LOS.5min$fraction_norm.y)/LOS.5min$fraction_norm.x

LOS.30min <- merge(control.LOS, Thirty.Min.LOS, by = "bait")
LOS.30min$LOS <- (LOS.30min$fraction_norm.x - LOS.30min$fraction_norm.y)/LOS.30min$fraction_norm.x

LOS.60min <- merge(control.LOS, Sixty.Min.LOS, by = "bait")
LOS.60min$LOS <- (LOS.60min$fraction_norm.x - LOS.60min$fraction_norm.y)/LOS.60min$fraction_norm.x

# HiChIP Enrichment

gg8.p2.hichip_df_chr2 <- hic_contacts %>%
  filter(geno == "gg8.p2.hichip") %>%
  group_by(bait) %>%
  summarise(norm = sum(norm)) %>%
  filter(bait %in% chr2_bins_bait$bait)

gg8.p2.ctrl_df_chr2 <- hic_contacts %>%
  filter(geno == "gg8.p2.ctrl") %>%
  group_by(bait) %>%
  summarise(norm = sum(norm)) %>%
  filter(bait %in% chr2_bins_bait$bait)

gg8.p2.norm <- left_join(gg8.p2.hichip_df_chr2, gg8.p2.ctrl_df_chr2, by = "bait") %>%
  mutate(foldchange = log(norm.x/norm.y))
#gg8.p2.norm <- gg8.p2.norm %>% tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) 

# EigenValues

eigen <- fread("/media/storageE/ariel/HiC/eigen/gg8tTA.tetOP2.het.ImmediateFix.Liquid.HiC.210224.AP29.50000.compartmentscore.cis.vecs.tsv",
               header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
eigen$bait <- paste(eigen$chrom, eigen$start, sep = "_")
eigen$eigen <- eigen$E1 %>% as.numeric() 
eigen$eigen <- eigen$eigen*-1
eigen <- eigen %>% select(bait, eigen) %>% as.data.frame()

OR_Clusters_50kb_bait$OR <- OR_Clusters_50kb_bait$type
OR_Clusters <- OR_Clusters_50kb_bait %>% select(bait, OR)

scatter.5min <- merge(LOS.5min, eigen, by = "bait")
scatter.5min$type <- "5-min"
scatter.5min <- scatter.5min %>% tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) 
scatter.5min <- merge(scatter.5min, OR_Clusters, by = "bait", all.x = TRUE)

scatter.30min <- merge(LOS.30min, eigen, by = "bait")
scatter.30min$type <- "30-min"
scatter.30min <- scatter.30min %>% tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) 
scatter.30min <- merge(scatter.30min, OR_Clusters, by = "bait", all.x = TRUE)

scatter.60min <- merge(LOS.60min, eigen, by = "bait")
scatter.60min$type <- "60-min"
scatter.60min <- scatter.60min %>% tidyr::separate('bait', into = c('bait_chr','bait_loc'), sep = '_', remove = F) 
scatter.60min <- merge(scatter.60min, OR_Clusters, by = "bait", all.x = TRUE)

scatter <- rbind(scatter.5min, scatter.30min, scatter.60min)

df <- left_join(scatter.5min, gg8.p2.norm, by = "bait")

#################################################
#################################################
### Correlation Plots ###
#################################################
#################################################

a <- ggplot(scatter.5min %>% sample_n(1000), aes(eigen, LOS)) + 
  ylim(-0.2, 0.6) + 
  ggtitle("5-min Liquid Hi-C") + 
  geom_point(alpha = 0.5) + 
  theme_classic()

b <- ggplot(scatter.30min%>% sample_n(1000), aes(eigen, LOS)) + 
  ylim(-0.2, 0.6) + 
  ggtitle("30-min Liquid Hi-C") + 
  geom_point(alpha = 0.5)+ 
  theme_classic()

c <- ggplot(scatter.60min%>% sample_n(1000), aes(eigen, LOS)) + 
  ylim(-0.2, 0.6) + 
  ggtitle("60-min Liquid Hi-C") + 
  geom_point(alpha = 0.5)+ 
  theme_classic()

d <- ggplot(df%>% sample_n(1000), aes(eigen, foldchange)) + 
  geom_point(alpha = 0.5)+ 
  ggtitle("H3K27ac HiChIP") + 
  theme_classic()

a + b + c + d + plot_layout(ncol = 4)


#################################################
#################################################
### Formatting ###
#################################################
#################################################

scatter_format.5min <- scatter.5min %>% sample_n(1000) %>% select(eigen, LOS) %>% mutate(condition = '5-min Pre-Dig')
scatter_format.30min <- scatter.30min %>% sample_n(1000) %>% select(eigen, LOS) %>% mutate(condition = '30-min Pre-Dig')
scatter_format.60min <- scatter.60min %>% sample_n(1000) %>% select(eigen, LOS) %>% mutate(condition = '60-min Pre-Dig')

scatter_format <- rbind(scatter_format.5min, scatter_format.30min, scatter_format.60min)
hichip_format <- df %>% sample_n(1000) %>% select(eigen, foldchange) %>% mutate(condition = 'hichip')

write.table(scatter_format, file = "/media/storageE/ariel/R/finalpaper_August2023/liquidhic/liquid.hic.LOS/SuppFig7b_LOS.txt", sep = '\t')
write.table(hichip_format, file = "/media/storageE/ariel/R/finalpaper_August2023/liquidhic/liquid.hic.LOS/SuppFig7c_HiChIPFoldChange.txt", sep = '\t')

