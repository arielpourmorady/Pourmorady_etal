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

# set directory to current directory
# this directory should contain all files in the GitHub folder

working_directory <- "/media/storageE/ariel/R/finalpaper_August2023/liquidhic/liquid.hic.heatmaps.replicates/"
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

dump_dir <- working_directory
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
gg8.p2.ctrl.r1.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.r0"
gg8.p2.ctrl.r1.trans_hic_path <-format_trans_path(gg8.p2.ctrl.r1.name)
gg8.p2.ctrl.r1.cis_hic_path <-format_cis_path(gg8.p2.ctrl.r1.name)
gg8.p2.ctrl.r1.total_count <- 99557522

gg8.p2.ctrl.r2.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.r1"
gg8.p2.ctrl.r2.trans_hic_path <-format_trans_path(gg8.p2.ctrl.r2.name)
gg8.p2.ctrl.r2.cis_hic_path <-format_cis_path(gg8.p2.ctrl.r2.name)
gg8.p2.ctrl.r2.total_count <- 26038914

gg8.p2.ctrl.r3.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.r2"
gg8.p2.ctrl.r3.trans_hic_path <-format_trans_path(gg8.p2.ctrl.r3.name)
gg8.p2.ctrl.r3.cis_hic_path <-format_cis_path(gg8.p2.ctrl.r3.name)
gg8.p2.ctrl.r3.total_count <- 237073823

gg8.p2.5m.r1.name <- "gg8tTA.tetOP2.het.control.5minPreDig.LiquidHiC.r0"
gg8.p2.5m.r1.trans_hic_path <-format_trans_path(gg8.p2.5m.r1.name)
gg8.p2.5m.r1.cis_hic_path <-format_cis_path(gg8.p2.5m.r1.name)
gg8.p2.5m.r1.total_count <- 22808819

gg8.p2.5m.r2.name <- "gg8tTA.tetOP2.het.control.5minPreDig.LiquidHiC.r1"
gg8.p2.5m.r2.trans_hic_path <-format_trans_path(gg8.p2.5m.r2.name)
gg8.p2.5m.r2.cis_hic_path <-format_cis_path(gg8.p2.5m.r2.name)
gg8.p2.5m.r2.total_count <- 133517472

gg8.p2.5m.r3.name <- "gg8tTA.tetOP2.het.control.5minPreDig.LiquidHiC.r2"
gg8.p2.5m.r3.trans_hic_path <-format_trans_path(gg8.p2.5m.r3.name)
gg8.p2.5m.r3.cis_hic_path <-format_cis_path(gg8.p2.5m.r3.name)
gg8.p2.5m.r3.total_count <- 110204448

gg8.p2.30m.r1.name <- "gg8tTA.tetOP2.het.control.30minPreDig.Liquid.HiC.210807.AP66"
gg8.p2.30m.r1.trans_hic_path <-format_trans_path(gg8.p2.30m.r1.name)
gg8.p2.30m.r1.cis_hic_path <-format_cis_path(gg8.p2.30m.r1.name)
gg8.p2.30m.r1.total_count <- 45376915

gg8.p2.30m.r2.name <- "gg8tTA.tetOP2.het.control.30minPreDig.Liquid.HiC.210807.AP67"
gg8.p2.30m.r2.trans_hic_path <-format_trans_path(gg8.p2.30m.r2.name)
gg8.p2.30m.r2.cis_hic_path <-format_cis_path(gg8.p2.30m.r2.name)
gg8.p2.30m.r2.total_count <- 66841694

gg8.p2.30m.r3.name <- "gg8tTA.tetOP2.het.control.30minPreDig.Liquid.HiC.210807.AP68"
gg8.p2.30m.r3.trans_hic_path <-format_trans_path(gg8.p2.30m.r3.name)
gg8.p2.30m.r3.cis_hic_path <-format_cis_path(gg8.p2.30m.r3.name)
gg8.p2.30m.r3.total_count <- 57828278

gg8.p2.60m.r1.name <- "gg8tTA.tetOP2.het.control.60minPreDig.Liquid.HiC.210807.AP69"
gg8.p2.60m.r1.trans_hic_path <-format_trans_path(gg8.p2.60m.r1.name)
gg8.p2.60m.r1.cis_hic_path <-format_cis_path(gg8.p2.60m.r1.name)
gg8.p2.60m.r1.total_count <- 53985594

gg8.p2.60m.r2.name <- "gg8tTA.tetOP2.het.control.60minPreDig.Liquid.HiC.210807.AP70"
gg8.p2.60m.r2.trans_hic_path <-format_trans_path(gg8.p2.60m.r2.name)
gg8.p2.60m.r2.cis_hic_path <-format_cis_path(gg8.p2.60m.r2.name)
gg8.p2.60m.r2.total_count <- 62692481

gg8.p2.60m.r3.name <- "gg8tTA.tetOP2.het.control.60minPreDig.Liquid.HiC.210807.AP71"
gg8.p2.60m.r3.trans_hic_path <-format_trans_path(gg8.p2.60m.r3.name)
gg8.p2.60m.r3.cis_hic_path <-format_cis_path(gg8.p2.60m.r3.name)
gg8.p2.60m.r3.total_count <- 83636927


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


libraries <- c("gg8.p2.ctrl.r1",
               "gg8.p2.ctrl.r2",
               "gg8.p2.ctrl.r3",
               "gg8.p2.5m.r1",
               "gg8.p2.5m.r2",
               "gg8.p2.5m.r3",
               "gg8.p2.30m.r1",
               "gg8.p2.30m.r2",
               "gg8.p2.30m.r3",
               "gg8.p2.60m.r1",
               "gg8.p2.60m.r2",
               "gg8.p2.60m.r3") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(gg8.p2.ctrl.r1.trans_hic_path,
                    gg8.p2.ctrl.r2.trans_hic_path,
                    gg8.p2.ctrl.r3.trans_hic_path,
                    gg8.p2.5m.r1.trans_hic_path,
                    gg8.p2.5m.r2.trans_hic_path,
                    gg8.p2.5m.r3.trans_hic_path,
                    gg8.p2.30m.r1.trans_hic_path,
                    gg8.p2.30m.r2.trans_hic_path,
                    gg8.p2.30m.r3.trans_hic_path,
                    gg8.p2.60m.r1.trans_hic_path,
                    gg8.p2.60m.r2.trans_hic_path,
                    gg8.p2.60m.r3.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path.trans_hic_path, , control.acid.r1.trans_hic_path.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.3.trans_hic_path0.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(gg8.p2.ctrl.r1.cis_hic_path,
                  gg8.p2.ctrl.r2.cis_hic_path,
                  gg8.p2.ctrl.r3.cis_hic_path,
                  gg8.p2.5m.r1.cis_hic_path,
                  gg8.p2.5m.r2.cis_hic_path,
                  gg8.p2.5m.r3.cis_hic_path,
                  gg8.p2.30m.r1.cis_hic_path,
                  gg8.p2.30m.r2.cis_hic_path,
                  gg8.p2.30m.r3.cis_hic_path,
                  gg8.p2.60m.r1.cis_hic_path,
                  gg8.p2.60m.r2.cis_hic_path,
                  gg8.p2.60m.r3.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(gg8.p2.ctrl.r1.total_count,
                 gg8.p2.ctrl.r2.total_count,
                 gg8.p2.ctrl.r3.total_count,
                 gg8.p2.5m.r1.total_count,
                 gg8.p2.5m.r2.total_count,
                 gg8.p2.5m.r3.total_count,
                 gg8.p2.30m.r1.total_count,
                 gg8.p2.30m.r2.total_count,
                 gg8.p2.30m.r3.total_count,
                 gg8.p2.60m.r1.total_count,
                 gg8.p2.60m.r2.total_count,
                 gg8.p2.60m.r3.total_count)

#############################
### MAKING 2 X 2 HEATMAPS ###
#############################
###################################
#### These Hi-C Files are Large and
#### therefore need to be filtered
#### to load the whole thing into R
###################################
# Two MB Island Filter

options(scipen = 999)
All_Bins_prey <- bed_to_prey("mm10_assembled.50kb.bed")

Islands_50kb <- read.table("Greek_Islands.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, type = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  select(prey, type) %>%
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique()

Bins.And.Islands <- left_join(All_Bins_prey, Islands_50kb)
Bins.And.Islands <- Bins.And.Islands %>% tidyr::separate('prey', into = c('prey_chr','prey_loc'), sep = '_')  %>% #split prey into chr and loc
  mutate_at(.vars = vars(matches('loc')), .funs = as.numeric) #convert "loc" column to numeric
Islands <- Bins.And.Islands[complete.cases(Bins.And.Islands),]

## Filter out Islands within 2 MB of an Island ##

TwoMB.Islands = data.frame()
for (i in seq_along(Islands[,1])) {
  df <- Bins.And.Islands %>% filter(Bins.And.Islands$prey_chr == Islands$prey_chr[i]) %>%
    filter((prey_loc > (Islands$prey_loc[i] - 1e6)) & (prey_loc < (Islands$prey_loc[i] + 1e6)))
  df$name <- Islands$type[i]
  df$order <- c(1:nrow(df))
  TwoMB.Islands <- rbind(TwoMB.Islands, df)
}

TwoMB.Islands$prey <- paste(TwoMB.Islands$prey_chr, TwoMB.Islands$prey_loc, sep = "_")
TwoMB.Islands <- TwoMB.Islands %>% dplyr::select(prey_chr, prey_loc, type, name, order, prey)

# Confirm that all Island 2-Mb window are within borders of genome
TwoMB.Islands %>% select(prey) %>% nrow()
TwoMB.Islands %>% dplyr::filter(prey %in% All_Bins_prey$prey) %>% select(prey) %>% nrow()


# ## ## ## ## ## ## 
## OR Insulation ##  
## ## ## ## ## ## # 

# Two MB OR Filter

OR_Clusters_50kb <- read.table("ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique()

Bins.And.ORs <- left_join(All_Bins_prey, OR_Clusters_50kb) 
Bins.And.ORs <- Bins.And.ORs %>% tidyr::separate('prey', into = c('prey_chr','prey_loc'), sep = '_')  %>% #split prey into chr and loc
  mutate_at(.vars = vars(matches('loc')), .funs = as.numeric) #convert "loc" column to numeric
ORs <- Bins.And.ORs[complete.cases(Bins.And.ORs),]

## Filter out ORs within 2 MB of an OR ##

TwoMB.ORs = data.frame()
for (i in seq_along(ORs[,1])) {
  df <- Bins.And.ORs %>% filter(Bins.And.ORs$prey_chr == ORs$prey_chr[i]) %>%
    filter((prey_loc > (ORs$prey_loc[i] - 1e6)) & (prey_loc < (ORs$prey_loc[i] + 1e6)))
  df$name <- ORs$type[i]
  df$order <- c(1:nrow(df))
  TwoMB.ORs <- rbind(TwoMB.ORs, df)
}


TwoMB.ORs$prey <- paste(TwoMB.ORs$prey_chr, TwoMB.ORs$prey_loc, sep = "_")
TwoMB.ORs <- TwoMB.ORs %>% dplyr::select(prey_chr, prey_loc, type, name, order, prey)

# Confirm that all Island 2-Mb window are within borders of genome
TwoMB.ORs %>% select(prey) %>% nrow()
TwoMB.ORs %>% dplyr::filter(prey %in% All_Bins_prey$prey) %>% select(prey) %>% nrow()

### The filter that I am creating includes all the OR Clusters as well as 2-MB around each Island

filter <- rbind(TwoMB.Islands, TwoMB.ORs) %>% dplyr::select(prey) %>% unique()


hic_contacts = data.frame()
for (i in c(1:12)) {
  trans_all = fread(trans_hic_path[i], col.names = c("bait", "prey", "contact"))
  trans_all$arch = 'trans'
  trans_all <- trans_all[trans_all$bait %in% filter$prey,] # You should activate this when you *ARE* FILTERING
  trans_all <- trans_all[trans_all$prey %in% filter$prey,] # You should activate this when you *ARE* FILTERING
  cis_all = fread(cis_hic_path[i], col.names = c("bait", "prey", "contact")) 
  cis_all <- cis_all[cis_all$bait %in% filter$prey,] # You should activate this when you *ARE* FILTERING
  cis_all <- cis_all[cis_all$prey %in% filter$prey,] # You should activate this when you *ARE* FILTERING
  cis_short <- filter(cis_all, cis_all$bait == cis_all$prey)
  cis_short$arch = 'cis_short'
  cis_long <- filter(cis_all, cis_all$bait != cis_all$prey)
  cis_long$arch = 'cis_long'
  hic_contacts_loop <- rbind(trans_all, cis_short, cis_long)
  hic_contacts_loop$geno <- print(libraries[i])
  hic_contacts_loop[is.na(hic_contacts_loop)] <- 0
  #totalcount <- hic_contacts_loop %>% group_by(geno) %>% summarise(sum = sum(contact)) # You should activate this when you ARE NOT FILTERING
  #hic_contacts_loop$norm <- hic_contacts_loop$contact/totalcount$sum # You should activate this when you ARE NOT FILTERING
  hic_contacts_loop$norm <- hic_contacts_loop$contact/total_count[i] # You should activate this when you *ARE* FILTERING
  hic_contacts_loop <- select(hic_contacts_loop, -contact)
  hic_contacts <- rbind(hic_contacts, hic_contacts_loop)
  rm(trans_all, cis_short, cis_long, cis_all, hic_contacts_loop)
}

##################
### DATAFRAMES ###
##################

# Here we go
df_ORs <- TwoMB.ORs %>% 
  mutate("loc" = prey) %>%
  select(loc, order, name)

df_islands <- TwoMB.Islands %>% 
  mutate("loc" = prey, "name" = sub("\\_.*", "", name)) %>%
  select(loc, order, name)


#########################################################
### PERCENT DECREASE IN CONTACT SPECIFICITY FUNCTIONS ###
#########################################################

order_tbl_df <- data.frame(bait_order = 1:39) 
order_tbl_df <- expand.grid(bait_order = order_tbl_df$bait_order, prey_order = order_tbl_df$bait_order) %>% as_tibble()

P2 <- "7_107050000"
P2_cluster <- "OR-Cluster-26_Olfr714,OR-Cluster-26_Olfr17"

activeP2_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
  df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
    select(bait, prey, geno, norm) # filter all contacts made in one cell
  df3 <- df_ORs %>% filter(name == P2_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df6 <- left_join(df5, df_islands, by = c("prey" = "loc")) %>%
    mutate("prey_OR" = name, "prey_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
  df7 <- df6[complete.cases(df6),]
  df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
    summarise(norm = sum(norm))
  df9 <- left_join(order_tbl_df, df8)
  df9 <- df9 %>% replace(is.na(.), 0)
  df9$norm <- df9$norm/sum(df9$norm)
  return(df9) 
}

inactiveP2_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
  df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
    select(bait, prey, geno, norm) # filter all contacts made in one cell
  df3 <- df_ORs %>% filter(name != P2_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df6 <- left_join(df5, df_islands, by = c("prey" = "loc")) %>%
    mutate("prey_OR" = name, "prey_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
  df7 <- df6[complete.cases(df6),]
  df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
    summarise(norm = sum(norm))
  df9 <- left_join(order_tbl_df, df8)
  df9 <- df9 %>% replace(is.na(.), 0)
  df9$norm <- df9$norm/sum(df9$norm)
  return(df9) 
}

#####################################
### CHANGE IN CONTACT SPECIFICITY ###
#####################################

# Rep 1

activeP2_to_enhancer.gg8.p2.ctrl.r1 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r1", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.5m.r1 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.30m.r1 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.60m.r1 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

activeP2_to_enhancer.gg8.p2.r1 <- as.data.frame(rbind(activeP2_to_enhancer.gg8.p2.ctrl.r1, 
                                                   activeP2_to_enhancer.gg8.p2.5m.r1,
                                                   activeP2_to_enhancer.gg8.p2.30m.r1,
                                                   activeP2_to_enhancer.gg8.p2.60m.r1)) %>%
  mutate(specificity = norm - activeP2_to_enhancer.gg8.p2.ctrl.r1$norm) %>%
  mutate(specificity = 100*specificity/activeP2_to_enhancer.gg8.p2.ctrl.r1$norm) %>% 
  select(norm, specificity)
activeP2_to_enhancer.gg8.p2.r1$type <- "active OR"
activeP2_to_enhancer.gg8.p2.r1$time <- c(0, 5, 30, 60)

inactiveP2_to_enhancer.gg8.p2.ctrl.r1 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r1", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.5m.r1 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.30m.r1 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.60m.r1 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r1", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

inactiveP2_to_enhancer.gg8.p2.r1 <- as.data.frame(rbind(inactiveP2_to_enhancer.gg8.p2.ctrl.r1, 
                                                     inactiveP2_to_enhancer.gg8.p2.5m.r1,
                                                     inactiveP2_to_enhancer.gg8.p2.30m.r1,
                                                     inactiveP2_to_enhancer.gg8.p2.60m.r1)) %>%
  mutate(specificity = norm - inactiveP2_to_enhancer.gg8.p2.ctrl.r1$norm) %>%
  mutate(specificity = 100*specificity/inactiveP2_to_enhancer.gg8.p2.ctrl.r1$norm) %>% 
  select(norm, specificity)
inactiveP2_to_enhancer.gg8.p2.r1$type <- "inactive OR"
inactiveP2_to_enhancer.gg8.p2.r1$time <- c(0, 5, 30, 60)


enhancer.gg8.p2.r1 <- rbind(activeP2_to_enhancer.gg8.p2.r1, inactiveP2_to_enhancer.gg8.p2.r1)

# Rep 2

activeP2_to_enhancer.gg8.p2.ctrl.r2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r2", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.5m.r2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.30m.r2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.60m.r2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

activeP2_to_enhancer.gg8.p2.r2 <- as.data.frame(rbind(activeP2_to_enhancer.gg8.p2.ctrl.r2, 
                                                      activeP2_to_enhancer.gg8.p2.5m.r2,
                                                      activeP2_to_enhancer.gg8.p2.30m.r2,
                                                      activeP2_to_enhancer.gg8.p2.60m.r2)) %>%
  mutate(specificity = norm - activeP2_to_enhancer.gg8.p2.ctrl.r2$norm) %>%
  mutate(specificity = 100*specificity/activeP2_to_enhancer.gg8.p2.ctrl.r2$norm) %>% 
  select(norm, specificity)
activeP2_to_enhancer.gg8.p2.r2$type <- "active OR"
activeP2_to_enhancer.gg8.p2.r2$time <- c(0, 5, 30, 60)

inactiveP2_to_enhancer.gg8.p2.ctrl.r2 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r2", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.5m.r2 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.30m.r2 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.60m.r2 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r2", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

inactiveP2_to_enhancer.gg8.p2.r2 <- as.data.frame(rbind(inactiveP2_to_enhancer.gg8.p2.ctrl.r2, 
                                                        inactiveP2_to_enhancer.gg8.p2.5m.r2,
                                                        inactiveP2_to_enhancer.gg8.p2.30m.r2,
                                                        inactiveP2_to_enhancer.gg8.p2.60m.r2)) %>%
  mutate(specificity = norm - inactiveP2_to_enhancer.gg8.p2.ctrl.r2$norm) %>%
  mutate(specificity = 100*specificity/inactiveP2_to_enhancer.gg8.p2.ctrl.r2$norm) %>% 
  select(norm, specificity)
inactiveP2_to_enhancer.gg8.p2.r2$type <- "inactive OR"
inactiveP2_to_enhancer.gg8.p2.r2$time <- c(0, 5, 30, 60)


enhancer.gg8.p2.r2 <- rbind( activeP2_to_enhancer.gg8.p2.r2, inactiveP2_to_enhancer.gg8.p2.r2)

# Rep 3

activeP2_to_enhancer.gg8.p2.ctrl.r3 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r3", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.5m.r3 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.30m.r3 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.60m.r3 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

activeP2_to_enhancer.gg8.p2.r3 <- as.data.frame(rbind(activeP2_to_enhancer.gg8.p2.ctrl.r3, 
                                                      activeP2_to_enhancer.gg8.p2.5m.r3,
                                                      activeP2_to_enhancer.gg8.p2.30m.r3,
                                                      activeP2_to_enhancer.gg8.p2.60m.r3)) %>%
  mutate(specificity = norm - activeP2_to_enhancer.gg8.p2.ctrl.r3$norm) %>%
  mutate(specificity = 100*specificity/activeP2_to_enhancer.gg8.p2.ctrl.r3$norm) %>% 
  select(norm, specificity)
activeP2_to_enhancer.gg8.p2.r3$type <- "active OR"
activeP2_to_enhancer.gg8.p2.r3$time <- c(0, 5, 30, 60)

inactiveP2_to_enhancer.gg8.p2.ctrl.r3 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.r3", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.5m.r3 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.30m.r3 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.60m.r3 <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m.r3", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

inactiveP2_to_enhancer.gg8.p2.r3 <- as.data.frame(rbind(inactiveP2_to_enhancer.gg8.p2.ctrl.r3, 
                                                        inactiveP2_to_enhancer.gg8.p2.5m.r3,
                                                        inactiveP2_to_enhancer.gg8.p2.30m.r3,
                                                        inactiveP2_to_enhancer.gg8.p2.60m.r3)) %>%
  mutate(specificity = norm - inactiveP2_to_enhancer.gg8.p2.ctrl.r3$norm) %>%
  mutate(specificity = 100*specificity/inactiveP2_to_enhancer.gg8.p2.ctrl.r3$norm) %>% 
  select(norm, specificity)
inactiveP2_to_enhancer.gg8.p2.r3$type <- "inactive OR"
inactiveP2_to_enhancer.gg8.p2.r3$time <- c(0, 5, 30, 60)


enhancer.gg8.p2.r3 <- rbind( activeP2_to_enhancer.gg8.p2.r3, inactiveP2_to_enhancer.gg8.p2.r3)




# Final Rep

enhancer.gg8.p2.r1$rep <- "r1"
enhancer.gg8.p2.r2$rep <- "r2"
enhancer.gg8.p2.r3$rep <- "r3"


enhancer.gg8.p2 <- rbind(enhancer.gg8.p2.r1, enhancer.gg8.p2.r2, enhancer.gg8.p2.r3)

enhancer.gg8.p2_mean_se_lineplot <- enhancer.gg8.p2 %>% 
  group_by(type, time) %>%
  filter(type != "islands") %>%
  summarise(mean = mean(specificity),
            se = sd(specificity)/sqrt(length(specificity)))

t_test_fxn <- function(df, num) {
  v1 <- df %>% filter(type == "active OR", time == num) %>% dplyr::select(specificity)
  v2 <- df %>% filter(type == "inactive OR", time == num) %>% dplyr::select(specificity)
  p <- t.test(v1, v2)$p.value
  return(p)
}

for (i in c(5, 30, 60)) {
  p <- t_test_fxn(enhancer.gg8.p2, i)
  print(p)
}

ggplot(enhancer.gg8.p2_mean_se_lineplot, aes(x = time, y = mean, colour = type)) + geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=4,  colour="black") + 
  geom_line( linewidth = 0.75) +
  geom_point( size = 2.5) + 
  scale_colour_manual(values=c("#409b73", "#b10101")) +
  ggtitle("% Decrease in trans Contacts \n square Norm") + 
  theme_classic()

####### Values

enhancer.gg8.p2 %>% group_by(time, type) %>%
  summarise(mean = mean(specificity),
            sd = sd(specificity))

write.table(enhancer.gg8.p2, file = "/media/storageE/ariel/R/finalpaper_August2023/liquidhic/liquid.hic.heatmaps.replicates/Fig3b_LiquidHiCChange.txt", sep = '\t')

