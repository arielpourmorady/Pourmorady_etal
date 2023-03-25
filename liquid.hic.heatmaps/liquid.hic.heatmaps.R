# Clear environment
list.of.packages <- c("vctrs", "tidyr", "tibble", "gridExtra", "populationPDXdesign", "multiHiCcompare")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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


## ESTABLISH BINS ##

bed_to_bait <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(bait = str_c(chrNo,"_",V2)) %>% select(bait)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

bed_to_prey <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(prey = str_c(chrNo,"_",V2)) %>% select(prey)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

All_Bins_bait <- bed_to_bait("/media/storageA/kevin/annotation/mm10_assembled.50kb.bed")
OR_Clusters_50kb_bait <- bed_to_bait("/media/storageA/kevin/annotation/ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% 
  mutate(type="OR_Clusters")
Islands_50kb_bait <- bed_to_bait("/media/storageA/kevin/annotation/Greek_Islands.50kb.mm10.bed") %>% 
  mutate(type="Islands")

All_Bins_prey <- bed_to_prey("/media/storageA/kevin/annotation/mm10_assembled.50kb.bed")
OR_Clusters_50kb_prey <- bed_to_prey("/media/storageA/kevin/annotation/ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% 
  mutate(type="OR_Clusters")
Islands_50kb_prey <- bed_to_prey("/media/storageA/kevin/annotation/Greek_Islands.50kb.mm10.bed") %>% 
  mutate(type="Islands")


# All Bins ... Juicer KR Normalized Merged Files


libraries <- c("gg8.p2.ctrl",
               "gg8.p2.5m",
               "gg8.p2.30m",
               "gg8.p2.60m") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(gg8.p2.ctrl.trans_hic_path,
                    gg8.p2.5m.trans_hic_path,
                    gg8.p2.30m.trans_hic_path,
                    gg8.p2.60m.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path, , control.acid.r1.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.30.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(gg8.p2.ctrl.cis_hic_path,
                  gg8.p2.5m.cis_hic_path,
                  gg8.p2.30m.cis_hic_path,
                  gg8.p2.60m.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(gg8.p2.ctrl.total_count,
                 gg8.p2.5m.total_count,
                 gg8.p2.30m.total_count,
                 gg8.p2.60m.total_count)


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
All_Bins_prey <- bed_to_prey("/media/storageA/kevin/annotation/mm10_assembled.50kb.bed")

Islands_50kb <- read.table("/media/storageE/ariel/annotations/Greek_Islands.bed") %>%
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

for (i in seq_along(Islands[,1])) {
  if (i == 1) {
    df <- Bins.And.Islands %>% filter(Bins.And.Islands$prey_chr == Islands$prey_chr[i]) %>%
      filter((prey_loc > (Islands$prey_loc[i] - 1e6)) & (prey_loc < (Islands$prey_loc[i] + 1e6)))
    df$name <- Islands$type[i]
    df$order <- c(1:nrow(df))
    TwoMB.Islands <- df
  } else {
    tryCatch({
      df.loop <- Bins.And.Islands %>% filter(Bins.And.Islands$prey_chr == Islands$prey_chr[i]) %>%
        filter((prey_loc > (Islands$prey_loc[i] - 1e6)) & (prey_loc < (Islands$prey_loc[i] + 1e6)))
      df.loop$name <- Islands$type[i]
      df.loop$order <- c(1:nrow(df.loop))
      TwoMB.Islands.loop <- df.loop
      TwoMB.Islands <- rbind(TwoMB.Islands, TwoMB.Islands.loop)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

TwoMB.Islands$prey <- paste(TwoMB.Islands$prey_chr, TwoMB.Islands$prey_loc, sep = "_")
TwoMB.Islands <- TwoMB.Islands %>% dplyr::select(prey_chr, prey_loc, type, name, order, prey)


# ## ## ## ## ## ## 
## OR Insulation ##  
## ## ## ## ## ## # 

# Two MB OR Filter

options(scipen = 999)
All_Bins_prey <- bed_to_prey("/media/storageA/kevin/annotation/mm10_assembled.50kb.bed")

OR_Clusters_50kb <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
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

for (i in seq_along(ORs[,1])) {
  if (i == 1) {
    df <- Bins.And.ORs %>% filter(Bins.And.ORs$prey_chr == ORs$prey_chr[i]) %>%
      filter((prey_loc > (ORs$prey_loc[i] - 1e6)) & (prey_loc < (ORs$prey_loc[i] + 1e6)))
    df$name <- ORs$type[i]
    df$order <- c(1:nrow(df))
    TwoMB.ORs <- df
  } else {
    tryCatch({
      df.loop <- Bins.And.ORs %>% filter(Bins.And.ORs$prey_chr == ORs$prey_chr[i]) %>%
        filter((prey_loc > (ORs$prey_loc[i] - 1e6)) & (prey_loc < (ORs$prey_loc[i] + 1e6)))
      df.loop$name <- ORs$type[i]
      df.loop$order <- c(1:nrow(df.loop))
      TwoMB.ORs.loop <- df.loop
      TwoMB.ORs <- rbind(TwoMB.ORs, TwoMB.ORs.loop)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

TwoMB.ORs$prey <- paste(TwoMB.ORs$prey_chr, TwoMB.ORs$prey_loc, sep = "_")
TwoMB.ORs <- TwoMB.ORs %>% dplyr::select(prey_chr, prey_loc, type, name, order, prey)

TwoMB.ORs$prey <- paste(TwoMB.ORs$prey_chr, TwoMB.ORs$prey_loc, sep = "_")
TwoMB.ORs <- TwoMB.ORs %>% dplyr::select(prey_chr, prey_loc, type, name, order, prey)

### The filter that I am creating includes all the OR Clusters as well as 2-MB around each Island

filter <- rbind(TwoMB.Islands, TwoMB.ORs) %>% dplyr::select(prey) %>% unique()

for (i in c(1:4)) {
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
  if (i == 1) {
    hic_contacts <- rbind(trans_all, cis_short, cis_long)
    hic_contacts$geno <- print(libraries[i])
    hic_contacts[is.na(hic_contacts)] <- 0
    #totalcount <- hic_contacts %>% group_by(geno) %>% summarise(sum = sum(contact)) # You should *DE-activate* this when you ARE FILTERING
    #hic_contacts$norm <- hic_contacts$contact/totalcount$sum # You should *DE-activate* this when you ARE FILTERING
    hic_contacts$norm <- hic_contacts$contact/total_count[i] # You should activate this when you *ARE* FILTERING
    hic_contacts <- select(hic_contacts, -contact)
    rm(trans_all, cis_short, cis_long, cis_all)
  } else {
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
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_ORs %>% filter(name == P2_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df6 <- left_join(df5, df_islands, by = c("prey" = "loc")) %>%
    mutate("prey_OR" = name, "prey_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
  df6 <- df6 %>% filter(bait_OR != prey_OR)
  df7 <- df6[complete.cases(df6),]
  df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
    summarise(norm = sum(norm))
  df9 <- left_join(order_tbl_df, df8)
  df9 <- df9 %>% replace(is.na(.), 0)
  if (normalization == "square"){
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9) 
  }
  if (normalization == "none"){
    return(df9) 
  }
  if (normalization == "row_bait"){
    df9 <- df9 %>%
      filter(bait_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}

inactiveP2_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans",  geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_ORs %>% filter(name != P2_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df6 <- left_join(df5, df_islands, by = c("prey" = "loc")) %>%
    mutate("prey_OR" = name, "prey_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
  df6 <- df6 %>% filter(bait_OR != prey_OR)
  df7 <- df6[complete.cases(df6),]
  df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
    summarise(norm = sum(norm))
  df9 <- left_join(order_tbl_df, df8)
  df9 <- df9 %>% replace(is.na(.), 0)
  if (normalization == "square"){
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9) 
  }
  if (normalization == "none"){
    return(df9) 
  }
  if (normalization == "row_bait"){
    df9 <- df9 %>%
      filter(bait_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}


enhancer_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_islands 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df6 <- left_join(df3, df5, by = c("loc" = "prey")) %>%
    dplyr::rename("prey" = "loc") %>%
    mutate("prey_OR" = name, "prey_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
  df6 <- df6 %>% filter(bait_OR != prey_OR)
  df7 <- df6[complete.cases(df6),]
  df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
    summarise(norm = sum(norm))
  df9 <- left_join(order_tbl_df, df8)
  df9 <- df9 %>% replace(is.na(.), 0)
  if (normalization == "square"){
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9) 
  }
  if (normalization == "none"){
    return(df9) 
  }
  if (normalization == "row_bait"){
    df9 <- df9 %>%
      filter(bait_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}


###########################
### MATRICES & HEATMAPS ###
###########################

# Generate Contact Specificity Score #
activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)

inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
##


enhancer_to_enhancer_matrix.gg8.p2.ctrl <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
enhancer_to_enhancer_matrix.gg8.p2.5m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
enhancer_to_enhancer_matrix.gg8.p2.30m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
enhancer_to_enhancer_matrix.gg8.p2.60m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(enhancer_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  max(enhancer_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  length.out=ceiling(paletteLength)))

a <- pheatmap(enhancer_to_enhancer_matrix.gg8.p2.ctrl[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

b <- pheatmap(enhancer_to_enhancer_matrix.gg8.p2.5m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

c <- pheatmap(enhancer_to_enhancer_matrix.gg8.p2.30m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

d <- pheatmap(enhancer_to_enhancer_matrix.gg8.p2.60m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

activeP2_to_enhancer_matrix.gg8.p2.ctrl <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
activeP2_to_enhancer_matrix.gg8.p2.5m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
activeP2_to_enhancer_matrix.gg8.p2.30m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
activeP2_to_enhancer_matrix.gg8.p2.60m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(activeP2_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  max(activeP2_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  length.out=ceiling(paletteLength)))

e <- pheatmap(activeP2_to_enhancer_matrix.gg8.p2.ctrl[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

f <- pheatmap(activeP2_to_enhancer_matrix.gg8.p2.5m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

g <- pheatmap(activeP2_to_enhancer_matrix.gg8.p2.30m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

h <- pheatmap(activeP2_to_enhancer_matrix.gg8.p2.60m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)



inactiveP2_to_enhancer_matrix.gg8.p2.ctrl <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
inactiveP2_to_enhancer_matrix.gg8.p2.5m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
inactiveP2_to_enhancer_matrix.gg8.p2.30m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
inactiveP2_to_enhancer_matrix.gg8.p2.60m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(inactiveP2_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  max(inactiveP2_to_enhancer_matrix.gg8.p2.ctrl[,-1]), 
                  length.out=ceiling(paletteLength)))

i <- pheatmap(inactiveP2_to_enhancer_matrix.gg8.p2.ctrl[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

j <- pheatmap(inactiveP2_to_enhancer_matrix.gg8.p2.5m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

k <- pheatmap(inactiveP2_to_enhancer_matrix.gg8.p2.30m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)

l <- pheatmap(inactiveP2_to_enhancer_matrix.gg8.p2.60m[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=2.5, cellwidth = 2.5,
              breaks = myBreaks)




grid.arrange(a[[4]], b[[4]], c[[4]], d[[4]], 
             e[[4]], f[[4]], g[[4]], h[[4]],
             i[[4]], j[[4]], k[[4]], l[[4]], ncol = 4, nrow = 3)

m <- arrangeGrob(a[[4]], b[[4]], c[[4]], d[[4]], 
                  e[[4]], f[[4]], g[[4]], h[[4]],
                  i[[4]], j[[4]], k[[4]], l[[4]], ncol = 4, nrow = 3)
m
ggsave("/media/storageE/ariel/R/finalpaper/liquid.hic.heatmaps/liquid.hic.heatmaps_1Mb.pdf", m)


#####################################
### CHANGE IN CONTACT SPECIFICITY ###
#####################################

enhancer_to_enhancer.gg8.p2.ctrl <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer.gg8.p2.5m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer.gg8.p2.30m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer.gg8.p2.60m <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

enhancer_to_enhancer.gg8.p2 <- as.data.frame(rbind(enhancer_to_enhancer.gg8.p2.ctrl, 
                                                   enhancer_to_enhancer.gg8.p2.5m,
                                                   enhancer_to_enhancer.gg8.p2.30m,
                                                   enhancer_to_enhancer.gg8.p2.60m)) %>%
  mutate(specificity = norm - enhancer_to_enhancer.gg8.p2.ctrl$norm) %>%
  mutate(specificity = 100*specificity/enhancer_to_enhancer.gg8.p2.ctrl$norm) %>% 
  select(norm, specificity)
enhancer_to_enhancer.gg8.p2$type <- "islands"
enhancer_to_enhancer.gg8.p2$time <- c(0, 5, 30, 60)


activeP2_to_enhancer.gg8.p2.ctrl <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.5m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.30m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer.gg8.p2.60m <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

activeP2_to_enhancer.gg8.p2 <- as.data.frame(rbind(activeP2_to_enhancer.gg8.p2.ctrl, 
                                                   activeP2_to_enhancer.gg8.p2.5m,
                                                   activeP2_to_enhancer.gg8.p2.30m,
                                                   activeP2_to_enhancer.gg8.p2.60m)) %>%
  mutate(specificity = norm - activeP2_to_enhancer.gg8.p2.ctrl$norm) %>%
  mutate(specificity = 100*specificity/activeP2_to_enhancer.gg8.p2.ctrl$norm) %>% 
  select(norm, specificity)
activeP2_to_enhancer.gg8.p2$type <- "active OR"
activeP2_to_enhancer.gg8.p2$time <- c(0, 5, 30, 60)

inactiveP2_to_enhancer.gg8.p2.ctrl <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl", "trans", "square") %>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.5m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.5m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.30m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.30m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer.gg8.p2.60m <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.60m", "trans", "square")%>%
  filter(bait_order == 20, prey_order == 20)

inactiveP2_to_enhancer.gg8.p2 <- as.data.frame(rbind(inactiveP2_to_enhancer.gg8.p2.ctrl, 
                                                   inactiveP2_to_enhancer.gg8.p2.5m,
                                                   inactiveP2_to_enhancer.gg8.p2.30m,
                                                   inactiveP2_to_enhancer.gg8.p2.60m)) %>%
  mutate(specificity = norm - inactiveP2_to_enhancer.gg8.p2.ctrl$norm) %>%
  mutate(specificity = 100*specificity/inactiveP2_to_enhancer.gg8.p2.ctrl$norm) %>% 
  select(norm, specificity)
inactiveP2_to_enhancer.gg8.p2$type <- "inactive OR"
inactiveP2_to_enhancer.gg8.p2$time <- c(0, 5, 30, 60)


enhancer.gg8.p2 <- rbind(enhancer_to_enhancer.gg8.p2, activeP2_to_enhancer.gg8.p2, inactiveP2_to_enhancer.gg8.p2)

ggplot(enhancer.gg8.p2, aes(x = time, y = specificity, colour = type)) + geom_line()




