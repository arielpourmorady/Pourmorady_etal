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
library(ggeasy)
library(patchwork)

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

# MERGED DATA

omp.gfp.name <- "OMPgfp.merge"
omp.gfp.trans_hic_path <-format_trans_path(omp.gfp.name)
omp.gfp.cis_hic_path <-format_cis_path(omp.gfp.name)
omp.gfp.total_count <- 564546231

red.name <- "OMPtTAtetOP2.mor28iCretdT.tdTpos.AP166.AP179.AP182"
red.trans_hic_path <-format_trans_path(red.name)
red.cis_hic_path <-format_cis_path(red.name)
red.total_count <- 517202034

yellow.name <- "OMPtTAtetOP2.mor28iCretdT.GFPpos.tdTpos.AP165.AP180.AP183"
yellow.trans_hic_path <-format_trans_path(yellow.name)
yellow.cis_hic_path <-format_cis_path(yellow.name)
yellow.total_count <- 448037007

green.name <- "OMPtTAtetOP2.mor28iCretdT.GFPpos.AP167.AP181.AP184"
green.trans_hic_path <-format_trans_path(green.name)
green.cis_hic_path <-format_cis_path(green.name)
green.total_count <- 511816455

# UNMERGED REPLICATES

AP165.name <- "AP165_r2"
AP165.trans_hic_path <-format_trans_path(AP165.name)
AP165.cis_hic_path <-format_cis_path(AP165.name)
AP165.total_count <- 182455307

AP166.name <- "AP166_r2"
AP166.trans_hic_path <-format_trans_path(AP166.name)
AP166.cis_hic_path <-format_cis_path(AP166.name)
AP166.total_count <- 202912243

AP167.name <- "AP167_r2"
AP167.trans_hic_path <-format_trans_path(AP167.name)
AP167.cis_hic_path <-format_cis_path(AP167.name)
AP167.total_count <- 204157250

AP179.name <- "AP179"
AP179.trans_hic_path <-format_trans_path(AP179.name)
AP179.cis_hic_path <-format_cis_path(AP179.name)
AP179.total_count <- 107579240

AP180.name <- "AP180"
AP180.trans_hic_path <-format_trans_path(AP180.name)
AP180.cis_hic_path <-format_cis_path(AP180.name)
AP180.total_count <- 116330862

AP181.name <- "AP181"
AP181.trans_hic_path <-format_trans_path(AP181.name)
AP181.cis_hic_path <-format_cis_path(AP181.name)
AP181.total_count <- 125114639

AP182.name <- "AP182"
AP182.trans_hic_path <-format_trans_path(AP182.name)
AP182.cis_hic_path <-format_cis_path(AP182.name)
AP182.total_count <- 107492505

AP183.name <- "AP183"
AP183.trans_hic_path <-format_trans_path(AP183.name)
AP183.cis_hic_path <-format_cis_path(AP183.name)
AP183.total_count <- 67198138

AP184.name <- "AP184"
AP184.trans_hic_path <-format_trans_path(AP184.name)
AP184.cis_hic_path <-format_cis_path(AP184.name)
AP184.total_count <- 94853375




## ESTABLISH BINS ##

bed_to_bait <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(bait = str_c(chrNo,"_",V2)) %>% dplyr::select(bait)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

bed_to_prey <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(prey = str_c(chrNo,"_",V2)) %>% dplyr::select(prey)
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

libraries <- c("omp.gfp_merge",
               "red_merge",
               "yellow_merge",
               "green_merge",
               "yellow_r1",
               "red_r1",
               "green_r1",
               "red_r2",
               "yellow_r2",
               "green_r2",
               "red_r3",
               "yellow_r3",
               "green_r3") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(omp.gfp.trans_hic_path,
                    red.trans_hic_path,
                    yellow.trans_hic_path,
                    green.trans_hic_path,
                    AP165.trans_hic_path,
                    AP166.trans_hic_path,
                    AP167.trans_hic_path,
                    AP179.trans_hic_path,
                    AP180.trans_hic_path,
                    AP181.trans_hic_path,
                    AP182.trans_hic_path,
                    AP183.trans_hic_path,
                    AP184.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path.trans_hic_path, .trans_hic_path, control.acid.r1.trans_hic_path.trans_hic_path, liquid.HiC.control.trans_hic_path.trans_hic_path, liquid.HiC.30.trans_hic_path.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(omp.gfp.cis_hic_path,
                  red.cis_hic_path,
                  yellow.cis_hic_path,
                  green.cis_hic_path,
                  AP165.cis_hic_path,
                  AP166.cis_hic_path,
                  AP167.cis_hic_path,
                  AP179.cis_hic_path,
                  AP180.cis_hic_path,
                  AP181.cis_hic_path,
                  AP182.cis_hic_path,
                  AP183.cis_hic_path,
                  AP184.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(omp.gfp.total_count,
                 red.total_count,
                 yellow.total_count,
                 green.total_count,
                 AP165.total_count,
                 AP166.total_count,
                 AP167.total_count,
                 AP179.total_count,
                 AP180.total_count,
                 AP181.total_count,
                 AP182.total_count,
                 AP183.total_count,
                 AP184.total_count)



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
###

### The filter that I am creating includes all the OR Clusters as well as 2-MB around each Island

filter <- rbind(TwoMB.Islands, TwoMB.ORs) %>% dplyr::select(prey) %>% unique()


for (i in c(1:13)) {
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
    hic_contacts <- dplyr::select(hic_contacts, -contact)
    rm(trans_all, cis_short, cis_long, cis_all)
  } else {
    hic_contacts_loop <- rbind(trans_all, cis_short, cis_long)
    hic_contacts_loop$geno <- print(libraries[i])
    hic_contacts_loop[is.na(hic_contacts_loop)] <- 0
    #totalcount <- hic_contacts_loop %>% group_by(geno) %>% summarise(sum = sum(contact)) # You should activate this when you ARE NOT FILTERING
    #hic_contacts_loop$norm <- hic_contacts_loop$contact/totalcount$sum # You should activate this when you ARE NOT FILTERING
    hic_contacts_loop$norm <- hic_contacts_loop$contact/total_count[i] # You should activate this when you *ARE* FILTERING
    hic_contacts_loop <- dplyr::select(hic_contacts_loop, -contact)
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
  dplyr::select(loc, order, name)

df_islands <- TwoMB.Islands %>% 
  mutate("loc" = prey, "name" = sub("\\_.*", "", name)) %>%
  dplyr::select(loc, order, name)


#################
### FUNCTIONS ###
#################

P2 <- "7_107050000"
P2_cluster <- "OR-Cluster-26_Olfr714,OR-Cluster-26_Olfr17"

mor28 <- "14_52450000"
mor28_cluster <- "OR-Cluster-55_Olfr1509,OR-Cluster-55_Olfr1508,OR-Cluster-55_Olfr1507"

order_tbl_df <- data.frame(bait_order = 1:39) 
order_tbl_df <- expand.grid(bait_order = order_tbl_df$bait_order, prey_order = order_tbl_df$bait_order) %>% as_tibble()

activemor28_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_ORs %>% filter(name == mor28_cluster) 
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
    #df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}


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
    #df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}

P2_to_mor28_matrix <- function(hic_contacts, genotype, arch, normalization) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans",  geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_ORs %>% filter(name == mor28_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df4 <- df_ORs %>% filter(name == P2_cluster)
  df6 <- left_join(df5, df4, by = c("prey" = "loc")) %>%
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


P2_to_mor28_matrix <- function(hic_contacts, genotype, arch, normalization) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans",  geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_ORs %>% filter(name == mor28_cluster) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df4 <- df_ORs %>% filter(name == P2_cluster)
  df6 <- left_join(df5, df4, by = c("prey" = "loc")) %>%
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
    #df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
  if (normalization == "row_prey"){
    df9 <- df9 %>%
      filter(prey_order == 20)
    #df9$norm <- df9$norm/sum(df9$norm)
    return(df9)
  }
}
#################
### LINEPLOTS ###
#################

activeP2_to_enhancer_matrix_red_merge <- activeP2_to_enhancer_matrix(hic_contacts, "red_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "red_merge") %>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix_yellow_merge <- activeP2_to_enhancer_matrix(hic_contacts, "yellow_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "yellow_merge")%>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix_green_merge <- activeP2_to_enhancer_matrix(hic_contacts, "green_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "green_merge")%>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix_omp.gfp_merge <- activeP2_to_enhancer_matrix(hic_contacts, "omp.gfp_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "omp.gfp_merge")%>% mutate(norm = norm*1e9) %>% filter(prey_order == 20) %>% select(norm)


activemor28_to_enhancer_matrix_red_merge <- activemor28_to_enhancer_matrix(hic_contacts, "red_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "red_merge")%>% mutate(norm = norm*1e9)
activemor28_to_enhancer_matrix_yellow_merge <- activemor28_to_enhancer_matrix(hic_contacts, "yellow_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "yellow_merge")%>% mutate(norm = norm*1e9)
activemor28_to_enhancer_matrix_green_merge <- activemor28_to_enhancer_matrix(hic_contacts, "green_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "green_merge")%>% mutate(norm = norm*1e9)
activemor28_to_enhancer_matrix_omp.gfp_merge <- activemor28_to_enhancer_matrix(hic_contacts, "omp.gfp_merge", "trans", "row_bait") %>%
  dplyr::mutate(geno = "omp.gfp_merge")%>% mutate(norm = norm*1e9) %>% filter(prey_order == 20) %>% select(norm)


b <- ggplot(activeP2_to_enhancer_matrix_red_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activeP2_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("firebrick3")) +
  scale_fill_manual(values=c("firebrick3")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("P2 to GI contacts")

a <- ggplot(activeP2_to_enhancer_matrix_green_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activeP2_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("green4")) +
  scale_fill_manual(values=c("green4")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("P2 to GI contacts")

c <- ggplot(activeP2_to_enhancer_matrix_yellow_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activeP2_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("gold1")) +
  scale_fill_manual(values=c("gold1")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("P2 to GI contacts")

e <- ggplot(activemor28_to_enhancer_matrix_red_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activemor28_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("firebrick3")) +
  scale_fill_manual(values=c("firebrick3")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("mor28 to GI contacts")

d <- ggplot(activemor28_to_enhancer_matrix_green_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activemor28_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("green4")) +
  scale_fill_manual(values=c("green4")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("mor28 to GI contacts")

f <- ggplot(activemor28_to_enhancer_matrix_yellow_merge, aes(x = prey_order, y = norm, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = activemor28_to_enhancer_matrix_omp.gfp_merge$norm, linetype = 'dashed') + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("gold1")) +
  scale_fill_manual(values=c("gold1")) +
  theme_classic() + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("mor28 to GI contacts")

a + b + c + d + e + f + plot_layout(ncol = 3, guides = "collect")

############################
### P2 to MOR28 CONTACTS ###
############################


P2_to_mor28_matrix.green_merge <- P2_to_mor28_matrix(hic_contacts, "green_merge", "trans", "square")  %>% 
  group_by(bait_order, prey_order) %>%
  summarise(norm = sum(norm)) %>%
  dcast(bait_order ~ prey_order, value.var = "norm")

P2_to_mor28_matrix.red_merge <- P2_to_mor28_matrix(hic_contacts, "red_merge", "trans", "square")  %>% 
  group_by(bait_order, prey_order) %>%
  summarise(norm = sum(norm)) %>%
  dcast(bait_order ~ prey_order, value.var = "norm")

P2_to_mor28_matrix.omp.gfp_merge <- P2_to_mor28_matrix(hic_contacts, "omp.gfp_merge", "trans", "square")  %>% 
  group_by(bait_order, prey_order) %>%
  summarise(norm = sum(norm)) %>%
  dcast(bait_order ~ prey_order, value.var = "norm")

P2_to_mor28_matrix(hic_contacts, "green_merge", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
P2_to_mor28_matrix(hic_contacts, "red_merge", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
P2_to_mor28_matrix(hic_contacts, "omp.gfp_merge", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)



paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(P2_to_mor28_matrix.omp.gfp_merge[,-1]), 
                  max(P2_to_mor28_matrix.omp.gfp_merge[,-1]), 
                  length.out=ceiling(paletteLength)))


q <- pheatmap(P2_to_mor28_matrix.green_merge[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=6, cellwidth =6,
              breaks=myBreaks,
              main = "GFP+ \n P2_to_mor28")


r <- pheatmap(P2_to_mor28_matrix.red_merge[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=6, cellwidth =6,
              breaks=myBreaks,
              main = "tdT+ \n P2_to_mor28")


s <- pheatmap(P2_to_mor28_matrix.omp.gfp_merge[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=6, cellwidth =6,
              breaks=myBreaks,
              main = "OMP-GFP \n P2_to_mor28")

grid.arrange(q[[4]], r[[4]], s[[4]], ncol = 3)

