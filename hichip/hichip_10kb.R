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
bin_size <- 10000

bed_to_prey <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(prey = str_c(chrNo,"_",V2)) %>% dplyr::select(prey)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA


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

gg8.p2.ctrl.IF.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.merge"
gg8.p2.ctrl.IF.trans_hic_path <-format_trans_path(gg8.p2.ctrl.IF.name)
gg8.p2.ctrl.IF.cis_hic_path <-format_cis_path(gg8.p2.ctrl.IF.name)
gg8.p2.ctrl.IF.total_count <- 362670259

omptta.tetop2.h27ac.HiChIP.name <- "AO.OMPtTAtetOP2.H3K27ac.HiChIP.220629"
omptta.tetop2.h27ac.HiChIP.trans_hic_path <-format_trans_path(omptta.tetop2.h27ac.HiChIP.name)
omptta.tetop2.h27ac.HiChIP.cis_hic_path <-format_cis_path(omptta.tetop2.h27ac.HiChIP.name)
omptta.tetop2.h27ac.HiChIP.total_count <- 10401866

gg8tta.tetop2.h27ac.HiChIP.IP.name <- "gg8tTAtetOP2.H3K27ac.HiChIP.AP151.AP152"
gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path <-format_trans_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path <-format_cis_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.total_count <- 332494243

gg8tta.tetop2.h27ac.HiChIP.input.name <- "gg8tTAtetOP2.H3K27ac.HiChIP.input.AP153"
gg8tta.tetop2.h27ac.HiChIP.input.trans_hic_path <-format_trans_path(gg8tta.tetop2.h27ac.HiChIP.input.name)
gg8tta.tetop2.h27ac.HiChIP.input.cis_hic_path <-format_cis_path(gg8tta.tetop2.h27ac.HiChIP.input.name)
gg8tta.tetop2.h27ac.HiChIP.input.total_count <- 66816284

omp.gfp.name <- "OMPgfp.merge"
omp.gfp.trans_hic_path <-format_trans_path(omp.gfp.name)
omp.gfp.cis_hic_path <-format_cis_path(omp.gfp.name)
omp.gfp.total_count <- 564546231

omp.gfp.h27ac.HiChIP.name <- "AO13"
omp.gfp.h27ac.HiChIP.trans_hic_path <-format_trans_path(omp.gfp.h27ac.HiChIP.name)
omp.gfp.h27ac.HiChIP.cis_hic_path <-format_cis_path(omp.gfp.h27ac.HiChIP.name)
omp.gfp.h27ac.HiChIP.total_count <- 117407088

# All Bins ... Juicer KR Normalized Merged Files


libraries <- c("gg8.p2.ctrl.IF",
               "omptta.tetop2.h27ac.HiChIP",
               "gg8tta.tetop2.h27ac.HiChIP.IP",
               "gg8tta.tetop2.h27ac.HiChIP.input",
               "omp.gfp",
               "omp.gfp.h27ac.HiChIP") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(gg8.p2.ctrl.IF.trans_hic_path,
                    omptta.tetop2.h27ac.HiChIP.trans_hic_path,
                    gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path,
                    gg8tta.tetop2.h27ac.HiChIP.input.trans_hic_path,
                    omp.gfp.trans_hic_path,
                    omp.gfp.h27ac.HiChIP.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path, , control.acid.r1.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.30.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(gg8.p2.ctrl.IF.cis_hic_path,
                  omptta.tetop2.h27ac.HiChIP.cis_hic_path,
                  gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path,
                  gg8tta.tetop2.h27ac.HiChIP.input.cis_hic_path,
                  omp.gfp.cis_hic_path,
                  omp.gfp.h27ac.HiChIP.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(gg8.p2.ctrl.IF.total_count,
                 omptta.tetop2.h27ac.HiChIP.total_count,
                 gg8tta.tetop2.h27ac.HiChIP.IP.total_count,
                 gg8tta.tetop2.h27ac.HiChIP.input.total_count,
                 omp.gfp.total_count,
                 omp.gfp.h27ac.HiChIP.total_count)


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
All_Bins_prey <- bed_to_prey("/media/storageA/kevin/annotation/mm10_assembled.10kb.bed")

Islands_10kb <- read.table("/media/storageE/ariel/annotations/Greek_Islands.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, type = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  select(prey, type) %>%
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique()

Bins.And.Islands <- left_join(All_Bins_prey, Islands_10kb)
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
All_Bins_prey <- bed_to_prey("/media/storageA/kevin/annotation/mm10_assembled.10kb.bed")

OR_Clusters_10kb <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique()

Bins.And.ORs <- left_join(All_Bins_prey, OR_Clusters_10kb) 
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

for (i in c(1,3)) {
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

order_tbl_df <- data.frame(bait_order = 1:199) 
order_tbl_df <- expand.grid(bait_order = order_tbl_df$bait_order, prey_order = order_tbl_df$bait_order) %>% as_tibble()

# Ten Kb
P2 <- "7_107090000"
P2_cluster <- "OR-Cluster-26_Olfr17"

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
# Must be done at Fifty-KB resolution


activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- activeP2_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[,-1]
activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input<- activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[,-1]

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  max(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  length.out=ceiling(paletteLength)))

a <- pheatmap(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

b <- pheatmap(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

## Log Fold Change in Contact Specificity

#activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20]/activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]
#2.221638 increase in contact specificity

#df <- log2(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP/activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input)
#df[sapply(df, is.infinite)] <- NA

df <- activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP - activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input
#df <- df[86:114, 86:114]
# Fifty Kb
ymin <- min(df)
ymax <- max(df)


paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(ymin, 0, length.out=ceiling(paletteLength/2)), 
              seq((ymax - ymin)/100, ymax, length.out=floor(paletteLength/2)))

c <- pheatmap(df, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)


enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- enhancer_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[,-1]
enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input<- enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[,-1]

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  max(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  length.out=ceiling(paletteLength)))

d <- pheatmap(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

e <- pheatmap(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

## Log Fold Change in Contact Specificity

#enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20]/enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]
#2.221638 increase in contact specificity

#df <- log2(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP/enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input)
#df[sapply(df, is.infinite)] <- NA

df <- enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP - enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input
#df <- df[86:114, 86:114]
# Fifty Kb
ymin <- min(df)
ymax <- max(df)

paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(ymin, 0, length.out=ceiling(paletteLength/2)), 
              seq((ymax - ymin)/100, ymax, length.out=floor(paletteLength/2)))

f <- pheatmap(df, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)





inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square") %>% dcast(bait_order ~ prey_order, value.var = "norm")
inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[,-1]
inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input <- inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square")%>% dcast(bait_order ~ prey_order, value.var = "norm")
inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input<- inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[,-1]

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  max(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP), 
                  length.out=ceiling(paletteLength)))

g <- pheatmap(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

h <- pheatmap(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

## Log Fold Change in Contact Specificity

#inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20]/inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]
#2.221638 increase in contact specificity

#df <- log2(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP/inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input)
#df[sapply(df, is.infinite)] <- NA

df <- inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP - inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input
#df <- df[85:115, 85:115]
# Fifty Kb
ymin <- min(df)
ymax <- max(df)

paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(ymin, 0, length.out=ceiling(paletteLength/2)), 
              seq((ymax - ymin)/100, ymax, length.out=floor(paletteLength/2)))

i <- pheatmap(df, 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

grid.arrange(b[[4]], a[[4]], c[[4]],
             e[[4]], d[[4]], f[[4]],
             h[[4]], g[[4]], i[[4]], ncol = 3)

grid.arrange(c[[4]],
             f[[4]],
             i[[4]], nrow = 3)


grid.arrange(b[[4]], a[[4]], h[[4]], g[[4]], nrow = 2)

z <- arrangeGrob(b[[4]], a[[4]], h[[4]], g[[4]], nrow = 2)

grid.arrange(e[[4]], d[[4]], f[[4]], ncol = 3)

y <- arrangeGrob(e[[4]], d[[4]], f[[4]], ncol = 3)

ggsave("/media/storageE/ariel/R/finalpaper/hichip/OR_hub_heatmaps.pdf", z)
ggsave("/media/storageE/ariel/R/finalpaper/hichip/island_hub_heatmaps.pdf", y)

# Must be done at Ten-KB resolution


P2_cluster_df <- df_ORs %>% filter(name == P2_cluster) %>% unique()
P2_cluster_df$order <- 105.1 + (P2_cluster_df$order - 1)*0.01
df <- hic_contacts %>% filter(bait %in% P2_cluster_df$loc) %>%
  filter(prey %in% Islands_10kb$prey) %>% 
  mutate(norm = norm*10000000) %>% 
  filter(arch == "trans")

df_hichip <- df %>% filter(geno == "gg8tta.tetop2.h27ac.HiChIP.IP")
df_hichip <- left_join(P2_cluster_df, df_hichip, by = c("loc" = "bait")) %>%
  group_by(order) %>% 
  summarise(norm = sum(norm))
df_hichip[is.na(df_hichip)] <- 0
df_hichip$geno <- "gg8tta.tetop2.h27ac.HiChIP.IP"

df_hic <- df %>% filter(geno == "gg8.p2.ctrl.IF")
df_hic <- left_join(P2_cluster_df, df_hic, by = c("loc" = "bait")) %>%
  group_by(order) %>% 
  summarise(norm = sum(norm))
df_hic[is.na(df_hic)] <- 0
df_hic$geno <- "gg8.p2.ctrl.IF"

df <- rbind(df_hic, df_hichip) 


a <- ggplot(df_hic, aes(x = order, y = norm)) + geom_line(colour = "darkblue") + theme_classic()+
  ylim(0, max(df_hichip$norm)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

b <- ggplot(df_hichip, aes(x = order, y = norm)) + geom_line(colour = "darkblue") + theme_classic() +
  ylim(0, max(df_hichip$norm)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

a + b

df_hic_narrow <- df_hic %>% filter(order < 107.09 + 1, order > 107.09 - 1)
df_hichip_narrow <- df_hichip %>% filter(order < 107.09 + 1, order > 107.09 - 1)


c <- ggplot(df_hic_narrow, aes(x = order, y = norm)) + geom_line(colour = "darkblue") + theme_classic()+
  ylim(0, max(df_hichip_narrow$norm)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

d <- ggplot(df_hichip_narrow, aes(x = order, y = norm)) + geom_line(colour = "darkblue") + theme_classic()+
  ylim(0, max(df_hichip_narrow$norm)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

a + b + c + d

# Single Island-to-Island Contacts #
# Must be done at Ten-KB resolution

island_to_island_hichip <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% Islands_10kb$prey, 
         arch == "trans", 
         bait != prey,
         geno == "gg8tta.tetop2.h27ac.HiChIP.IP") %>%
  mutate(norm = norm*1e9) %>%
  dcast(bait ~ prey, value.var = "norm")
island_to_island_hichip[is.na(island_to_island_hichip)] <- 0

island_to_island_hic <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% Islands_10kb$prey, 
         arch == "trans", 
         bait != prey,
         geno == "gg8.p2.ctrl.IF") %>%
  mutate(norm = norm*1e9) %>%
  dcast(bait ~ prey, value.var = "norm")
island_to_island_hic[is.na(island_to_island_hic)] <- 0


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(island_to_island_hic[,-1]), 
                  max(island_to_island_hic[,-1]), 
                  length.out=ceiling(paletteLength)))

a <- pheatmap(island_to_island_hic[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              treeheight_row = 0,
              treeheight_col = 0,
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

b <- pheatmap(island_to_island_hichip[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE,
              treeheight_row = 0,
              treeheight_col = 0,
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks
)


grid.arrange(a[[4]], b[[4]],ncol = 2)

####################
### Violin Plots ###
####################

df1 <- Islands_10kb$prey %>% as.data.frame() %>% rename("bait" = ".")
df2 <- Islands_10kb$prey %>% as.data.frame() %>% rename("bait" = ".")

df1$geno <- "gg8.p2.ctrl.IF"
df2$geno <- "gg8tta.tetop2.h27ac.HiChIP.IP"

df <- rbind(df1, df2)


island_to_island <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% Islands_10kb$prey, 
         arch == "trans", 
         bait != prey) %>%
  mutate(norm = norm*1e9) %>%
  group_by(bait, geno) %>%
  summarise(sum_norm = sum(norm))
island_to_island <- left_join(df, island_to_island, by = c("bait", "geno"))
island_to_island[is.na(island_to_island)] <- 0

island_to_inactiveORs <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% OR_Clusters_10kb$prey, 
         arch == "trans", 
         bait != prey) %>%
  mutate(norm = norm*1e9) %>%
  group_by(bait, geno) %>%
  summarise(sum_norm = sum(norm)) 
island_to_inactiveORs <- left_join(df, island_to_inactiveORs, by = c("bait", "geno"))
island_to_inactiveORs[is.na(island_to_inactiveORs)] <- 0

island_to_activeOR <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey == P2, 
         arch == "trans", 
         bait != prey) %>%
  mutate(norm = norm*1e9) %>%
  group_by(bait, geno) %>%
  summarise(sum_norm = sum(norm)) 
island_to_activeOR <- left_join(df, island_to_activeOR, by = c("bait", "geno"))
island_to_activeOR[is.na(island_to_activeOR)] <- 0


a <- ggplot(island_to_activeOR, aes(x = geno, y = sum_norm)) + 
 # geom_violin(fill = "grey", alpha ) +
#  geom_boxplot(width = 0.2) + 
  geom_dotplot(binaxis = "y",
       stackdir = "center", 
       dotsize = 0.5, 
       binwidth = 0.2,
       alpha = 0.5) + 
  theme_classic()
a


b <- ggplot(island_to_island, aes(x = geno, y = sum_norm)) + 
  #geom_violin(fill = "grey") +
 # geom_boxplot(width = 0.2) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               dotsize = 0.5, 
               binwidth = 0.2,
               alpha = 0.5) + 
  theme_classic()

c <- ggplot(island_to_inactiveORs, aes(x = geno, y = sum_norm)) + 
 # geom_violin(fill = "grey") +
 # geom_boxplot(width = 0.2, colour = "#b10101") + 
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               dotsize = 0.5, 
               binwidth = 0.2,
               alpha = 0.5) + 
  theme_classic()

a + b + c

t_test_fxn <- function(df) {
  v1 <- df %>% filter(geno == "gg8.p2.ctrl.IF") %>% as.data.frame() %>% dplyr::select(sum_norm)
  v2 <- df %>% filter(geno == "gg8tta.tetop2.h27ac.HiChIP.IP") %>% as.data.frame() %>% dplyr::select(sum_norm)
  p <- t.test(v1, v2)$p.value
  return(p)
}

t_test_fxn(island_to_activeOR)
t_test_fxn(island_to_island)
t_test_fxn(island_to_inactiveORs)

island_to_activeOR %>% 
  group_by(geno) %>%
  summarise(mean = mean(sum_norm),
            sd = sd(sum_norm))

island_to_inactiveORs%>% 
  group_by(geno) %>%
  summarise(mean = mean(sum_norm),
            sd = sd(sum_norm))

island_to_island%>% 
  group_by(geno) %>%
  summarise(mean = mean(sum_norm),
            sd = sd(sum_norm))






#####################################
### One Enhancer to All Enhancers ###
#####################################


oneenhancer_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization, one_enhancer) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_islands %>% filter(name == one_enhancer) 
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

oneenhancer_to_P2_matrix <- function(hic_contacts, genotype, arch, normalization, one_enhancer) {
  if (arch == "trans"){
    df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
      select(bait, prey, geno, norm) # filter all contacts made in one cell
  } else {
    df1 <- hic_contacts %>% filter(geno == genotype) %>%
      select(bait, prey, geno, norm)
  }
  df3 <- df_islands %>% filter(name == one_enhancer) 
  df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
    dplyr::rename("bait" = "loc") %>%
    mutate("bait_OR" = name, "bait_order" = order) %>% 
    select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
  df_P2 <- df_ORs %>% filter(name == P2_cluster) 
  df6 <- left_join(df5, df_P2, by = c("prey" = "loc")) %>%
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








datalist = list()
for (i in Islands$type) {
  print(i)
  df <-oneenhancer_to_P2_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square", i ) %>% 
    filter(bait_order == 20) 
  df$island <- i
  datalist[[i]] <- df
}
oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP = do.call(rbind, datalist)

datalist = list()
for (i in Islands$type) {
  print(i)
  df <-oneenhancer_to_P2_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square", i ) %>% 
    filter(bait_order == 20) 
  df$island <- i
  datalist[[i]] <- df
}
oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF = do.call(rbind, datalist)

oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP <- oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP %>%
  rename("norm_hichip" = "norm")

oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF <- oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF %>%
  rename("norm_input" = "norm")

oneenhancer_to_P2_matrix_log2fc <- left_join(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP, oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF, by = c("bait_order", "prey_order", "island")) %>%
  mutate("FC" = norm_hichip/norm_input) %>%
  mutate("log2fc" = log2(FC)) %>%
  dcast(island ~ prey_order, value.var = "log2fc") %>%
  select(!island)
oneenhancer_to_P2_matrix_log2fc <- oneenhancer_to_P2_matrix_log2fc[is.finite(rowSums(oneenhancer_to_P2_matrix_log2fc)),]

paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths

ymin <- min(oneenhancer_to_P2_matrix_log2fc)
ymax <- max(oneenhancer_to_P2_matrix_log2fc)

myBreaks <- c(seq(ymin, 0, length.out=ceiling(paletteLength/2)), 
              seq((ymax - ymin)/100, ymax, length.out=floor(paletteLength/2)))


pheatmap(oneenhancer_to_P2_matrix_log2fc, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = myColor,
         cellheight=3, cellwidth = 3,
         breaks = myBreaks)

oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP <- oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP %>%
  dcast(island ~ prey_order, value.var = "norm_hichip") %>%
  select(!island)

oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF <- oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF %>%
  dcast(island ~ prey_order, value.var = "norm_input") %>%
  select(!island)

oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP <- oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP[order(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP[,20],decreasing=TRUE),]
oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP <- na.omit(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP)

oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF <- oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF[order(oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF[,20],decreasing=TRUE),]
oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF <- na.omit(oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF)


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP), 
                  max(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP), 
                  length.out=ceiling(paletteLength)))

g <- pheatmap(oneenhancer_to_P2_matrix_gg8tta.tetop2.h27ac.HiChIP.IP[,11:29], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=6, cellwidth = 6,
              breaks = myBreaks)

h <- pheatmap(oneenhancer_to_P2_matrix_gg8.p2.ctrl.IF[,11:29], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=6, cellwidth = 6,
              breaks = myBreaks)

grid.arrange(g[[4]], h[[4]], ncol = 2)



# ### 10 kb Row ###
# 
# enhancer_to_enhancer_matrix <- function(hic_contacts, genotype, arch, normalization) {
#   if (arch == "trans"){
#     df1 <- hic_contacts %>% filter(arch == "trans", geno == genotype) %>%
#       dplyr::select(bait, prey, geno, norm) # filter all contacts made in one cell
#   } else {
#     df1 <- hic_contacts %>% filter(geno == genotype) %>%
#       dplyr::select(bait, prey, geno, norm)
#   }
#   df3 <- df_islands 
#   df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
#     dplyr::rename("bait" = "loc") %>%
#     mutate("bait_OR" = name, "bait_order" = order) %>% 
#     dplyr::select(bait, bait_order, bait_OR, prey, norm) # prey is 2 MB region around all islands
#   df6 <- left_join(df3, df5, by = c("loc" = "prey")) %>%
#     dplyr::rename("prey" = "loc") %>%
#     mutate("prey_OR" = name, "prey_order" = order) %>% 
#     dplyr::select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, norm)# bait is just the islands in the hub
#   df6 <- df6 %>% filter(bait_OR != prey_OR)
#   df7 <- df6[complete.cases(df6),]
#   df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
#     summarise(norm = sum(norm))
#   matrix <- dcast(df8,  bait_order ~ prey_order, value.var = "norm")
#   matrix[is.na(matrix)] <- min(matrix, na.rm = T)
#   if (normalization == "square"){
#     matrix <- matrix/(sum(df8$norm))
#     return(matrix[,102:300]) 
#   }
#   if (normalization == "row"){
#     matrix <- matrix[,102:300]
#     matrix<-t(apply(matrix,1, function(x) x/sum(x)))
#     return(matrix)
#   }
#   if (normalization == "none"){
#     matrix <- matrix[,102:300]
#     return(matrix)
#   }
# }
# 
# hic_contacts %>% filter(bait %in% Islands_10kb$prey) %>%
#   filter(prey )
# 
# 
# enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- enhancer_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP","trans", "row")[150,] %>% melt()
# enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP <- enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP %>%  mutate(order = row.names(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP))
# enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP$geno <- "hichip"
# 
# enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF <- enhancer_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF","trans", "row")[150,] %>% melt()
# enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF <- enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF %>%  mutate(order = row.names(enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF))
# enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF$geno <- "hic"
# 
# df <- rbind(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP, enhancer_to_enhancer_matrix.gg8.p2.ctrl.IF)
# df$value <- as.numeric(df$value)
# df$order <- as.numeric(df$order)
# 
# 
# ggplot(df, aes(x = order, y = value)) + geom_line(aes(colour = geno))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Heatmaps #
# 
# df_hichip <- df %>% filter(geno == "gg8tta.tetop2.h27ac.HiChIP.IP") 
# df_hichip <- left_join(P2_cluster_df, df_hichip, by = c("loc" = "bait")) 
# df_hichip <-  df_hichip %>% dcast(prey ~ loc, value.var = "norm")
# df_hichip[is.na(df_hichip)] <- 0
# 
# paletteLength <- 100
# myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
# myBreaks <- c(seq(min(df_hichip[,-1]), 
#                   0.5, 
#                   length.out=ceiling(paletteLength)))
# 
# d <- pheatmap(df_hichip[,-1], 
#               cluster_rows=FALSE, 
#               cluster_cols=FALSE, 
#               border_color = NA,
#               show_rownames = FALSE,
#               show_colnames = FALSE,
#               col = myColor,
#               cellheight=1, cellwidth = 1,
#               breaks = myBreaks)
# 
# e <- pheatmap(enhancer_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input, 
#               cluster_rows=FALSE, 
#               cluster_cols=FALSE, 
#               border_color = NA,
#               show_rownames = FALSE,
#               show_colnames = FALSE,
#               col = myColor,
#               cellheight=3, cellwidth = 3,
#               breaks = myBreaks)
# 
# 
# df_hic <- df %>% filter(geno == "gg8.p2.ctrl.IF")
# df_hic <- left_join(P2_cluster_df, df_hic, by = c("loc" = "bait")) %>%
#   group_by(order) %>% 
#   summarise(norm = sum(norm))
# df_hic[is.na(df_hic)] <- 0
# df_hic$geno <- "gg8.p2.ctrl.IF"
# 
# # c <- ggplot(df, aes(x = order, y = norm, colour = geno)) + 
# #   geom_line() + 
# #   theme_classic() +
# #   scale_color_manual(values = c("red", "blue")) +
# #   ylim(0, max(df_hichip$norm))
# # c
# 
# 
# 
# island_to_island <- hic_contacts %>%
#   filter(bait %in% Islands_10kb$prey, 
#          prey %in% Islands_10kb$prey, 
#          arch == "trans", 
#          bait != prey) %>%
#   mutate(norm = norm*1e9)
# 
# ggplot(island_to_island, aes(x = geno, y = norm)) + geom_boxplot() +
#   scale_y_log10() + ggtitle("cis single GI-to-GI \n contacts - 10kb")
# 
# island_to_island_hic <- hic_contacts %>%
#   filter(bait %in% Islands_10kb$prey, 
#          prey %in% Islands_10kb$prey, 
#          arch == "trans", 
#          bait != prey,
#          geno == "gg8.p2.ctrl.IF") %>%
#   mutate(norm = norm*1e9) %>%
#   dcast(bait ~ prey, value.var = "norm")
# island_to_island_hic[is.na(island_to_island_hic)] <- 0
# 
# 
# p2_to_island <- hic_contacts %>%
#   filter(bait %in% "7_107090000", 
#          prey %in% Islands_10kb$prey, 
#          arch == "trans") %>%
#   mutate(norm = norm*1e9) 
# 
# c <- ggplot(p2_to_island, aes(x = prey, y = norm)) + 
#   geom_point(aes(colour = geno)) + 
#   theme_classic()
# 
# c
# grid.arrange(a[[4]], b[[4]], c,
#              ncol = 3)
# 
# 
# 
# ggplot(island_to_island, aes(x = bait, y = norm)) + geom_point(aes(colour = geno))
# 
# 
# 
# 
# 
