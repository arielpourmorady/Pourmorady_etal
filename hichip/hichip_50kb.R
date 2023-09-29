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

# set directory to current directory
# this directory should contain all files in the GitHub folder

working_directory <- "/media/storageE/ariel/R/finalpaper_August2023/hichip/"
setwd(working_directory) 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

## Juicer Merged Files with KR Normalization ##
options(scipen = 999)

dump_dir <- '/media/storageD/distiller/dump/'
bin_size <- 50000

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

gg8tta.tetop2.h27ac.HiChIP.IP.name <- "gg8tTAtetOP2.H3K27ac.HiChIP.AP151.AP152_230921"
gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path <-format_trans_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path <-format_cis_path(gg8tta.tetop2.h27ac.HiChIP.IP.name)
gg8tta.tetop2.h27ac.HiChIP.IP.total_count <- 298853841


# All Bins ... Juicer KR Normalized Merged Files


libraries <- c("gg8.p2.ctrl.IF",
               "gg8tta.tetop2.h27ac.HiChIP.IP") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(gg8.p2.ctrl.IF.trans_hic_path,
                    gg8tta.tetop2.h27ac.HiChIP.IP.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path, , control.acid.r1.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.30.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(gg8.p2.ctrl.IF.cis_hic_path,
                  gg8tta.tetop2.h27ac.HiChIP.IP.cis_hic_path)#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(gg8.p2.ctrl.IF.total_count,
                 gg8tta.tetop2.h27ac.HiChIP.IP.total_count)


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
for (i in c(1:2)) {
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

# Fifty KB
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


###########################
### MATRICES & HEATMAPS ###
###########################
# Must be done at Fifty-KB resolution

# Generate Contact Specificity Score
activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
activeP2_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ctrl.IF", "trans", "square")%>% filter(bait_order == 20, prey_order == 20)
inactiveP2_to_enhancer_matrix(hic_contacts, "gg8tta.tetop2.h27ac.HiChIP.IP", "trans", "square") %>% filter(bait_order == 20, prey_order == 20)


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

100*(activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20] - activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20])/activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]
#
activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20] - activeP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]



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
 
c <- pheatmap(df, #df[16:24,16:24] for Zoom In
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               cellheight=3, cellwidth = 3,
              breaks = myBreaks)

grid.arrange(a[[4]], b[[4]], c[[4]])


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
100*(inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20] - inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20])/inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]
inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.IP[20,20] - inactiveP2_to_enhancer_matrix.gg8tta.tetop2.h27ac.HiChIP.input[20,20]

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

i <- pheatmap(df, #df[16:24,16:24]
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=3, cellwidth = 3,
              breaks = myBreaks)

grid.arrange(g[[4]], h[[4]], i[[4]])

grid.arrange(b[[4]], a[[4]], c[[4]],
             h[[4]], g[[4]], i[[4]], ncol = 3)
