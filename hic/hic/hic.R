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

working_directory <- "/media/storageE/ariel/R/finalpaper_August2023/hic/hic/"
setwd(working_directory) 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

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

omp.gfp.name <- "OMPgfp.merge"
omp.gfp.trans_hic_path <-format_trans_path(omp.gfp.name)
omp.gfp.cis_hic_path <-format_cis_path(omp.gfp.name)
omp.gfp.total_count <- 564546231

p2.name <- "P2_merged"
p2.trans_hic_path <-format_trans_path(p2.name)
p2.cis_hic_path <-format_cis_path(p2.name)
p2.total_count <- 495587086

omp.p2.name <- "OMPtTA.tetOP2.NormalFood.ImmediateFix.Liquid.HiC.AP60.AP94.merge"
omp.p2.trans_hic_path <-format_trans_path(omp.p2.name)
omp.p2.cis_hic_path <-format_cis_path(omp.p2.name)
omp.p2.total_count <- 228623115

omp.p2.ko.name <- "OMPtTAtetOP2KO.AP187.AP188_230921"
omp.p2.ko.trans_hic_path <-format_trans_path(omp.p2.ko.name)
omp.p2.ko.cis_hic_path <-format_cis_path(omp.p2.ko.name)
omp.p2.ko.total_count <- 215358174

gg8.p2.name <- "gg8tTA.tetOP2.het.control.ImmediateFix.LiquidHiC.merge"
gg8.p2.trans_hic_path <-format_trans_path(gg8.p2.name)
gg8.p2.cis_hic_path <-format_cis_path(gg8.p2.name)
gg8.p2.total_count <- 362670259

gg8.p2.ko.name <- "gg8tTAtetOP2KO.GFPPos.AP147.AP149"
gg8.p2.ko.trans_hic_path <-format_trans_path(gg8.p2.ko.name)
gg8.p2.ko.cis_hic_path <-format_cis_path(gg8.p2.ko.name)
gg8.p2.ko.total_count <- 99715072

omp.m71.ko.name <- "OMPtTAtetOM71KO.AP145.AP146"
omp.m71.ko.trans_hic_path <-format_trans_path(omp.m71.ko.name)
omp.m71.ko.cis_hic_path <-format_cis_path(omp.m71.ko.name)
omp.m71.ko.total_count <- 101318384

atf5.name <- "HiC.Atf5-TR-WT.iRFP"
atf5.trans_hic_path <-format_trans_path(atf5.name)
atf5.cis_hic_path <-format_cis_path(atf5.name)
atf5.total_count <- 183989588

ngn.name <- "NGN1gfp.merge"
ngn.trans_hic_path <-format_trans_path(ngn.name)
ngn.cis_hic_path <-format_cis_path(ngn.name)
ngn.total_count <- 621909519

hbc.name <- "HBC.merge"
hbc.trans_hic_path <-format_trans_path(hbc.name)
hbc.cis_hic_path <-format_cis_path(hbc.name)
hbc.total_count <- 439746366

## ESTABLISH BINS ##

bed_to_bait <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(bait = str_c(chrNo,"_",V2)) %>% dplyr::select(bait)
  return(y)
} # will remove intervals from non numeric chromosomes with a warning about NA

bed_to_prey <- function(x) {
  y <- read.delim(x,header=FALSE) %>% mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>% filter(!is.na(chrNo)) %>% arrange(chrNo,V2) %>% mutate(prey = str_c(chrNo,"_",V2)) %>% dplyr::select(prey)
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


libraries <- c("omp.gfp",
               "p2",
               "omp.p2",
               "omp.p2.ko",
               "gg8.p2",
               "gg8.p2.ko",
               "omp.m71.ko",
               "atf5",
               "ngn",
               "hbc") # Here is the name of each library "P2","OMP", , "liquid.HiC.control", "liquid.HiC.30", "liquid.HiC.60"
trans_hic_path <- c(omp.gfp.trans_hic_path,
                    p2.trans_hic_path,
                    omp.p2.trans_hic_path,
                    omp.p2.ko.trans_hic_path,
                    gg8.p2.trans_hic_path,
                    gg8.p2.ko.trans_hic_path,
                    omp.m71.ko.trans_hic_path,
                    atf5.trans_hic_path,
                    ngn.trans_hic_path,
                    hbc.trans_hic_path) # Here are the trans conctacts  control.acid.r1.trans_hic_path, , control.acid.r1.trans_hic_path, liquid.HiC.control.trans_hic_path, liquid.HiC.30.trans_hic_path, liquid.HiC.60.trans_hic_path 
cis_hic_path <- c(omp.gfp.cis_hic_path,
                  p2.cis_hic_path,
                  omp.p2.cis_hic_path,
                  omp.p2.ko.cis_hic_path,
                  gg8.p2.cis_hic_path,
                  gg8.p2.ko.cis_hic_path,
                  omp.m71.ko.cis_hic_path,
                  atf5.cis_hic_path,
                  ngn.cis_hic_path,
                  hbc.cis_hic_path
                  )#Here are the cis contacts  , P2.cis_hic_path, OMP.cis_hic_path, , control.acid.r1.cis_hic_path, liquid.HiC.control.cis_hic_path,liquid.HiC.30.cis_hic_path,liquid.HiC.60.cis_hic_path
total_count <- c(omp.gfp.total_count,
                 p2.total_count,
                 omp.p2.total_count,
                 omp.p2.ko.total_count,
                 gg8.p2.total_count,
                 gg8.p2.ko.total_count,
                 omp.m71.ko.total_count,
                 atf5.total_count,
                 ngn.total_count,
                 hbc.total_count
                 )


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
for (i in c(1:10)) {
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
  return(df9) 
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

###########################################
## CUMMULATIVE GI CONTACTS OVER P2 LOCUS ##
###########################################

activeP2_to_enhancer_matrix.omp.gfp <- activeP2_to_enhancer_matrix(hic_contacts, "omp.gfp", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "omp.gfp")

activeP2_to_enhancer_matrix.atf5 <- activeP2_to_enhancer_matrix(hic_contacts, "atf5", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "atf5")

activeP2_to_enhancer_matrix.p2 <- activeP2_to_enhancer_matrix(hic_contacts, "p2", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "p2")

activeP2_to_enhancer_matrix.gg8.p2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "gg8.p2")

activeP2_to_enhancer_matrix.gg8.p2.ko <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ko", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "gg8.p2.ko")

activeP2_to_enhancer_matrix.omp.p2 <- activeP2_to_enhancer_matrix(hic_contacts, "omp.p2", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "omp.p2")

activeP2_to_enhancer_matrix.omp.p2.ko <- activeP2_to_enhancer_matrix(hic_contacts, "omp.p2.ko", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "omp.p2.ko")

activeP2_to_enhancer_matrix.omp.m71.ko <- activeP2_to_enhancer_matrix(hic_contacts, "omp.m71.ko", "trans", "none") %>% 
  filter(bait_order == 20) %>% mutate(geno = "omp.m71.ko")

df <- rbind(activeP2_to_enhancer_matrix.omp.gfp,
            activeP2_to_enhancer_matrix.atf5,
            activeP2_to_enhancer_matrix.p2,
            activeP2_to_enhancer_matrix.gg8.p2,
            activeP2_to_enhancer_matrix.gg8.p2,
            activeP2_to_enhancer_matrix.gg8.p2.ko,
            activeP2_to_enhancer_matrix.omp.p2,
            activeP2_to_enhancer_matrix.omp.p2.ko,
            activeP2_to_enhancer_matrix.omp.m71.ko) %>%
  mutate(prey_order = as.numeric(prey_order),
         norm = norm*1000000000)
  

m <- ggplot(df, aes(x = prey_order, y = norm, colour = geno)) + 
  geom_line(size = 0.75) + 
  geom_area(aes(fill = geno), position = "identity", alpha = 0.2) + 
  ylab("contacts per billion") + 
  ggtitle("trans Greek Island \n contacts to P2 locus") + 
  theme_classic()
m

df %>% filter(prey_order == 20)

write.table(df, file = "/media/storageE/ariel/R/finalpaper_August2023/hic/hic/SuppFig10d_P2toGIcontacts.txt")



#################################
#################################

activeP2_to_enhancer_matrix_gg8.p2 <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2", "trans", "none") %>% filter(bait_order == 20) %>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix.omp.p2 <- activeP2_to_enhancer_matrix(hic_contacts, "omp.p2", "trans", "none") %>% filter(bait_order == 20) %>% mutate(norm = norm*1e9)

activeP2_to_enhancer_matrix_gg8.p2.ko <- activeP2_to_enhancer_matrix(hic_contacts, "gg8.p2.ko", "trans", "none") %>% filter(bait_order == 20) %>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix.omp.p2.ko <- activeP2_to_enhancer_matrix(hic_contacts, "omp.p2.ko", "trans", "none") %>% filter(bait_order == 20) %>% mutate(norm = norm*1e9)
activeP2_to_enhancer_matrix.omp.gfp <- activeP2_to_enhancer_matrix(hic_contacts, "omp.gfp", "trans", "none") %>% filter(bait_order == 20) %>% mutate(norm = norm*1e9)

a <- ggplot(activeP2_to_enhancer_matrix_gg8.p2.ko, aes(x = prey_order, y = norm)) + 
  geom_line(color = 'green4') + 
  geom_point(color = 'green4') + 
  geom_area(position = "identity", alpha = 0.05, fill = 'green4') + 
  theme_classic() + 
  geom_hline(yintercept = as.numeric(activeP2_to_enhancer_matrix.omp.gfp[20,3]), linetype = 'dashed') + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("P2 to GI contacts \n gg8.p2.ko")
a

b <- ggplot(activeP2_to_enhancer_matrix.omp.p2.ko, aes(x = prey_order, y = norm)) + 
  geom_line(color = 'green4') + 
  geom_point(color = 'green4') + 
  geom_area(position = "identity", alpha = 0.05, fill = 'green4') + 
  theme_classic() + 
  geom_hline(yintercept = as.numeric(activeP2_to_enhancer_matrix.omp.gfp[20,3]), linetype = 'dashed') + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,500) + 
  ggtitle("P2 to GI contacts \n omp.p2.ko")
b
a + b + plot_layout(ncol = 2)


c <- ggplot(activeP2_to_enhancer_matrix_gg8.p2, aes(x = prey_order, y = norm)) + 
  geom_line(color = 'green4') + 
  geom_point(color = 'green4') + 
  geom_area(position = "identity", alpha = 0.05, fill = 'green4') + 
  theme_classic() + 
  geom_hline(yintercept = as.numeric(activeP2_to_enhancer_matrix.omp.gfp[20,3]), linetype = 'dashed') + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,700) + 
  ggtitle("P2 to GI contacts \n gg8.p2")
c

d <- ggplot(activeP2_to_enhancer_matrix.omp.p2, aes(x = prey_order, y = norm)) + 
  geom_line(color = 'green4') + 
  geom_point(color = 'green4') + 
  geom_area(position = "identity", alpha = 0.05, fill = 'green4') + 
  theme_classic() + 
  geom_hline(yintercept = as.numeric(activeP2_to_enhancer_matrix.omp.gfp[20,3]), linetype = 'dashed') + 
  ylab("contacts per billion") + 
  xlab("Average GI locus ± 1Mb") + 
  ylim(0,700) + 
  ggtitle("P2 to GI contacts \n omp.p2")
d
c + d + plot_layout(ncol = 2)

#################################################################
#### GREEK ISLAND TO GREEK ISLAND HEATMAPS DURING DEVELOPMENT ###
#################################################################

enhancer_to_enhancer_matrix.hbc <- enhancer_to_enhancer_matrix(hic_contacts, "hbc", "trans", "none") %>% mutate(norm = norm*1e6)
enhancer_to_enhancer_matrix.hbc %>% filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer_matrix.hbc <- enhancer_to_enhancer_matrix.hbc %>% dcast(bait_order ~ prey_order, value.var = "norm")

enhancer_to_enhancer_matrix.ngn <- enhancer_to_enhancer_matrix(hic_contacts, "ngn", "trans", "none") %>% mutate(norm = norm*1e6)
enhancer_to_enhancer_matrix.ngn %>% filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer_matrix.ngn <- enhancer_to_enhancer_matrix.ngn %>% dcast(bait_order ~ prey_order, value.var = "norm")

enhancer_to_enhancer_matrix.atf5 <- enhancer_to_enhancer_matrix(hic_contacts, "atf5", "trans", "none") %>% mutate(norm = norm*1e6)
enhancer_to_enhancer_matrix.atf5 %>% filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer_matrix.atf5 <- enhancer_to_enhancer_matrix.atf5 %>% dcast(bait_order ~ prey_order, value.var = "norm")

enhancer_to_enhancer_matrix.omp.gfp <- enhancer_to_enhancer_matrix(hic_contacts, "omp.gfp", "trans", "none") %>% mutate(norm = norm*1e6)
enhancer_to_enhancer_matrix.omp.gfp %>% filter(bait_order == 20, prey_order == 20)
enhancer_to_enhancer_matrix.omp.gfp <- enhancer_to_enhancer_matrix.omp.gfp %>% dcast(bait_order ~ prey_order, value.var = "norm")

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(enhancer_to_enhancer_matrix.omp.gfp[,-1]), 
                  max(enhancer_to_enhancer_matrix.omp.gfp[,-1]), 
                  length.out=ceiling(paletteLength)))

h <- pheatmap(enhancer_to_enhancer_matrix.omp.gfp[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=4, cellwidth = 4,
              breaks=myBreaks,
              main = "omp.gfp \n GI to GI contacts")

i <- pheatmap(enhancer_to_enhancer_matrix.hbc[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=4, cellwidth = 4,
              breaks=myBreaks,
              main = "hbc \n GI to GI contacts")

j <- pheatmap(enhancer_to_enhancer_matrix.atf5[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=4, cellwidth = 4,
              breaks=myBreaks,
              main = "atf5 \n GI to GI contacts")

k <- pheatmap(enhancer_to_enhancer_matrix.ngn[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              cellheight=4, cellwidth = 4,
              breaks=myBreaks,
              main = "ngn \n GI to GI contacts")

grid.arrange(i[[4]], k[[4]], j[[4]], h[[4]], ncol = 4)

