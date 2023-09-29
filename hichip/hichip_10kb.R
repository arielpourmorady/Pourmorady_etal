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

working_directory <- "/media/storageE/ariel/R/finalpaper_August2023/hichip"
setwd(working_directory) 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################


## Juicer Merged Files with KR Normalization ##
options(scipen = 999)

dump_dir <- '/media/storageD/distiller/dump/'
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
All_Bins_prey <- bed_to_prey("mm10_assembled.10kb.bed")

Islands_10kb <- read.table("Greek_Islands.bed") %>%
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

OR_Clusters_10kb <- read.table("ORs-by-cluster.bed") %>%
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

##############################################################################
##############################################################################
##############################################################################

#####################
# FIGURE GENERATION #
#####################

# Single Island-to-Island Contacts #
# Must be done at Ten-KB resolution

Islands_10kb_trans_combinations <- expand.grid(Islands_10kb$prey, Islands_10kb$prey) %>%
  dplyr::rename("bait" = "Var1", "prey" = "Var2") %>%
  separate(bait, sep = "_", into = c("bait_chr", "bait_loc")) %>%
  separate(prey, sep = "_", into = c("prey_chr", "prey_loc")) %>%
  dplyr::filter(bait_chr != "X",
                prey_chr != "X") %>%
  mutate(bait_chr = as.numeric(bait_chr),
         prey_chr = as.numeric(prey_chr)) %>%
  mutate(arch = case_when(bait_chr == prey_chr ~ 'cis',
                          bait_chr != prey_chr ~ 'trans')) %>%
  dplyr::filter(arch == 'trans') %>%
  mutate(bait = paste(bait_chr, bait_loc, sep = "_"),
         prey = paste(prey_chr, prey_loc, sep = "_")) %>%
  dplyr::select(bait, prey, arch)


island_to_island_hichip <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% Islands_10kb$prey, 
         arch == "trans", 
         bait != prey,
         geno == "gg8tta.tetop2.h27ac.HiChIP.IP") %>%
  mutate(norm = norm*1e9) %>%
  dplyr::select(!geno)

island_to_island_hichip <- left_join(Islands_10kb_trans_combinations, 
                                     island_to_island_hichip, 
                                     by = c("bait", "prey", "arch"))
island_to_island_hichip[is.na(island_to_island_hichip)] <- 0
island_to_island_hichip <- island_to_island_hichip %>% dcast(bait ~ prey, value.var = "norm")
island_to_island_hichip[is.na(island_to_island_hichip)] <- 0



island_to_island_hic <- hic_contacts %>%
  filter(bait %in% Islands_10kb$prey, 
         prey %in% Islands_10kb$prey, 
         arch == "trans", 
         bait != prey,
         geno == "gg8.p2.ctrl.IF") %>%
  mutate(norm = norm*1e9) %>%
  dplyr::select(!geno)

island_to_island_hic <- left_join(Islands_10kb_trans_combinations, 
                                     island_to_island_hic, 
                                     by = c("bait", "prey", "arch"))
island_to_island_hic[is.na(island_to_island_hic)] <- 0
island_to_island_hic <- island_to_island_hic %>% dcast(bait ~ prey, value.var = "norm")
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
              breaks = myBreaks)

sum(island_to_island_hic[,-1])
sum(island_to_island_hichip[,-1])

grid.arrange(a[[4]], b[[4]],ncol = 2)

####################
### Violin Plots ###
####################

P2 <- "7_107090000"

Islands_10kb_filtered <- Islands_10kb %>% filter(!grepl("7_", prey))

df1 <- Islands_10kb_filtered$prey %>% as.data.frame() %>% rename("bait" = ".")
df2 <- Islands_10kb_filtered$prey %>% as.data.frame() %>% rename("bait" = ".")

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
  geom_violin(fill = "darkgreen", alpha = 0.5) +
  geom_boxplot(width = 0.2) + 
  theme_classic()
a


b <- ggplot(island_to_island, aes(x = geno, y = sum_norm)) + 
    geom_violin(fill = "grey", alpha = 0.5) +
  geom_boxplot(width = 0.2) + 
  theme_classic()
b
c <- ggplot(island_to_inactiveORs, aes(x = geno, y = sum_norm)) + 
  geom_violin(fill = "firebrick3", alpha = 0.5) +
  geom_boxplot(width = 0.2) + 
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
            sd = sd(sum_norm),
            n = n())

island_to_inactiveORs%>% 
  group_by(geno) %>%
  summarise(mean = mean(sum_norm),
            sd = sd(sum_norm),
            n = n())

island_to_island%>% 
  group_by(geno) %>%
  summarise(mean = mean(sum_norm),
            sd = sd(sum_norm),
            n = n())


########## Formatting ########## 

island_to_activeOR$category <- "active OR to GIs"
island_to_inactiveORs$category <- "inactive OR to GIs"
island_to_island$category <- "GIs to GIs"

islandSupplement <- rbind(island_to_activeOR,island_to_inactiveORs,  island_to_island)

write.table(islandSupplement, file = "/media/storageE/ariel/R/finalpaper_August2023/hichip/SuppFig7efg_HiChIP10kb.txt", sep = "\t")






