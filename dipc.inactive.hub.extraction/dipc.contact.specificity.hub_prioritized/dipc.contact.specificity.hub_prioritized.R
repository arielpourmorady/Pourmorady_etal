# Clear environment

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
options(scipen=999)
'%out%' <- function(x,y)!('%in%'(x,y))
bin_size <- 50000


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

### Haplotyping Bait and Prey Files

# Function creates haplotyped dataframes
haplotype <- function(df) {
  df_mat <- df
  df_mat$prey_chr <- sprintf('%s,1', df_mat$prey_chr)
  df_mat$name <- sprintf('%s,1', df_mat$name)
  df_pat <- df
  df_pat$prey_chr <- sprintf('%s,0', df_pat$prey_chr)
  df_pat$name <- sprintf('%s,0', df_pat$name)
  df2 <- rbind (df_mat, df_pat)
  df2$prey <- paste(df2$prey_chr, "_", df2$prey_loc, sep ="")
  return(df2)
}

TwoMB.Islands_hap <- haplotype(TwoMB.Islands)
TwoMB.ORs_hap <- haplotype(TwoMB.ORs)

# Loading P2 Data
P2_ImmFix_dipc_male <- fread(file = '/media/storageE/ariel_dipc/P2-ImmFix-dipc/male_pairs/P2-ImmFix-dipc.male.50kb.contacts.haplotyped.txt') %>% as.data.frame()
P2_ImmFix_dipc_female <- fread(file = '/media/storageE/ariel_dipc/P2-ImmFix-dipc/female_pairs/P2-ImmFix-dipc.female.50kb.contacts.haplotyped.txt') %>% as.data.frame()
P2_ImmFix_dipc <- rbind(P2_ImmFix_dipc_male, P2_ImmFix_dipc_female)
P2_ImmFix_dipc$geno <- "ctrl"
rm(P2_ImmFix_dipc_male, P2_ImmFix_dipc_female)
P2_dipc <- P2_ImmFix_dipc %>% as.data.table()

# Loading MOR28 Data
MOR28_ImmFix_dipc_male <- fread(file = '/media/storageE/ariel_dipc/APM28DIP/male_pairs/APM28DIP.male.50kb.contacts.haplotyped.txt') %>% as.data.frame()
MOR28_ImmFix_dipc_female <- fread(file = '/media/storageE/ariel_dipc/APM28DIP/female_pairs/APM28DIP.female.50kb.contacts.haplotyped.txt') %>% as.data.frame()
MOR28_ImmFix_dipc <- rbind(MOR28_ImmFix_dipc_male, MOR28_ImmFix_dipc_female)
MOR28_ImmFix_dipc$geno <- "ctrl"
rm(MOR28_ImmFix_dipc_male, MOR28_ImmFix_dipc_female)
MOR28_dipc <- MOR28_ImmFix_dipc %>% as.data.table()

## extracting contacts

contact_extraction <- function(df, OR) {
  require(data.table)
  df <- df[, c("bait_chr", "bait_loc") := tstrsplit(bait, "_", fixed = TRUE)]
  df <- df[, c("prey_chr", "prey_loc") := tstrsplit(prey, "_", fixed = TRUE)]
  df_cis <- df %>% filter(bait_chr == prey_chr, bait_haplo == prey_haplo)
  df_cis$arch <- "cis"
  df_trans <- df %>% filter((bait_chr != prey_chr) | (bait_chr == prey_chr & bait_haplo != prey_haplo))
  df_trans$arch <- "trans"
  df <- rbind(df_cis, df_trans) 
  df$bait <- paste(df$bait_chr, ",", df$bait_haplo, "_", df$bait_loc, sep ="")
  df$prey <- paste(df$prey_chr, ",", df$prey_haplo, "_", df$prey_loc, sep ="")
  df <- df %>% select(bait, prey, size, cell, arch, geno)
  df$OR <- OR
  return(df)
}

P2_dipc <- contact_extraction(P2_dipc, "P2")
MOR28_dipc <- contact_extraction(MOR28_dipc, "MOR28")


####### ####### DATAFRAMES ####### #######

df_ORs_h <- TwoMB.ORs_hap %>% 
  mutate("loc" = prey) %>%
  select(loc, order, name)

df_Islands_h <- TwoMB.Islands_hap %>% 
  mutate("loc" = prey) %>%
  select(loc, order, name)

islands_percell_dataframe_P2 <- read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/islands_percell_dataframe_P2.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell)

islands_percell_dataframe_MOR28 <- read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/islands_percell_dataframe_MOR28.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell)

islands_percell_dataframe_MOR28_plot <- islands_percell_dataframe_MOR28 %>%
  group_by(state, cell) %>%
  summarise(total = n())

islands_percell_dataframe_P2_plot <- islands_percell_dataframe_P2  %>%
  group_by(state, cell) %>%
  summarise(total = n())

islands_percell_dataframe_plot <- rbind(islands_percell_dataframe_MOR28_plot, islands_percell_dataframe_P2_plot)

ggplot(islands_percell_dataframe_plot, aes(x = state, y = total)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) +
  ylab("GIs per hub") + 
  ggtitle("Equal Sampling of GIHs") + 
  theme_classic()

active_island_count <- islands_percell_dataframe_plot %>% filter(state == "active") %>% as.data.frame() %>% dplyr::select(total)
inactive_island_count <- islands_percell_dataframe_plot %>% filter(state == "inactive") %>% as.data.frame() %>% dplyr::select(total)

t.test(active_island_count$total , inactive_island_count$total)

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>%
  separate(col = "type", sep = "_", into  =c("cluster", "name")) %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name) 

inactiveORs_percell_dataframe_P2 <-  read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/inactiveORs_percell_dataframe_P2.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., df_test, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique() %>%
  dplyr::rename("name" = "type") %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name)

inactiveORs_percell_dataframe_P2 <- inactiveORs_percell_dataframe_P2 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., df_test, by = "prey") %>%
  dplyr::select(cell, prey, name)





df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>%
  separate(col = "type", sep = "_", into  =c("cluster", "name")) %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name) 

inactiveORs_percell_dataframe_MOR28 <-  read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/inactiveORs_percell_dataframe_MOR28.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., df_test, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique() %>%
  dplyr::rename("name" = "type") %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name)

inactiveORs_percell_dataframe_MOR28 <- inactiveORs_percell_dataframe_MOR28 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., df_test, by = "prey") %>%
  dplyr::select(cell, prey, name)

#### Renaming DF of Inactive ORs in active GIH ####

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>%
  separate(col = "type", sep = "_", into  =c("cluster", "name")) %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name) 

activeGIH_ORs_percell_datalist_P2 <-  read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/activeGIH_ORs_percell_datalist_P2.txt",sep = "\t") %>%
  dplyr::filter(element != "active") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., df_test, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique() 

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique() %>%
  dplyr::rename("name" = "type") %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name)

activeGIH_ORs_percell_datalist_P2 <- activeGIH_ORs_percell_datalist_P2 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., df_test, by = "prey") %>%
  dplyr::select(cell, prey, name)





df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>%
  separate(col = "type", sep = "_", into  =c("cluster", "name")) %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name) 

activeGIH_ORs_percell_datalist_MOR28 <-  read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/activeGIH_ORs_percell_datalist_MOR28.txt",sep = "\t") %>%
  dplyr::filter(element != "active") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., df_test, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

df_test <- read.table("/media/storageA/kevin/annotation/ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>% 
  group_by(prey) %>%
  mutate(type = paste(type, collapse = ",")) %>% 
  unique() %>%
  dplyr::rename("name" = "type") %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name)

activeGIH_ORs_percell_datalist_MOR28 <- activeGIH_ORs_percell_datalist_MOR28 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., df_test, by = "prey") %>%
  dplyr::select(cell, prey, name)

####### ####### FUNCTIONS ####### #######

# This first function establishes the OR compartment background
# which serves as the "expected" that the Island data is 
# normalized to

order_tbl_df <- data.frame(bait_order = 1:39) 
order_tbl_df <- expand.grid(bait_order = order_tbl_df$bait_order, prey_order = order_tbl_df$bait_order) %>% as_tibble()

active_islands_datalist <- function(dipc_contacts, df_activeislands, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_activeislands %>% 
      filter(cell == i, state == "active")
    print("Islands active in cell")
    print(df2$name)
    df4 <- inner_join(df_Islands_h, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name")
    print(df5)
    df6 <- df5 %>%
      dplyr::filter(prey_name %in% df2$name) %>%
      dplyr::filter(bait_name %in% df2$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    df7 <- left_join(order_tbl_df, df6)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    if (normalization == "square"){
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "none"){
      datalist[[i]] <- df7
    }
    if (normalization == "row_bait"){
      df7 <- df7 %>%
        filter(bait_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "row_prey"){
      df7 <- df7 %>%
        filter(prey_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
  }
  return(datalist)
}  


inactive_islands_datalist <- function(dipc_contacts, df_activeislands, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_activeislands %>% 
      filter(cell == i, state == "inactive")
    print("Islands active in cell")
    print(df2$name)
    df4 <- inner_join(df_Islands_h, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name")
    print(df5)
    df6 <- df5 %>%
      dplyr::filter(prey_name %in% df2$name) %>%
      dplyr::filter(bait_name %in% df2$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    df7 <- left_join(order_tbl_df, df6)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    if (normalization == "square"){
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "none"){
      datalist[[i]] <- df7
    }
    if (normalization == "row_bait"){
      df7 <- df7 %>%
        filter(bait_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "row_prey"){
      df7 <- df7 %>%
        filter(prey_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
  }
  return(datalist)
}  

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

active_islands_datalist_P2 <- active_islands_datalist(P2_dipc, islands_percell_dataframe_P2, "none")
inactive_islands_datalist_P2 <- inactive_islands_datalist(P2_dipc, islands_percell_dataframe_P2, "none")

active_islands_datalist_MOR28 <- active_islands_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "none")
inactive_islands_datalist_MOR28 <- inactive_islands_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "none")

big_data = do.call(rbind, active_islands_datalist_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data = do.call(rbind, inactive_islands_datalist_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

a <- pheatmap(active_matrix[,-1], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = myColor,
      #   breaks = myBreaks,
         cellheight=3, cellwidth = 3)

b <- pheatmap(inactive_matrix[,-1], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = myColor,
        # breaks = myBreaks,
         cellheight=3, cellwidth = 3)


big_data = do.call(rbind, active_islands_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data = do.call(rbind, inactive_islands_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

c <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              #   breaks = myBreaks,
              cellheight=3, cellwidth = 3)

d <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              # breaks = myBreaks,
              cellheight=3, cellwidth = 3)

grid.arrange(a[[4]], b[[4]], c[[4]],d[[4]],  ncol = 2)



big_data_P2 = do.call(rbind, active_islands_datalist_P2)
big_data_MOR28 = do.call(rbind, active_islands_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, inactive_islands_datalist_P2)
big_data_MOR28 = do.call(rbind, inactive_islands_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

e <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active GIHs",
              breaks = myBreaks,
              cellheight=3, cellwidth = 3)

f <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive GIHs",
              breaks = myBreaks,
              cellheight=3, cellwidth = 3)

g <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active GIHs",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)

h <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive GIHs",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)


grid.arrange(a[[4]], b[[4]], 
             c[[4]], d[[4]],
              e[[4]], f[[4]],  ncol = 2)

grid.arrange(e[[4]], f[[4]], g[[4]], h[[4]],  ncol = 2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


activeOR_hub_datalist <- function(dipc_contacts,df_activeislands, activeOR, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df3 <- df_activeislands %>% 
      filter(cell == i, state == "active")
    print("Islands active in cell")
    print(df3$name)
    df2 <- df_ORs_h %>% filter(name == activeOR)
    df4 <- inner_join(df2, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    print(df4)
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    print(df5)
     df7 <- left_join(order_tbl_df, df5)
     df7 <- df7 %>% replace(is.na(.), 0)
     print(df7)
    if (normalization == "square"){
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "none"){
      datalist[[i]] <- df7
    }
    if (normalization == "row_bait"){
      df7 <- df7 %>%
        filter(bait_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
    if (normalization == "row_prey"){
      df7 <- df7 %>%
        filter(prey_order == 20)
      df7$size <- df7$size/sum(df7$size)
      datalist[[i]] <- df7
    }
  }
  return(datalist)
} 

activeOR_hub_datalist_P2 <- activeOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, "OR-Cluster-26_Olfr714,OR-Cluster-26_Olfr17,0"  ,"none")
activeOR_hub_datalist_MOR28 <- activeOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "OR-Cluster-55_Olfr1509,OR-Cluster-55_Olfr1508,OR-Cluster-55_Olfr1507,0"  ,"none")



inactiveOR_hub_datalist <- function(dipc_contacts, df_activeislands, df_inactiveORs, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df3 <- df_activeislands %>% 
      filter(cell == i, state == "inactive")
    print("Islands active in cell")
    print(df3$name)
    df_inactiveORs <- df_inactiveORs %>% filter(cell == i)
    df2 <- df_ORs_h %>% filter(name %in% df_inactiveORs$name)
    df4 <- inner_join(df2, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    print(df4)
    df5_a <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      group_by(bait_name) %>%
      summarise(size = sum(size))
    df5_a <- df5_a[which.max(df5_a$size),]
    print(df5_a)
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      dplyr::filter(bait_name %in% df5_a$bait_name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    print(df5)
    df7 <- left_join(order_tbl_df, df5)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    datalist[[i]] <- df7
  }
  return(datalist)
} 

notmax_inactiveOR_hub_datalist <- function(dipc_contacts, df_activeislands, df_inactiveORs, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df3 <- df_activeislands %>% 
      filter(cell == i, state == "inactive")
    print("Islands active in cell")
    print(df3$name)
    df_inactiveORs <- df_inactiveORs %>% filter(cell == i)
    print(df_inactiveORs)
    df2 <- df_ORs_h %>% filter(name %in% df_inactiveORs$name)
    df4 <- inner_join(df2, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    print(df4)
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    print(df5)
    df7 <- left_join(order_tbl_df, df5)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    datalist[[i]] <- df7
  }
  return(datalist)
} 

all_ORs_inactiveOR_hub_datalist <- function(dipc_contacts, df_activeislands, df_inactiveORs, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df3 <- df_activeislands %>% 
      filter(cell == i, state == "inactive")
    print("Islands active in cell")
    print(df3$name)
    df2 <- df_ORs_h
    df4 <- inner_join(df2, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    print(df4)
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    print(df5)
    df7 <- left_join(order_tbl_df, df5)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    datalist[[i]] <- df7
  }
  return(datalist)
} 



inactiveOR_hub_datalist_P2 <- inactiveOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, inactiveORs_percell_dataframe_P2)
inactiveOR_hub_datalist_MOR28 <- inactiveOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, inactiveORs_percell_dataframe_MOR28)

notmax_inactiveOR_hub_datalist_P2 <- notmax_inactiveOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, inactiveORs_percell_dataframe_P2)
notmax_inactiveOR_hub_datalist_MOR28 <- notmax_inactiveOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, inactiveORs_percell_dataframe_MOR28)

all_ORs_inactiveOR_hub_datalist_P2 <- all_ORs_inactiveOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, inactiveORs_percell_dataframe_P2)
all_ORs_inactiveOR_hub_datalist_MOR28 <- all_ORs_inactiveOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, inactiveORs_percell_dataframe_MOR28)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

big_data = do.call(rbind, activeOR_hub_datalist_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))

big_data = do.call(rbind, inactiveOR_hub_datalist_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

a <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active OR-GIH P2",
              #   breaks = myBreaks,
              cellheight=3, cellwidth = 3)

b <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive OR-GIH P2",
              # breaks = myBreaks,
              cellheight=3, cellwidth = 3)


big_data = do.call(rbind, activeOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))

big_data = do.call(rbind, inactiveOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

c <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active OR-GIH MOR28",
              #   breaks = myBreaks,
              cellheight=3, cellwidth = 3)

d <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive OR-GIH MOR28",
               breaks = myBreaks,
              cellheight=3, cellwidth = 3)

grid.arrange(a[[4]], b[[4]], c[[4]],d[[4]],  ncol = 2)



big_data_P2 = do.call(rbind, activeOR_hub_datalist_P2)
big_data_MOR28 = do.call(rbind, activeOR_hub_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, inactiveOR_hub_datalist_P2)
big_data_MOR28 = do.call(rbind, inactiveOR_hub_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, notmax_inactiveOR_hub_datalist_P2)
big_data_MOR28 = do.call(rbind, notmax_inactiveOR_hub_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
notmax_inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, all_ORs_inactiveOR_hub_datalist_P2)
big_data_MOR28 = do.call(rbind, all_ORs_inactiveOR_hub_datalist_MOR28)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
all_ORs_inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)#/(sum(big_data_sum$sum))


paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  max(active_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

e <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
             # main = "Active GIHs",
              breaks = myBreaks,
              cellheight=3, cellwidth = 3)

f <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
             # main = "Inactive GIHs",
              breaks = myBreaks,
              cellheight=3, cellwidth = 3)

g <- pheatmap(notmax_inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active OR GIH",
              breaks = myBreaks,
              cellheight=3, cellwidth = 3)

h <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive OR GIH",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)


grid.arrange(a[[4]], b[[4]], 
             c[[4]], d[[4]],
             e[[4]], f[[4]],  ncol = 2)

k <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active OR GIH",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)

l <- pheatmap(notmax_inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "notmax_inactive OR GIH",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)




o <- pheatmap(active_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Active OR GIH",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)

p <- pheatmap(all_ORs_inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "all_ORs_inactive OR GIH",
              #breaks = myBreaks,
              cellheight=3, cellwidth = 3)


grid.arrange(m[[4]], n[[4]], o[[4]], p[[4]],  ncol = 2)

grid.arrange( g[[4]], h[[4]], k[[4]], l[[4]], o[[4]], p[[4]],  ncol = 2)

 ############## ############## ############## 
############## OR-to-OR Contacts ############## 
 ############## ############## ############## 


ORs_in_GIH_datalist <- function(dipc_contacts, df_ORs_in_GIH){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_ORs_in_GIH %>% 
      filter(cell == i)
    print("ORs in  GIH for cell")
    print(df2$name)
    df4 <- inner_join(df_ORs_h, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    df5 <- inner_join(df_ORs_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name")
    print(df5)
    df6 <- df5 %>%
      dplyr::filter(prey_name %in% df2$name) %>%
      dplyr::filter(bait_name %in% df2$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    df7 <- left_join(order_tbl_df, df6)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    datalist[[i]] <- df7
  }
  return(datalist)
}  

ORs_in_GIH_datalist_active_P2 <- ORs_in_GIH_datalist(P2_dipc, activeGIH_ORs_percell_datalist_P2)
ORs_in_GIH_datalist_inactive_P2 <- ORs_in_GIH_datalist(P2_dipc, inactiveORs_percell_dataframe_P2)
ORs_in_GIH_datalist_active_MOR28 <- ORs_in_GIH_datalist(MOR28_dipc, activeGIH_ORs_percell_datalist_MOR28)
ORs_in_GIH_datalist_inactive_MOR28 <- ORs_in_GIH_datalist(MOR28_dipc, inactiveORs_percell_dataframe_MOR28)

big_data_P2 = do.call(rbind, ORs_in_GIH_datalist_active_P2)  %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data_MOR28 = do.call(rbind, ORs_in_GIH_datalist_active_MOR28) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30)%>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_OR_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, ORs_in_GIH_datalist_inactive_P2) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data_MOR28 = do.call(rbind, ORs_in_GIH_datalist_inactive_MOR28) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_OR_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, active_islands_datalist_P2) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data_MOR28 = do.call(rbind, active_islands_datalist_MOR28) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
active_island_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

big_data_P2 = do.call(rbind, inactive_islands_datalist_P2) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data_MOR28 = do.call(rbind, inactive_islands_datalist_MOR28) %>% 
  #filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30) %>%
  mutate("bait_order" = ceiling(bait_order/2)*2) %>%
  mutate("prey_order" = ceiling(prey_order/2)*2)
big_data <- rbind(big_data_P2, big_data_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
inactive_island_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_island_matrix[,-1]), 
                  max(active_island_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

aa <- pheatmap(active_island_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "GIs active GIH",
              breaks = myBreaks,
              cellheight=6, cellwidth =6)

bb <- pheatmap(inactive_island_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "GIs inactive GIH",
              breaks = myBreaks,
              cellheight=6, cellwidth =6)

cc <- pheatmap(active_OR_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "ORs active GIH",
               breaks = myBreaks,
               cellheight=6, cellwidth =6)

dd <- pheatmap(inactive_OR_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "ORs inactive GIH",
               breaks = myBreaks,
               cellheight=6, cellwidth =6)




grid.arrange(aa[[4]], bb[[4]],cc[[4]], dd[[4]],  ncol = 4)

############## ############## ############## ##############
############## active OR-to-OR Contacts ############## 
############## ############## ############################

activeOR_ORcomp_datalist <- function(dipc_contacts, df_activeORs, activeOR){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df3 <- df_activeORs %>% 
      filter(cell == i)
    print("ORs active in cell")
    print(df3$name)
    df2 <- df_ORs_h %>% filter(name == activeOR)
    df4 <- inner_join(df2, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    print(df4)
    df5 <- inner_join(df_ORs_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") %>%
      dplyr::filter(prey_name %in% df3$name) %>%
      group_by(bait_order, prey_order) %>%
      summarise(size = sum(size))
    print(df5)
    df7 <- left_join(order_tbl_df, df5)
    df7 <- df7 %>% replace(is.na(.), 0)
    print(df7)
    datalist[[i]] <- df7
  }
  return(datalist)
} 

oneOR_ORcomp_datalist <- function(dipc_contacts, df_ORs_in_GIH){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_ORs_in_GIH %>% 
      filter(cell == i) 
    if (nrow(df2) >= 1){
      df2 <- sample_n(df2, size = 1)
      print("ORs in  GIH for cell")
      print(df2$name)
      df4 <- inner_join(df_ORs_h, df1, by = c("loc" = "bait")) %>%
        dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
      df5 <- inner_join(df_ORs_h, df4, by = c("loc" = "prey")) %>%
        dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name")
      print(df5)
      df6 <- df5 %>%
        dplyr::filter(prey_name %in% df2$name) %>%
        dplyr::filter(bait_name %in% df2$name) %>%
        group_by(bait_order, prey_order) %>%
        summarise(size = sum(size))
      df7 <- left_join(order_tbl_df, df6)
      df7 <- df7 %>% replace(is.na(.), 0)
      print(df7)
      datalist[[i]] <- df7
    }
  }
  return(datalist)
}  


activeOR_ORcomp_datalist_P2 <- activeOR_ORcomp_datalist(P2_dipc, activeGIH_ORs_percell_datalist_P2, "OR-Cluster-26_Olfr714,OR-Cluster-26_Olfr17,0")
activeOR_ORcomp_datalist_MOR28 <- activeOR_ORcomp_datalist(MOR28_dipc, activeGIH_ORs_percell_datalist_MOR28, "OR-Cluster-55_Olfr1509,OR-Cluster-55_Olfr1508,OR-Cluster-55_Olfr1507,0")

big_data_P2 = do.call(rbind, activeOR_ORcomp_datalist_P2) #%>% filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30)
big_data_MOR28 = do.call(rbind, activeOR_ORcomp_datalist_MOR28) #%>% filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30)
big_data <- rbind(big_data_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
activeOR_ORcomp_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

oneOR_ORcomp_datalist_active_P2 <- oneOR_ORcomp_datalist(P2_dipc, activeGIH_ORs_percell_datalist_P2)
oneOR_ORcomp_datalist_active_MOR28 <- oneOR_ORcomp_datalist(MOR28_dipc, activeGIH_ORs_percell_datalist_MOR28)

big_data_P2 = do.call(rbind, oneOR_ORcomp_datalist_active_P2) #%>% filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30)
big_data_MOR28 = do.call(rbind, oneOR_ORcomp_datalist_active_MOR28) #%>% filter(bait_order > 10, bait_order < 30, prey_order > 10, prey_order < 30)
big_data <- rbind(big_data_P2)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
oneOR_ORcomp_matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))



paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(activeOR_ORcomp_matrix[,-1]), 
                  max(activeOR_ORcomp_matrix[,-1]), 
                  length.out=ceiling(paletteLength)))

ee <- pheatmap(activeOR_ORcomp_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "GIs active GIH",
              # breaks = myBreaks,
               cellheight=3, cellwidth =3)

ff <- pheatmap(oneOR_ORcomp_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "GIs inactive GIH",
              3breaks = myBreaks,
               cellheight=3, cellwidth =3)

grid.arrange(ee[[4]], ff[[4]], ncol = 2)
























# or_compartment_datalist <- function(dipc_contacts, normalization) {
#   datalist = list()
#   all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
#   for (i in all_cells$cell) {
#     print("Hello, it's cell ...")
#     print(i)
#     df1 <- dipc_contacts %>%
#       filter(cell == i, arch == "trans") %>%
#       select(bait, prey, size) # filter all trans contacts made in one cell
#     df3 <- df_ORs_h 
#     df5 <- left_join(df3, df1, by = c("loc" = "bait")) %>%
#       dplyr::rename("bait" = "loc") %>%
#       mutate("bait_OR" = name, "bait_order" = order) %>% 
#       select(bait, bait_order, bait_OR, prey, size) # prey is 2 MB region around all islands
#     df6 <- left_join(df3, df5, by = c("loc" = "prey")) %>%
#       dplyr::rename("prey" = "loc") %>%
#       mutate("prey_OR" = name, "prey_order" = order) %>% 
#       select(bait, bait_order, bait_OR, prey, prey_order, prey_OR, size)# bait is just the islands in the hub
#     df6 <- df6 %>% filter(bait_OR != prey_OR)
#     df7 <- df6[complete.cases(df6),]
#     df8 <- df7 %>% group_by(bait_order, prey_order) %>% 
#       summarise(size = sum(size))
#     df9 <- left_join(order_tbl_df, df8)
#     df9 <- df9 %>% replace(is.na(.), 0)
#     if (normalization == "square"){
#       df9$size <- df9$size/sum(df9$size)
#       datalist[[i]] <- df9
#     }
#     if (normalization == "none"){
#       datalist[[i]] <- df9
#     }
#     if (normalization == "row_bait"){
#       df9 <- df9 %>%
#         filter(bait_order == 20)
#       df9$size <- df9$size/sum(df9$size)
#       datalist[[i]] <- df9
#     }
#     if (normalization == "row_prey"){
#       df9 <- df9 %>%
#         filter(prey_order == 20)
#       df9$size <- df9$size/sum(df9$size)
#       datalist[[i]] <- df9
#     }
#   }
#   return(datalist)
# }  



inactiveOR_hub_datalist <- function(dipc_contacts, activeOR, df_activeislands, normalization){
  datalist = list()
  all_cells <- dipc_contacts %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_activeislands %>% filter(cell == i)
    print(df2)
    df4 <- inner_join(df_ORs_h, df1, by = c("loc" = "bait")) %>%
      dplyr::rename("bait" = "loc", "bait_order" = "order", "bait_name" = "name")
    df5 <- inner_join(df_Islands_h, df4, by = c("loc" = "prey")) %>%
      dplyr::rename("prey" = "loc", "prey_order" = "order", "prey_name" = "name") 
    active_islands_contacting <- df5 %>%
      dplyr::filter(bait_name %in% df2$name) #%>%
    #dplyr::filter(bait_name == activeOR) #%>%
    #dplyr::filter(bait_order == 20) %>%
    #dplyr::filter(prey_order == 20)
    print(active_islands_contacting)
    inactive_islands_contacting <- df5 %>%
      dplyr::filter(prey_name %out% df2$name) %>%
      dplyr::filter(bait_name != activeOR) %>%
      dplyr::filter(bait_order == 20) %>%
      dplyr::filter(prey_order == 20)
    print(inactive_islands_contacting)
    # set.seed(99)
    # if (nrow(inactive_islands_contacting) > nrow(active_islands_contacting)) {inactive_islands_contacting <- sample_n(inactive_islands_contacting, size = nrow(active_islands_contacting))}
    # inactive_contacts <- paste(inactive_islands_contacting$prey_name, inactive_islands_contacting$bait_name)
    # print(inactive_contacts)
    # df7 <- left_join(order_tbl_df, df5)
    # df7 <- df7 %>% replace(is.na(.), 0)
    # print(df7)
    # if (normalization == "square"){
    #   df7$size <- df7$size/sum(df7$size)
    #   datalist[[i]] <- df7
    # }
    # if (normalization == "none"){
    #   datalist[[i]] <- df7
    # }
    # if (normalization == "row_bait"){
    #   df7 <- df7 %>%
    #     filter(bait_order == 20)
    #   df7$size <- df7$size/sum(df7$size)
    #   datalist[[i]] <- df7
    # }
    # if (normalization == "row_prey"){
    #   df7 <- df7 %>%
    #     filter(prey_order == 20)
    #   df7$size <- df7$size/sum(df7$size)
    #   datalist[[i]] <- df7
    # }
  }
  return(datalist)
}

inactiveOR_hub_datalist(MOR28_dipc, "OR-Cluster-55_Olfr1509,OR-Cluster-55_Olfr1508,OR-Cluster-55_Olfr1507,1", df_activeislands_MOR28)
############## ############## ############## #####
### Establishing the OR Compartment Background ###
############## ############## ############## #####

or_compartment_datalist_MOR28 <- or_compartment_datalist(MOR28_dipc, "square")
or_compartment_datalist_P2 <- or_compartment_datalist(P2_dipc)
or_compartment_datalist_combined <- c(or_compartment_datalist_P2, or_compartment_datalist_MOR28)

big_data = do.call(rbind, or_compartment_datalist_combined)
big_data_sum <- big_data %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
matrix_OR <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))

############## ############## ############## ###
### ### Establishing Greek Island Signal ### ### 
############## ############## ############## ###

active_islands_datalist_MOR28 <- active_islands_datalist(MOR28_dipc, df_activeislands_MOR28, "square")
active_islands_datalist_P2 <- active_islands_datalist(P2_dipc, df_activeislands_P2)
active_islands_datalist_combined <- c(active_islands_datalist_MOR28, active_islands_datalist_P2)

big_data = do.call(rbind, active_islands_datalist_combined)
big_data_sum <- big_data %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))
matrix[is.na(matrix)] <- min(matrix, na.rm = T)
matrix <- log(matrix/matrix_OR)

paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2)), 
              seq(0.06, 4, length.out=floor(paletteLength/2)))

pheatmap(matrix[11:29, 12:30], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = myColor,
         cellheight=12, cellwidth = 12,
         breaks = myBreaks)

big_data_percell <- big_data %>% filter(bait_order == 20) %>% select(cell, prey_order, size)
matrix_percell <- dcast(big_data_percell, cell ~ prey_order)
matrix_percell[is.na(matrix_percell)] <- 0
matrix_percell[2:29] <- matrix_percell[2:29]/rowSums(matrix_percell[2:29])
matrix_percell <- matrix_percell[order(-matrix_percell$`20`),]
nrow(matrix_percell)

breaksList_2 = seq(0, 1, by = 0.02)

pheatmap(matrix_percell[7:25], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = colorRampPalette(c("white", "firebrick3"))(50),
         cellheight=4.5, cellwidth = 7.5,
         breaks = breaksList_2)


#### Inactive Hub

inactive_islands_datalist_MOR28 <- inactive_islands_datalist(MOR28_dipc, df_activeislands_MOR28)
big_data_MOR28 = do.call(rbind, inactive_islands_datalist_MOR28)
big_data_MOR28$OR <- "MOR28"

inactive_islands_datalist_P2 <- inactive_islands_datalist(P2_dipc, df_activeislands_P2)
big_data_P2 = do.call(rbind, inactive_islands_datalist_P2)
big_data_P2$OR <- "P2"

big_data = rbind(big_data_MOR28, big_data_P2)
big_data_sum <- big_data %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
matrix <- dcast(big_data_sum, bait_order ~ prey_order)/(sum(big_data_sum$sum))
matrix[is.na(matrix)] <- min(matrix, na.rm = T)
matrix <- log(matrix/matrix_OR)


paletteLength <- 100
myColor <- colorRampPalette(c("blue3", "white", "firebrick3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2)), 
              seq(0.06, 4, length.out=floor(paletteLength/2)))

pheatmap(matrix[11:29, 12:30], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = myColor,
         cellheight=12, cellwidth = 12,
         breaks = myBreaks)


big_data_percell <- big_data %>% filter(bait_order == 20) %>% select(cell, prey_order, size, OR)
matrix_percell <- dcast(big_data_percell, cell + OR ~ prey_order, value.var = "size")
matrix_percell[is.na(matrix_percell)] <- 0
matrix_percell[3:39] <- matrix_percell[3:39]/rowSums(matrix_percell[3:39])
matrix_percell <- matrix_percell[order(-matrix_percell$`20`),]
nrow(matrix_percell)

breaksList_2 = seq(0, 1, by = 0.02)

pheatmap(matrix_percell[13:30], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         col = colorRampPalette(c("white", "firebrick3"))(50),
         cellheight=4.5*(56/110), cellwidth = 7.5,
         breaks = breaksList_2)


#######################
## Island Insulation ##
#######################

hic_contacts_OneIslandtoOneIsland <- merge(dipc, Islands_50kb_bait, by = c("bait", "bait_haplo")) %>%
  merge(., TwoMB.Islands, by = c("prey","prey_haplo"), allow.cartesian=TRUE) 

## Intersection
Island.Island.Matrix <- dcast(hic_contacts_OneIslandtoOneIsland, bait + bait_haplo + name + prey_haplo + geno + cell + OR ~ order, value.var = "size") # this matrix will be separated into all 6 libraries in order to increase the amount you are sampling from P2
Island.Island.Matrix[is.na(Island.Island.Matrix)] <- 0

# Active vs. Inactive Islands
Island.Island.Matrix$prey_identity <- paste(Island.Island.Matrix$name, Island.Island.Matrix$prey_haplo, Island.Island.Matrix$cell, sep = "_")
Island.Island.Matrix$bait_identity <- paste(Island.Island.Matrix$bait, Island.Island.Matrix$bait_haplo, Island.Island.Matrix$cell, sep = "_")

Island.Island.Matrix.MOR28 <- Island.Island.Matrix %>% filter(OR == "mor28")
Island.Island.Matrix.P2 <- Island.Island.Matrix %>% filter(OR == "P2")

Island.Island.Matrix.MOR28.inactive <- Island.Island.Matrix.MOR28 %>% filter(Island.Island.Matrix.MOR28$bait_identity %out% islands_mat_MOR28_ImmFix_dipc$bait_identity)
Island.Island.Matrix.MOR28.active <- Island.Island.Matrix.MOR28 %>% filter(Island.Island.Matrix.MOR28$bait_identity %in% islands_mat_MOR28_ImmFix_dipc$bait_identity)

Island.Island.Matrix.MOR28.inactive <- Island.Island.Matrix.MOR28.inactive %>% filter(Island.Island.Matrix.MOR28.inactive$prey_identity %out% islands_mat_MOR28_ImmFix_dipc$prey_identity)
Island.Island.Matrix.MOR28.active <- Island.Island.Matrix.MOR28.active %>% filter(Island.Island.Matrix.MOR28.active$prey_identity %in% islands_mat_MOR28_ImmFix_dipc$prey_identity)

Island.Island.Matrix.P2.inactive <- Island.Island.Matrix.P2 %>% filter(Island.Island.Matrix.P2$bait_identity %out% islands_mat_P2_ImmFix_dipc$bait_identity)
Island.Island.Matrix.P2.active <- Island.Island.Matrix.P2 %>% filter(Island.Island.Matrix.P2$bait_identity %in% islands_mat_P2_ImmFix_dipc$bait_identity)

Island.Island.Matrix.P2.inactive <- Island.Island.Matrix.P2.inactive %>% filter(Island.Island.Matrix.P2.inactive$prey_identity %out% islands_mat_P2_ImmFix_dipc$prey_identity)
Island.Island.Matrix.P2.active <- Island.Island.Matrix.P2.active %>% filter(Island.Island.Matrix.P2.active$prey_identity %in% islands_mat_P2_ImmFix_dipc$prey_identity)



## Active Island - to - Active Island ##
MOR28 <- Island.Island.Matrix.MOR28.active
MOR28 <- MOR28[,c(6,8:46)] 
MOR28 <- MOR28 %>% group_by(cell) %>% summarise_each(funs(sum)) 
MOR28 <- MOR28[,2:40] %>% as.matrix()
MOR28.colMeans <- MOR28 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
MOR28.colMeans$order <- c(1:39)
MOR28.colMeans$contact <- MOR28.colMeans$contact/sum(MOR28.colMeans$contact) # normalize neighborhood specificity by all
MOR28.colMeans$olfr <- "mor28"

P2 <- Island.Island.Matrix.P2.active %>% filter(OR == "P2")
P2 <- P2[,c(6,8:46)] 
P2 <- P2 %>% group_by(cell) %>% summarise_each(funs(sum)) 
P2 <- P2[,2:40] %>% as.matrix()
P2.colMeans <- P2 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
P2.colMeans$order <- c(1:39)
P2.colMeans$contact <- P2.colMeans$contact/sum(P2.colMeans$contact) # normalize neighborhood specificity by all
P2.colMeans$olfr <- "P2"

summary <- rbind(MOR28.colMeans, P2.colMeans)
geno.levels <- c("mor28", "P2")
summary$geno <- factor(summary$geno,levels=geno.levels, ordered=TRUE)

a <- ggplot(summary, aes(x = order, y = contact, fill = olfr)) + 
  geom_point() + 
  geom_line() + 
  geom_area(position = "identity", alpha = 0.1) +
  theme_classic() +
  theme(legend.position="bottom") +
  ylim(0, 0.32) +
  ggtitle("Active Islands") +
  ylab("Normalized Contacts") +
  xlab("Average Greek Island") 

## inactive Island - to - inactive Island ##
MOR28 <- Island.Island.Matrix.MOR28.inactive
MOR28 <- MOR28[,c(6,8:46)] 
MOR28 <- MOR28 %>% group_by(cell) %>% summarise_each(funs(sum)) 
MOR28 <- MOR28[,2:40] %>% as.matrix()
MOR28.colMeans <- MOR28 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
MOR28.colMeans$order <- c(1:39)
MOR28.colMeans$contact <- MOR28.colMeans$contact/sum(MOR28.colMeans$contact) # normalize neighborhood specificity by all
MOR28.colMeans$olfr <- "mor28"

P2 <- Island.Island.Matrix.P2.inactive %>% filter(OR == "P2")
P2 <- P2[,c(6,8:46)] 
P2 <- P2 %>% group_by(cell) %>% summarise_each(funs(sum)) 
P2 <- P2[,2:40] %>% as.matrix()
P2.colMeans <- P2 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
P2.colMeans$order <- c(1:39)
P2.colMeans$contact <- P2.colMeans$contact/sum(P2.colMeans$contact) # normalize neighborhood specificity by all
P2.colMeans$olfr <- "P2"

summary <- rbind(MOR28.colMeans, P2.colMeans)
geno.levels <- c("mor28", "P2")
summary$geno <- factor(summary$geno,levels=geno.levels, ordered=TRUE)

b <- ggplot(summary, aes(x = order, y = contact, fill = olfr)) + 
  geom_point() + 
  geom_line() + 
  geom_area(position = "identity", alpha = 0.1) +
  theme_classic() +
  theme(legend.position="bottom") +
  ylim(0, 0.32) +
  ggtitle("Inactive Islands") +
  ylab("Normalized Contacts") +
  xlab("Average Greek Island") 

# ## ## ## ## ## ## 
## OR Insulation ##  
## ## ## ## ## ## # 

# Two MB OR Filter

options(scipen = 999)

OR.Names <- read.table("/media/storageA/kevin/annotation/ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% select(V1, V2, V4) %>%
  mutate(chr = str_sub(V1, 4, -1))
OR.Names$prey <- paste(OR.Names$chr, OR.Names$V2, sep = "_")
OR.Names$type <- OR.Names$V4
OR.Names.bait<- OR.Names %>% rename("bait" = "prey") #not working for some reason

Bins.And.ORs <- left_join(All_Bins_prey, OR.Names)
Bins.And.ORs <- Bins.And.ORs %>% tidyr::separate('prey', into = c('prey_chr','prey_loc'), sep = '_')  %>% #split prey into chr and loc
  mutate_at(.vars = vars(matches('loc')), .funs = as.numeric) #convert "loc" column to numeric
ORs <- Bins.And.ORs[complete.cases(Bins.And.ORs),]

## Filter out ORs within 2 MB of an OR ##

for (i in seq_along(ORs[,1])) {
  if (i == 1) {
    TwoMB.ORs <- Bins.And.ORs %>% filter(Bins.And.ORs$prey_chr == ORs$prey_chr[i]) %>%
      filter((prey_loc > (ORs$prey_loc[i] - 1e6)) & (prey_loc < (ORs$prey_loc[i] + 1e6)))
    TwoMB.ORs$name <- ORs$type[i]
    TwoMB.ORs$order <- c(1:39)
  } else {
    tryCatch({
      TwoMB.ORs.loop <- Bins.And.ORs %>% filter(Bins.And.ORs$prey_chr == ORs$prey_chr[i]) %>%
        filter((prey_loc > (ORs$prey_loc[i] - 1e6)) & (prey_loc < (ORs$prey_loc[i] + 1e6)))
      TwoMB.ORs.loop$name <- ORs$type[i]
      TwoMB.ORs.loop$order <- c(1:39)
      TwoMB.ORs <- rbind(TwoMB.ORs, TwoMB.ORs.loop)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

TwoMB.ORs$prey <- paste(TwoMB.ORs$prey_chr, TwoMB.ORs$prey_loc, sep = "_")
TwoMB.ORs <- TwoMB.ORs %>% select(prey, type, name, order)

TwoMB.ORs_mat <- TwoMB.ORs
TwoMB.ORs_mat$prey_haplo <- 1
TwoMB.ORs_pat <- TwoMB.ORs
TwoMB.ORs_pat$prey_haplo <- 0
TwoMB.ORs <- rbind (TwoMB.ORs_mat, TwoMB.ORs_pat)



hic_contacts_OneORtoOneOR <- merge(dipc, OR_Clusters_50kb_bait, by = c("bait", "bait_haplo")) %>%
  merge(., TwoMB.ORs, by = c("prey","prey_haplo"), allow.cartesian=TRUE) 

## Intersection
OR.OR.Matrix <- dcast(hic_contacts_OneORtoOneOR, bait + bait_haplo + name + prey_haplo + geno + cell + OR ~ order, value.var = "size") # this matrix will be separated into all 6 libraries in order to increase the amount you are sampling from P2
OR.OR.Matrix[is.na(OR.OR.Matrix)] <- 0


## gg8tta TETOP2 HET ##
MOR28 <- OR.OR.Matrix %>% filter(OR == "mor28")
MOR28 <- MOR28[,c(6,8:46)] 
MOR28 <- MOR28 %>% group_by(cell) %>% summarise_each(funs(sum)) 
MOR28 <- MOR28[,2:40] %>% as.matrix()
MOR28.colMeans <- MOR28 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
MOR28.colMeans$order <- c(1:39)
MOR28.colMeans$contact <- MOR28.colMeans$contact/sum(MOR28.colMeans$contact) # normalize neighborhood specificity by all
MOR28.colMeans$olfr <- "mor28"

P2 <- OR.OR.Matrix %>% filter(OR == "P2")
P2 <- P2[,c(6,8:46)] 
P2 <- P2 %>% group_by(cell) %>% summarise_each(funs(sum)) 
P2 <- P2[,2:40] %>% as.matrix()
P2.colMeans <- P2 %>% colMeans() %>% as.data.frame() %>% dplyr::rename("contact" = ".")
P2.colMeans$order <- c(1:39)
P2.colMeans$contact <- P2.colMeans$contact/sum(P2.colMeans$contact) # normalize neighborhood specificity by all
P2.colMeans$olfr <- "P2"

summary <- rbind(MOR28.colMeans, P2.colMeans)
geno.levels <- c("mor28", "P2")
summary$geno <- factor(summary$geno,levels=geno.levels, ordered=TRUE)

c <- ggplot(summary, aes(x = order, y = contact, fill = olfr)) + 
  geom_point() + 
  geom_line() + 
  geom_area(position = "identity", alpha = 0.1) +
  theme_classic() +
  theme(legend.position="bottom") +
  ylim(0, 0.32) +
  ggtitle("ORs") +
  ylab("Normalized Contacts") +
  xlab("Average OR") 

grid.arrange(a, b, c, ncol = 3, top = "Average Island-to-Island or OR-to-OR Interactions \n mor28-IRES-GFP and gg8-tTA > tetO-P2 Dip-C")


