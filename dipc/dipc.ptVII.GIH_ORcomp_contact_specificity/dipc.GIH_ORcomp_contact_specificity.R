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

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity") 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

##############################################################################
# Dataframe Generation
# The initial part of the script is used to generate dataframe that are used
# for examining contacts around OR genes and GIs
##############################################################################

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
# 


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

#########################################
# Loading Dip-C contacts Binned at 50-kb Resolution
#########################################

# Loading P2 Data
P2_ImmFix_dipc_male <- fread(file = 'P2-ImmFix-dipc.male.50kb.contacts.haplotyped.txt') %>% as.data.frame()
P2_ImmFix_dipc_female <- fread(file = 'P2-ImmFix-dipc.female.50kb.contacts.haplotyped.txt') %>% as.data.frame()
P2_ImmFix_dipc <- rbind(P2_ImmFix_dipc_male, P2_ImmFix_dipc_female)
P2_ImmFix_dipc$geno <- "ctrl"
rm(P2_ImmFix_dipc_male, P2_ImmFix_dipc_female)
P2_dipc <- P2_ImmFix_dipc %>% as.data.table()
rm(P2_ImmFix_dipc)

# Loading MOR28 Data
MOR28_ImmFix_dipc_male <- fread(file = 'APM28DIP.male.50kb.contacts.haplotyped.txt') %>% as.data.frame()
MOR28_ImmFix_dipc_female <- fread(file = 'APM28DIP.female.50kb.contacts.haplotyped.txt') %>% as.data.frame()
MOR28_ImmFix_dipc <- rbind(MOR28_ImmFix_dipc_male, MOR28_ImmFix_dipc_female)
MOR28_ImmFix_dipc$geno <- "ctrl"
rm(MOR28_ImmFix_dipc_male, MOR28_ImmFix_dipc_female)
MOR28_dipc <- MOR28_ImmFix_dipc %>% as.data.table()
rm(MOR28_ImmFix_dipc)

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


# Generation of dataframes of Island and OR gene composition
# of GIHs in both Dip-C Datasets

df_ORs_h <- TwoMB.ORs_hap %>% 
  mutate("loc" = prey) %>%
  select(loc, order, name)

df_Islands_h <- TwoMB.Islands_hap %>% 
  mutate("loc" = prey) %>%
  select(loc, order, name)


############################################################################
############################################################################

        # Formatting of DFs of Greek Islands in Greek Island Hubs #

############################################################################
############################################################################


islands_percell_dataframe_P2 <- read.table("islands_percell_dataframe_P2_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell)

islands_percell_dataframe_MOR28 <- read.table("islands_percell_dataframe_MOR28_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell)

############################################################################
############################################################################

        # Formatting of DFs of ORs in Inactive Greek Island Hubs #

############################################################################
############################################################################

OR_formatting_df_ptI <- read.table("ORs-by-cluster.bed") %>%
  mutate(chr = V1, loc = (floor(V3/bin_size))*bin_size, cluster = V5, or = V4) %>%
  separate(chr, sep = "r", into = c("chr_1", "chr_2")) %>%
  mutate(prey = paste(chr_2, loc, sep = "_")) %>%
  mutate(type = paste(cluster, or, sep = "_")) %>%
  select(prey, type) %>%
  separate(col = "type", sep = "_", into  =c("cluster", "name")) %>%
  separate(col = "prey", sep = "_", into = c("prey_chr", "prey_loc")) %>%
  haplotype() %>% select(prey, name) 

OR_formatting_df_ptII <- read.table("ORs-by-cluster.bed") %>%
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

############################################################################
############################################################################

inactiveORs_percell_dataframe_P2 <-  read.table("inactiveORs_percell_dataframe_P2_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "OR", "haplo")) %>%
  mutate("name" = paste(OR, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., OR_formatting_df_ptI, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

inactiveORs_percell_dataframe_MOR28 <-  read.table("inactiveORs_percell_dataframe_MOR28_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., OR_formatting_df_ptI, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

inactiveORs_percell_dataframe_P2 <- inactiveORs_percell_dataframe_P2 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., OR_formatting_df_ptII, by = "prey") %>%
  dplyr::select(cell, prey, name)

inactiveORs_percell_dataframe_MOR28 <- inactiveORs_percell_dataframe_MOR28 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., OR_formatting_df_ptII, by = "prey") %>%
  dplyr::select(cell, prey, name)


############################################################################
############################################################################

          # Formatting of DFs of ORs in Active Greek Island Hubs #

############################################################################
############################################################################

activeORs_percell_dataframe_P2 <-  read.table("activeORs_percell_dataframe_P2_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "OR", "haplo")) %>%
  mutate("name" = paste(OR, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., OR_formatting_df_ptI, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

activeORs_percell_dataframe_MOR28 <-  read.table("activeORs_percell_dataframe_MOR28_230927.txt",sep = "\t") %>%
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, cell) %>%
  left_join(., OR_formatting_df_ptI, by = "name") %>% 
  group_by(cell, prey) %>%
  mutate(name = paste(name, collapse = ",")) %>% 
  unique()

activeORs_percell_dataframe_P2 <- activeORs_percell_dataframe_P2 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., OR_formatting_df_ptII, by = "prey") %>%
  dplyr::select(cell, prey, name)

activeORs_percell_dataframe_MOR28 <- activeORs_percell_dataframe_MOR28 %>%
  mutate("bad_names" = name) %>%
  dplyr::select(cell, prey, bad_names) %>%
  left_join(., OR_formatting_df_ptII, by = "prey") %>%
  dplyr::select(cell, prey, name)



############################################################################
############################################################################
############################################################################
############################################################################

order_tbl_df <- data.frame(bait_order = 1:39) 
order_tbl_df <- expand.grid(bait_order = order_tbl_df$bait_order, prey_order = order_tbl_df$bait_order) %>% as_tibble()

############## ############## ############## 
############## GI-to-GI Contacts ########### 
############## ############## ############## 

islands_datalist <- function(dipc_contacts, df_activeislands, activity){
  master = data.frame()
  master_QC = data.frame()
  all_cells <- P2_dipc %>% dplyr::select(cell) %>% unique()
  for (i in all_cells$cell) {
    print(c("Hello, it's cell ...",i))
    df1 <- dipc_contacts %>%
      filter(cell == i, arch == "trans") %>%
      select(bait, prey, size, arch) # filter all trans contacts made in one cell
    df2 <- df_activeislands %>% 
      filter(cell == i, state == activity)
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
    master <- rbind(master, df7)
    df_QC <- data.frame(cell = c(i), numberOfIslands = c(length(df2$name)), contactsPerMillion =c(sum(df6$size)*1e6) )
    master_QC <- rbind(master_QC, df_QC)
  }
  return(master)
}  


############## ############## ############## 
############## OR-to-OR Contacts ########### 
############## ############## ############## 


ORs_in_GIH_datalist <- function(dipc_contacts, df_ORs_in_GIH){
  master = data.frame()
  master_QC = data.frame()
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
    master <- rbind(master, df7)
    df_QC <- data.frame(cell = c(i), numberOfIslands = c(length(df2$name)), contactsPerMillion =c(sum(df6$size)*1e6) )
    master_QC <- rbind(master_QC, df_QC)
  }
  return(master)
}  


# Data Analysis

active_islands_datalist_P2 <- islands_datalist(P2_dipc, islands_percell_dataframe_P2, "active")
inactive_islands_datalist_P2 <- islands_datalist(P2_dipc, islands_percell_dataframe_P2, "inactive")
active_islands_datalist_MOR28 <- islands_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "active")
inactive_islands_datalist_MOR28 <- islands_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "inactive")

ORs_in_GIH_datalist_active_P2 <- ORs_in_GIH_datalist(P2_dipc, activeORs_percell_dataframe_P2)
ORs_in_GIH_datalist_inactive_P2 <- ORs_in_GIH_datalist(P2_dipc, inactiveORs_percell_dataframe_P2)
ORs_in_GIH_datalist_active_MOR28 <- ORs_in_GIH_datalist(MOR28_dipc, activeORs_percell_dataframe_MOR28)
ORs_in_GIH_datalist_inactive_MOR28 <- ORs_in_GIH_datalist(MOR28_dipc, inactiveORs_percell_dataframe_MOR28)

### Collecting Data ###

big_data <- rbind(ORs_in_GIH_datalist_active_P2, ORs_in_GIH_datalist_active_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
active_OR_matrix <- dcast(big_data_sum, bait_order ~ prey_order)

big_data <- rbind(ORs_in_GIH_datalist_inactive_P2, ORs_in_GIH_datalist_inactive_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
inactive_OR_matrix <- dcast(big_data_sum, bait_order ~ prey_order)

big_data <- rbind(active_islands_datalist_P2, active_islands_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
active_island_matrix <- dcast(big_data_sum, bait_order ~ prey_order)

big_data <- rbind(inactive_islands_datalist_P2, inactive_islands_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>% summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
inactive_island_matrix <- dcast(big_data_sum, bait_order ~ prey_order)

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_island_matrix[,-1]), 
                  active_island_matrix[20,21], 
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
               cellheight=3, cellwidth =3)

bb <- pheatmap(inactive_island_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "GIs inactive GIH",
               breaks = myBreaks,
               cellheight=3, cellwidth =3)

cc <- pheatmap(active_OR_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "ORs active GIH",
               breaks = myBreaks,
               cellheight=3, cellwidth =3)

dd <- pheatmap(inactive_OR_matrix[,-1], 
               cluster_rows=FALSE, 
               cluster_cols=FALSE, 
               border_color = NA,
               show_rownames = FALSE,
               show_colnames = FALSE,
               col = myColor,
               main = "ORs inactive GIH",
               breaks = myBreaks,
               cellheight=3, cellwidth =3)




# Figure 2F
grid.arrange(aa[[4]], bb[[4]],cc[[4]], dd[[4]],  ncol = 4)

write.table(active_island_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/Fig2e_ActiveGItoGIContacts.txt",
            sep = '\t', col.names = FALSE)
write.table(inactive_island_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/Fig2f_InactiveGItoGIContacts.txt",
            sep = '\t', col.names = FALSE)
write.table(active_OR_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/SuppFig6b_ORtoORContactsActiveGIH.txt",
            sep = '\t', col.names = FALSE)
write.table(inactive_OR_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/SuppFig6c_ORtoORContactsInactiveGIH.txt",
            sep = '\t', col.names = FALSE)

library(plotly)

#############################################################################
#############################################################################

            ## Extracting OR and GI contacts within the same hubs ##

#############################################################################
#############################################################################


####### ####### FUNCTIONS ####### #######

# examine active OR - GI contacts in the active GIH

activeOR_hub_datalist <- function(dipc_contacts,df_activeislands, activeOR){
  master = data.frame()
  master_QC = data.frame()
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
    df7$cell <- i
    print(df7)
    master <- rbind(master, df7)
    df_QC <- data.frame(cell = c(i), numberOfIslands = c(length(df3$name)), contactsPerMillion =c(sum(df5$size)*1e6) )
    master_QC <- rbind(master_QC, df_QC)
  }
  return(master)
} 

# examine inactive OR inactive GI contacts in inactive GIH of same size as active GIH

inactiveOR_hub_datalist <- function(dipc_contacts, df_activeislands, df_inactiveORs_original, normalization){
  master = data.frame()
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
    df_inactiveORs <- df_inactiveORs_original %>% filter(cell == i)
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
    df7$cell <- i
    print(df7)
    master <- rbind(master, df7)
  }
  return(master)
} 

activeOR_hub_datalist_P2 <- activeOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, "OR-Cluster-26_Olfr714,OR-Cluster-26_Olfr17,0")
activeOR_hub_datalist_MOR28 <- activeOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, "OR-Cluster-55_Olfr1509,OR-Cluster-55_Olfr1508,OR-Cluster-55_Olfr1507,0")

inactiveOR_hub_datalist_P2 <- inactiveOR_hub_datalist(P2_dipc, islands_percell_dataframe_P2, inactiveORs_percell_dataframe_P2)
inactiveOR_hub_datalist_MOR28 <- inactiveOR_hub_datalist(MOR28_dipc, islands_percell_dataframe_MOR28, inactiveORs_percell_dataframe_MOR28)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


big_data <- rbind(activeOR_hub_datalist_P2, activeOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>%
  summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
active_matrix <- dcast(big_data_sum, bait_order ~ prey_order)


big_data <- rbind(inactiveOR_hub_datalist_P2, inactiveOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>%
  summarise(sum = sum(size))
big_data_sum$sum <- big_data_sum$sum/sum(big_data_sum$sum)
big_data_sum %>% filter(bait_order == 20, prey_order == 20)
inactive_matrix <- dcast(big_data_sum, bait_order ~ prey_order)

paletteLength <- 100
myColor <- colorRampPalette(c("white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(active_matrix[,-1]), 
                  active_matrix[20,21], 
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
              cellheight=3, cellwidth =3)

f <- pheatmap(inactive_matrix[,-1], 
              cluster_rows=FALSE, 
              cluster_cols=FALSE, 
              border_color = NA,
              show_rownames = FALSE,
              show_colnames = FALSE,
              col = myColor,
              main = "Inactive GIHs",
              breaks = myBreaks,
              cellheight=3, cellwidth =3)


#grid.arrange(e[[4]], f[[4]], ncol = 2)

grid.arrange(aa[[4]], bb[[4]], e[[4]], f[[4]], ncol = 4)

grid.arrange(cc[[4]], dd[[4]], ncol = 2)

write.table(active_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/Fig2g_ActiveORtoActiveGIContacts.txt",
            sep = '\t', col.names = FALSE)
write.table(inactive_matrix[,-1], file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/Fig2h_InactiveORtoInactiveGIContacts.txt",
            sep = '\t', col.names = FALSE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

big_data <- rbind(activeOR_hub_datalist_P2, activeOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>%
  summarise(sum = sum(size))
big_data_sum_column <- big_data_sum %>% dplyr::filter(prey_order == 20)
big_data_sum_column$sum <- big_data_sum_column$sum/sum(big_data_sum_column$sum)
big_data_sum_column %>% filter(bait_order == 20, prey_order == 20)
activeGIH_sum_column <- big_data_sum_column
activeGIH_sum_column$geno <- 'active'


big_data <- rbind(inactiveOR_hub_datalist_P2, inactiveOR_hub_datalist_MOR28)
big_data_sum <- big_data[is.finite(big_data$size),] %>% group_by(bait_order, prey_order) %>%
  summarise(sum = sum(size))
big_data_sum_column <- big_data_sum %>% dplyr::filter(prey_order == 20)
big_data_sum_column$sum <- big_data_sum_column$sum/sum(big_data_sum_column$sum)
big_data_sum_column %>% filter(bait_order == 20, prey_order == 20)
inactiveGIH_sum_column <- big_data_sum_column
inactiveGIH_sum_column$geno <- 'inactive'

z <- ggplot(activeGIH_sum_column, aes(x = bait_order, y = sum, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("green4")) +
  scale_fill_manual(values=c("green4")) +
  theme_classic() + 
  ylim(0, max(activeGIH_sum_column$sum)) + 
  theme_classic() + 
  ylab("Normalized Contacts") + 
  xlab("P2 locus ± 1Mb") + 
  ggtitle("Active GI to Active OR contacts")

y <- ggplot(inactiveGIH_sum_column, aes(x = bait_order, y = sum, color = geno, fill = geno)) + 
  geom_line() + 
  geom_point() + 
  geom_area(position = "identity", alpha = 0.05) + 
  scale_color_manual(values=c("firebrick3")) +
  scale_fill_manual(values=c("firebrick3")) +
  theme_classic() + 
  ylim(0, max(activeGIH_sum_column$sum)) + 
  theme_classic() + 
  ylab("Normalized Contacts") + 
  xlab("Max Inactive OR locus ± 1Mb") + 
  ggtitle("Inactive GI to Inactive OR contacts")

z + y

GIH_sum_column <- rbind(activeGIH_sum_column, inactiveGIH_sum_column)
write.table(GIH_sum_column, file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.GIH_ORcomp_contact_specificity/SuppFig6_LinePlots.txt")

