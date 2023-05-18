library(pheatmap)
library(dendextend)
library(ggplot2)
library(tidyr)
library(dplyr)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_finalized/dipc/dipc.hubs.per.cell_2D") 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

#####################################################################################################
# The code below is used to extract the average/stdev, and various ranges of the # GIs per active GIH.
#####################################################################################################

names <- read.table("/media/storageE/ariel_dipc/annotations/Greek_Islands.nonX.names")

islands_percell_dataframe_P2_mean <- read.table("islands_percell_dataframe_P2.txt",sep = "\t") %>% # generated using dipc.inactive.hub.extraction.R
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell) %>%
  dplyr::filter(state == "active") %>%
  mutate(geno = "P2")

islands_percell_dataframe_MOR28_mean <- read.table("islands_percell_dataframe_MOR28.txt",sep = "\t") %>% # generated using dipc.inactive.hub.extraction.R
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell) %>%
  dplyr::filter(state == "active") %>%
  mutate(geno = "MOR28")

islands_percell_dataframe <- rbind(islands_percell_dataframe_P2_mean, islands_percell_dataframe_MOR28_mean) %>%
  group_by(cell, geno, state) %>%
  summarise(sum = n()) %>%
  group_by(state) %>%
  summarise(mean = mean(sum),
            sd = sd(sum),
            lower_1sd = mean - sd,
            upper_1sd = mean + sd,
            lower_2sd = mean - 2*sd,
            upper_2sd = mean + 2*sd)
islands_percell_dataframe

###################################################################################
# The code below is used to extract the # of GIHs per cell that fall within
# specific ranges of potential GIH size
###################################################################################

#################
### Functions ###
#################

cluster_count_mean <- function(paired_distance_path){
  df <- read.table(paired_distance_path,col.names = names$V1)
  rownames(df) <- names$V1
  df_sc <- as.data.frame(scale(df))
  dist_mat <- dist(df_sc, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  plot(hclust_avg)
  cut_avg <- cutree(hclust_avg, h = 2.5)
  cut_avg_df <- cut_avg %>% as.data.frame()
  cut_avg_df <- cut_avg_df %>% dplyr::rename("cluster" = ".")
  valid_clusters <- cut_avg_df %>% group_by(cluster) %>% 
    summarise(total = n()) %>% 
    filter(total >= islands_percell_dataframe$mean)
  a <- nrow(valid_clusters)
  return(a)
}

cluster_count_1sd <- function(paired_distance_path){
  df <- read.table(paired_distance_path,col.names = names$V1)
  rownames(df) <- names$V1
  df_sc <- as.data.frame(scale(df))
  dist_mat <- dist(df_sc, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  plot(hclust_avg)
  cut_avg <- cutree(hclust_avg, h = 2.5)
  cut_avg_df <- cut_avg %>% as.data.frame()
  cut_avg_df <- cut_avg_df %>% dplyr::rename("cluster" = ".")
  valid_clusters <- cut_avg_df %>% group_by(cluster) %>% 
    summarise(total = n()) %>% 
    filter(total >= islands_percell_dataframe$lower_1sd,
           total <= islands_percell_dataframe$upper_1sd)
  a <- nrow(valid_clusters)
  return(a)
}

cluster_count_2sd <- function(paired_distance_path){
  df <- read.table(paired_distance_path,col.names = names$V1)
  rownames(df) <- names$V1
  df_sc <- as.data.frame(scale(df))
  dist_mat <- dist(df_sc, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  plot(hclust_avg)
  cut_avg <- cutree(hclust_avg, h = 2.5)
  cut_avg_df <- cut_avg %>% as.data.frame()
  cut_avg_df <- cut_avg_df %>% dplyr::rename("cluster" = ".")
  valid_clusters <- cut_avg_df %>% group_by(cluster) %>% 
    summarise(total = n()) %>% 
    filter()
  a <- nrow(valid_clusters)
  return(a)
}

############################################
############################################
############################################
############################################

###########################
# gg8-tTA > tetO-P2 Dip-C #
###########################

names <- read.table("Greek_Islands.nonX.names")
`%notin%` <- Negate(`%in%`)

male_files <- list.files(path = "matrix.P2-ImmFix-dipc.male.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                         "*.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.P2-ImmFix-dipc.female.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                           "*.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

male_cluster_count_mean <- lapply(male_files, cluster_count_mean) 
female_cluster_count_mean <- lapply(female_files, cluster_count_mean)
male_cluster_count_mean <- as.data.frame(do.call(rbind, male_cluster_count_mean))
female_cluster_count_mean <- as.data.frame(do.call(rbind, female_cluster_count_mean))
cluster_count_mean_p2 <- rbind(male_cluster_count_mean, female_cluster_count_mean)
cluster_count_mean_p2$geno <- "p2"

male_cluster_count_1sd <- lapply(male_files, cluster_count_1sd) 
female_cluster_count_1sd <- lapply(female_files, cluster_count_1sd)
male_cluster_count_1sd <- as.data.frame(do.call(rbind, male_cluster_count_1sd))
female_cluster_count_1sd <- as.data.frame(do.call(rbind, female_cluster_count_1sd))
cluster_count_1sd_p2 <- rbind(male_cluster_count_1sd, female_cluster_count_1sd)
cluster_count_1sd_p2$geno <- "p2"

male_cluster_count_2sd <- lapply(male_files, cluster_count_2sd) 
female_cluster_count_2sd <- lapply(female_files, cluster_count_2sd)
male_cluster_count_2sd <- as.data.frame(do.call(rbind, male_cluster_count_2sd))
female_cluster_count_2sd <- as.data.frame(do.call(rbind, female_cluster_count_2sd))
cluster_count_2sd_p2 <- rbind(male_cluster_count_2sd, female_cluster_count_2sd)
cluster_count_2sd_p2$geno <- "p2"

###################
# mor28iGFP Dip-C #
###################

female_files <- list.files(path = "matrix.mor28.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                           "APM28DIP.*.female.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

male_files <- list.files(path = "matrix.mor28.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                         "APM28DIP.*.male.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)

male_files <- male_files[male_files %notin% female_files]


male_cluster_count_mean <- lapply(male_files, cluster_count_mean) 
female_cluster_count_mean <- lapply(female_files, cluster_count_mean)
male_cluster_count_mean <- as.data.frame(do.call(rbind, male_cluster_count_mean))
female_cluster_count_mean <- as.data.frame(do.call(rbind, female_cluster_count_mean))
cluster_count_mean_MOR28 <- rbind(male_cluster_count_mean, female_cluster_count_mean)
cluster_count_mean_MOR28$geno <- "MOR28"

male_cluster_count_1sd <- lapply(male_files, cluster_count_1sd) 
female_cluster_count_1sd <- lapply(female_files, cluster_count_1sd)
male_cluster_count_1sd <- as.data.frame(do.call(rbind, male_cluster_count_1sd))
female_cluster_count_1sd <- as.data.frame(do.call(rbind, female_cluster_count_1sd))
cluster_count_1sd_MOR28 <- rbind(male_cluster_count_1sd, female_cluster_count_1sd)
cluster_count_1sd_MOR28$geno <- "MOR28"

male_cluster_count_2sd <- lapply(male_files, cluster_count_2sd) 
female_cluster_count_2sd <- lapply(female_files, cluster_count_2sd)
male_cluster_count_2sd <- as.data.frame(do.call(rbind, male_cluster_count_2sd))
female_cluster_count_2sd <- as.data.frame(do.call(rbind, female_cluster_count_2sd))
cluster_count_2sd_MOR28 <- rbind(male_cluster_count_2sd, female_cluster_count_2sd)
cluster_count_2sd_MOR28$geno <- "MOR28"

##########################
# Generating Final Plots #
##########################

cluster_count_mean_df <- rbind(cluster_count_mean_MOR28, cluster_count_mean_p2) %>%
  mutate(type = "> mean size of active GIH")
cluster_count_1sd_df <- rbind(cluster_count_1sd_MOR28, cluster_count_1sd_p2) %>%
  mutate(type = "mean ± 1 sd size of active GIH")
cluster_count_2sd_df <- rbind(cluster_count_2sd_MOR28, cluster_count_2sd_p2) %>%
  mutate(type = "mean ± 2 sd size of active GIH")

cluster_count_df <- rbind(cluster_count_mean_df, cluster_count_1sd_df)

cluster_count_df %>% group_by(type) %>%
  summarise(mean = mean(V1),
            sd = sd(V1))

ggplot(cluster_count_mean_df, aes(x = type, y = V1)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylim(0, 25) +
  ylab("enhancer clusters per cell") +
  theme_classic()

ggplot(cluster_count_1sd_df, aes(x = type, y = V1)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylim(0, 25) +
  ylab("enhancer clusters per cell") +
  theme_classic()

