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
  filter(state == 'active')
  
ggplot(islands_percell_dataframe, aes(y=sum)) +
  geom_boxplot() + 
  #geom_area() + 
  theme_classic()


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
    summarise(total = n()) 
  return(valid_clusters)
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

##########################
# Generating Final Plots #
##########################

cluster_count_mean_df <- rbind(cluster_count_mean_MOR28, cluster_count_mean_p2) %>%
  mutate(hub = 'all_hubs') %>%
  select(total,hub)

islands_percell_dataframe <- islands_percell_dataframe %>%
  as.data.frame() %>%
  select(total, hub)
  
df <- rbind(cluster_count_mean_df, islands_percell_dataframe)
ggplot(df, aes(total, fill = hub)) + 
  geom_density(adjust = 2, alpha = 0.4) +
  scale_fill_manual(values = c('forestgreen', 'black')) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$total)) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$total) -sd(islands_percell_dataframe$total), linetype = 'dashed' ) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$total) +sd(islands_percell_dataframe$total), linetype = 'dashed' ) + 
  theme_classic()
 
##########################
# Generating Final Plots #
##########################

islands_percell_dataframe %>%
  as.data.frame() %>%
  select(total, hub) %>%
  group_by(hub) %>%
  summarise(mean = mean(total),
            sd = sd(total),
            count = n())


rbind(cluster_count_mean_MOR28, cluster_count_mean_p2) %>%
  mutate(hub = 'all_hubs') %>%
  select(total,hub) %>%
  group_by(hub) %>%
  summarise(mean = mean(total),
            sd = sd(total),
            count = n())

rbind(cluster_count_mean_MOR28, cluster_count_mean_p2) %>%
  mutate(hub = 'all_hubs') %>%
  select(total,hub) %>%
  dplyr::filter(total < (5.59 + 2.88), total > (5.59 - 2.88)) %>%
  group_by(hub) %>%
  summarise(mean = mean(total),
            sd = sd(total),
            count = n())

