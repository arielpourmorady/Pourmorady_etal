library(pheatmap)
library(dendextend)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(data.table)

`%notin%` <- Negate(`%in%`)
set.seed(99)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVIII.inactive.hub.extraction.tan/") 

###################################################################################
###################################################################################
###################################################################################
############# Extracting all GIHs in the 2019 Tan et al. dataset ##################
###################################################################################
###################################################################################
###################################################################################

names <- read.table("ORs.and.GIs.nonX.names")

adult_files <- list.files(path = "matrix.MOE_adult.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
newborn_files <- list.files(path = "matrix.MOE_newborn.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)


#############
# MOE ADULT #
#############

islands_percell_dataframe = data.frame()

for (i in c(1:length(adult_files))){
  # Extact the data
  name <- adult_files[i]
  name <- sub(".20k.*", "", name)
  name <- sub(".*pd/MOE_adult.","",name)
  print(name)
  df <- fread(adult_files[i]) %>% as.data.frame()
  rownames(df) <- names$V1
  colnames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame()
  df_sc <- as.data.frame(scale(df1)) #scale data
  dist_mat <- dist(df_sc, method = 'euclidean') #distance matrix
  hclust_avg <- hclust(dist_mat, method = 'average')
  #cut tree
  cut_avg_df <- cutree(hclust_avg, h = 2.5) %>% #change to 2.5 to look at 5 p.r. diameter hubs
    as.data.frame() %>%
    dplyr::rename("cluster" = ".")
  cut_avg_df_dup <- cut_avg_df
  
  #find valid clusters, extract number of enhancers in active hub, extract a similarly sized inactive hub
  clusters <- cut_avg_df %>% # hubs that have at least 3 or more GIs
    group_by(cluster) %>% 
    summarise(total = n()) %>%
    filter(total > 2)
  hubs_islands <- cut_avg_df %>%
    dplyr::filter(cluster %in% clusters$cluster) %>%
    mutate(cell = name) %>% as.data.frame()
  hubs_islands$island = row.names(hubs_islands)
  islands_percell_dataframe <- rbind(islands_percell_dataframe, hubs_islands)
}
  
islands_percell_dataframe_adult <- islands_percell_dataframe

###############
# MOE NEWBORN #
###############

islands_percell_dataframe = data.frame()

for (i in c(1:length(newborn_files))){
  # Extact the data
  name <- newborn_files[i]
  name <- sub(".20k.*", "", name)
  name <- sub(".*pd/MOE_newborn.","",name)
  print(name)
  df <- fread(newborn_files[i]) %>% as.data.frame()
  rownames(df) <- names$V1
  colnames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame()
  df_sc <- as.data.frame(scale(df1)) #scale data
  dist_mat <- dist(df_sc, method = 'euclidean') #distance matrix
  hclust_avg <- hclust(dist_mat, method = 'average')
  #cut tree
  cut_avg_df <- cutree(hclust_avg, h = 2.5) %>% #change to 2.5 to look at 5 p.r. diameter hubs
    as.data.frame() %>%
    dplyr::rename("cluster" = ".")
  cut_avg_df_dup <- cut_avg_df
  
  #find valid clusters, extract number of enhancers in active hub, extract a similarly sized inactive hub
  clusters <- cut_avg_df %>% # hubs that have at least 3 or more GIs
    group_by(cluster) %>% 
    summarise(total = n()) %>%
    filter(total > 2)
  hubs_islands <- cut_avg_df %>%
    dplyr::filter(cluster %in% clusters$cluster) %>%
    mutate(cell = name) %>% as.data.frame()
  hubs_islands$island = row.names(hubs_islands)
  islands_percell_dataframe <- rbind(islands_percell_dataframe, hubs_islands)
}

islands_percell_dataframe_newborn <- islands_percell_dataframe

islands_percell_dataframe <- rbind(islands_percell_dataframe_adult, islands_percell_dataframe_newborn)


write.table(islands_percell_dataframe,
            "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVIII.inactive.hub.extraction.tan/islands_percell_dataframe.txt",
            sep = "\t",
            col.names = TRUE)
