library(pheatmap)
library(dendextend)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

`%notin%` <- Negate(`%in%`)
set.seed(99)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.inactive.hub.extraction.tan/") 

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
inactiveORs_percell_dataframe = data.frame()
for (i in c(1:length(adult_files))){
  # Extact the data
  df <- read.table(adult_files[i],col.names = names$V1)
  name <- adult_files[i]
  name <- sub(".20k.*", "", name)
  name <- sub(".*pd/MOE_adult.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame()
  df_sc <- as.data.frame(scale(df1)) #scale data
  dist_mat <- dist(df_sc, method = 'euclidean') #distance matrix
  hclust_avg <- hclust(dist_mat, method = 'average') #clustering to make dendrogram
  plot(hclust_avg)
  
  #cut tree
  cut_avg_df <- cutree(hclust_avg, h = 2.5) %>% #change to 2.5 to look at 5 p.r. diameter hubs
    as.data.frame() %>%
    dplyr::rename("cluster" = ".")
  cut_avg_df_dup <- cut_avg_df
  
  #find valid clusters, extract number of enhancers in active hub, extract a similarly sized inactive hub
  clusters <- cut_avg_df %>% 
    group_by(cluster) %>% 
    summarise(total = n()) %>%
    filter(total > 2)
  hubs_islands <- cut_avg_df %>%
    dplyr::filter(cluster %in% clusters$cluster) %>%
    mutate(cell = name) %>% as.data.frame()
  hubs_islands$island = row.names(hubs_islands)
  islands_percell_dataframe <- rbind(islands_percell_dataframe, hubs_islands)
  
  for (i in unique(hubs_islands$cluster)){
    print(i)
    inactive_islands <- hubs_islands %>%
      dplyr::filter(cluster == i)
    inactiveORs_in_hub <- df %>% 
      dplyr::select(gsub(",", ".", inactive_islands$island)) %>% t() %>%
      as.data.frame()
    inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t()
    inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
    if (nrow(inactive_islands) == 1){
      inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
      inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>%
        rownames() %>% as.data.frame() %>%
        dplyr::rename("element" = ".") %>%
        mutate("cell" = name) %>%
        mutate(cluster = i)
      inactiveORs_percell_dataframe <- rbind(inactiveORs_percell_dataframe, inactiveORs_in_hub)
    } else {
      inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
      inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
        rownames() %>% as.data.frame() %>%
        dplyr::rename("element" = ".") %>%
        mutate("cell" = name) %>%
        mutate(cluster = i)
      inactiveORs_percell_dataframe <- rbind(inactiveORs_percell_dataframe, inactiveORs_in_hub)
    }
  }
}
  
inactiveORs_percell_dataframe_adult <- inactiveORs_percell_dataframe
islands_percell_dataframe_adult <- islands_percell_dataframe

###############
# MOE NEWBORN #
###############

islands_percell_dataframe = data.frame()
inactiveORs_percell_dataframe = data.frame()
for (i in c(1:length(newborn_files))){
  # Extact the data
  df <- read.table(newborn_files[i],col.names = names$V1)
  name <- newborn_files[i]
  name <- sub(".20k.*", "", name)
  name <- sub(".*pd/MOE_newborn.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame()
  df_sc <- as.data.frame(scale(df1)) #scale data
  dist_mat <- dist(df_sc, method = 'euclidean') #distance matrix
  hclust_avg <- hclust(dist_mat, method = 'average') #clustering to make dendrogram
  plot(hclust_avg)
  
  #cut tree
  cut_avg_df <- cutree(hclust_avg, h = 2.5) %>% #change to 2.5 to look at 5 p.r. diameter hubs
    as.data.frame() %>%
    dplyr::rename("cluster" = ".")
  cut_avg_df_dup <- cut_avg_df
  
  #find valid clusters, extract number of enhancers in active hub, extract a similarly sized inactive hub
  clusters <- cut_avg_df %>% 
    group_by(cluster) %>% 
    summarise(total = n()) %>%
    filter(total > 2)
  hubs_islands <- cut_avg_df %>%
    dplyr::filter(cluster %in% clusters$cluster) %>%
    mutate(cell = name) %>% as.data.frame()
  hubs_islands$island = row.names(hubs_islands)
  islands_percell_dataframe <- rbind(islands_percell_dataframe, hubs_islands)
  
  for (i in unique(hubs_islands$cluster)){
    print(i)
    inactive_islands <- hubs_islands %>%
      dplyr::filter(cluster == i)
    inactiveORs_in_hub <- df %>% 
      dplyr::select(gsub(",", ".", inactive_islands$island)) %>% t() %>%
      as.data.frame()
    inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t()
    inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
    if (nrow(inactive_islands) == 1){
      inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
      inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>%
        rownames() %>% as.data.frame() %>%
        dplyr::rename("element" = ".") %>%
        mutate("cell" = name) %>%
        mutate(cluster = i)
      inactiveORs_percell_dataframe <- rbind(inactiveORs_percell_dataframe, inactiveORs_in_hub)
    } else {
      inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
      inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
        rownames() %>% as.data.frame() %>%
        dplyr::rename("element" = ".") %>%
        mutate("cell" = name) %>%
        mutate(cluster = i)
      inactiveORs_percell_dataframe <- rbind(inactiveORs_percell_dataframe, inactiveORs_in_hub)
    }
  }
}

inactiveORs_percell_dataframe_newborn <- inactiveORs_percell_dataframe
islands_percell_dataframe_newborn <- islands_percell_dataframe


inactiveORs_percell_dataframe <- rbind(inactiveORs_percell_dataframe_adult, inactiveORs_percell_dataframe_newborn)
islands_percell_dataframe <- rbind(islands_percell_dataframe_adult, islands_percell_dataframe_newborn)


# write.table(islands_percell_dataframe,
#             "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.inactive.hub.extraction.tan/islands_percell_dataframe.txt",
#             sep = "\t",
#             col.names = TRUE)
# 
# write.table(inactiveORs_percell_dataframe,
#             "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVII.inactive.hub.extraction.tan/inactiveORs_percell_dataframe.txt",
#             sep = "\t",
#             col.names = TRUE)


