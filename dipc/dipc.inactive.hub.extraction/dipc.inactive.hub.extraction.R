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

setwd("/media/storageE/ariel/R/finalpaper_finalized/dipc/dipc.inactive.hub.extraction") 

###################################################################################
###################################################################################
###################################################################################
############# Extracting Active GIH, Inacitve GIH, and ORs in GIHs ################
###################################################################################
###################################################################################
###################################################################################

###################################################################################
# The code below is used to extract the GIs in the active GIH, the GIs in an
# inactive hub of the same size as the active GIH in each cell, and all of the OR
# genes within each of these hubs
###################################################################################

###########################
# gg8-tTA > tetO-P2 Dip-C #
###########################

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active'


male_files <- list.files(path = "matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.P2-ImmFix-dipc.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

islands_percell_datalist = list()
inactiveORs_percell_datalist = list()
clusters_datalist = list()

for (i in length(male_files)){
  # Extact the data
  df <- read.table(male_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  islands_percell$cell <- name
  islands_percell_datalist[[i]] <- islands_percell
  clusters$cell <- name
  clusters_datalist[[i]] <- clusters
  
  
  #extract inactive ORs participating in the hub
  inactiveORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", inactive_islands$element)) %>% t() %>%
    as.data.frame()
  inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t() 
  inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
  if (nrow(inactive_islands) == 1){
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub
  } else {
  inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
  inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("cell" = name)
  inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub}
  
}

inactiveORs_percell_dataframe_male_P2 = do.call(rbind, inactiveORs_percell_datalist)
islands_percell_dataframe_male_P2 = do.call(rbind, islands_percell_datalist)
clusters_dataframe_male_P2 = do.call(rbind, clusters_datalist)


islands_percell_datalist = list()
inactiveORs_percell_datalist = list()
clusters_datalist = list()

for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  islands_percell$cell <- name
  islands_percell_datalist[[i]] <- islands_percell
  clusters$cell <- name
  clusters_datalist[[i]] <- clusters
  
  #extract inactive ORs participating in the hub
  inactiveORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", inactive_islands$element)) %>% t() %>%
    as.data.frame()
  inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t() 
  inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
  if (nrow(inactive_islands) == 1){
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub
  } else {
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub}
  
}

inactiveORs_percell_dataframe_female_P2 = do.call(rbind, inactiveORs_percell_datalist)
islands_percell_dataframe_female_P2 = do.call(rbind, islands_percell_datalist)
clusters_dataframe_female_P2 = do.call(rbind, clusters_datalist)


islands_percell_dataframe_P2 <- rbind(islands_percell_dataframe_male_P2, islands_percell_dataframe_female_P2)
inactiveORs_percell_dataframe_P2 <- rbind(inactiveORs_percell_dataframe_male_P2, inactiveORs_percell_dataframe_female_P2)
clusters_dataframe_P2 <- rbind(clusters_dataframe_male_P2, clusters_dataframe_female_P2)

#write.table(islands_percell_dataframe_P2,
#            "/media/storageE/ariel/R/finalpaper_finalized/dipc.inactive.hub.extraction/islands_percell_dataframe_P2.txt",
#            sep = "\t",
#            col.names = TRUE)

#write.table(inactiveORs_percell_dataframe_P2,
#            "/media/storageE/ariel/R/finalpaper_finalized/dipc.inactive.hub.extraction/inactiveORs_percell_dataframe_P2.txt",
#            sep = "\t",
#            col.names = TRUE)


######################
# Heatmap Generation #
######################

set.seed(99)
names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active'

heatmap <- function(paired_distance_path){
  df <- read.table(paired_distance_path,col.names = names$V1)
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
  clusters <- cut_avg_df %>% 
    group_by(cluster) %>% 
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  cut_avg_df$island <- rownames(cut_avg_df)
  cut_avg_df_active <- cut_avg_df %>% filter(cluster == 1)
  cut_avg_df_inactive <- cut_avg_df %>% filter(cluster == random_inactive_hub)
  cut_avg_df_non_annotated <- cut_avg_df %>% filter(cluster != 1, cluster != random_inactive_hub) %>%
    mutate("cluster" = 0)
  cut_avg_df <- rbind(cut_avg_df_active, cut_avg_df_inactive, cut_avg_df_non_annotated) 
  cut_avg_df <- cut_avg_df %>%
    dplyr::slice(match(names$V1, island))
  annotation <- data.frame(as.factor(cut_avg_df$cluster))
  rownames(annotation) <- colnames(df1)
  paletteLength <- 100
  myColor <- colorRampPalette(c("firebrick3", "white"))(paletteLength)
  myBreaks <- c(seq(0, 2.5, length.out=ceiling(paletteLength/2)), 
                seq(2.55, 5, length.out=floor(paletteLength/2)))
  a <- pheatmap(df1, 
                cluster_rows=TRUE, 
                cluster_cols=TRUE, 
                border_color = NA,
                show_rownames = FALSE,
                show_colnames = FALSE,
                col = myColor,
                cellheight=2, cellwidth = 2, breaks = myBreaks,
                annotation_col = annotation)
  return(a)
}


# Figure 2E generation
a <- heatmap("matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd/P2-ImmFix-dipc.10.male.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt")




###################################################################################
###################################################################################
###################################################################################
############# Extracting Active GIH, Inacitve GIH, and ORs in GIHs ################
###################################################################################
###################################################################################
###################################################################################

###################################################################################
# The code below is used to extract the GIs in the active GIH, the GIs in an
# inactive hub of the same size as the active GIH in each cell, and all of the OR
# genes within each of these hubs
###################################################################################

########################
# mor28-IRES-GFP Dip-C #
########################

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr14,52491331,Olfr1507,0",] <- 'active'


male_files <- list.files(path = "matrix.APM28DIP.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

islands_percell_datalist = list()
inactiveORs_percell_datalist = list()
clusters_datalist = list()
for (i in c(1:length(male_files))){
  # Extact the data
  df <- read.table(male_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/APM28DIP.","",name)
  print(name)
  islands_percell$cell <- name
  islands_percell_datalist[[i]] <- islands_percell
  clusters$cell <- name
  clusters_datalist[[i]] <- clusters
  
  
  #extract inactive ORs participating in the hub
  inactiveORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", inactive_islands$element)) %>% t() %>%
    as.data.frame()
  inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t() 
  inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
  if (nrow(inactive_islands) == 1){
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub
  } else {
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub}
}

inactiveORs_percell_dataframe_male_MOR28 = do.call(rbind, inactiveORs_percell_datalist)
islands_percell_dataframe_male_MOR28 = do.call(rbind, islands_percell_datalist)
clusters_dataframe_male_MOR28 = do.call(rbind, clusters_datalist)


islands_percell_datalist = list()
inactiveORs_percell_datalist = list()
clusters_datalist = list()

for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/APM28DIP.","",name)
  print(name)
  islands_percell$cell <- name
  islands_percell_datalist[[i]] <- islands_percell
  clusters$cell <- name
  clusters_datalist[[i]] <- clusters
  
  #extract inactive ORs participating in the hub
  inactiveORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", inactive_islands$element)) %>% t() %>%
    as.data.frame()
  inactiveORs_in_hub <- select(inactiveORs_in_hub, contains("Olfr")) %>% t() 
  inactiveORs_in_hub <- inactiveORs_in_hub %>% as.data.frame()
  if (nrow(inactive_islands) == 1){
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[inactiveORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub
  } else {
    inactiveORs_in_hub$olfr <- rownames(inactiveORs_in_hub)
    inactiveORs_in_hub <- inactiveORs_in_hub[!rowSums(inactiveORs_in_hub[,-(ncol(inactiveORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    inactiveORs_percell_datalist[[i]] <- inactiveORs_in_hub}
  
}

inactiveORs_percell_dataframe_female_MOR28 = do.call(rbind, inactiveORs_percell_datalist)
islands_percell_dataframe_female_MOR28 = do.call(rbind, islands_percell_datalist)
clusters_dataframe_female_MOR28 = do.call(rbind, clusters_datalist)


islands_percell_dataframe_MOR28 <- rbind(islands_percell_dataframe_male_MOR28, islands_percell_dataframe_female_MOR28)
inactiveORs_percell_dataframe_MOR28 <- rbind(inactiveORs_percell_dataframe_male_MOR28, inactiveORs_percell_dataframe_female_MOR28)
clusters_dataframe_MOR28 <- rbind(clusters_dataframe_male_MOR28, clusters_dataframe_female_MOR28)

#write.table(islands_percell_dataframe_MOR28,
#            "/media/storageE/ariel/R/finalpaper_finalized/dipc.inactive.hub.extraction/islands_percell_dataframe_MOR28.txt",
#            sep = "\t",
#            col.names = TRUE)

#write.table(inactiveORs_percell_dataframe_MOR28,
#            "/media/storageE/ariel/R/finalpaper_finalized/dipc.inactive.hub.extraction/inactiveORs_percell_dataframe_MOR28.txt",
#            sep = "\t",
#            col.names = TRUE)


######################
# Heatmap Generation #
######################

set.seed(99)
names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr14,52491331,Olfr1507,0",] <- 'active'

heatmap <- function(paired_distance_path){
  df <- read.table(paired_distance_path,col.names = names$V1)
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
  clusters <- cut_avg_df %>% 
    group_by(cluster) %>% 
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  cut_avg_df$island <- rownames(cut_avg_df)
  cut_avg_df_active <- cut_avg_df %>% filter(cluster == 1)
  cut_avg_df_inactive <- cut_avg_df %>% filter(cluster == random_inactive_hub)
  cut_avg_df_non_annotated <- cut_avg_df %>% filter(cluster != 1, cluster != random_inactive_hub) %>%
    mutate("cluster" = 0)
  cut_avg_df <- rbind(cut_avg_df_active, cut_avg_df_inactive, cut_avg_df_non_annotated) 
  cut_avg_df <- cut_avg_df %>%
    dplyr::slice(match(names$V1, island))
  annotation <- data.frame(as.factor(cut_avg_df$cluster))
  rownames(annotation) <- colnames(df1)
  paletteLength <- 100
  myColor <- colorRampPalette(c("firebrick3", "white"))(paletteLength)
  myBreaks <- c(seq(0, 2.5, length.out=ceiling(paletteLength/2)), 
                seq(2.55, 5, length.out=floor(paletteLength/2)))
  a <- pheatmap(df1, 
                cluster_rows=TRUE, 
                cluster_cols=TRUE, 
                border_color = NA,
                show_rownames = FALSE,
                show_colnames = FALSE,
                col = myColor,
                cellheight=2, cellwidth = 2, breaks = myBreaks,
                annotation_col = annotation)
  return(a)
}


heatmap("matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd/APM28DIP.34.female.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt")


###################################################################################
# The code below is a proof of principle that I am indeed extracting enhancers hubs
# where enhancers participating in the same hub are all on average ~2.5 p.r. away
# from each other.
###################################################################################

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr14,52491331,Olfr1507,0",] <- 'active'


male_files <- list.files(path = "matrix.APM28DIP.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

df_GI_distances_datalist = list()

for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  print(name)
  df_melt <- melt(as.matrix(df1), value.name = "dist") %>%
    filter(Var1 != "active") %>%
    filter(Var2 != "active") %>% 
    filter(dist != 0) %>%
    as.data.frame() %>%
    mutate("Var2" = gsub(".",",",as.character(Var2), fixed = TRUE))
  df_melt_active <- df_melt %>%
    filter(Var1 %in% active_islands$element) %>%
    filter(Var2 %in% active_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "active")
  df_melt_inactive <- df_melt %>%
    filter(Var1 %in% inactive_islands$element) %>%
    filter(Var2 %in% inactive_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "inactive")
  df_GI_distances <- rbind(df_melt_active, df_melt_inactive)
  df_GI_distances_datalist[[i]] <- df_GI_distances
}

df_GI_distances_datalist_female_MOR28 = do.call(rbind, df_GI_distances_datalist)


df_GI_distances_datalist = list()

for (i in c(1:length(male_files))){
  # Extact the data
  df <- read.table(male_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  print(name)
  df_melt <- melt(as.matrix(df1), value.name = "dist") %>%
    filter(Var1 != "active") %>%
    filter(Var2 != "active") %>% 
    filter(dist != 0) %>%
    as.data.frame() %>%
    mutate("Var2" = gsub(".",",",as.character(Var2), fixed = TRUE))
  df_melt_active <- df_melt %>%
    filter(Var1 %in% active_islands$element) %>%
    filter(Var2 %in% active_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "active")
  df_melt_inactive <- df_melt %>%
    filter(Var1 %in% inactive_islands$element) %>%
    filter(Var2 %in% inactive_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "inactive")
  df_GI_distances <- rbind(df_melt_active, df_melt_inactive)
  df_GI_distances_datalist[[i]] <- df_GI_distances
}

df_GI_distances_datalist_male_MOR28 = do.call(rbind, df_GI_distances_datalist)

# P2

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active'


male_files <- list.files(path = "matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.P2-ImmFix-dipc.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

df_GI_distances_datalist = list()

for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  df_melt <- melt(as.matrix(df1), value.name = "dist") %>%
    filter(Var1 != "active") %>%
    filter(Var2 != "active") %>% 
    filter(dist != 0) %>%
    as.data.frame() %>%
    mutate("Var2" = gsub(".",",",as.character(Var2), fixed = TRUE))
  df_melt_active <- df_melt %>%
    filter(Var1 %in% active_islands$element) %>%
    filter(Var2 %in% active_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "active")
  df_melt_inactive <- df_melt %>%
    filter(Var1 %in% inactive_islands$element) %>%
    filter(Var2 %in% inactive_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "inactive")
  df_GI_distances <- rbind(df_melt_active, df_melt_inactive)
  df_GI_distances_datalist[[i]] <- df_GI_distances
}

df_GI_distances_datalist_female_P2 = do.call(rbind, df_GI_distances_datalist)


df_GI_distances_datalist = list()

for (i in c(1:length(male_files))){
  # Extact the data
  df <- read.table(male_files[i],col.names = names$V1)
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
    summarise(total = n()) 
  count_active_hub_enhancers <- clusters$total[1] - 1
  random_inactive_hub <- clusters %>% filter(cluster != 1)
  random_inactive_hub <- random_inactive_hub[sample(nrow(random_inactive_hub)),] # randomly reorder rows
  random_inactive_hub <- random_inactive_hub[which.min(abs(count_active_hub_enhancers-random_inactive_hub$total)),] %>% # this function will find the cluster that has the nearest number of enhancers as the active hub
    select(cluster) %>% as.numeric()
  inactive_islands <- cut_avg_df %>% 
    filter(cluster == random_inactive_hub) %>% 
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "inactive")
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active")
  islands_percell <- rbind(inactive_islands, active_islands)
  islands_percell
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  df_melt <- melt(as.matrix(df1), value.name = "dist") %>%
    filter(Var1 != "active") %>%
    filter(Var2 != "active") %>% 
    filter(dist != 0) %>%
    as.data.frame() %>%
    mutate("Var2" = gsub(".",",",as.character(Var2), fixed = TRUE))
  df_melt_active <- df_melt %>%
    filter(Var1 %in% active_islands$element) %>%
    filter(Var2 %in% active_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "active")
  df_melt_inactive <- df_melt %>%
    filter(Var1 %in% inactive_islands$element) %>%
    filter(Var2 %in% inactive_islands$element) %>%
    mutate(cell = name) %>%
    mutate(state = "inactive")
  df_GI_distances <- rbind(df_melt_active, df_melt_inactive)
  df_GI_distances_datalist[[i]] <- df_GI_distances
}

df_GI_distances_datalist_male_P2 = do.call(rbind, df_GI_distances_datalist)


df_GI_distances_datalist_aggregated <- rbind(df_GI_distances_datalist_female_MOR28,
                                             df_GI_distances_datalist_male_MOR28,
                                             df_GI_distances_datalist_female_P2,
                                             df_GI_distances_datalist_male_P2)

ggplot(df_GI_distances_datalist_aggregated, aes(x = state, y = dist)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylab("intra-hub GI distance") + 
  theme_classic()

ggplot(df_GI_distances_datalist_aggregated, aes(x = dist, fill = state)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("#00FF00", "black")) + 
  xlab("intra-hub GI distance") + 
  theme_classic()

df_GI_distances_datalist_aggregated %>%
  group_by(state) %>%
  summarise(
    mean = mean(dist),
    sd = sd(dist)
  )

###################################################################################
# The code below is a proof of principle that I am indeed extracting enhancers hubs
# that contain similar numbers of enhancer elements
###################################################################################

islands_percell_dataframe_P2$geno <- "P2"
islands_percell_dataframe_MOR28$geno <- "MOR28"

islands_percell_dataframe <- rbind(islands_percell_dataframe_P2, islands_percell_dataframe_MOR28) %>%
  group_by(geno, cell, state) %>%
  summarise(total = n())

ggplot(islands_percell_dataframe, aes(x = state, y = total)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylab("GIs per GIH") + 
  theme_classic()

ggplot(islands_percell_dataframe, aes(x = total, fill = state)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c("#00FF00", "black")) + 
  ylab("GIs per GIH") + 
  theme_classic()

islands_percell_dataframe %>%
  group_by(state) %>%
  summarise(
    mean = mean(total),
    sd = sd(total)
  )


