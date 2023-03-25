library(pheatmap)
library(dendextend)
library(ggplot2)
library(dplyr)

`%notin%` <- Negate(`%in%`)
set.seed(99)
# P2 Analysis

# P2 Analysis
names <- read.table("/media/storageE/ariel_dipc/annotations/ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active'


male_files <- list.files(path = "/media/storageE/ariel_dipc/P2-ImmFix-dipc/matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "/media/storageE/ariel_dipc/P2-ImmFix-dipc/matrix.P2-ImmFix-dipc.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)


activeORs_percell_datalist = list()
for (i in c(1)){
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
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active") %>%
    mutate("cell" = name)
  print(name)
  
  #extract inactive ORs participating in the hub
  activeORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", active_islands$element)) %>% t() %>%
    as.data.frame()
  activeORs_in_hub <- select(activeORs_in_hub, c(contains("Olfr"), active)) %>% t() 
  activeORs_in_hub <- activeORs_in_hub %>% as.data.frame()
  if (nrow(active_islands) == 0){
  }
  if (nrow(active_islands) == 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[activeORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub
  } 
  if (nrow(active_islands) > 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[!rowSums(activeORs_in_hub[,-(ncol(activeORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub}
}

activeORs_percell_datalist_male_P2 = do.call(rbind, activeORs_percell_datalist)

activeORs_percell_datalist = list()
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
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active") %>%
    mutate("cell" = name)
  print(name)
  
  #extract inactive ORs participating in the hub
  activeORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", active_islands$element)) %>% t() %>%
    as.data.frame()
  activeORs_in_hub <- select(activeORs_in_hub, c(contains("Olfr"), active)) %>% t() 
  activeORs_in_hub <- activeORs_in_hub %>% as.data.frame()
  if (nrow(active_islands) == 0){
  }
  if (nrow(active_islands) == 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[activeORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub
  } 
  if (nrow(active_islands) > 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[!rowSums(activeORs_in_hub[,-(ncol(activeORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub}
}

activeORs_percell_datalist_female_P2 = do.call(rbind, activeORs_percell_datalist)

activeORs_percell_datalist_P2 <- rbind(activeORs_percell_datalist_male_P2, activeORs_percell_datalist_female_P2)

write.table(activeORs_percell_datalist,
            "/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/activeGIH_ORs_percell_datalist_P2.txt",
            sep = "\t",
            col.names = TRUE)




#### Mor28 Analysis


names <- read.table("/media/storageE/ariel_dipc/annotations/ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr14,52491331,Olfr1507,0",] <- 'active'


male_files <- list.files(path = "/media/storageE/AP_LB_dipc/matrix.APM28DIP.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "/media/storageE/AP_LB_dipc/matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

activeORs_percell_datalist = list()
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
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active") %>%
    mutate("cell" = name)
  print(name)
  
  #extract inactive ORs participating in the hub
  activeORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", active_islands$element)) %>% t() %>%
    as.data.frame()
  activeORs_in_hub <- select(activeORs_in_hub, c(contains("Olfr"), active)) %>% t() 
  activeORs_in_hub <- activeORs_in_hub %>% as.data.frame()
  if (nrow(active_islands) == 0){
  }
  if (nrow(active_islands) == 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[activeORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub
  } 
  if (nrow(active_islands) > 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[!rowSums(activeORs_in_hub[,-(ncol(activeORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub}
}

activeORs_percell_datalist_male_MOR28 = do.call(rbind, activeORs_percell_datalist)

activeORs_percell_datalist = list()
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
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  active_islands <- cut_avg_df %>%
    dplyr::filter(cluster == 1) %>%
    rownames() %>% as.data.frame() %>%
    dplyr::rename("element" = ".") %>%
    mutate("state" = "active") %>%
    dplyr::filter(element != "active") %>%
    mutate("cell" = name)
  print(name)
  
  #extract inactive ORs participating in the hub
  activeORs_in_hub <- df %>% 
    dplyr::select(gsub(",", ".", active_islands$element)) %>% t() %>%
    as.data.frame()
  activeORs_in_hub <- select(activeORs_in_hub, c(contains("Olfr"), active)) %>% t() 
  activeORs_in_hub <- activeORs_in_hub %>% as.data.frame()
  if (nrow(active_islands) == 0){
  }
  if (nrow(active_islands) == 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[activeORs_in_hub[,1] < 5,] %>% 
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub
  } 
  if (nrow(active_islands) > 1){
    activeORs_in_hub$olfr <- rownames(activeORs_in_hub)
    activeORs_in_hub <- activeORs_in_hub[!rowSums(activeORs_in_hub[,-(ncol(activeORs_in_hub))] > 5),] %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename("element" = ".") %>%
      mutate("cell" = name)
    activeORs_percell_datalist[[i]] <- activeORs_in_hub}
}

activeORs_percell_datalist_female_MOR28 = do.call(rbind, activeORs_percell_datalist)

activeORs_percell_datalist_MOR28 <- rbind(activeORs_percell_datalist_male_MOR28, activeORs_percell_datalist_female_MOR28)

write.table(activeORs_percell_datalist_MOR28,
            "/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/activeGIH_ORs_percell_datalist_MOR28.txt",
            sep = "\t",
            col.names = TRUE)
