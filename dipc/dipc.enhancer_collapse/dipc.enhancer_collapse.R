library(pheatmap)
library(dendextend)

# set directory to current directory
# this directory should contain all files in the GitHub folder

working_directory <- "/media/storageE/ariel/R/finalpaper_finalized/dipc/dipc.enhancer_collapse/"
setwd(working_directory) 

`%notin%` <- Negate(`%in%`)

##############################################################################
############################## READING IN FILES ##############################
##############################################################################

# Greek Island DF Name
names <- read.table("Greek_Islands.nonX.names")
`%notin%` <- Negate(`%in%`)

pca_percell <- read.table("pca_percell.txt", header = TRUE, sep = "\t") %>%
  mutate("cell" = sub(".cpg_b1m.color2.txt.gz*", "", X0)) 
pca_percell <- pca_percell %>%separate(col = cell, sep = "_cell_", into = c("trash", "cell")) %>%
  dplyr::select(principal.component.1, principal.component.2, target, cell)


adult_files <- list.files(path = "matrix.MOE_adult.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                          "*.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                          full.names = TRUE,
                          recursive = FALSE)

newborn_files <- list.files(path = "matrix.MOE_newborn.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd",
                            "*.20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt",
                            full.names = TRUE,
                            recursive = FALSE)

##################################################################
######################### FUNCTIONS ##############################
##################################################################

heatmap <- function(paired_distance_path){
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
    filter(total > median(enhancers_perhub_df$V1))
  cut_avg_df_dup <- cut_avg_df
  cut_avg_df_dup$island <- rownames(cut_avg_df)
  cut_avg_df_dup_valid <- cut_avg_df_dup %>% filter(cut_avg_df_dup$cluster %in% valid_clusters$cluster)
  cut_avg_df_dup_invalid <- cut_avg_df_dup %>% filter(cut_avg_df_dup$cluster %notin% valid_clusters$cluster)
  cut_avg_df_dup_invalid$cluster <- 0
  cut_avg_df <- rbind(cut_avg_df_dup_valid, cut_avg_df_dup_invalid) 
  cut_avg_df <- cut_avg_df %>%
    dplyr::slice(match(names$V1, island))
  annotation <- data.frame(as.factor(cut_avg_df$cluster))
  rownames(annotation) <- colnames(df)
  paletteLength <- 100
  myColor <- colorRampPalette(c("firebrick3", "gold1", "white"))(paletteLength)
  myBreaks <- c(seq(0, 5, length.out=ceiling(paletteLength/2)), 
                seq(5.2, 15, length.out=floor(paletteLength/2)))
  a <- pheatmap(df, 
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

cluster_count_newborn <- function(paired_distance_path){
  cell_name <- paired_distance_path 
  cell_name <- sub(".*MOE_newborn.cell_", "", cell_name)  
  cell_name <- sub(".20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt.*", "", cell_name)
  print(cell_name)
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
    summarise(total = n()) #%>% 
  #filter(total > 3)
  a <- nrow(valid_clusters) %>% as.data.frame() %>%
    rename("clusters" = ".")
  score <- pca_percell %>% filter(cell == cell_name)
  a <- cbind(a, score)
  return(a)
}

cluster_count_adult <- function(paired_distance_path){
  cell_name <- paired_distance_path 
  cell_name <- sub(".*MOE_adult.cell_", "", cell_name)  
  cell_name <- sub(".20k.1.clean.Greek_Islands.nonX.L1.Greek_Islands.nonX.L2.pd.txt.*", "", cell_name)
  print(cell_name)
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
    summarise(total = n()) #%>% 
  #filter(total > 3)
  a <- nrow(valid_clusters) %>% as.data.frame() %>%
    rename("clusters" = ".")
  print(a)
  score <- pca_percell %>% filter(cell == cell_name)
  a <- cbind(a, score)
  return(a)
}

######################################################################
############################## ANALYSIS ##############################
######################################################################

newborn_cluster_count <- lapply(newborn_files, cluster_count_newborn) 
adult_cluster_count <- lapply(adult_files, cluster_count_adult)

newborn_cluster_count <- as.data.frame(do.call(rbind, newborn_cluster_count))
newborn_cluster_count$geno <- "newborn"
adult_cluster_count <- as.data.frame(do.call(rbind, adult_cluster_count))
adult_cluster_count$geno <- "adult"
cluster_count_tan <- rbind(newborn_cluster_count, adult_cluster_count)

cluster_count_tan$principal.component.1 %>% range()

aaa <- ggplot(cluster_count_tan, aes(x = principal.component.1, y = principal.component.2)) + geom_point(aes(colour = target)) + 
  scale_color_manual(values = c("black", "red")) + 
  scale_x_reverse() + 
  theme_classic()

cluster_count_tan_one <-cluster_count_tan %>% filter( principal.component.1 < -24.57934 + (45.16282 - -24.57934)/4)
cluster_count_tan_one$geno <- "one"
cluster_count_tan_two <-cluster_count_tan %>% filter(-24.57934 + 2*(45.16282 - -24.57934)/4 > principal.component.1, principal.component.1> -24.57934 + (45.16282 - -24.57934)/4)
cluster_count_tan_two$geno <- "two"
cluster_count_tan_three <-cluster_count_tan %>% filter(-24.57934 + 3*(45.16282 - -24.57934)/4 > principal.component.1, principal.component.1> -24.57934 + 2*(45.16282 - -24.57934)/4)
cluster_count_tan_three$geno <- "three"
cluster_count_tan_four <-cluster_count_tan %>% filter(-24.57934 + 3*(45.16282 - -24.57934)/4 < principal.component.1)
cluster_count_tan_four$geno <- "four"

cluster_count_tan <- rbind(cluster_count_tan_one, cluster_count_tan_two, cluster_count_tan_three, cluster_count_tan_four)
cluster_count_tan$geno <- factor(cluster_count_tan$geno, levels = c("four", "three", "two", "one"))


bbb <- ggplot(cluster_count_tan, aes(x = geno, y = clusters)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylab("enhancer clusters per cell") + 
  theme_classic()

bbb / aaa +  plot_layout(heights = c(1, 2))

# GIH Average parameters
cluster_count_tan %>% group_by(geno) %>%
  summarise(mean = mean(clusters),
            sd = sd(clusters),
            n = n())
