library(factoextra)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggalt)
library(ggforce)
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

`%notin%` <- Negate(`%in%`)
set.seed(99)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.tanPCAandGIHsize/") 


#############################################################
######################### FILE LIST #########################
#############################################################

adult_files <- list.files(path = "MOE_adult",
                          "*.cpg_b1m.color2.txt.gz",
                          full.names = TRUE,
                          recursive = FALSE)

newborn_files <- list.files(path = "MOE_newborn",
                            "*.cpg_b1m.color2.txt.gz",
                            full.names = TRUE,
                            recursive = FALSE)

lomvardas_files <- list.files(path = "Lomvardas_downsample",
                              "*.cpg_b1m.color2",
                              full.names = TRUE,
                              recursive = FALSE)


#############################################################
###################### CPG EXCTRACTION ######################
#############################################################

extract_cpg_fromfile <- function(files, dataset_name){
  df = data.frame()
  for (i in files){
    df1 <- read.table(i) %>%
      dplyr::rename("chr" = V1,
                    "loc" = V2,
                    "cpg" = V3) %>%
      mutate(cell = i,
             dataset = dataset_name)
    df <- rbind(df, df1)
  }
  return(df)
}

tan_adult_cpg <- extract_cpg_fromfile(adult_files, "tan_adult")
tan_newborn_cpg <- extract_cpg_fromfile(newborn_files, "tan_newborn")
lomvardas_cpg <- extract_cpg_fromfile(lomvardas_files, "lomvardas")

cpg <- rbind(tan_adult_cpg, tan_newborn_cpg,
             lomvardas_cpg) %>%
  dplyr::filter(chr != "chrX") %>%
  dplyr::filter(chr != "chrY") %>%
  mutate(loc = paste(chr, loc, sep = "_"),
         cpg = as.numeric(cpg)) %>%
  dplyr::select(cell, loc, cpg) %>% 
  dcast(cell ~ loc, value.var = "cpg")

cpg <- cpg %>% dplyr::select_if(~ !any(is.na(.)))


##########################################
########### Project onto PCA  ############
##########################################

cpg_tan <- rbind(tan_adult_cpg, tan_newborn_cpg) %>%
  dplyr::filter(chr != "chrX") %>%
  dplyr::filter(chr != "chrY") %>%
  mutate(loc = paste(chr, loc, sep = "_"),
         cpg = as.numeric(cpg)) %>%
  dplyr::select(cell, loc, cpg) %>% 
  dcast(cell ~ loc, value.var = "cpg")

cpg_tan <- cpg_tan %>% dplyr::select_if(~ !any(is.na(.)))


cpg_else <- rbind(lomvardas_cpg) %>%
  dplyr::filter(chr != "chrX") %>%
  dplyr::filter(chr != "chrY") %>%
  mutate(loc = paste(chr, loc, sep = "_"),
         cpg = as.numeric(cpg)) %>%
  dplyr::select(cell, loc, cpg) %>% 
  dcast(cell ~ loc, value.var = "cpg")

cpg_else <- cpg_else %>% dplyr::select_if(~ !any(is.na(.)))

common_colnames <- intersect(colnames(cpg_else), colnames(cpg_tan))
cpg_else <- cpg_else %>%dplyr::select(common_colnames)

#################################################
###################### PCA ######################
#################################################

res.pca_tan <- prcomp(cpg_tan[,-1], scale = TRUE)
res.ind <- get_pca_ind(res.pca_tan)

df_coordinates <- res.ind$coord %>% as.data.frame() %>%
  dplyr::select(Dim.1, Dim.2)

df_coordinates$cell <- cpg_tan$cell

df_coordinates <- df_coordinates %>%
  separate(cell, sep = "/", into=c("age", "cell"))

df_coordinates$cell <- sub(".cpg.*", "", df_coordinates$cell)
df_coordinates <- df_coordinates %>% separate(cell, sep = "cell_", into = c(NA, "cell"))

####
res.pca_else <- predict(res.pca_tan, newdata = cpg_else[,-1]) %>% 
  as.data.frame() %>% 
  dplyr::select(PC1, PC2) %>%
  dplyr::rename("Dim.1" = PC1, "Dim.2" = PC2)
res.pca_else$cell <- cpg_else[,1]

res.pca_else <- res.pca_else %>%
  separate(cell, sep = "/", into=c("age", "cell"))

res.pca_else_lomvardas <- res.pca_else %>%
  filter(., grepl("25M", cell)) %>%
  mutate(age = c("Atf5", "Atf5", "Krt5", "Krt5", "Mash1Ngn", "Omp", "Omp")) 
res.pca_else_lomvardas <- res.pca_else_lomvardas[-10,]

#################################################
################### PLOTTING ####################
#################################################

for_plot <- res.pca_else_lomvardas %>%
  group_by(age) %>%
  summarise(Dim.1 = mean(Dim.1),
            Dim.2 = mean(Dim.2))

a <- ggplot(df_coordinates, aes(x = -Dim.1, y = -Dim.2, color = age)) +
  geom_point() +
  scale_color_manual(values = c( 'black', 'firebrick2')) + 
  #scale_color_manual(values = c( 'blue3', 'chocolate4', 'darkorange2', 'black', 'firebrick2','darkmagenta')) + 
  #geom_point(data = for_plot, aes(x = -Dim.1, y = -Dim.2, color = age), size = 2) + 
  geom_vline(xintercept = -for_plot$Dim.1, linetype = 'dashed') + 
  theme_classic()
a

##########################################################################
################### HUB ANALYSIS OVER DIFFERENTIATION ####################
##########################################################################


islands_percell_dataframe <- read.table("/media/storageE/ariel/R/dipc.inactive.hub.extraction.tan/islands_percell_dataframe.txt",sep = "\t") 

islands_percell_dataframe_grouped <- islands_percell_dataframe %>%
  group_by(cell, cluster) %>%
  summarise(island_count = n()) 

summarized_df_grouped <- islands_percell_dataframe_grouped %>%
  group_by(cell) %>%
  summarise(cluster_count = n(),
            cluster_size = mean(island_count))

df_coordinates <- df_coordinates %>%
  mutate(cell = paste('cell_', cell, sep=""))

df_coordinates_grouped_sd <- df_coordinates %>%
  left_join(summarized_df_grouped, by = "cell") 

b <- ggplot(df_coordinates_grouped_sd, aes(x = -Dim.1, y = cluster_size, color =  cluster_count)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  xlab("Differentiation") + 
  ylab("Average # GIs per GIH") + 
  labs(color = "# GIHs per cell") + 
  geom_vline(xintercept = -for_plot$Dim.1, linetype = 'dashed') + 
  theme_classic()
b

a + b

