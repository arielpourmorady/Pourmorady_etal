library(pheatmap)
library(dendextend)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

`%notin%` <- Negate(`%in%`)
set.seed(99)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptII.enhancer.preferences/") 

###############################################################
######################### P2 ANALYSIS #########################
###############################################################

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active' # Location of P2 allele


male_files <- list.files(path = "matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.P2-ImmFix-dipc.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

distance_from_active_allele_datalist = list()
for (i in c(1:length(male_files))){
  # Extact the data
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  df <- fread(male_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame() %>% select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_male_P2 = do.call(rbind, distance_from_active_allele_datalist)

distance_from_active_allele_datalist = list()
for (i in c(1:length(female_files))){
  # Extact the data
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  df <- fread(female_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame() %>% select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_female_P2 = do.call(rbind, distance_from_active_allele_datalist)
distance_from_active_allele_dataframe_P2 <- rbind(distance_from_active_allele_dataframe_male_P2, 
                                                  distance_from_active_allele_dataframe_female_P2) %>%
  mutate("OR" = "P2")

##################################################################
####################### ZONAL OMP ANALYSIS #######################
##################################################################

names <- read.table("ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr7,107095985,Olfr17,0",] <- 'active' # Location of P2 allele


files <- list.files(path = "/media/storageE/ariel_dipc/zonal_omp/matrix.zonal_omp.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)

distance_from_active_allele_datalist = list()
for (i in c(1:length(files))){
  # Extact the data
  name <- files[i]
  name <- sub(".20k.*", "", name)
  name <- sub(".*z","",name)
  print(name)
  df <- fread(files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))  
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% t() %>% as.data.frame() %>% select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_zonal_OMP = do.call(rbind, distance_from_active_allele_datalist) %>%
  mutate("OR" = "P2")


###################################################
#################### TILE PLOT ####################
###################################################

distance_from_active_allele_dataframe_P2_mut <- distance_from_active_allele_dataframe_P2 %>%
  filter(island != 'active') %>%
  mutate(in_hub = case_when(active < 2.5 ~ 1, active > 2.5 ~ 0)) %>%
  separate(island, sep = ',', into = c('chr', 'loc', 'island', 'haplo')) %>%
  separate(chr, sep = 'r', into = c(NA, 'chr')) %>%
  mutate("chr" = as.numeric(chr)) %>%
  mutate("loc" = as.numeric(loc)) %>%
  mutate(geno = 'P2')

distance_from_active_allele_dataframe_zonal_OMP_mut <- distance_from_active_allele_dataframe_zonal_OMP %>%
  filter(island != 'active') %>%
  mutate(in_hub = case_when(active < 2.5 ~ 1, active > 2.5 ~ 0)) %>%
  separate(island, sep = ',', into = c('chr', 'loc', 'island', 'haplo')) %>%
  separate(chr, sep = 'r', into = c(NA, 'chr')) %>%
  mutate("chr" = as.numeric(chr)) %>%
  mutate("loc" = as.numeric(loc)) %>%
  mutate(geno = 'zonal_OMP')

distance_from_active_allele_dataframe_mut <- rbind(distance_from_active_allele_dataframe_P2_mut, 
                                                   distance_from_active_allele_dataframe_zonal_OMP_mut) 


distance_from_active_allele_dataframe_mut <- distance_from_active_allele_dataframe_mut[
  with(distance_from_active_allele_dataframe_mut, order(haplo, chr, loc)),
  ] %>%
  mutate(ordering = paste(haplo, chr, loc, sep = ","))

distance_from_active_allele_dataframe_mut$ordering <- factor(distance_from_active_allele_dataframe_mut$ordering, 
                                                             levels=unique(distance_from_active_allele_dataframe_mut$ordering), 
                                                             ordered=TRUE)

distance_from_active_allele_dataframe_mut <- distance_from_active_allele_dataframe_mut %>%
  mutate(geno_cell = paste(geno, cell, sep = "-"))

distance_from_active_allele_dataframe_mut$geno_cell <- factor(distance_from_active_allele_dataframe_mut$geno_cell,
                                                              levels =rev(unique(distance_from_active_allele_dataframe_mut$geno_cell)), ordered = TRUE )


distance_from_active_allele_dataframe_mut_P2 <- distance_from_active_allele_dataframe_mut %>% filter(geno == 'P2')
distance_from_active_allele_dataframe_mut_zonal_OMP <- distance_from_active_allele_dataframe_mut %>% filter(geno == 'zonal_OMP')
zonal_OMP_cells <- distance_from_active_allele_dataframe_mut_zonal_OMP$cell %>% unique() %>% sample(40)
distance_from_active_allele_dataframe_mut_zonal_OMP <- distance_from_active_allele_dataframe_mut_zonal_OMP %>%
  filter(cell %in% zonal_OMP_cells)

distance_from_active_allele_dataframe_mut <- rbind(distance_from_active_allele_dataframe_mut_P2, distance_from_active_allele_dataframe_mut_zonal_OMP)

a <- ggplot(distance_from_active_allele_dataframe_mut, 
            aes(x = ordering, y = geno_cell, fill = in_hub)) + 
  scale_fill_gradientn(colors = c("white", "firebrick3")) +
  geom_tile() + 
  theme_bw() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
a

b <- ggplot(distance_from_active_allele_dataframe_mut , aes(x = ordering, y = 0, fill = as.character(chr))) + 
  geom_tile() +
  theme_void()
b

a /b  + plot_layout(heights = c(10, 1), guides = "collect")


#### Supplementary File

supplementary_file <- distance_from_active_allele_dataframe_mut %>%
  dplyr::select(geno, cell, ordering, island, haplo, in_hub)

write.table(supplementary_file, file = '/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptII.enhancer.preferences/Fig1b_P2_GIsInHub.txt', 
            sep = '\t', col.names = TRUE)
