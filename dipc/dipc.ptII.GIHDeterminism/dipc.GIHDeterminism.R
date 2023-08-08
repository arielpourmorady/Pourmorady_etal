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

setwd("/media/storageE/ariel/R/finalpaper_finalized/dipc/dipc.enhancer_preferences") 

##################################################################
######################### MOR28 ANALYSIS #########################
##################################################################

names <- read.table("/media/storageE/ariel_dipc/annotations/ORs.and.GIs.nonX.names")
names$V1 <- as.character(names$V1)
names[names$V1 == "chr14,52491331,Olfr1507,0",] <- 'active' # Location of MOR28 allele


male_files <- list.files(path = "matrix.APM28DIP.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

distance_from_active_allele_datalist = list()
for (i in c(1:length(male_files))){
  # Extact the data
  df <- read.table(male_files[i],col.names = names$V1)
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% 
    t() %>% as.data.frame() %>%
    select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_male_MOR28 = do.call(rbind, distance_from_active_allele_datalist)

distance_from_active_allele_datalist = list()
for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*APM28DIP.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% 
    t() %>% as.data.frame() %>%
    select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_female_MOR28 = do.call(rbind, distance_from_active_allele_datalist)
distance_from_active_allele_dataframe_MOR28 <- rbind(distance_from_active_allele_dataframe_male_MOR28, 
                                                     distance_from_active_allele_dataframe_female_MOR28) %>%
  mutate("OR" = "MOR28")

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
  df <- read.table(male_files[i],col.names = names$V1)
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% 
    t() %>% as.data.frame() %>%
    select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_male_P2 = do.call(rbind, distance_from_active_allele_datalist)

distance_from_active_allele_datalist = list()
for (i in c(1:length(female_files))){
  # Extact the data
  df <- read.table(female_files[i],col.names = names$V1)
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  rownames(df) <- names$V1
  df1 <- select(df, !contains("Olfr"))
  df1 <- df1 %>% t() %>% as.data.frame()
  df1 <- select(df1, !contains("Olfr")) %>% 
    t() %>% as.data.frame() %>%
    select(active)
  df1 <- df1 %>% mutate("island" = rownames(df1)) %>%
    mutate("cell" = name)
  distance_from_active_allele_datalist[[i]] <- df1
}

distance_from_active_allele_dataframe_female_P2 = do.call(rbind, distance_from_active_allele_datalist)
distance_from_active_allele_dataframe_P2 <- rbind(distance_from_active_allele_dataframe_male_P2, 
                                                  distance_from_active_allele_dataframe_female_P2) %>%
  mutate("OR" = "P2")


##################################################################
######################### DETERMINISM P2 #########################
##################################################################

distance_from_active_allele_dataframe_P2_filter <- distance_from_active_allele_dataframe_P2 %>%
  filter(island != 'active', active < 2.5) %>% 
  filter(!grepl('chr7', island)) %>%
  filter(!grepl('chr14', island))

P2_df2 <- split(distance_from_active_allele_dataframe_P2_filter$island, distance_from_active_allele_dataframe_P2_filter$cell)
P2_df3 <- P2_df2[lapply(P2_df2,length)>0]

start <- Sys.time()
P2_df_clean <- data.frame()
for (i in c(1:length(P2_df3))){
  #print(i)
  name_i <- names(P2_df3[i])
  for (j in c(1:length(P2_df3))){
    name_j <- names(P2_df3[j])
   # print(j)
    if (i != j){
      a <- intersect(P2_df3[[i]], P2_df3[[j]]) %>% length()
      P2_df_iterate <- data.frame(name_i, name_j, a)
      P2_df_clean <- rbind(P2_df_clean, P2_df_iterate)
    }
  }
}
print( Sys.time() - start)

P2_grouped_GI_freq <- P2_df_clean %>%
  group_by(a) %>%
  summarise(count = n()) %>%
  mutate(count = count/sum(count)*100)

a_p2 <- ggplot(P2_grouped_GI_freq, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line() +
  ylim(0,100) + 
  ylab('percent of cell pairs with GIs in common') + 
  xlab('number of GIs in common') + 
  theme_classic()
a_p2

#####################################################################
######################### DETERMINISM MOR28 #########################
#####################################################################

distance_from_active_allele_dataframe_MOR28_filter <- distance_from_active_allele_dataframe_MOR28 %>%
  filter(island != 'active', active < 2.5) %>% 
  filter(!grepl('chr7', island)) %>%
  filter(!grepl('chr14', island))

MOR28_df2 <- split(distance_from_active_allele_dataframe_MOR28_filter$island, distance_from_active_allele_dataframe_MOR28_filter$cell)
MOR28_df3 <- MOR28_df2[lapply(MOR28_df2,length)>0]

start <- Sys.time()
MOR28_df_clean <- data.frame()
for (i in c(1:length(MOR28_df3))){
  #print(i)
  name_i <- names(MOR28_df3[i])
  for (j in c(1:length(MOR28_df3))){
    name_j <- names(MOR28_df3[j])
    # print(j)
    if (i != j){
      a <- intersect(MOR28_df3[[i]], MOR28_df3[[j]]) %>% length()
      MOR28_df_iterate <- data.frame(name_i, name_j, a)
      MOR28_df_clean <- rbind(MOR28_df_clean, MOR28_df_iterate)
    }
  }
}
print( Sys.time() - start)

MOR28_grouped_GI_freq <- MOR28_df_clean %>%
  group_by(a) %>%
  summarise(count = n()) %>%
  mutate(count = count/sum(count)*100)

a_MOR28 <- ggplot(MOR28_grouped_GI_freq, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line() +
  ylim(0,100) + 
  ylab('percent of cell pairs with GIs in common') + 
  xlab('number of GIs in common') + 
  theme_classic()
a_MOR28

##########################################################################
######################### DETERMINISM MOR28 & P2 #########################
##########################################################################

start <- Sys.time()
combined_df_clean <- data.frame()
for (i in c(1:length(MOR28_df3))){
  name_i <- names(MOR28_df3[i])
  for (j in c(1:length(P2_df3))){
    name_j <- names(P2_df3[j])
    a <- intersect(MOR28_df3[[i]], P2_df3[[j]]) %>% length()
    combined_df_iterate <- data.frame(name_i, name_j, a)
    combined_df_clean <- rbind(combined_df_clean, combined_df_iterate)
  }
}
print( Sys.time() - start)

combined_grouped_GI_freq <- combined_df_clean %>%
  group_by(a) %>%
  summarise(count = n()) %>%
  mutate(count = count/sum(count)*100)

a_combined <- ggplot(combined_grouped_GI_freq, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line() +
  ylim(0,100) + 
  ylab('percent of cell pairs with GIs in common') + 
  xlab('number of GIs in common') + 
  theme_classic()
a_combined

P2_grouped_GI_freq$comp <- 'P2vP2'
MOR28_grouped_GI_freq$comp <- 'mor28vmor28'
combined_grouped_GI_freq$comp <- 'P2vmor28'
final_data <- rbind(P2_grouped_GI_freq, MOR28_grouped_GI_freq, combined_grouped_GI_freq)

ggplot(final_data, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line(aes(linetype = comp)) +
  ylim(0,100) + 
  ylab('percent of cell pairs with GIs in common') + 
  xlab('number of GIs in common') + 
  theme_classic()
