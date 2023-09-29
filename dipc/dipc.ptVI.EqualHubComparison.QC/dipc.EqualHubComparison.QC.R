library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(patchwork)

`%notin%` <- Negate(`%in%`)
set.seed(99)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptVI.EqualHubComparison.QC") 

###################################################################################
###################################################################################
###################################################################################
          ############# Extracting Intra GIH QC Parameters ################
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

male_files <- list.files(path = "matrix.P2-ImmFix-dipc.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.P2-ImmFix-dipc.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

islands_percell_dataframe_P2_230927 <- read.table("islands_percell_dataframe_P2_230927.txt")

intraIslandDistance = data.frame()
for (i in c(1:length(male_files))){
  # Extract Islands
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  activeIslands <- islands_percell_dataframe_P2_230927 %>% dplyr::filter(cell == name, state == "active") %>% select(element)
  inactiveIslands <- islands_percell_dataframe_P2_230927 %>% dplyr::filter(cell == name, state == "inactive") %>% select(element)
  # Extact the data
  df <- fread(male_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  #Intra-active Island Distance 
  df1_active <- select(df, activeIslands$element)  
  df1_active <- df1_active %>% t() %>% as.data.frame()
  df2_active <- select(df1_active, activeIslands$element) %>% t() %>% as.data.frame()
  df2_active <- df2_active %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "active", cell = name)
  #Intra-inactive Island Distance 
  df1_inactive <- select(df, inactiveIslands$element)  
  df1_inactive <- df1_inactive %>% t() %>% as.data.frame()
  df2_inactive <- select(df1_inactive, inactiveIslands$element) %>% t() %>% as.data.frame()
  df2_inactive <- df2_inactive %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "inactive", cell = name)
  intraIslandDistance <- rbind(intraIslandDistance, df2_active, df2_inactive)
}

intraIslandDistance_P2_male <- intraIslandDistance


for (i in c(1:length(female_files))){
  # Extract Islands
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/P2-ImmFix-dipc.","",name)
  print(name)
  activeIslands <- islands_percell_dataframe_P2_230927 %>% dplyr::filter(cell == name, state == "active") %>% select(element)
  inactiveIslands <- islands_percell_dataframe_P2_230927 %>% dplyr::filter(cell == name, state == "inactive") %>% select(element)
  # Extact the data
  df <- fread(female_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  #Intra-active Island Distance 
  df1_active <- select(df, activeIslands$element)  
  df1_active <- df1_active %>% t() %>% as.data.frame()
  df2_active <- select(df1_active, activeIslands$element) %>% t() %>% as.data.frame()
  df2_active <- df2_active %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "active", cell = name)
  #Intra-inactive Island Distance 
  df1_inactive <- select(df, inactiveIslands$element)  
  df1_inactive <- df1_inactive %>% t() %>% as.data.frame()
  df2_inactive <- select(df1_inactive, inactiveIslands$element) %>% t() %>% as.data.frame()
  df2_inactive <- df2_inactive %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "inactive", cell = name)
  intraIslandDistance <- rbind(intraIslandDistance, df2_active, df2_inactive)
}

intraIslandDistance_P2_female <- intraIslandDistance



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

male_files <- list.files(path = "matrix.APM28DIP.male.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                         "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)
female_files <- list.files(path = "matrix.APM28DIP.female.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd",
                           "*.20k.1.clean.ORs.and.GIs.nonX.L1.ORs.and.GIs.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

islands_percell_dataframe_MOR28_230927 <- read.table("islands_percell_dataframe_MOR28_230927.txt")


intraIslandDistance = data.frame()
for (i in c(1:length(male_files))){
  # Extract Islands
  name <- male_files[i]
  name <- sub(".male.20k.*", "", name)
  name <- sub(".*pd/APM28DIP.","",name)
  print(name)
  activeIslands <- islands_percell_dataframe_MOR28_230927 %>% dplyr::filter(cell == name, state == "active") %>% select(element)
  inactiveIslands <- islands_percell_dataframe_MOR28_230927 %>% dplyr::filter(cell == name, state == "inactive") %>% select(element)
  # Extact the data
  df <- fread(male_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  #Intra-active Island Distance 
  df1_active <- select(df, activeIslands$element)  
  df1_active <- df1_active %>% t() %>% as.data.frame()
  df2_active <- select(df1_active, activeIslands$element) %>% t() %>% as.data.frame()
  df2_active <- df2_active %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "active", cell = name)
  #Intra-inactive Island Distance 
  df1_inactive <- select(df, inactiveIslands$element)  
  df1_inactive <- df1_inactive %>% t() %>% as.data.frame()
  df2_inactive <- select(df1_inactive, inactiveIslands$element) %>% t() %>% as.data.frame()
  df2_inactive <- df2_inactive %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "inactive", cell = name)
  intraIslandDistance <- rbind(intraIslandDistance, df2_active, df2_inactive)
}

intraIslandDistance_MOR28_male <- intraIslandDistance

intraIslandDistance = data.frame()
for (i in c(1:length(female_files))){
  # Extract Islands
  name <- female_files[i]
  name <- sub(".female.20k.*", "", name)
  name <- sub(".*pd/APM28DIP.","",name)
  print(name)
  activeIslands <- islands_percell_dataframe_MOR28_230927 %>% dplyr::filter(cell == name, state == "active") %>% select(element)
  inactiveIslands <- islands_percell_dataframe_MOR28_230927 %>% dplyr::filter(cell == name, state == "inactive") %>% select(element)
  # Extact the data
  df <- fread(female_files[i],col.names = names$V1) %>% as.data.frame()
  rownames(df) <- names$V1
  #Intra-active Island Distance 
  df1_active <- select(df, activeIslands$element)  
  df1_active <- df1_active %>% t() %>% as.data.frame()
  df2_active <- select(df1_active, activeIslands$element) %>% t() %>% as.data.frame()
  df2_active <- df2_active %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "active", cell = name)
  #Intra-inactive Island Distance 
  df1_inactive <- select(df, inactiveIslands$element)  
  df1_inactive <- df1_inactive %>% t() %>% as.data.frame()
  df2_inactive <- select(df1_inactive, inactiveIslands$element) %>% t() %>% as.data.frame()
  df2_inactive <- df2_inactive %>% as.matrix() %>% reshape2::melt() %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::rename(islandLeft = Var1, islandRight = Var2) %>%
    mutate(state = "inactive", cell = name)
  intraIslandDistance <- rbind(intraIslandDistance, df2_active, df2_inactive)
}

intraIslandDistance_MOR28_female <- intraIslandDistance

#### Combinding Data Number of Islands

islands_percell_dataframe_P2_230927$geno <- "P2"
islands_percell_dataframe_MOR28_230927$geno <- "MOR28"

islands_percell_dataframe <- rbind(islands_percell_dataframe_P2_230927, islands_percell_dataframe_MOR28_230927)%>%
  group_by(geno, cell, state) %>%
  summarise(total = n())

a <- ggplot(islands_percell_dataframe, aes(x = state, y = total)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylab("GIs per GIH") + 
  theme_classic()

islands_percell_dataframe %>%
  group_by(state) %>%
  summarise(
    mean = mean(total),
    sd = sd(total)
  )

#### Combinding Data Distance of Islands

intraIslandDistance <- rbind(intraIslandDistance_MOR28_female, intraIslandDistance_MOR28_male,
                             intraIslandDistance_P2_female, intraIslandDistance_P2_male)

b <- ggplot(intraIslandDistance, aes(x = state, y = value)) + 
  geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.2) + 
  ylab("GI Distance") + 
  theme_classic()

a + b

intraIslandDistance %>% group_by(state) %>%
  summarise(mean = mean(value),
            sd = sd(value))

islands_percell_dataframe %>% group_by(state) %>%
  summarise(mean = mean(total),
            sd = sd(total))
islands_percell_dataframe %>% group_by(geno, cell) %>% summarise(n()) %>% nrow()

write.table(intraIslandDistance, "SuppFig6d_GIDistance.txt", sep = '\t')
write.table(islands_percell_dataframe, "SuppFig6e_GINumber.txt", sep = '\t')

