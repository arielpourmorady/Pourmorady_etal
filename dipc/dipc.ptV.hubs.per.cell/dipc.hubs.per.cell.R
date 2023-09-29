library(pheatmap)
library(dendextend)
library(ggplot2)
library(tidyr)
library(dplyr)

# set directory to current directory
# this directory should contain all files in the GitHub folder

setwd("/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptV.hubs.per.cell") 

`%notin%` <- Negate(`%in%`)

##############################################################################
##############################################################################
##############################################################################

# Cell Names
P2_male_cells <- read.table("/media/storageE/ariel_dipc/P2-ImmFix-dipc/male.txt")
P2_female_cells <- read.table("/media/storageE/ariel_dipc/P2-ImmFix-dipc/female.txt")
P2_cells <- rbind(P2_male_cells, P2_female_cells) %>%
  dplyr::rename(cell = V1) %>%
  mutate(geno = "P2")

MOR28_male_cells <- read.table("/media/storageE/ariel_dipc/APM28DIP/male.txt")
MOR28_female_cells <- read.table("/media/storageE/ariel_dipc/APM28DIP/female.txt")
MOR28_cells <- rbind(MOR28_male_cells, MOR28_female_cells) %>%
  dplyr::rename(cell = V1) %>%
  mutate(geno = "MOR28")

all_cells <- rbind(P2_cells, MOR28_cells)

#####################################################################################################
# The code below is used to extract the average/stdev, and various ranges of the # GIs per active GIH.
#####################################################################################################

names <- read.table("/media/storageE/ariel_dipc/annotations/Greek_Islands.nonX.names")

islands_percell_dataframe_P2 <- read.table("islands_percell_dataframe_P2_230927.txt",sep = "\t") %>% # generated using dipc.inactive.hub.extraction.R
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell) %>%
  dplyr::filter(state == "active") %>%
  mutate(geno = "P2")

islands_percell_dataframe_MOR28 <- read.table("islands_percell_dataframe_MOR28_230927.txt",sep = "\t") %>% # generated using dipc.inactive.hub.extraction.R
  separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("name" = paste(island, haplo, sep = ",")) %>%
  dplyr::select(name, state, cell) %>%
  dplyr::filter(state == "active") %>%
  mutate(geno = "MOR28")

islands_percell_dataframe <- rbind(islands_percell_dataframe_P2, islands_percell_dataframe_MOR28) %>%
  group_by(cell, geno, state) %>%
  summarise(sum = n()) %>%
  filter(state == 'active')

islands_percell_dataframe <- left_join(all_cells, islands_percell_dataframe, by = c("geno", "cell")) 
islands_percell_dataframe$state[is.na(islands_percell_dataframe$state)] <- 'active'
islands_percell_dataframe$sum[is.na(islands_percell_dataframe$sum)] <- 0

########################################################################################################
# The code below is used to extract the average/stdev, and various ranges of the # GIs per inactive GIH.
########################################################################################################

hubs_percell_dataframe_P2 <- read.table("hubs_percell_dataframe_P2_230927.txt",sep = "\t") %>%
  dplyr::filter(cluster != 1) %>% # remove the active hub
  mutate(geno = "P2", state = 'inactive', sum = total) %>%
  select(!c(cluster,total))

hubs_percell_dataframe_MOR28 <- read.table("hubs_percell_dataframe_MOR28_230927.txt",sep = "\t") %>%
  dplyr::filter(cluster != 1) %>% # remove the active hub
  mutate(geno = "MOR28", state = 'inactive', sum = total) %>%
  select(!c(cluster,total))

hubs_percell_dataframe <- rbind(hubs_percell_dataframe_P2, hubs_percell_dataframe_MOR28)

########################################################################################################
# PLOT
########################################################################################################

All_Islands <- rbind(islands_percell_dataframe, hubs_percell_dataframe)

ggplot(All_Islands, aes(sum, fill = state)) + 
  geom_density(adjust = 2, alpha = 0.4) +
  scale_fill_manual(values = c('forestgreen', 'black')) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$sum)) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$sum) -sd(islands_percell_dataframe$sum), linetype = 'dashed' ) + 
  geom_vline(xintercept = mean(islands_percell_dataframe$sum) +sd(islands_percell_dataframe$sum), linetype = 'dashed' ) + 
  theme_classic()

All_Islands %>% group_by(state) %>%
  summarise(mean = mean(sum),
            sd = sd(sum),
            count = n())

write.table(All_Islands, file = "/media/storageE/ariel/R/finalpaper_August2023/dipc/dipc.ptV.hubs.per.cell/SuppFig5d_HubsPerCell.txt",
            sep = '\t')

