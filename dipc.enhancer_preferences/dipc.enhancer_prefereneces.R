library(pheatmap)
library(dendextend)
library(ggplot2)
library(dplyr)
library(patchwork)

`%notin%` <- Negate(`%in%`)
set.seed(99)

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

#### MOR 28 ####

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

distance_from_active_allele_dataframe_P2_point <- distance_from_active_allele_dataframe_P2 %>%
  group_by(island, OR) %>%
  summarise(mean = mean(active),
            se = sd(active)/sqrt(n())) %>%
  dplyr::filter(island != "active") %>%
  separate(island, sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("island" = paste(island, haplo, sep = ",")) %>%
  mutate("zero" = 0)
#distance_from_active_allele_dataframe_P2_point$island <- factor(distance_from_active_allele_dataframe_P2_point$island,levels=distance_from_active_allele_dataframe_P2_point$island, ordered=TRUE)

distance_from_active_allele_dataframe_MOR28_point <- distance_from_active_allele_dataframe_MOR28 %>%
  group_by(island, OR) %>%
  summarise(mean = mean(active),
            se = sd(active)/sqrt(n())) %>%
  dplyr::filter(island != "active") %>%
  separate(island, sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  mutate("island" = paste(island, haplo, sep = ",")) %>%
  mutate("zero" = 0)
#distance_from_active_allele_dataframe_MOR28_point$island <- factor(distance_from_active_allele_dataframe_MOR28_point$island,levels=distance_from_active_allele_dataframe_MOR28_point$island, ordered=TRUE)

distance_from_active_allele_dataframe_point <- rbind(distance_from_active_allele_dataframe_P2_point,
                                                     distance_from_active_allele_dataframe_MOR28_point) %>%
  separate(chr, sep = "r", into = c(NA, "chr")) %>%
  mutate("chr" = as.numeric(chr),
         "loc" = as.numeric(loc),
         "haplo" = as.numeric(haplo))
  

distance_from_active_allele_dataframe_point <- distance_from_active_allele_dataframe_point[
  with(distance_from_active_allele_dataframe_point, order(haplo, chr, loc, OR)),
  ] %>%
  mutate(ordering = paste(chr, loc, haplo, OR, sep = ","))

distance_from_active_allele_dataframe_point$ordering <- factor(distance_from_active_allele_dataframe_point$ordering,levels=distance_from_active_allele_dataframe_point$ordering, ordered=TRUE)

OR_colors <- c("P2" = "black", "MOR28" = "white")
a <- ggplot(distance_from_active_allele_dataframe_point, aes(x = ordering, y = mean)) + 
  geom_point(aes(fill = OR), pch = 21) + 
  scale_fill_manual(values = OR_colors) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
    
b <- ggplot(distance_from_active_allele_dataframe_point, aes(x = ordering, y = zero, fill = as.character(chr))) + 
    geom_tile()+
    theme_void()

a / b + plot_layout(heights = c(10, 1), guides = "collect") 

# Making Heatmaps with geom_tile

distance_from_active_allele_dataframe_P2_df <- distance_from_active_allele_dataframe_P2 %>%
  dplyr::filter(island != "active") %>%
  separate(island, sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  separate(chr, sep = "r", into = c(NA, "chr")) %>%
  mutate("chr" = as.numeric(chr)) %>%
  mutate("loc" = as.numeric(loc))

distance_from_active_allele_dataframe_MOR28_df <- distance_from_active_allele_dataframe_MOR28 %>%
  dplyr::filter(island != "active") %>%
  separate(island, sep = ",", into = c("chr", "loc", "island", "haplo")) %>%
  separate(chr, sep = "r", into = c(NA, "chr")) %>%
  mutate("chr" = as.numeric(chr)) %>%
  mutate("loc" = as.numeric(loc))

distance_from_active_allele_dataframe_df <- rbind(distance_from_active_allele_dataframe_P2_df, distance_from_active_allele_dataframe_MOR28_df)

distance_from_active_allele_dataframe_df <- distance_from_active_allele_dataframe_df[
  with(distance_from_active_allele_dataframe_df, order(haplo, chr, loc)),
  ] %>%
  mutate(name_cell = paste(cell, OR, sep = "_")) %>%
  mutate(ordering = paste(OR, chr, loc, haplo, sep = ","))

distance_from_active_allele_dataframe_df$ordering <- factor(distance_from_active_allele_dataframe_df$ordering, 
                                                               levels=unique(distance_from_active_allele_dataframe_df$ordering), 
                                                               ordered=TRUE)
distance_from_active_allele_dataframe_df$active[distance_from_active_allele_dataframe_df$active < 2.5] <- 0
distance_from_active_allele_dataframe_df$active[distance_from_active_allele_dataframe_df$active >= 2.5] <- 1



c <- ggplot(distance_from_active_allele_dataframe_df %>% filter(OR == "P2"), aes(x = ordering, y = cell, fill = as.character(active))) + 
  geom_tile() +
  scale_fill_manual(values = c("firebrick3", "white")) + 
  theme_void()
c

d <- ggplot(distance_from_active_allele_dataframe_df %>% filter(OR == "MOR28"), aes(x = ordering, y = cell, fill = as.character(active))) + 
  geom_tile() +
  scale_fill_manual(values = c("firebrick3", "white")) + 
  theme_void()
d

e <- ggplot(distance_from_active_allele_dataframe_point, aes(x = ordering, y = zero, fill = as.character(haplo))) + 
  geom_tile()+
  theme_void()

a /c / d / b /e + plot_layout(heights = c(10, 10, 10, 1, 1),guides = "collect")

 #####################################
##### GI GIH Probability Score ##### 
#####################################
# 
# islands_percell_dataframe_P2 <- read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/islands_percell_dataframe_P2.txt",sep = "\t") %>%
#   separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) 
# 
# islands_percell_dataframe_P2 <- islands_percell_dataframe_P2 %>%
#   group_by(chr, loc, island, haplo, state) %>%
#   summarise(prob = n()/(islands_percell_dataframe_P2 %>% group_by(cell) %>% count() %>% nrow())) %>%
#   mutate("OR" = "P2")
# 
# islands_percell_dataframe_MOR28 <- read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/islands_percell_dataframe_MOR28.txt",sep = "\t") %>%
#   separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) 
# 
# islands_percell_dataframe_MOR28 <- islands_percell_dataframe_MOR28 %>%
#   group_by(chr, loc, island, haplo, state) %>%
#   summarise(prob = n()/(islands_percell_dataframe_MOR28 %>% group_by(cell) %>% count() %>% nrow())) %>%
#   mutate("OR" = "MOR28")
# 
# islands_percell_dataframe_active <- rbind(islands_percell_dataframe_P2, islands_percell_dataframe_MOR28) %>% filter(state == "active")
# islands_percell_dataframe_inactive <- rbind(islands_percell_dataframe_P2, islands_percell_dataframe_MOR28) %>% filter(state == "inactive")
# 
# islands_percell_dataframe <- merge(islands_percell_dataframe_active, islands_percell_dataframe_inactive, by = c("chr", "loc", "island", "haplo", "OR")) %>%
#   mutate("prob" = prob.x/prob.y) %>%
#   select(chr, loc, island, haplo, OR, prob)  %>%
#   separate(chr, sep = "r", into = c(NA, "chr")) %>%
#   mutate("chr" = as.numeric(chr)) %>%
#   mutate("loc" = as.numeric(loc)) %>%
#   mutate("haplo" = as.numeric(haplo))
# 
# islands_percell_dataframe <- islands_percell_dataframe[
#   with(islands_percell_dataframe, order(haplo, chr, loc)),
#   ] %>%
#   mutate(ordering = paste(chr, loc, haplo, sep = ","))
# 
# islands_percell_dataframe$ordering <- factor(islands_percell_dataframe$ordering, 
#                                                             levels=unique(islands_percell_dataframe$ordering), 
#                                                             ordered=TRUE)
# 
# f <- ggplot(islands_percell_dataframe, aes(x = ordering, y = prob, colour = OR)) + 
#   geom_point() + theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# f/ b /e + plot_layout(heights = c( 10, 1, 1),guides = "collect")
# 
# 
# 
# islands_percell_dataframe_P2 <- read.table("/media/storageE/ariel/R/finalpaper/dipc.inactive.hub.extraction/islands_percell_dataframe_P2.txt",sep = "\t") %>%
#   separate(col = "element", sep = ",", into = c("chr", "loc", "island", "haplo")) 
# 
# islands_percell_dataframe_P2 <- islands_percell_dataframe_P2 %>%
#   group_by(chr, loc, island, haplo, state) %>%
#   summarise(prob = n()) %>%
#   mutate("OR" = "P2")
# 
# islands_percell_dataframe_P2_active <- islands_percell_dataframe_P2 %>% filter(state == "active")
# islands_percell_dataframe_P2_inactive <- islands_percell_dataframe_P2 %>% filter(state == "inactive")
# islands_percell_dataframe_P2 <- merge(islands_percell_dataframe_P2_active, islands_percell_dataframe_P2_inactive, by = c("chr", "loc", "island", "haplo")) %>%
#   mutate(ordering = paste(haplo, chr, loc, island, sep = ","))
# 
# ggplot(islands_percell_dataframe_P2 %>% filter(), aes(x = prob.x, y = prob.y)) + geom_point()
# 
# 
