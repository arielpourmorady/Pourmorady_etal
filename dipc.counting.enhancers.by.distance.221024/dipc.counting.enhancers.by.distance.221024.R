library(pheatmap)
library(dendextend)
library(reshape2)
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)

# P2 Analysis

# Enhancers per Active Hub
names_col <- read.table("/media/storageE/ariel_dipc/annotations/OR.no_GI_overlap_50kb.nonX.names") %>%
  separate(V1, into = c("chr", "loc", "olfr", "allele"), sep = ",") %>%
  mutate(name = paste(olfr, allele, sep = "_"))
names_row <- read.table("/media/storageE/ariel_dipc/annotations/Greek_Islands.nonX.names")


### FUNCTIONS ###
enhancer_heatmap <- function(enhancers_distance, olfr, chr) {
  df <- enhancers_distance %>% 
    filter(OR == olfr, !grepl(chr,island), value < 12.5) %>% 
    mutate(distance= cut(value, breaks=seq(0, 12.5, 2.5))) %>%
    group_by(cell, distance) %>%
    summarise(count = n()) %>%
    mutate(distance = as.numeric(distance)*2.5) %>%
    mutate(count = count/(3.14*2*distance)) %>%
    dcast(distance ~ cell, value.var = "count")
  df[is.na(df)] <- 0
  rownames(df) <- df$distance
  df <- df[,-1]
  df <- melt(as.matrix(df))
  return(df)
}

enhancer_norm_distance <- function(enhancers_distance, olfr, chr) {
  df <- enhancers_distance %>% 
    filter(OR == olfr, !grepl(chr,island), value < 12.5) %>% 
    mutate(distance= cut(value, breaks=seq(0, 12.5, 2.5))) %>%
    group_by(cell, distance) %>%
    summarise(count = n()) %>%
    mutate(distance = as.numeric(distance)*2.5) %>%
    mutate(count = count/(3.14*2*distance)) %>%
    dcast(distance ~ cell, value.var = "count")
  df[is.na(df)] <- 0
  rownames(df) <- df$distance
  df <- df[,-1]
  df <- melt(as.matrix(df))
  df <- df %>% group_by(Var1) %>%
    summarise(
      mean = mean(value),
      se = sd(value)/sqrt(length(value))
    )
  return(df)
}

t_test_fxn <- function(df1, df2, num) {
  v1 <- df1 %>% filter(Var1 == num) %>% dplyr::select(value)
  v2 <- df2 %>% filter(Var1 == num) %>% dplyr::select(value)
  p <- t.test(v1, v2)$p.value
  return(p)
}

enhancers_bydist <- function(path){
  df <- read.table(path, col.names = names_col$name)
  df$island <- names_row$V1
  df <-  melt(df) %>% dplyr::rename("OR" = "variable")
  df$cell <- path
  return(df)
}

# Extracting P2 Files #

male_P2_files <- list.files(path = "/media/storageE/ariel_dipc/P2-ImmFix-dipc/matrix.P2-ImmFix-dipc.male.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd",
                         "*.20k.1.clean.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd.txt",
                         full.names = TRUE,
                         recursive = FALSE)

male_P2_enhancers_bydist <- lapply(male_P2_files, enhancers_bydist)
male_P2_enhancers_bydist <- as.data.frame(do.call(rbind, male_P2_enhancers_bydist))

female_P2_files <- list.files(path = "/media/storageE/ariel_dipc/P2-ImmFix-dipc/matrix.P2-ImmFix-dipc.female.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd",
                           "*.20k.1.clean.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd.txt",
                           full.names = TRUE,
                           recursive = FALSE)

female_P2_enhancers_bydist <- lapply(female_P2_files, enhancers_bydist)
female_P2_enhancers_bydist <- as.data.frame(do.call(rbind, female_P2_enhancers_bydist))

enhancers_distance_P2 <- rbind(male_P2_enhancers_bydist, female_P2_enhancers_bydist)

# Extracting mor28 Files #

male_mor28_files <- list.files(path = "/media/storageE/AP_LB_dipc/matrix.APM28DIP.male.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd",
                               "*.20k.1.clean.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd.txt",
                               full.names = TRUE,
                               recursive = FALSE)

male_mor28_enhancers_bydist <- lapply(male_mor28_files, enhancers_bydist)
male_mor28_enhancers_bydist <- as.data.frame(do.call(rbind, male_mor28_enhancers_bydist))

female_mor28_files <- list.files(path = "/media/storageE/AP_LB_dipc/matrix.APM28DIP.female.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd",
                                 "*.20k.1.clean.Greek_Islands.nonX.L1.OR.no_GI_overlap_50kb.nonX.L2.pd.txt",
                                 full.names = TRUE,
                                 recursive = FALSE)

female_mor28_enhancers_bydist <- lapply(female_mor28_files, enhancers_bydist)
female_mor28_enhancers_bydist <- as.data.frame(do.call(rbind, female_mor28_enhancers_bydist))

enhancers_distance_mor28 <- rbind(male_mor28_enhancers_bydist, female_mor28_enhancers_bydist)

### Single Cell Heatmaps in P2 Data ###

p2_bl6 <- enhancer_heatmap(enhancers_distance_P2, "Olfr17_0", "chr7")
p2_CAST <- enhancer_heatmap(enhancers_distance_P2, "Olfr17_1", "chr7")
p2_bl6_dist <- enhancer_norm_distance(enhancers_distance_P2, "Olfr17_0", "chr7")
p2_CAST_dist <- enhancer_norm_distance(enhancers_distance_P2, "Olfr17_1", "chr7")

a <- ggplot(p2_bl6, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, max(p2_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

b <- ggplot(p2_CAST, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white",  high = "firebrick3", limits = c(0, max(p2_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

c <- ggplot(p2_bl6_dist, aes(x = Var1, y = mean)) +
  geom_line(group = 1, colour = "darkgreen") + 
  geom_point(colour = "darkgreen") + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), colour = "darkgreen", width = 0.5) + 
  geom_line(data = p2_CAST_dist, aes(x = Var1, y = mean), group = 1, colour = "red") + 
  geom_point(data = p2_CAST_dist, aes(x = Var1, y = mean), colour = "red") + 
  geom_errorbar(data = p2_CAST_dist, aes(ymax = mean + se, ymin = mean - se), colour = "red", width = 0.5) +
  theme_classic() + 
  ylim(0, 0.09) 

z <- c | (a + b + plot_layout(guides = "collect"))
z
for (i in seq(2.5, 12.5, 2.5)) {
  p <- t_test_fxn(p2_bl6, p2_CAST, i)
  print(p)
}

# analyzing mor28 in the p2 dataset

mor28_bl6 <- enhancer_heatmap(enhancers_distance_P2, "Olfr1507_0", "chr14")
mor28_CAST <- enhancer_heatmap(enhancers_distance_P2, "Olfr1507_1", "chr14")
mor28_bl6_dist <- enhancer_norm_distance(enhancers_distance_P2, "Olfr1507_0", "chr14")
mor28_CAST_dist <- enhancer_norm_distance(enhancers_distance_P2, "Olfr1507_1", "chr14")

d <- ggplot(mor28_bl6, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, max(mor28_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

e <- ggplot(mor28_CAST, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white",  high = "firebrick3", limits = c(0, max(mor28_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

f <- ggplot(mor28_bl6_dist, aes(x = Var1, y = mean)) +
  geom_line(group = 1, colour = "darkgreen") + 
  geom_point(colour = "darkgreen") + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), colour = "darkgreen", width = 0.5) + 
  geom_line(data = mor28_CAST_dist, aes(x = Var1, y = mean), group = 1, colour = "red") + 
  geom_point(data = mor28_CAST_dist, aes(x = Var1, y = mean), colour = "red") + 
  geom_errorbar(data = mor28_CAST_dist, aes(ymax = mean + se, ymin = mean - se), colour = "red", width = 0.5) +
  ylim(0, 0.09) + 
  theme_classic() 

y <- f | (d + e + plot_layout(guides = "collect"))
y
for (i in seq(2.5, 12.5, 2.5)) {
  p <- t_test_fxn(mor28_bl6, mor28_CAST, i)
  print(p)
}

z / y 

### Single Cell Heatmaps in mor28 Data ###
# analyzing P2 in the mor28 data
p2_bl6 <- enhancer_heatmap(enhancers_distance_mor28, "Olfr17_0", "chr7")
p2_CAST <- enhancer_heatmap(enhancers_distance_mor28, "Olfr17_1", "chr7")
p2_bl6_dist <- enhancer_norm_distance(enhancers_distance_mor28, "Olfr17_0", "chr7")
p2_CAST_dist <- enhancer_norm_distance(enhancers_distance_mor28, "Olfr17_1", "chr7")

g <- ggplot(p2_bl6, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, max(p2_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

h <- ggplot(p2_CAST, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white",  high = "firebrick3", limits = c(0, max(p2_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

i <- ggplot(p2_bl6_dist, aes(x = Var1, y = mean)) +
  geom_line(group = 1, colour = "darkgreen") + 
  geom_point(colour = "darkgreen") + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), colour = "darkgreen", width = 0.5) + 
  geom_line(data = p2_CAST_dist, aes(x = Var1, y = mean), group = 1, colour = "red") + 
  geom_point(data = p2_CAST_dist, aes(x = Var1, y = mean), colour = "red") + 
  geom_errorbar(data = p2_CAST_dist, aes(ymax = mean + se, ymin = mean - se), colour = "red", width = 0.5) +
  theme_classic() +
  ylim(0, 0.09)

x <- i | (g + h + plot_layout(guides = "collect"))
x
for (i in seq(2.5, 12.5, 2.5)) {
  p <- t_test_fxn(p2_bl6, p2_CAST, i)
  print(p)
}

# analyzing mor28 in the mor28 dataset

mor28_bl6 <- enhancer_heatmap(enhancers_distance_mor28, "Olfr1507_0", "chr14")
mor28_CAST <- enhancer_heatmap(enhancers_distance_mor28, "Olfr1507_1", "chr14")
mor28_bl6_dist <- enhancer_norm_distance(enhancers_distance_mor28, "Olfr1507_0", "chr14")
mor28_CAST_dist <- enhancer_norm_distance(enhancers_distance_mor28, "Olfr1507_1", "chr14")

j <- ggplot(mor28_bl6, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, max(mor28_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

k <- ggplot(mor28_CAST, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white",  high = "firebrick3", limits = c(0, max(mor28_bl6$value))) + 
  theme_void() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(),
        axis.ticks.x = element_line()) 

l <- ggplot(mor28_bl6_dist, aes(x = Var1, y = mean)) +
  geom_line(group = 1, colour = "darkgreen") + 
  geom_point(colour = "darkgreen") + 
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), colour = "darkgreen", width = 0.5) + 
  geom_line(data = mor28_CAST_dist, aes(x = Var1, y = mean), group = 1, colour = "red") + 
  geom_point(data = mor28_CAST_dist, aes(x = Var1, y = mean), colour = "red") + 
  geom_errorbar(data = mor28_CAST_dist, aes(ymax = mean + se, ymin = mean - se), colour = "red", width = 0.5) +
  theme_classic() + 
  ylim(0,0.09)

w <- l | (j + k + plot_layout(guides = "collect"))
w
for (i in seq(2.5, 12.5, 2.5)) {
  p <- t_test_fxn(mor28_bl6, mor28_CAST, i)
  print(p)
}

x / w 

(z / y) | (x / w)

