#### Trajectory Analysis ####
# Here I am using monocle to
# do trajectory analysis

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyverse)
library(cowplot)

set.seed(1234)
options(scipen=999)

moe <- readRDS("/data/finalpaper_August2023/multiome_partI.input_processing/moe.rds")

######################################################
######################################################
############## Pseudotime Projection ################# 
######################################################
######################################################

neurons <- subset(moe, idents = c(0, 4, 10))
DefaultAssay(moe) <- "RNA"


moe_rna.cds <- as.cell_data_set(moe)
moe_rna.cds <- cluster_cells(cds = moe_rna.cds, reduction_method = "UMAP")
moe_rna.cds <- learn_graph(moe_rna.cds, use_partition = TRUE)

# ORDERING CELLS #

#####################
### VERY IMPORTANT###
#####################
## Make sure to pick the BEGINNING most cell in the GBC cluster
## all the way to the left
## this will ensure accurate pseudotiming
##

moe_rna.cds <- order_cells(moe_rna.cds, reduction_method = "UMAP")

DimPlot(moe, reduction = "umap", label = TRUE)

plot_cells(
  cds = moe_rna.cds,
  color_cells_by = "pseudotime",
  norm_method = c("size_only"),
  show_trajectory_graph = FALSE, 
  label_branch_points = FALSE, 
  label_leaves = FALSE, 
  label_roots = FALSE,
  scale_to_range = TRUE
)

moe_rna.cds@principal_graph_aux@listData$UMAP$pseudotime[is.infinite(moe_rna.cds@principal_graph_aux@listData$UMAP$pseudotime)] <- NA 

moe <- AddMetaData(
  object = moe,
  metadata = moe_rna.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

# Visualising

neurons <- subset(moe, idents = c(0, 4, 10))

plot1 <- DimPlot(neurons, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend() + xlim (c(-10, 5)) + ylim (c(-13, -1))

plot2 <- FeaturePlot(neurons, 
                     reduction = "umap",
                     features = "pseudotime", 
                     order = TRUE) & scale_color_viridis_c() & xlim (c(-10, 5)) & ylim (c(-13, -1))

plot1 + plot2


###########################################################
###########################################################
# Confirmation of Pseudotime with Gene Expression Markers # 
###########################################################
###########################################################

metadata <- data.frame(neurons@meta.data)

rna <- GetAssayData(object = neurons, 
                    assay = "SCT",
                    slot = "data") 

cells <- Cells(rna)
Ascl1 <- rna[rna@Dimnames[[1]] == "Ascl1"]
Neurod1 <- rna[rna@Dimnames[[1]] == "Neurod1"]
Tex15 <- rna[rna@Dimnames[[1]] == "Tex15"]
Gap43 <- rna[rna@Dimnames[[1]] == "Gap43"]
Carns1 <- rna[rna@Dimnames[[1]] == "Carns1"]

df2 <- data.frame(cells, Ascl1, Neurod1, Tex15, Gap43, Carns1)
df2 <- df2 %>% remove_rownames %>% column_to_rownames(var="cells")

time <- metadata %>% dplyr::select(pseudotime)
time$pseudotime <- round(time$pseudotime)
time[is.na(time)] <- 0
time <- merge(time, df2, by = "row.names")

df <- time %>% 
  group_by(pseudotime) %>% 
  summarise(
    Ascl1_mean = mean(Ascl1),
    Ascl1_se = sqrt(var(Ascl1) / length(Ascl1)),
    Neurod1_mean = mean(Neurod1),
    Neurod1_se = sqrt(var(Neurod1) / length(Neurod1)),
    Tex15_mean = mean(Tex15),
    Tex15_se = sqrt(var(Tex15) / length(Tex15)),
    Gap43_mean = mean(Gap43),
    Gap43_se = sqrt(var(Gap43) / length(Gap43)),
    Carns1_mean = mean(Carns1),
    Carns1_se = sqrt(var(Carns1) / length(Carns1))
  )

df$Ascl1_mean <- df$Ascl1_mean/(max(df$Ascl1_mean))
df$Ascl1_se <- df$Ascl1_se/(max(df$Ascl1_mean))

df$Neurod1_mean <- df$Neurod1_mean/(max(df$Neurod1_mean))
df$Neurod1_se <- df$Neurod1_se/(max(df$Neurod1_mean))

df$Tex15_mean <- df$Tex15_mean/(max(df$Tex15_mean))
df$Tex15_se <- df$Tex15_se/(max(df$Tex15_mean))

df$Gap43_mean <- df$Gap43_mean/(max(df$Gap43_mean))
df$Gap43_se <- df$Gap43_se/(max(df$Gap43_mean))

df$Carns1_mean <- df$Carns1_mean/(max(df$Carns1_mean))
df$Carns1_se <- df$Carns1_se/(max(df$Carns1_mean))

write.table(df, file = '/data/finalpaper_August2023/multiome_partII.pseudotime_projection/Fig1c_MarkerGenesOverPseudotime.txt', 
            sep = '\t', col.names = TRUE)

ggplot(df, aes(x = pseudotime)) + 
  geom_line(aes(y = Ascl1_mean), colour = "red1") +
  geom_errorbar(aes(ymin=Ascl1_mean-Ascl1_se, ymax=Ascl1_mean+Ascl1_se), colour = "red1") +
  geom_line(aes(y = Neurod1_mean), colour = "darkorange") +
  geom_errorbar(aes(ymin=Neurod1_mean-Neurod1_se, ymax=Neurod1_mean+Neurod1_se), colour = "darkorange") +
  geom_line(aes(y = Tex15_mean), colour = "yellow2") +
  geom_errorbar(aes(ymin=Tex15_mean-Tex15_se, ymax=Tex15_mean+Tex15_se), colour = "yellow2") +
  geom_line(aes(y = Gap43_mean), colour = "seagreen") +
  geom_errorbar(aes(ymin=Gap43_mean-Gap43_se, ymax=Gap43_mean+Gap43_se), colour = "seagreen") +
  geom_line(aes(y = Carns1_mean), colour = "royalblue3") +
  geom_errorbar(aes(ymin=Carns1_mean-Carns1_se, ymax=Carns1_mean+Carns1_se), colour = "royalblue3") + 
  theme_light() + theme(element_blank()) + ylab("Norm Exp")


###########################################################
###########################################################
######## Expression of Single Genes over Pseudotime ####### 
###########################################################
###########################################################

# cells <- Cells(rna)
# Lhx2 <- rna[rna@Dimnames[[1]] == "Lhx2"]
# 
# df2 <- data.frame(cells, Lhx2)
# df2 <- df2 %>% remove_rownames %>% column_to_rownames(var="cells")
# 
# time <- metadata %>% dplyr::select(pseudotime)
# time$pseudotime <- round(time$pseudotime)
# time[is.na(time)] <- 0
# time <- merge(time, df2, by = "row.names")
# 
# df <- time %>% 
#   group_by(pseudotime) %>% 
#   summarise(
#     Lhx2_mean = mean(Lhx2),
#     Lhx2_se = sqrt(var(Lhx2) / length(Lhx2))
#   )
# 
# df$Lhx2_mean <- df$Lhx2_mean/(max(df$Lhx2_mean))
# df$Lhx2_se <- df$Lhx2_se/(max(df$Lhx2_mean))
# 
# ggplot(df, aes(x = pseudotime)) + 
#   geom_line(aes(y = Lhx2_mean), colour = "red1") +
#   geom_errorbar(aes(ymin=Lhx2_mean-Lhx2_se, ymax=Lhx2_mean+Lhx2_se), colour = "red1") +
#   theme_light() + theme(element_blank()) + ylab("Norm Exp")

#saveRDS(neurons, file = "/data/finalpaper/multiome_partII.pseudotime_projection/neurons.rds")
