library(dplyr)
library(Signac)
library(Seurat)
library(SeuratDisk)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(tidyr)
library(patchwork)
library(reshape2)

set.seed(1234)
options(scipen=999)

fragpath <- "/data/outs/atac_fragments.tsv.gz"
fragments <- CreateFragmentObject(fragpath)

neurons <- readRDS("/data/finalpaper_August2023/multiome_partIII.pseudotime_plots/neurons.rds")
mOSNs <- subset(neurons, idents = 0)


##################################################################################
##################################################################################  
############## Finding Cells Expressing an OR and identifying which OR ########### 
##################################################################################
##################################################################################

# You are writing code to add top expressed OR to metadata

olfr_anno <- read.table("/data/annotations/Olfr-CDS-ucsc.bed") %>%
  as.data.frame() %>%
  dplyr::select(V4, V1) %>%
  dplyr::rename("olfr" = "V4", "chr" = "V1")

rna <- GetAssayData(object = mOSNs, 
                    assay = "SCT",
                    slot = "data") 
olfrs <- rna %>% 
  row.names() %>% 
  as.data.frame() %>% 
  dplyr::filter(grepl("Olfr", .))  

# Extracting the most highly expressed Olfr per Cell
df <- rna %>% as.data.frame() 
maxOR_exp_percell_unique_cells <- df %>% mutate(genes = row.names(df)) %>% melt() %>% 
  dplyr::filter(genes %in% olfrs$.,
                value > 0) %>%
  group_by(variable) %>%
  dplyr::filter(value == max(value)) %>%
  group_by(variable) %>%
  summarise(count = n()) %>%
  dplyr::filter(count == 1)

maxOR_exp_percell<- df %>% mutate(genes = row.names(df)) %>% melt() %>% 
  dplyr::filter(genes %in% olfrs$.,
                value > 0) %>%
  group_by(variable) %>%
  dplyr::filter(value == max(value),
                variable %in% maxOR_exp_percell_unique_cells$variable)  %>%
  as.data.frame()

row.names(maxOR_exp_percell) <- maxOR_exp_percell$variable
maxOR_exp_percell <- maxOR_exp_percell %>% dplyr::select(genes)

mOSNs <- AddMetaData(object = mOSNs,
                     metadata = maxOR_exp_percell,
                     col.name = "OR_exp")

list_of_cells_expressing_OR <- row.names(maxOR_exp_percell)

#################################################################################
#################################################################################
#################################################################################
#################################################################################

gi_bed <- as.data.frame(read.table("/data/multiome/Greek_Islands.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

df <- FeatureMatrix(
  fragments,
  features = gi_bed[,1],
  cells = Cells(mOSNs),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE) %>% as.data.frame()
df$loc <- rownames(df)

df1 <- melt(df) %>% 
  dplyr::filter(value > 0) %>%
  dplyr::rename('cell' = 'variable', 
                'island' = 'loc') %>%
  dplyr::select(cell, island)

df2 <- df1 %>% dplyr::filter(cell %in% list_of_cells_expressing_OR)

df3 <- split(df2$island, df2$cell)
df3 <- df3[lapply(df3,length)>0]

df_clean <- data.frame()
for (i in c(1:length(df3))){
  #print(i)
  name_i <- names(df3[i])
  saveRDS(df_clean, file = '/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df_clean.rds')
  for (j in c(1:length(df3))){
    name_j <- names(df3[j])
    # print(j)
    if (i != j){
      a <- intersect(df3[[i]], df3[[j]]) %>% length()
      df_iterate <- data.frame(name_i, name_j, a)
      df_clean <- rbind(df_clean, df_iterate)
    }
  }
}

grouped_GI_freq <- df_clean %>%
  group_by(a) %>%
  summarise(count = n())

ggplot(grouped_GI_freq, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line() +
  scale_y_log10() + 
  theme_classic()


df4 <- data.frame()
df4_a <- data.frame()
for (k in grouped_GI_freq$a){
  num <- k
  df_clean_num_left <- df_clean %>%
    dplyr::filter(a == num) %>%
    dplyr::select(name_i) %>% 
    unique() %>% as.list()
  df_clean_num_right <- df_clean %>%
    dplyr::filter(a == num) %>%
    dplyr::select(name_j) %>% unique() %>% as.list()
  df_clean_union <- union(df_clean_num_left$name_i, df_clean_num_right$name_j)
  ORsPerCombo <- maxOR_exp_percell[df_clean_union,] %>% 
    as.data.frame() %>%
    dplyr::rename('olfr' = '.') %>%
    mutate(common_islands = num)
  ORsPerCombo_a <- ORsPerCombo %>%
    group_by(olfr) %>%
    summarise(count = n()/nrow(ORsPerCombo)) %>%
    mutate(common_islands = num)
  df4 <- rbind(df4, ORsPerCombo)
  df4_a <- rbind(df4_a, ORsPerCombo_a)
}

df4_a_mean_se <- df4_a %>%
  group_by(common_islands) %>% 
  summarise(
    mean = mean(count),
    se = sqrt(var(count) / length(count)))

df4_freq <- df4 %>% group_by(common_islands) %>%
  summarise(count = n(),
            freq = 1/count)

saveRDS(df4, file = '/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df4.rds')
saveRDS(df4_a, file = '/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df4_a.rds')
saveRDS(df4_a_mean_se, file = '/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df4_a_mean_se.rds')
saveRDS(df4_freq, file = '/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df4_freq.rds')

