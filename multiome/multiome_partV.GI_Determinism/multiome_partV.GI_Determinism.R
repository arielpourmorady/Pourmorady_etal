## 
##
##

# Before running the code below you have to run the 'cli_multiome_partVII.GI_Determinism.R' 
# within 'cli_program' directory



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

neurons <- readRDS("/data/finalpaper/multiome_partIII.pseudotime_plots/neurons.rds")
mOSNs <- subset(neurons, idents = 0)


#######################################################################
#######################################################################
############## Finding Cells and the ORs they Express #################
#######################################################################
#######################################################################

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
 
# #################################################################################
# ######################## Calculating Levels of Determinism ######################
# #################################################################################

df_clean <- readRDS('/data/finalpaper_revisions/multiome_partVII.GI_Determinism/cli_program/df_clean.rds')

grouped_GI_freq <- df_clean %>%
  group_by(a) %>%
  summarise(count = n()) %>%
  mutate(count = count/sum(count)*100)

a <- ggplot(grouped_GI_freq, aes(x = a, y = count)) + 
  geom_point() + 
  geom_line() +
  ylim(0,100) + 
  ylab('percent of cell pairs with GIs in common') + 
  xlab('number of GIs in common') + 
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
    mutate(common_islands = num,
           unique_cells = length(df_clean_union))
  df4 <- rbind(df4, ORsPerCombo)
  df4_a <- rbind(df4_a, ORsPerCombo_a)
}

df5 <- df4 %>%
  group_by(common_islands, olfr) %>%
  summarise(count = n())

ggplot(df4_a, aes(x = common_islands, y = count, group = common_islands)) + 
  geom_boxplot() +
  theme_classic()

df4_a_mean_se <- df4_a %>%
  group_by(common_islands) %>% 
  summarise(
    mean = mean(count),
    se = sqrt(var(count) / length(count)))

df4_freq <- df4 %>% group_by(common_islands) %>%
  summarise(count = n(),
            freq = 1/count)

c <- ggplot(df4_a_mean_se, aes(x = common_islands)) + 
  geom_point(aes(y = mean), colour = "black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour = "black") +
  geom_line(data = df4_freq, aes(x = common_islands, y = freq, color = 'firebrickred4', group = 1)) + 
  theme_classic() +
  scale_y_log10() + 
  theme(legend.position = "none") + 
  xlab('number of GIs in common') + 
  ylab('OR frequency')

cols = rainbow(length(unique(df4_a$olfr)), s=.6, v=.9)[sample(1:length(unique(df4_a$olfr)),length(unique(df4_a$olfr)))]
cols = magma(length(unique(df4_a$olfr)))[sample(1:length(unique(df4_a$olfr)),length(unique(df4_a$olfr)))]
d <- ggplot(df4_a, aes(x = common_islands, y = count, fill = olfr)) + 
  geom_bar(stat = 'identity')  + 
  scale_fill_manual(values=cols) + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab('number of GIs in common') + 
  ylab('OR frequency')

a + c + d