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

#################################################################################
#################################################################################
#################################################################################
#################################################################################

##################################################################################
##################################################################################  
############## Cumulative GI Accessibility per GI per Pseudotime ################# 
##################################################################################
##################################################################################

gi_bed <- as.data.frame(read.table("/data/multiome/Greek_Islands.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE)) %>%
  mutate(loc = paste(V1, V2, V3, sep = "-")) %>%
  dplyr::select(loc)

df <- FeatureMatrix(
  fragments,
  features = gi_bed[,1],
  cells = Cells(neurons),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE) %>% as.data.frame()
df$loc <- rownames(df)
df <- melt(df) #%>% dplyr::filter(value > 0) 

df <- df %>% mutate(cell = variable) %>% dplyr::select(-variable)

pseudotime_df <- neurons@meta.data %>% dplyr::select(pseudotime, nCount_ATAC)
pseudotime_df <- pseudotime_df %>% mutate(cell = row.names(pseudotime_df))

df1 <- left_join(df, pseudotime_df, by = "cell")
df1 <- df1 %>% mutate(value = value/nCount_ATAC) %>%
  mutate(pseudotime = round(pseudotime))

df1 <- df1[order(df1$pseudotime),]
df1$cell <- factor(df1$cell, levels = unique(df1$cell))

df2 <- df1 %>% 
  group_by(pseudotime,loc) %>% 
  summarise(mean = mean(value)) %>%
  separate(loc, sep = "-", into = c("chr", "beg", "end")) %>%
  separate(chr, sep = 'r', into =c(NA, "chr")) %>%
  mutate(ordering = paste(chr, beg, end, sep = "-")) #%>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(zero = 0)

df2 <- df2[order(df2$chr),]
df2$ordering <- factor(df2$ordering, levels = unique(df2$ordering))

a <- ggplot(df2, aes(x = pseudotime, y = ordering, fill = mean)) + geom_tile() +
  theme_classic() + 
  scale_fill_gradientn(colors = c("white", "firebrick4")) +
  coord_fixed() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank())

b <- ggplot(df2, aes(y = ordering, x = 0, fill = as.character(chr))) + 
  geom_tile()+
  theme_void()

b + a + plot_layout(guides = "collect", widths = c(1, 10))


