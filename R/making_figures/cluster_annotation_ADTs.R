library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
library(scales)
theme_set(theme_classic())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "cDC", "CD14 Mono", "B", "pDC", "CD4 T", "CD8 T", "NK", "Plasma")
rna.obj$CellType <- factor(rna.obj$cluster_celltype, levels = celltype.order)
table(rna.obj$CellType, useNA = "always")

###################################### Make the ADT average expression per cluster plot ################################################

current.assay <- "ADT_DSB"

## Create an seurat object with only cells with non-zero ADT expression
adt.counts.per.bc <- colSums(rna.obj@assays[[current.assay]]@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.adt <- subset(rna.obj, cells = bcs.keep)

## Scale the DSB assay
rna.obj.adt <- ScaleData(rna.obj.adt, assay = current.assay)
rna.obj.adt@assays$ADT_DSB@data[1:5, 1:5]
rna.obj.adt@assays$ADT_DSB@scale.data[1:5, 1:5]

## Compile a data frame with the average expression per cluster
avg.exp <- lapply(levels(rna.obj.adt$CellType), function(x) {
  print(x)
  bc.list <- rna.obj.adt@meta.data %>% filter(CellType == x) %>% rownames()
  df <- data.frame(expression_mean = rowMeans(rna.obj.adt[[current.assay]]@scale.data[, bc.list], na.rm = T))
  colnames(df) <- x
  return(df)
})
avg.exp.matrix <- do.call(cbind, avg.exp)
avg.exp.df <- avg.exp.matrix %>% rownames_to_column("ADT") %>% pivot_longer(cols = colnames(avg.exp.matrix), names_to = "cluster", values_to = "average_expression")


## Manually selecting tags)
tags.keep <- c("Hu.CD90", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD36", "Hu.CD41", "Hu.CD123", "Hu.CD11c", "Hu.CD45RA",  "Hu.CD14-M5E2",  "Hu.CD20-2H7", "Hu.CD19", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD3-UCHT1", "Hu.CD56", "Hu.CD16", "Hu.CD7")
dev.off()
pdf("figures/current_figure_drafts/RNA_ADT_heatmap.pdf", width = 4, height = 4)
pheatmap::pheatmap(t(avg.exp.matrix[tags.keep, celltype.order]), 
                   # scale = "column",
                   cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-2, 2, length.out = 50))
dev.off()



######################################## Extra ###########################################

## Pull ADT markers for each cluster
Idents(rna.obj.adt) <- "CellType"
cluster.markers.adt <- FindAllMarkers(rna.obj.adt, assay = current.assay, only.pos = T, features = rna.obj.adt@assays$ADT_DSB@data %>% rownames(), group.by = "CellType", mean.fxn = rowMeans)

## Pull the top N unique tags for each cluster
num.markers.to.pull = 3
top.tags <- cluster.markers.adt %>%
  mutate(metric = -log10(p_val) * avg_log2FC) %>% 
  arrange(-metric) %>% 
  distinct(gene, .keep_all = T) %>%
  mutate(cluster = factor(cluster, levels = celltype.order)) %>% 
  group_by(cluster) %>% 
  slice_max(n = num.markers.to.pull, order_by = metric) %>% 
  arrange(cluster) %>% 
  pull(gene)

avg.exp.df %>% 
  filter(ADT %in% top.tags) %>% 
  ggplot(aes(x = ADT, y = cluster, color = average_expression)) +
  geom_point() +
  # scale_color_gradient2(low = "blue", mid = "white", high = "red", limits=c(-2, 2), oob=squish) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pheatmap::pheatmap(t(avg.exp.matrix[top.tags, celltype.order]), cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-1, 1, length.out = 50))



pct.exp <- lapply(levels(rna.obj.adt$CellType), function(x) {
  print(x)
  bc.list <- rna.obj.adt@meta.data %>% filter(CellType == x) %>% rownames()
  df <- data.frame(
    pct.exp = rowSums(rna.obj.adt@assays$ADT_CLR@counts[, bc.list] == 0, na.rm = T) / length(bc.list)
  )
  colnames(df) <- x
  return(df)
})
pct.exp.df <- do.call(cbind, pct.exp)
pct.exp.df <- pct.exp.df %>% rownames_to_column("ADT") %>% pivot_longer(cols = colnames(pct.exp.df), names_to = "cluster", values_to = "pct_expression")

m <- rna.obj.adt@assays$ADT_DSB@data %>% as.matrix() %>%  t()
df <- merge(rna.obj.adt@meta.data, m, by = "row.names")
df %>% 
  filter(seurat_clusters %in% c(3, 14, 15, 8, 27)) %>% 
  mutate(seurat_clusters = factor(seurat_clusters, levels = c(3, 14, 15, 8, 27))) %>% 
  ggplot(aes(x = seurat_clusters, y = Hu.CD90, fill = seurat_clusters)) + 
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("3", "14"), c("3", "15"), c("3", "8"), c("3", "27"))) +
  theme(legend.position = "none")

