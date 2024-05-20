library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
theme_set(theme_classic())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "cDC", "CD14 Mono", "B", "pDC", "CD4 T", "CD8 T", "NK", "Plasma")
rna.obj$CellType <- factor(rna.obj$cluster_celltype, levels = celltype.order)
table(rna.obj$CellType, useNA = "always")

###################################### Gene expression heatmap to annotate cell types ################################################

cluster.markers <- read_csv("data/cluster_markers/VEXAS_final_cluster_markers_20231108.csv")

## Pull the top N unique tags for each cluster
num.markers.to.pull = 5
top.genes <- cluster.markers %>%
  filter(!is.na(cluster)) %>%
  arrange(p_val, -avg_log2FC) %>%
  distinct(gene, .keep_all = T) %>%
  mutate(cluster = factor(cluster, levels = celltype.order)) %>%
  group_by(cluster) %>%
  arrange(p_val) %>%
  slice_min(n = num.markers.to.pull, order_by = p_val) %>%
  slice_max(n = num.markers.to.pull, order_by = avg_log2FC) %>%
  pull(gene)

## Scale the RNA assay for just these genes
rna.obj <- ScaleData(rna.obj, features = top.genes, assay = "RNA")

## Downsample to 50 BCs per cluster
downsampled.bc.list <- rna.obj@meta.data %>%
  rownames_to_column("Barcode") %>%
  group_by(CellType) %>%
  slice_sample(n = 50, replace = T) %>%
  pull(Barcode)

## Pull metadata for annotation row
annotation.df <- rna.obj@meta.data[downsampled.bc.list, ] %>%
  select(CellType)

## Heatmap
cluster.colors <- cell.type.palette[levels(rna.obj$CellType)]
dev.off()
pdf("figures/current_figure_drafts/RNA_cluster_annotation_marker_gene_heatmap_colors_updated.pdf", width = 8, height = 10)
pheatmap::pheatmap(rna.obj@assays$RNA@scale.data[top.genes, downsampled.bc.list],
                   annotation_col = annotation.df,
                   annotation_colors = list(`CellType` = cluster.colors),
                   show_colnames = F, cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c("navy", "white", "red"))(50),  breaks = seq(-2, 2, length.out = 50))
dev.off()


