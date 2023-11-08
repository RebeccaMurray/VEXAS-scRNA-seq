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

## Pulling cluster markers
DefaultAssay(rna.obj) <- "RNA"
Idents(rna.obj) <- "CellType"
cluster.markers <- FindAllMarkers(rna.obj, assay = "RNA", group.by = "CellType", only.pos = T)
cluster.markers %>% 
  write_csv("data/cluster_markers/VEXAS_final_cluster_markers_20231187.csv")

