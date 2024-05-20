library(Seurat)
library(tidyverse)
library(pheatmap)
library(patchwork)
set.seed(1234)

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_doublets_filtered_NormalizeData.rds")

## Pull marker genes for each assigned cell type
Idents(rna.obj) <- "cluster_celltype"
DefaultAssay(rna.obj) <- "ADT_DSB"
markers.all.clusters <- FindAllMarkers(rna.obj, logfc.threshold = 0.2, only.pos = T, assay = "ADT_DSB")

write_csv(markers.all.clusters, file = "data/cluster_markers/landau_and_NIH_08_20_23_resolution_1_cluster_celltype_ADT_DSB.csv")