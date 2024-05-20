library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratDisk)
library(Matrix)
library(ggpubr)
library(ggrepel)

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_doublets_filtered_NormalizeData_cell_types.rds")

kallisto.output.dictionary <- list(
  `VEXAS_BM2_CD34_plus` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_plus/counts_unfiltered/",
  `VEXAS_BM2_CD34_negative` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_minus/counts_unfiltered/",
  `VX_BM14_Fpos` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fpos/counts_unfiltered/",
  `VX_BM14_Fneg` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fneg/counts_unfiltered/",
  `VX_BM12_A` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM12_A/counts_unfiltered/",
  `VX_BM13_B` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM13_B/counts_unfiltered/",
  `VX_BM15a` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15a/counts_unfiltered/",
  `VX_BM15b` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15b/counts_unfiltered/"
)

## Pull suffixes and prefixes for sample
suffix.df <- rna.obj@meta.data %>% rownames_to_column("Barcode") %>% group_by(orig.ident) %>% summarize(suffix = gsub("landau_.*-1", "", Barcode) %>% unique())
suffix.dictionary <- suffix.df$suffix
names(suffix.dictionary) <- suffix.df$orig.ident


prefix.df <- rna.obj@meta.data %>% rownames_to_column("Barcode") %>% group_by(orig.ident) %>% summarize(prefix = gsub("[ACGT]*-1_[0-9]{1,2}", "", Barcode) %>% unique())
prefix.dictionary <- prefix.df$prefix
names(prefix.dictionary) <- prefix.df$orig.ident


pull.adt.counts <- function(sample.name) {
  print(sample.name)
  base.asap.path <- kallisto.output.dictionary[[sample.name]]
  adt_data <- Matrix::readMM(paste0(base.asap.path, "cells_x_features.mtx"))
  adt_bcs <- read.csv(paste0(base.asap.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(prefix.dictionary[[sample.name]], V1, "-1", suffix.dictionary[[sample.name]])) %>% pull(V1) 
  adt_tag_names <- read.csv(paste0(base.asap.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
  rownames(adt_data) <- adt_bcs
  colnames(adt_data) <- adt_tag_names
  adt_data <- as.matrix(adt_data)
  return(adt_data)
}


## Sample list to pull counts for
sample.list <- rna.obj$orig.ident %>% table() %>% names
names(sample.list) <- sample.list
sample.list <- sample.list[!grepl("BM1_|UPN", names(sample.list))]

## Pull ADT counts and combine
adt.counts.list <- lapply(sample.list, pull.adt.counts)
adt.counts <- do.call(rbind, adt.counts.list)

## Add NA counts for remaining samples
# If there are any atac barcodes not present in adt, add them to the adt output
rna.obj.bcs <- rna.obj@meta.data %>% rownames()
bc_missing_from_adt <- rna.obj.bcs[!(rna.obj.bcs %in% rownames(adt.counts))]
print("Number of barcodes missing from adt.counts data for:")
print(length(bc_missing_from_adt))
missing_matrix <- matrix(NA_real_, ncol=ncol(adt.counts), nrow=length(bc_missing_from_adt), dimnames=list(bc_missing_from_adt, colnames(adt.counts)))
adt.counts <- rbind(adt.counts, missing_matrix)

## Add as a new assay, normalize, and scale
# rna.obj@assays$ADT_CLR <- NULL
adt.counts <- adt.counts[rna.obj@meta.data %>% rownames, ]
rna.obj[["ADT_CLR"]] <- CreateAssayObject(counts = t(adt.counts))
DefaultAssay(rna.obj) <- "ADT_CLR"
rna.obj <- NormalizeData(rna.obj, normalization.method = "CLR", margin = 2, assay = "ADT_CLR")
rna.obj <- ScaleData(rna.obj)

## Correct CD90 for BM2
BM2.barcodes <- rna.obj@meta.data %>% filter(orig.ident == "VEXAS_BM2_CD34_plus" | orig.ident == "VEXAS_BM2_CD34_negative") %>% rownames()
rna.obj@assays$ADT_CLR@counts["Hu.CD90", BM2.barcodes] <- NA_real_
rna.obj@assays$ADT_CLR@data["Hu.CD90", BM2.barcodes] <- NA_real_
rna.obj@assays$ADT_CLR@scale.data["Hu.CD90", BM2.barcodes] <- NA_real_



saveRDS(rna.obj, "~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_doublets_filtered_NormalizeData_cell_types.rds")

# 
# ################
# 
# adt.counts.total <- data.frame()
# for (sample in names(sample.list)) {
#   adt.dir.name = sample.list[[sample]]
#   adt.counts <- pull.adt.counts(sample, adt.dir.name)
#   print(dim(adt.counts))
#   adt.counts.total <- rbind(adt.counts.total, adt.counts)
# }
# 
# 
# 
# 
# # If there are any atac barcodes not present in adt, add them to the protein data
# atac.obj.bcs <- atac.obj@meta.data %>% rownames()
# bc_missing_from_adt <- atac.obj.bcs[!(atac.obj.bcs %in% rownames(adt.counts.total))]
# print("Number of barcodes missing from adt data:")
# print(length(bc_missing_from_adt))
# missing_matrix <- matrix(NA_real_, nrow=length(bc_missing_from_adt), ncol=ncol(adt.counts.total), dimnames=list(bc_missing_from_adt, colnames(adt.counts.total)))
# adt.counts.total <- rbind(adt.counts.total, missing_matrix)
# 
# ## Also filter for just ATAC barcodes in the protein data
# adt.counts.total <- adt.counts.total[atac.obj.bcs, ]
# 
# ## Create assay object
# atac.obj[["ADT_CLR"]] <-  CreateAssayObject(counts = t(adt.counts.total))
# 
# # Normalize the ADT data
# DefaultAssay(atac.obj) <- "ADT_CLR"
# atac.obj <- NormalizeData(atac.obj, normalization.method = "CLR", margin = 2, assay = "ADT_CLR")
# atac.obj <- ScaleData(atac.obj, assay = "ADT_CLR")
# 
# ## Done!
# 

# 
# 
# ####################### Plot some specific genes ##################################
# 
# 
# View(rownames(atac.obj[["ADT"]]) %>% sort())
# 
# p1 <- VlnPlot(atac.obj, features = "Hu.CD19", group.by = "predicted.l2") + NoLegend()
# p2 <- VlnPlot(atac.obj, features = "Hu.CD14-M5E2", group.by = "predicted.l2") + NoLegend()
# p3 <- VlnPlot(atac.obj, features = "Hu.CD16", group.by = "predicted.l2") + NoLegend()
# 
# p1 / p2 / p3
# 
# p1 <- VlnPlot(atac.obj, features = "Hu.CD3-UCHT1", group.by = "predicted.l2") + NoLegend()
# p2 <- VlnPlot(atac.obj, features = "Hu.CD4-RPA.T4", group.by = "predicted.l2") + NoLegend()
# p3 <- VlnPlot(atac.obj, features = "Hu.CD8", group.by = "predicted.l2") + NoLegend()
# 
# p1 / p2 / p3
# 
# Assays(atac.obj)
# DefaultAssay(atac.obj) <- "RNA"
# FeaturePlot(atac.obj, features = "IL4")
# 
# FeaturePlot(atac.obj, features = c("adt_Hu.CD4-RPA.T4", "Hu.CD3-UCHT1", "Hu.CD8", "Hu.CD19"),  min.cutoff = "q10", max.cutoff = "q90")
# 
# to_plot <- atac.obj@meta.data %>% 
#   group_by(predicted.l2) %>% 
#   filter(n() > 15) %>% 
#   pull(predicted.l2)
# 
# Idents(atac.obj) <- "predicted.l2"
# DotPlot(atac.obj, assay = "ADT", idents = to_plot, features = c("Hu.CD19", "Hu.CD3-UCHT1", "adt_Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD14-M5E2",  "Hu.CD16", "Hu.CD56", "Isotype-RTK2071", "Isotype-RTK2758", "Isotype-RTK4174",  "Isotype-G0114F7", "Isotype-HTK888"), group.by = "predicted.l2") +
#   labs(title = "ADT tag average expression", subtitle = "BM5, BM6, BM7 CD34+") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_x_discrete(labels = c("CD19", "CD3", "CD4", "CD8", "CD14",  "CD16", "CD56", "Isotype-RTK2071", "Isotype-RTK2758", "Isotype-RTK4174",  "Isotype-G0114F7", "Isotype-HTK888"))
# 
# DotPlot(atac.obj, assay = "ADT", idents = to_plot, features = c("Isotype-RTK2071", "Isotype-RTK2758", "Isotype-RTK4174",  "Isotype-G0114F7", "Isotype-HTK888"), group.by = "predicted.l2") +
#   labs(title = "ADT tag average expression, controls", subtitle = "BM5, BM6, BM7 CD34+") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# 
# # Scatter plot - T cells
# bcs_feature_scatter <- atac.obj@meta.data %>% 
#   filter(predicted.l1 %in% c("CD4 T", "CD8 T")) %>% 
#   rownames()
# Idents(atac.obj) <- "predicted.l1"
# p1 <- FeatureScatter(atac.obj, feature1 = "adt_Hu.CD4-RPA.T4", feature2 = "Hu.CD8", cells = bcs_feature_scatter) + 
#   labs(title = "ADT tag expression", subtitle = "BM5, BM6, BM7 CD34+")
# p1$layers[[1]]$aes_params$alpha = 0.3
# p1
# 
# # Scatter plot - Monocytes
# bcs_feature_scatter <- atac.obj@meta.data %>% 
#   filter(predicted.l2 %in% c("CD14 Mono", "CD16 Mono")) %>% 
#   rownames()
# Idents(atac.obj) <- "predicted.l2"
# p1 <- FeatureScatter(atac.obj, feature1 = "Hu.CD14-M5E2", feature2 = "Hu.CD16", cells = bcs_feature_scatter) + 
#   labs(title = "ADT tag expression", subtitle = "BM5, BM6, BM7 CD34+")
# p1$layers[[1]]$aes_params$alpha = 0.3
# 
# # ViolinPlot
# bcs_vln <- atac.obj@meta.data %>% 
#   filter(predicted.l1 %in% c("CD4 T", "CD8 T", "B")) %>% 
#   rownames()
# atac.obj.subset <- subset(atac.obj, subset = predicted.l2 == "CD4 T")
# 
# Idents(atac.obj) <- "predicted.l1"
# VlnPlot(atac.obj, features = "Hu.CD16", group.by = "predicted.l1", idents = c("CD4 T", "CD8 T", "B"), split.by = "best_genotyping_call", cols = c("green", "red", "blue"))
# 
# VlnPlot(atac.obj, features = c("Hu.CD16", "adt_Hu.CD4-RPA.T4", "Hu.CD8"), group.by = "predicted.l1", idents = c("CD4 T", "CD8 T", "B"), split.by = "best_genotyping_call", cols = c("green", "red", "blue"))
# 
# CoveragePlot(atac.obj, region = "IL4", assay = "ATAC", group.by = "predicted.l1",)
# 
# ################### Make a nice cell counts plot ####################
# 
# 
# atac.obj@meta.data %>% 
#   group_by(predicted.l1) %>% 
#   summarize(count = n()) %>% 
#   ggplot(aes(x = reorder(predicted.l1, -count), y = log10(count), fill = predicted.l1, label = count)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'none', axis.title.x = element_blank()) +
#   geom_text(vjust = -0.5) +
#   labs(title = "Cell type counts", subtitle = "BM5, BM6, BM7 CD34+\nAnnotated via bridge integration, Azimuth BMMC dataset")
# 
# atac.obj@meta.data %>% 
#   group_by(predicted.l2) %>% 
#   summarize(count = n()) %>% 
#   ggplot(aes(x = reorder(predicted.l2, -count), y = log10(count), fill = predicted.l2, label = count)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'none', axis.title.x = element_blank()) +
#   geom_text(vjust = -0.5) +
#   labs(title = "Cell type counts", subtitle = "BM5, BM6, BM7 CD34+\nAnnotated via bridge integration, Azimuth BMMC dataset")
# 
# 
# #################### Plot average expression per cluster #########################
# 
# Idents(atac.obj) <- "predicted.l2"
# cluster.adt.averages <- AverageExpression(atac.obj, assays = c('ADT', 'ATAC', 'RNA'), slot="data")
# cluster.adt.averages$ADT[1:5, 1:5]
# cluster.adt.averages$ATAC[1:5, 1:5]
# cluster.adt.averages$RNA[1:5, 1:5]
# 
# gex %>% subset(predicted.celltype.l2 == "CD4 Naive") %>% 
#   VlnPlot(., features = c("adt_CD3D.1", "adt_CD4.1", "adt_CD8A.1"))
# 
# plots <- list()
# 
# gene_lists <- list (
#   b_cell_markers = c('CD19', 'MS4A1'),
#   t_cell_markers =  c('CD3D', 'CD4', 'CD8A'),
#   nk_markers = c('ITGAL', 'ITGAM', 'ITGAX', 'NCAM1', 'FCGR3A', 'IL7R'),
#   monocyte_markers = c('CD14', 'FCGR3A', 'NCAM1'),
#   granulocyte_markers = c('CD44', 'LAMP1'),
#   dc_markers = c('CR2', 'FCER2'),
#   megakaryocyte_markers = c("CD36", "ITGA2B", "GP1BB"),
#   extra = c('CD27', 'CD38', 'CD5')
# )
# 
# plist <- list()
# 
# for (gene_set in names(gene_lists)) {
#   for (gene in gene_lists[[gene_set]]) {
#     adt_data <- cluster.adt.averages$ADT[paste0("Hu.", gene), ]
#     rna_data <- cluster.adt.averages$RNA[gene, ]
#     
#     combined <- as.data.frame(t(rbind(adt_data, rna_data)))
#     
#     combined$cluster <- rownames(combined)
#     
#     p <- ggscatter(combined, 
#                    x = "rna_data", 
#                    y = "adt_data",
#                    # label = "cluster",
#                    add = "reg.line",  # Add regression line
#                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                    conf.int = TRUE, # Add confidence interval
#                    title = gene
#     ) + 
#       stat_cor(
#         method = "pearson",
#         label.y = 0.2,
#         label.x.npc = "center",
#         position = position_nudge(y = 0.015 * 0.5)
#       ) +
#       geom_text_repel(label = combined$cluster)
#     
#     plist[[gene]] <- p
#     print(p)
#   }
# }
# 
# 
# # fetch_and_format_adt_counts <- function(sample.name, atac.obj) {
# #     base.asap.path <- paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/output/", sample.name, "_CITE_reformatted/counts_unfiltered/")
# #     
# #     adt_data <- Matrix::readMM(paste0(base.asap.path, "cells_x_features.mtx"))
# #     adt_bcs <- read.csv(paste0(base.asap.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% pull(V1) %>% paste0(sample.name, "_", ., "-1")
# #     adt_tag_names <- read.csv(paste0(base.asap.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
# #     
# #     rownames(adt_data) <- adt_bcs
# #     colnames(adt_data) <- adt_tag_names
# #     
# #     # Transpose matrix so the ADT tags are rows and the barcodes are columns
# #     adt_data <- t(adt_data)
# #     
# #     # Check the level of mismatch between ADT and ATAC barcodes
# #     atac_bcs <- atac.obj@meta.data %>% filter(orig.ident == sample.name) %>% rownames()
# #     
# #     print(adt_bcs[1:5])
# #     print(atac_bcs[1:5])
# #     print(length(adt_bcs))
# #     print(length(atac_bcs))
# #     print(length(atac_bcs[!(atac_bcs %in% adt_bcs)]))
# #     print(length(adt_bcs[!(adt_bcs %in% atac_bcs)]))
# #     
# #     # Make a list of the barcodes present in atac, but missing in adt
# #     bc_missing_from_adt <- atac_bcs[!(atac_bcs %in% adt_bcs)]
# #     
# #     missing_matrix <- matrix(0, nrow=nrow(adt_data), ncol=length(bc_missing_from_adt), dimnames=list(rownames(adt_data), bc_missing_from_adt))
# #     adt_expanded <- cbind(adt_data, missing_matrix)
# #     return(adt_expanded[ , atac_bcs])
# # }
# 
# 
# adt_assay_obj <- CreateAssayObject(
#   counts = adt_expanded[ , atac_bcs]
# )
# 
# atac.obj.subset <- subset(atac.obj, orig.ident == "BM7_LLL_CD34_plus")
# atac.obj.subset[["ADT"]] <- adt_assay_obj
# 
# # Sanity check a tag
# rownames(atac.obj.subset@assays$ADT)
# VlnPlot(atac.obj.subset, features = c("Hu.CD19"), group.by = "predicted.l1")
# 
# # adt_data <- Read10X("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/output/BM6_LLL_CD34_plus_CITE_reformatted/counts_unfiltered/")
# # adt_files <- lapply(unique(atac.obj$orig.ident), function(sample.name) {
# #   print(sample.name)
# #   adt.path <- paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/", sample.name, "_CITE_reformatted/outs/filtered_feature_bc_matrix/")
# #   adt <- Read10X(adt.path)
# # })
# 
# adt_files <- list()
# adt_files[["BM6_LLL_CD34_plus"]] <- Read10X("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/BM6_LLL_CD34_plus_CITE_reformatted/outs/filtered_feature_bc_matrix/")
# colnames(adt_files[["BM6_LLL_CD34_plus"]]) <- paste0("BM6_LLL_CD34_plus", "_", colnames(adt_files[["BM6_LLL_CD34_plus"]]))
# 
# adt_files[["BM5_LLL_CD34_plus"]] <- Read10X("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/BM5_LLL_CD34_plus_cite_reformatted/outs/filtered_feature_bc_matrix/")
# colnames(adt_files[["BM5_LLL_CD34_plus"]]) <- paste0("BM5_LLL_CD34_plus", "_", colnames(adt_files[["BM5_LLL_CD34_plus"]]))
# 
# adt_files[["BM7_LLL_CD34_plus"]] <- Read10X("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/antibodies/BM7_LLL_CD34_plus_CITE_reformatted/outs/filtered_feature_bc_matrix/")
# colnames(adt_files[["BM7_LLL_CD34_plus"]]) <- paste0("BM7_LLL_CD34_plus", "_", colnames(adt_files[["BM7_LLL_CD34_plus"]]))
# 
# adt <- do.call(cbind, adt_files)
# 
# adt_bcs <- colnames(adt)
# atac_bcs <- rownames(atac.obj@meta.data)
# 
# print(length(atac_bcs[!(atac_bcs %in% adt_bcs)]))
# print(length(adt_bcs[!(adt_bcs %in% atac_bcs)]))
# 
# bc_missing_from_adt <- atac_bcs[!(atac_bcs %in% adt_bcs)]
# 
# missing_matrix <- matrix(0, nrow=nrow(adt), ncol=length(bc_missing_from_adt), dimnames=list(rownames(adt), bc_missing_from_adt))
# adt_expanded <- cbind(adt, missing_matrix)
# 
# adt_assay_obj <- CreateAssayObject(
#   counts = adt_expanded[ , atac_bcs]
# )
# 
# atac.obj[["ADT"]] <- adt_assay_obj
# rownames(atac.obj$ADT)
# 
# atac.obj@assays$ADT["Hu.CD14-M5E2", ][2000:2015]
# gex <- NormalizeData(gex, normalization.method = "CLR", margin = 2, assay = "ADT")
# gex <- ScaleData(gex, assay = "ADT")
# 
# FeaturePlot(atac.obj, "adt_Hu.CD4-RPA.T4")
# 
# # adt_assay <- CreateAssayObject(counts = tenx_data$`Antibody Capture`[ ,rownames(filter(atac.obj@meta.data, orig.ident == "")])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############## Example ########################
# 
# #Convert to Seurat object
# obj <- CreateSeuratObject(counts = BM2_cd34_negative_10x$`Gene Expression`, min.cells = 3, min.features = 200, project = "BM2_cd34_negative_CITE")
# BM2_cd34_positive <- CreateSeuratObject(counts = BM2_cd34_positive_10x$`Gene Expression`, min.cells = 3, min.features = 200, project = "BM2_cd34_positive_CITE")
# # mds_p605 <- CreateSeuratObject(counts = mds605$`Gene Expression`, min.cells = 3, min.features = 200, project = "MDS_P605")
# 
# #Create assay object to add ADT
# BM2_cd34_negative[["ADT"]] <- CreateAssayObject(counts = BM2_cd34_negative_10x$`Antibody Capture`[ ,rownames(BM2_cd34_negative@meta.data)])
# BM2_cd34_positive[["ADT"]] <- CreateAssayObject(counts = BM2_cd34_positive_10x$`Antibody Capture`[ ,rownames(BM2_cd34_positive@meta.data)])
# # mds_p605[["ADT"]] <- CreateAssayObject(counts = mds605$`Antibody Capture`[ ,rownames(mds_p605@meta.data)])
# 
# #Calculate percent mitochondrial genes ####
# BM2_cd34_negative <- PercentageFeatureSet(BM2_cd34_negative, pattern = "^MT-", col.name = "percent.mito")
# BM2_cd34_positive <- PercentageFeatureSet(BM2_cd34_positive, pattern = "^MT-", col.name = "percent.mito")
# 
# # mds_p605 <- PercentageFeatureSet(mds_p605, pattern = "^MT-", col.name = "percent.mito")
# 
# #Calculate percent ribosomal genes ####
# BM2_cd34_negative <- PercentageFeatureSet(BM2_cd34_negative, pattern = "^RPS|^RPL", col.name = "percent.ribo")
# BM2_cd34_positive <- PercentageFeatureSet(BM2_cd34_positive, pattern = "^RPS|^RPL", col.name = "percent.ribo")
# # mds_p605 <- PercentageFeatureSet(mds_p605, pattern = "^RPS|^RPL", col.name = "percent.ribo")
# 
# #Calculate the percent cell cycle genes #### 
# BM2_cd34_negative <- CellCycleScoring(BM2_cd34_negative, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
# BM2_cd34_positive <- CellCycleScoring(BM2_cd34_positive, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
# # mds_p605 <- CellCycleScoring(mds_p605, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
# 
# #Create list of Seurat objects to be integrated ####
# # mds.sample.list <- list("mds_p545"= mds_p545, "mds_p605"= mds_p605)
# BM2_sample_list <- list("BM2_cd34_negative"=BM2_cd34_negative, "BM2_cd34_positive"=BM2_cd34_positive)
