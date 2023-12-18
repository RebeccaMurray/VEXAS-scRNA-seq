library(Seurat)
library(tidyverse)
library(patchwork)
theme_set(theme_bw())
set.seed(1)

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
# 
# 
# rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data.rds")
# 
# ########################################### Make QC plots ########################################
# 
# ## nCount_RNA
# sample.order <- rna.obj@meta.data %>% group_by(orig.ident) %>% summarize(mean_nCount_RNA = mean(nCount_RNA)) %>% arrange(-mean_nCount_RNA) %>% pull(orig.ident)
# rna.obj@meta.data %>% 
#   mutate(orig.ident = factor(orig.ident, levels = sample.order)) %>% 
#   ggplot(aes(x = orig.ident, y = log10(nCount_RNA), fill = orig.ident)) +
#   geom_boxplot() +
#   xlab(element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
# ggsave("figures/current_figure_drafts/qc_plots_nCount_RNA.pdf", device = "pdf", dpi = 300, width = 5, height = 4)
# 
# ## nFeature_RNA
# sample.order <- rna.obj@meta.data %>% group_by(orig.ident) %>% summarize(mean_nFeature_RNA = mean(nFeature_RNA)) %>% arrange(-mean_nFeature_RNA) %>% pull(orig.ident)
# rna.obj@meta.data %>% 
#   mutate(orig.ident = factor(orig.ident, levels = sample.order)) %>% 
#   ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) +
#   geom_boxplot() +
#   xlab(element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
# ggsave("figures/current_figure_drafts/qc_plots_nFeature_RNA.pdf", device = "pdf", dpi = 300, width = 5, height = 4)
# 
# ## percent.mito
# sample.order <- rna.obj@meta.data %>% group_by(orig.ident) %>% summarize(mean_pct_mito = mean(percent.mito)) %>% arrange(-mean_pct_mito) %>% pull(orig.ident)
# rna.obj@meta.data %>% 
#   mutate(orig.ident = factor(orig.ident, levels = sample.order)) %>% 
#   ggplot(aes(x = orig.ident, y = percent.mito, fill = orig.ident)) +
#   geom_violin(scale = "width") +
#   xlab(element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
# ggsave("figures/current_figure_drafts/qc_plots_percent_mito_RNA.pdf", device = "pdf", dpi = 300, width = 5, height = 4)
# 
# ## Sample of origin - split into separate panes
# p.sample.origin.split <- DimPlot(rna.obj, split.by = "orig.ident", group.by = "orig.ident", ncol = 9, pt.size = 0.1) + NoLegend() + 
#   coord_fixed() +
#   ggtitle(element_blank()) + theme(plot.title = element_text(size=2))
# p.sample.origin.split
# ggsave("figures/current_figure_drafts/qc_plots_rna_UMAP_orig_ident_split.pdf", plot = p.sample.origin.split, device = "pdf", dpi = 300, width = 18, height = 4, unit = "in")
# 
# 
# ########################################### UMAPs with cell type labels ########################################
# 
# p2 <- DimPlot(rna.obj, group.by = "predicted.celltype.l2", label = T, label.box = T, repel = T) + NoLegend()
# p3 <- DimPlot(rna.obj, group.by = "seurat_clusters", label = T, label.box = T, repel = T) + NoLegend()
# p4 <- DimPlot(rna.obj, group.by = "predicted_CellType", label = T, label.box = T, repel = T) + NoLegend()
# 
# p2 | p3
# p4 | p3
# 
# ########################################### Adjust clustering if necessary ##################################
# 
# ## Redo clustering if necessary
# DefaultAssay(rna.obj) <- "integrated"
# rna.obj <- FindClusters(rna.obj, resolution = 2, random.seed = 1)
# 
# ################################## Automatic assignment ######################################################
# 
# cell.type.mapping <- rna.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   count(predicted.celltype.l2) %>% 
#   slice_max(n) %>% 
#   select(seurat_clusters, predicted.celltype.l2) %>% 
#   mutate(predicted.celltype.l2 = gsub(".*CD8.*", "CD8 T", predicted.celltype.l2)) %>% 
#   mutate(predicted.celltype.l2 = gsub(".*CD4.*", "CD4 T", predicted.celltype.l2)) %>% 
#   mutate(predicted.celltype.l2 = gsub("Memory ", "", predicted.celltype.l2)) %>% 
#   mutate(predicted.celltype.l2 = gsub("Prog Mk", "MkP", predicted.celltype.l2)) %>% 
#   mutate(predicted.celltype.l2 = gsub("cDC.*", "cDC", predicted.celltype.l2))
# 
# ## Based on resolution == 2 - manual correction for these HSC clusters based on lack of CD90 expression
# cell.type.mapping$predicted.celltype.l2[cell.type.mapping$seurat_clusters == 35] <- "EMP"
# cell.type.mapping$predicted.celltype.l2[cell.type.mapping$seurat_clusters == 44] <- "EMP"
# 
# rna.obj$cluster_celltype <- plyr::mapvalues(rna.obj$seurat_clusters, from = cell.type.mapping$seurat_clusters, to = cell.type.mapping$predicted.celltype.l2)
# rna.obj$CellType <- rna.obj$cluster_celltype
# 
# ## Identify clusters we want to remove based on high mitochondrial content
# FeaturePlot(rna.obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "S.Score", "G2M.Score", "CC.Difference"), ncol = 3)
# VlnPlot(rna.obj, features = "percent.mito", group.by = "seurat_clusters", pt.size = 0) + NoLegend()
# clusters.remove <- c(10, 35, 36, 40, 46, 47, 48)
# 
# ## Make a nice plot documenting the average mito expression, to show why we're removing these clusters
# cluster.order <- rna.obj@meta.data %>% group_by(seurat_clusters) %>%  summarize(mean_percent_mito = mean(percent.mito)) %>% arrange(-mean_percent_mito) %>% pull(seurat_clusters)
# p.percent.mito.violins <- rna.obj@meta.data %>% 
#   mutate(`Keep?` = !(seurat_clusters %in% clusters.remove)) %>% 
#   mutate(seurat_clusters = factor(seurat_clusters, levels = cluster.order)) %>% 
#   ggplot(aes(x = seurat_clusters, y = percent.mito, fill = `Keep?`)) +
#   geom_violin(scale = "width", draw_quantiles = 0.5) +
#   scale_fill_manual(values = c(`FALSE` = "salmon", `TRUE` = "grey")) +
#   xlab("Seurat cluster")
# p.percent.mito.violins
# ggsave("figures/current_figure_drafts/percent_mito_violins_annotated.pdf", plot = p.percent.mito.violins, device = "pdf", dpi = 300, width = 10.5, height = 3, unit = "in")
# p.percent.mito.hist <- rna.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% summarize(mean_percent_mito = mean(percent.mito)) %>% 
#   ggplot(aes(x = mean_percent_mito)) +
#   geom_histogram(binwidth = 0.25) +
#   geom_vline(xintercept = 5, linetype = "longdash", color = "darkred") + 
#   xlab("Mean percent.mito per cluster")
# p.percent.mito.hist
# ggsave("figures/current_figure_drafts/percent_mito_hist.pdf", plot = p.percent.mito.hist, device = "pdf", dpi = 300, width = 3, height = 3, unit = "in")
# p.mito_rna.count.scatter <- rna.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize(mean_pct_mito = mean(percent.mito), mean_num_umis = mean(nCount_RNA)) %>% 
#   mutate(`Keep?` = !(seurat_clusters %in% clusters.remove)) %>% 
#   ggplot(aes(x = mean_pct_mito, y = mean_num_umis, label = seurat_clusters, fill = `Keep?`)) +
#   scale_fill_manual(values = c(`FALSE` = "salmon", `TRUE` = "grey")) +
#   geom_point() +
#   geom_label_repel() +
#   xlab("Mean percent mito reads") +
#   ylab("Mean UMI count")
# p.mito_rna.count.scatter
# ggsave("figures/current_figure_drafts/percent_mito_ncount_umi_scatter.pdf", plot = p.mito_rna.count.scatter, device = "pdf", dpi = 300, width = 5, height = 3.5, unit = "in")
# p.cluster.10 <- DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(seurat_clusters == 10) %>% rownames()) +
#   ggtitle("Cluster 10 highlighted")
# p.cluster.10
# 
# 
# ## Remove the cells from these clusters from the seurat object
# cell.barcodes.keep <- rna.obj@meta.data %>% filter(!(seurat_clusters %in% clusters.remove)) %>% rownames()
# rna.obj <- subset(rna.obj, cells = cell.barcodes.keep)
# 
# p2 <- DimPlot(rna.obj, group.by = "predicted.celltype.l2", label = T, label.box = T, repel = T) + NoLegend()
# p3 <- DimPlot(rna.obj, group.by = "cluster_celltype", label = T, label.box = T, repel = T, cols = cell.type.palette) + NoLegend()
# 
# p2 | p3
# 
# saveRDS(rna.obj, "~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")

#################################### Regenerate UMAP ##########################################

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")

## Regenerate the UMAP so that it is reproducible and can be used for the control projection
DefaultAssay(rna.obj) <- "integrated"
rna.obj <- RunUMAP(rna.obj, reduction = "pca", dims = 1:50, seed.use = 1, return.model = TRUE)

pdf("updated_UMAP.pdf")
DimPlot(rna.obj, group.by = "cluster_celltype", label = T, label.box = T, repel = T, cols = cell.type.palette) + NoLegend()
dev.off()

saveRDS(rna.obj, "~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")

# 
# ################################# Extra analysis #######################################
# 
# rna.obj@meta.data %>% filter(seurat_clusters == 7 & predicted.celltype.l1.score > 0.7) %>% count(predicted.celltype.l1) %>% arrange(-n) %>% head
# rna.obj@meta.data %>% filter(seurat_clusters == 35 & predicted.celltype.l2.score > 0) %>% count(predicted.celltype.l2) %>% arrange(-n) %>% head
# rna.obj@meta.data %>% filter(seurat_clusters == 27 & predicted.celltype.l2.score > 0) %>% count(predicted_CellType) %>% arrange(-n) %>% head
# rna.obj@meta.data %>% filter(seurat_clusters == 24) %>% count(predicted.celltype.l2, Genotype) %>% arrange(-n)
# rna.obj@meta.data %>% filter(seurat_clusters == 43) %>% count(Genotype)
# rna.obj@meta.data %>% filter(seurat_clusters == 8) %>% count(orig.ident) %>% arrange(-n) %>% head
# rna.obj@meta.data %>% filter(predicted.celltype.l2 == "Prog Mk") %>% count(seurat_clusters) %>% arrange(-n) %>% head
# 
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(seurat_clusters == 22) %>% rownames())
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(predicted.celltype.l2 == "HSC" & predicted.celltype.l2.score > 0.7) %>% rownames()) 
# 
# DimHeatmap(rna.obj, dims = 45:60, cells = 1000, balanced = TRUE)
# 
# p2 <- VlnPlot(rna.obj, features = "AVP", group.by = "seurat_clusters", pt.size = 0, idents = c(6, 19, 47, 7, 13, 27, 44, 15, 26, 25, 35, 36), assay = "RNA") + NoLegend()
# p1 <- VlnPlot(rna.obj, features = "HLF", group.by = "seurat_clusters", pt.size = 0, idents = c(6, 19, 47, 7, 13, 27, 44, 15, 26, 25, 35, 36), assay = "RNA") + NoLegend()
# p4 <- VlnPlot(rna.obj, features = "KIT", group.by = "seurat_clusters", pt.size = 0, idents = c(6, 19, 47, 7, 13, 27, 44, 15, 26, 25, 35, 36), assay = "RNA") + NoLegend()
# p2 / p3 / p4
# 
# ## format as heatmap
# m <- AverageExpression(rna.obj, assays = "RNA", features = c("AVP", "HLF", "SLAMF1", "KLF1", "HBD", "CA1", "HOXA9", "KIT", "LTB", "GATA1", "GATA2", "MPO", "CD38", "PF4"))
# m$RNA
# pheatmap::pheatmap(m$RNA[, c("6", "19", "47", "7", "13", "27", "44", "15", "26", "25", "35", "36", "14")], color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-1, 1, length.out = 50), scale = "row")
# 
# 
# ########################## Refining assignment- ADTs #########################
# adt.counts.per.bc <- colSums(rna.obj@assays$ADT_DSB@data)
# bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
# rna.obj.adt <- subset(rna.obj, cells = bcs.keep)
# p3 <- VlnPlot(rna.obj.adt, features = "Hu.CD90", group.by = "seurat_clusters", pt.size = 0, idents = c(6, 19, 47, 7, 13, 27, 44, 15, 26, 25, 35, 36)) + NoLegend() + coord_cartesian(ylim=c(0, 25))
# p4 <- VlnPlot(rna.obj.adt, features = "Hu.CD71", group.by = "seurat_clusters", pt.size = 0, idents = c(6, 19, 47, 7, 13, 27, 44, 15, 26, 25, 35, 36)) + NoLegend() + coord_cartesian(ylim=c(0, 25))
# 
# 
# current.assay <- "ADT_DSB"
# 
# ## Scale the DSB assay
# rna.obj.adt <- ScaleData(rna.obj.adt, assay = current.assay)
# rna.obj.adt@assays$ADT_DSB@data[1:5, 1:5]
# rna.obj.adt@assays$ADT_DSB@scale.data[1:5, 1:5]
# 
# ## Compile a data frame with the average expression per cluster
# avg.exp <- lapply(levels(rna.obj.adt$seurat_clusters), function(x) {
#   print(x)
#   if (x == 47) {
#     return(NA)
#   }
#   bc.list <- rna.obj.adt@meta.data %>% filter(seurat_clusters == x) %>% rownames()
#   df <- data.frame(expression_mean = rowMeans(rna.obj.adt[[current.assay]]@scale.data[, bc.list], na.rm = T))
#   colnames(df) <- x
#   return(df)
# })
# avg.exp.matrix <- do.call(cbind, avg.exp)
# avg.exp.df <- avg.exp.matrix %>% rownames_to_column("ADT") %>% pivot_longer(cols = colnames(avg.exp.matrix), names_to = "cluster", values_to = "average_expression")
# 
# 
# ## Manually selecting tags)
# tags.keep = c("Hu.CD90", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD36", "Hu.CD41", "Hu.CD45RA")
# pheatmap::pheatmap(avg.exp.matrix[tags.keep, c("6", "19", "7", "13", "27", "44", "15", "26", "25", "35", "36", "14")], 
#                    # scale = "column",
#                    cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-1, 1, length.out = 50))
# 
# 
# 
# 
# ######################### Extra #############################################
# 
# de.genes <- FindMarkers(rna.obj, ident.1 = "16", ident.2 = "3", group.by = "seurat_clusters", features = rna.obj@assays$RNA@counts %>% rownames())
# 
# Idents(rna.obj) <- "predicted.celltype.l2"
# se <- subset(rna.obj, predicted.celltype.l2 %in% c("HSC", "CD4 Memory", "CD14 Mono", "Late Eryth", "Early Eryth"))
# DimPlot(se, reduction = "pca", group.by = "predicted.celltype.l2", label = T, label.box = T)
# 
# ########################################### Make sure cell cycle, read depth, etc. isn't affecting integration ###########################
# 
# FeaturePlot(rna.obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", "S.Score", "G2M.Score"), ncol = 3, raster = F)
# FeaturePlot(rna.obj, features = "S.Score")
# FeaturePlot(rna.obj, features = "nCount_RNA", max.cutoff = 20000)
# VlnPlot(rna.obj, features = "nCount_RNA", group.by = "seurat_clusters", pt.size = 0) + NoLegend() + coord_cartesian(ylim = c(0, 10000))
# VlnPlot(rna.obj, features = "percent.mito", group.by = "seurat_clusters", pt.size = 0) + NoLegend()
# 
# rna.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize(mean_pct_mito = median(percent.mito)) %>% 
#   arrange(-mean_pct_mito) 
# 
# rna.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize(mean_pct_mito = median(percent.mito)) %>% 
#   pull(mean_pct_mito) %>% hist()
# 
# ########################################### Analyze the cluster cell types - misc. extra functions ########################################
# 
# ## Make a bar plot of the number of cell types per cluster - Azimuth
# md %>%
#   filter(predicted.celltype.l2.score > 0.5) %>%
#   group_by(predicted.celltype.l2, seurat_clusters) %>%
#   count() %>% 
#   group_by(seurat_clusters) %>%
#   top_n(5) %>% 
#   ggplot(aes(x = reorder(predicted.celltype.l2, -n), y = n, fill = predicted.celltype.l2)) +
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   facet_wrap(~ seurat_clusters, scales = "free", ncol = 9) +
#   # scale_fill_manual(values = cell.type.palette) +
#   theme(legend.position = "none") +
#   labs(x = "Azimuth prediction", y = "Cell count") +
#   ggtitle("Prediction score > 0.5")
# 
# ## Make a bar plot of the number of cell types per cluster - Bone Marrow Map
# md %>%
#   group_by(predicted_CellType, seurat_clusters) %>%
#   # filter(predicted_CellType_prob < 0.75) %>% 
#   count() %>% 
#   group_by(seurat_clusters) %>%
#   top_n(5) %>% 
#   ggplot(aes(x = reorder(predicted_CellType, -n), y = n, fill = predicted_CellType)) +
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   facet_wrap(~ seurat_clusters, scales = "free", ncol = 11) +
#   # scale_fill_manual(values = cell.type.palette) +
#   theme(legend.position = "none") +
#   labs(x = "BoneMarrowMap prediction", y = "Cell count") +
#   ggtitle("Celltype prob > 0.75")
# 
# ## Make a plot splitting by cluster
# DimPlot(rna.obj, group.by = "seurat_clusters", split.by = "seurat_clusters", shuffle = T, ncol = 8) + NoLegend()
# 
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(seurat_clusters == 16) %>% rownames()) + NoLegend()
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(seurat_clusters == 4) %>% rownames()) + NoLegend()
# 
# 
# ## Plot splitting by donor, cluster
# rna.obj$Donor <- gsub("(_CD34_minus|_CD34_plus|_plus|_minus|_A|_B|-T.*)$", "", rna.obj$Sample)
# DimPlot(rna.obj, group.by = "orig.ident", shuffle = T)
# DimPlot(rna.obj, group.by = "Donor", shuffle = T)
# DimPlot(rna.obj, group.by = "Donor", split.by = "Donor", shuffle = T, ncol = 3)
# 
# ## Plot by QC metrics
# VlnPlot(rna.obj, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0) + NoLegend() + ylim(c(0, 2000))
# 
# md %>% filter(seurat_clusters == 4) %>% count(predicted_CellType) %>% arrange(-n) %>% head
# 
# 
# ## Exploring clusters of interest
# md %>% filter(seurat_clusters == 31) %>% count(orig.ident) %>% arrange(-n)
# 
# ## Looking at smaller HSC clusters
# markers.31 <- FindMarkers(rna.obj, ident.1 = 31, ident.2 = c(5, 4, 18, 26, 14), group.by = "seurat_clusters", logfc.threshold = 0.2)
# plot_volcano(markers.31, stat_column = "p_val_adj")
# fgsea.outs <- run_fgsea_for_list(markers.31, pval_col = "p_val")
# 
# rna.obj.hscs <- subset(rna.obj, seurat_clusters %in% c(31, 26, 5, 18, 4, 14))
# markers.all.hspc.clusters <- FindAllMarkers(rna.obj.hscs, group.by = "seurat_clusters", logfc.threshold = 0.2, only.pos = T)
# 
# markers.31.vs.5 <- FindMarkers(rna.obj, ident.1 = 31, ident.2 = 5, group.by = "seurat_clusters", logfc.threshold = 0.2)
# plot_volcano(markers.31.vs.5, stat_column = "p_val_adj")
# fgsea.outs <- run_fgsea_for_list(markers.31.vs.5, pval_col = "p_val")
# 
# mut.vs.wt.cluster.5 <- FindMarkers(rna.obj, ident.1 = "MUT", ident.2 = "WT", group.by = "Genotype", subset.ident = 5, logfc.threshold = 0.2)
# plot_volcano(mut.vs.wt.cluster.5, stat_column = "p_val")
# fgsea.outs <- run_fgsea_for_list(mut.vs.wt.cluster.5, pval_col = "p_val")
# 
# ## Comparing two GMP clusters
# DefaultAssay(rna.obj) <- "RNA"
# gmp.14.vs.13 <- FindMarkers(rna.obj, ident.1 = 14, ident.2 = 13, group.by = "seurat_clusters")
# plot_volcano(gmp.14.vs.13, stat_column = "p_val_adj", effect_line = 1)
# 
# ## Comparing EMP/Eryth clusters
# DefaultAssay(rna.obj) <- "RNA"
# emp.subset <- subset(rna.obj, seurat_clusters %in% c(5, 12, 17, 19, 12, 22, 25, 24, 27))
# EMP.markers <- FindAllMarkers(emp.subset, group.by = "seurat_clusters")
# EMP.markers %>% 
#   group_by(cluster) %>% 
#   slice_min(order_by = p_val) -> top10
# 
# DoHeatmap(emp.subset, features = EMP.markers$gene, group.by = "seurat_clusters")
# 
# rna.obj$genotype <- if_else(is.na(rna.obj$Genotype_Approx_Match_No_Gene_Thresh), "NA", rna.obj$Genotype_Approx_Match_No_Gene_Thresh) 
# DimPlot(rna.obj, cols = c(HET="green", MUT="red", WT="blue"), group.by = "Genotype", order = c("HET", "MUT", "WT"), na.value = "grey80")
# 
# p1 <- DimPlot(rna.obj, cols = c(HET="green", MUT="red", WT="blue"), group.by = "Genotype", order = c("HET", "MUT", "WT"), na.value = "grey80")
# rna.obj.v1$Genotype <- if_else(is.na(rna.obj.v1$Genotype_Approx_Match_No_Gene_Thresh), "NA", rna.obj.v1$Genotype_Approx_Match_No_Gene_Thresh) 
# p0 <- DimPlot(rna.obj.v1, cols = c(HET="green", MUT="red", WT="blue"), group.by = "Genotype", order = c("HET", "MUT", "WT"), na.value = "grey80")
# p1 | p0
# 
# ## Highlighting specific clusters
# seurat_cluster_highlight = 41
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(seurat_clusters == seurat_cluster_highlight) %>% rownames()) + ggtitle(paste0("Seurat cluster = ", seurat_cluster_highlight))
# 
# ## Looking at breakdown of celltypes per cluster
# rna.obj@meta.data %>% filter(seurat_clusters == 36) %>% 
#   filter(predicted.celltype.l2.score > 0.3) %>%
#   count(predicted.celltype.l2) %>% arrange(-n) %>% head
# 
# ## Looking at breakdown of celltypes per cluster
# rna.obj@meta.data %>% filter(predicted.celltype.l2 == "CD16 Mono") %>% 
#   # filter(predicted.celltype.l2.score > 0.5) %>% 
#   count(seurat_clusters) %>% arrange(-n) %>% head
# 
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(predicted.celltype.l2 == "EMP" & predicted.celltype.l2.score > 0.7) %>% rownames()) + ggtitle("EMP")
# DimPlot(rna.obj, cells.highlight = rna.obj@meta.data %>% filter(predicted.celltype.l2 == "HSC" & predicted.celltype.l2.score > 0.7) %>% rownames()) + ggtitle("HSC")
# 
# ## Looking at some markers
# ## From Nir: HLA-DR, CD163 and AIF1 are Macrophage markers
# DefaultAssay(rna.obj) <- "RNA"
# FeaturePlot(rna.obj, features = c("AIF1", "CD163", "HLA-DRA", "HLA-DRB1", "HLA-DRB5"), min.cutoff ="q01", max.cutoff = "q99", ncol = 2)
# FeaturePlot(rna.obj, features = c("CD4", "CD8A"))
# DefaultAssay(rna.obj) <- "ADT"
# FeaturePlot(rna.obj, features = c("CD163"))
# 
# 
# ## Looking at progenitors
# highlights <- list()
# for (cell.type in c("EMP")) {
#   highlights[[cell.type]] <- rna.obj@meta.data %>% filter(predicted.celltype.l2 == cell.type & predicted.celltype.l2.score > 0.6) %>% rownames()
# }
# 
# DimPlot(rna.obj, cells.highlight = highlights, cols.highlight = c("red", "blue", "orange", "green", "purple", "gray"), shuffle = T)
# 
# 
# DimPlot(rna.obj, group.by = "seurat_clusters", split.by = "seurat_clusters", ncol = 7) + NoLegend()
# 
# rna.obj@meta.data %>% 
#   filter(predicted.celltype.l2 == "EMP") %>% 
#   filter(predicted.celltype.l2.score > 0.3) %>% 
#   count(seurat_clusters) %>% 
#   arrange(-n)
# 
# rna.obj@meta.data %>% 
#   filter(seurat_clusters == 27) %>% 
#   filter(predicted.celltype.l2.score > 0.7) %>%
#   count(predicted.celltype.l2) %>% 
#   arrange(-n)
# 
# ## Looking at CD90 expression
# rna.obj.subset <- subset(rna.obj, Donor %in% c("UPN2", "UPN12", "UPN13", "UPN14", "UPN15"))
# df <- rna.obj.subset
# 
# 
# ########################################### Assign a cell type to each cluster ########################################
# 
# ## Based on 2 resolution clustering, NormalizeData, with NIH 5' samples included, percent.mito regressed
# ## percent.mito = 10, 30 pcs used for integration and clustering
# # celltype.conversions <- list(
# #   `0` = "CD8 T",
# #   `1` = "HSC",
# #   `2` = "CD8 T",
# #   `3` = "LMPP", ## Alternatively early GMP?
# #   `4` = "EMP", ## Also has a lot of LMPP, and BoneMarrowMap says LMPP
# #   `5` = "CD14 Mono",
# #   `6` = "CD14 Mono",
# #   `7` = "Late Eryth",
# #   `8` = "HSC",
# #   `9` = "CD4 T",
# #   `10` = "HSC", # MPP according to bonemarrowmap
# #   `11` = "CD8 T",
# #   `12` = "CD4 T",
# #   `13` = "GMP-Neu",
# #   `14` = "CD8 T",
# #   `15` = "MkP", ## Also equal ish number Early Eryth
# #   `16` = "Early Eryth",
# #   `17` = "GMP-Mono",
# #   `18` = "cDC",
# #   `19` = "B",
# #   `20` = "CD8 T",
# #   `21` = "HSC",
# #   `22` = "Late Eryth",
# #   `23` = "CD8 T",
# #   `24` = "NK",
# #   `25` = "CD8 T",
# #   `26` = "pDC",
# #   `27` = "EMP", ## Also equal number Early Eryth
# #   `28` = "Early Eryth",
# #   `29` = "Late Eryth",
# #   `30` = "Early Eryth", ## Also a lot of LMPP, cycling
# #   `31` = "CD14 Mono",
# #   `32` = "CD8 T",
# #   `33` = "Late Eryth",
# #   `34` = "BaEoMa",
# #   `35` = "CD8 T",
# #   `36` = "HSC",
# #   `37` = "CD8 T",
# #   `38` = "CD14 Mono",
# #   `39` = "Plasma",
# #   `40` = "CD8 T",
# #   `41` = "CD14 Mono",
# #   `42` = "CD4 T",
# #   `43` = "CD14 Mono",
# #   `44` = "CLP",
# #   `45` = "Late Eryth",
# #   `46` = "Platelet",
# #   `47` = "GMP-Mono",
# #   `48` = "Late Eryth"
# # )
# 
# ## Based on 2 resolution clustering, NormalizeData, with NIH 5' samples included, nCount_RNA, S.Score, G2M.Score, percent.mito regressed
# ## percent.mito = 10, 50 pcs used for integration and clustering
# celltype.conversions <- list(
#   `0` = "HSC",
#   `1` = "CD8 T",
#   `2` = "CD8 T",
#   `3` = "CD14 Mono",
#   `4` = "CD4 T",
#   `5` = "HSC",
#   `6` = "LMPP",
#   `7` = "CD8 T",
#   `8` = "HSC",
#   `9` = "Late Eryth",
#   `10` = "GMP-Neu",
#   `11` = "Late Eryth",
#   `12` = "Late Eryth",
#   `13` = "HSC",
#   `14` = "CD8 T",
#   `15` = "MkP", ## Also lots of early eryth
#   `16` = "Early Eryth",
#   `17` = "Early Eryth",
#   `18` = "B",
#   `19` = "cDC",
#   `20` = "CD8 T",
#   `21` = "EMP",
#   `22` = "CD14 Mono",
#   `23` = "CD8 T",
#   `24` = "CD8 T",
#   `25` = "NK",
#   `26` = "pDC",
#   `27` = "CD4 T",
#   `28` = "CD14 Mono",
#   `29` = "Early Eryth", ## also lots of lmpp?
#   `30` = "GMP-Mono",
#   `31` = "HSC",
#   `32` = "GMP-Mono",
#   `33` = "BaEoMa",
#   `34` = "CD14 Mono",
#   `35` = "Plasma",
#   `36` = "CD14 Mono",
#   `37` = "CD14 Mono",
#   `38` = "CD8 T",
#   `39` = "CLP",
#   `40` = "Early Eryth",
#   `41` = "CD8 T",
#   `42` = "MkP",
#   `43` = "Late Eryth",
#   `44` = "Stromal"
# )
# 
# celltype.conversions <- list(
#   `4` = "HSC",
#   `0` = "HSC",
#   `21` = "HSC",
#   `11` = "EMP",
#   `8` = "LMPP",
#   `6` = "MkP"
# )
# 
# rna.obj$cluster_celltype <- plyr::mapvalues(rna.obj$seurat_clusters, from = names(celltype.conversions), to = unlist(celltype.conversions))
# rna.obj$CellType <- rna.obj$cluster_celltype
# 
# p2 <- DimPlot(rna.obj, group.by = "predicted.celltype.l2", label = T, label.box = T, repel = T) + NoLegend()
# p3 <- DimPlot(rna.obj, group.by = "cluster_celltype", label = T, label.box = T, repel = T) + NoLegend()
# 
# p2 | p3
# 
# 
# 
