library(Seurat)
library(ggh4x)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
theme_set(theme_classic())
set.seed(1)

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
current.plots.path <- c("figures/current_figure_drafts/")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")


######################################## Making a metadata data frame with coordinates ####################

## Adding pseudotime
vexas_pseudotime = readRDS('/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pseudotime/vexas_pseudotime_cells.Rds')
rna.obj <- AddMetaData(rna.obj, vexas_pseudotime@meta.data %>% select(monocle3_pseudotime))
rm(vexas_pseudotime)

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "cDC", "CD14 Mono", "B", "pDC", "CD4 T", "CD8 T", "NK", "Plasma")
rna.obj$CellType <- factor(rna.obj$cluster_celltype, levels = celltype.order)
table(rna.obj$CellType, useNA = "always")

## Get metadata plot with coordinates
md <- rna.obj[[]]
coords <- Embeddings(rna.obj[["umap"]])
md <- cbind(md, coords) %>% dplyr::rename(`UMAP1` = umap_1) %>% dplyr::rename(`UMAP2` = umap_2)
md <- md[sample(1:nrow(md)), ] ## Shuffle the rows in the dataframe, so plotting order is random

## To make small arrows in UMAP corner only
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "in")
)


###################################### Making plots - cell type UMAPS ################################################

## Seurat clusters
p.clusters <- DimPlot(rna.obj, group.by = "seurat_clusters", label = T, label.box = T, repel = T) + NoLegend() + 
  coord_fixed() +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) 
p.clusters
ggsave(paste0(current.plots.path, "/UMAP_RNA_seurat_clusters.pdf"), plot = p.clusters + theme(legend.position = "none"), device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")

## Sample of origin
p.sample.origin <- ggplot(md, aes(x = UMAP1, y = UMAP2, color = Donor )) +
  geom_point(size = 1) + theme_classic() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) +
  scale_color_manual(values = scales::hue_pal()(length(unique(md$Donor))))
p.sample.origin
ggsave(paste0(current.plots.path, "/UMAP_RNA_Donor.pdf"), plot = p.sample.origin, device = "pdf", dpi = 300, width = 7, height = 5, unit = "in")
ggsave(paste0(current.plots.path, "/UMAP_RNA_Donor_no_legend.pdf"), plot = p.sample.origin + theme(legend.position = "none"), device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")


## Sample of origin - split
p.sample.origin.split <- DimPlot(rna.obj, split.by = "Donor", group.by = "Donor", ncol = 5) + NoLegend()
p.sample.origin.split
ggsave("figures/current_figure_drafts/RNA_UMAP_sample_origin_split.tiff", plot = p.sample.origin.split, device = "tiff", dpi = 300, width = 12, height = 5.5)


## Bar plot of cell type proportions per donor
donor.order <- md %>% 
  group_by(Donor) %>% 
  mutate(total_per_sample = n()) %>% 
  group_by(Donor, total_per_sample) %>% 
  count(CellType) %>% 
  mutate(fraction = n / total_per_sample) %>% 
  filter(CellType == "CD14 Mono") %>% 
  arrange(-fraction) %>% 
  pull(Donor)
p.celltype.fractions.per.donor <- md %>% 
  mutate(Donor = factor(Donor, levels = donor.order)) %>% 
  group_by(Donor) %>% 
  mutate(total_per_sample = n()) %>% 
  group_by(Donor, total_per_sample) %>% 
  count(CellType) %>% 
  mutate(fraction = n / total_per_sample) %>% 
  ggplot(aes(x = Donor, y = fraction, fill = CellType)) +
  geom_col() +
  scale_fill_manual(values = cell.type.palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(element_blank()) +
  ylab("Percent of total cells") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "right") + guides(fill=guide_legend(ncol=3,byrow=TRUE))
p.celltype.fractions.per.donor
ggsave(paste0(current.plots.path, "/celltype_proportions_per_donor.pdf"), plot = p.celltype.fractions.per.donor, device = "pdf", dpi = 300, width = 6, height = 2.5, unit = "in")

## Azimuth labels
to.plot.tags <- rna.obj@meta.data %>% group_by(predicted.celltype.l2) %>% count() %>% filter(n > 50) %>% pull(predicted.celltype.l2)
rna.obj$predicted.celltype.l2.plot <- rna.obj$predicted.celltype.l2
rna.obj$predicted.celltype.l2.plot[!(rna.obj$predicted.celltype.l2.plot %in% to.plot.tags)] <- "Other"
p.azimuth <- DimPlot(rna.obj, group.by = "predicted.celltype.l2.plot", label = T, label.box = T, repel = T) + NoLegend() + 
  coord_fixed() +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) 
ggsave(paste0(current.plots.path, "/UMAP_RNA_azimuth_no_legend.pdf"), plot = p.azimuth + theme(legend.position = "none"), device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Plot the azimuth label/predicted celltype label overlap
predicted.celltype.l2.order <- c("HSC" ,  "LMPP"   ,   "EMP" , "Early Eryth", "Late Eryth",  "MkP", "BaEoMa",   "CLP", "GMP", "pre-mDC", "cDC", "CD14 Mono", "B", "pre-pDC", "pDC" , "CD4 T", "CD8 T", "NK", "Plasma", "CD16 Mono", "Macrophage", "ASDC", "MAIT", "NK CD56+", "Stromal")
p.azimuth.comparison <- md %>% 
  filter(!is.na(predicted.celltype.l2)) %>% ## Excludes 12 cells
  mutate(predicted.celltype.l2 = gsub(".*CD8.*", "CD8 T", predicted.celltype.l2)) %>%
  mutate(predicted.celltype.l2 = gsub(".*CD4.*", "CD4 T", predicted.celltype.l2)) %>%
  mutate(predicted.celltype.l2 = gsub(".* B", "B", predicted.celltype.l2)) %>% 
  mutate(predicted.celltype.l2 = gsub("Memory ", "", predicted.celltype.l2)) %>%
  mutate(predicted.celltype.l2 = gsub("Prog Mk", "MkP", predicted.celltype.l2)) %>%
  mutate(predicted.celltype.l2 = gsub("cDC.*", "cDC", predicted.celltype.l2)) %>%
  group_by(CellType) %>% 
  count(predicted.celltype.l2) %>% 
  mutate(percent = (n / sum(n))*100) %>% select(-n) %>% 
  pivot_wider(names_from = CellType, values_from = percent, values_fill = 0) %>% 
  column_to_rownames("predicted.celltype.l2")
pdf("figures/current_figure_drafts/RNA_cell_type_labels_heatmap.pdf", width = 5, height = 5)
pheatmap::pheatmap(p.azimuth.comparison[predicted.celltype.l2.order, ], cluster_rows = F, cluster_cols = F, main = "Assigned cell type labels", color = inferno(50))
dev.off()

## Cell type labels, no labels
p.celltypes <- ggplot(md, aes(x = UMAP1, y = UMAP2, fill = CellType, label = CellType)) +
  geom_point(pch = 21, size = 1, stroke = 0.1, color = "black") + ## size controls width of point, stroke controls width of border, color is border color
  scale_fill_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.celltypes
ggsave(paste0(current.plots.path, "/UMAP_RNA_celltypes_no_labels.pdf"), plot = p.celltypes, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Cell type labels, no labels, no dot outline
p.celltypes <- ggplot(md, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(pch = 19, size = 1, stroke = 0.1) + ## size controls width of point, stroke controls width of border, color is border color
  scale_color_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
p.celltypes
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
ggsave(paste0(current.plots.path, "/UMAP_RNA_celltypes_no_labels_no_outlines.pdf"), plot = p.celltypes, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")


## Predicted cell types, with labels (using DimPlot)
p.celltypes <- DimPlot(rna.obj, group.by = "CellType", cols = cell.type.palette, label = T, label.box = T, repel = T) + NoLegend() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + ggtitle(element_blank())
p.celltypes
ggsave(paste0(current.plots.path, "/UMAP_RNA_celltypes_with_labels.pdf"), plot = p.celltypes, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Predicted cell types, with labels (using DimPlot), VEXAS and Control cells
p.celltypes <- DimPlot(rna.obj, group.by = "CellType", cols = cell.type.palette, label = T, label.box = T, repel = T) + NoLegend() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + ggtitle(element_blank())
p.celltypes
ggsave(paste0(current.plots.path, "/UMAP_RNA_celltypes_with_labels_controls_included.pdf"), plot = p.celltypes, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Highlight just the monocytes
md.mono <- md %>% filter(cluster_celltype == "CD14 Mono")
md.na <- md
p.highlight.mono <- ggplot(md.na, aes(x = UMAP1, y = UMAP2)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.mono, pch = 21, size = 1, stroke = 0.1, color = "black", fill = "darkgreen") +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.highlight.mono
ggsave(paste0(current.plots.path, "/UMAP_mono_highlight.tiff"), plot = p.highlight.mono, device = "tiff", dpi = 300, width = 5, height = 4, unit = "in")

## Highlight just the T cells/NK cells
md.t.nk <- md %>% filter(cluster_celltype %in% c("CD4 T", "CD8 T", "NK"))
md.na <- md
p.highlighting.mature.lymphocytes <- ggplot(md.na, aes(x = UMAP1, y = UMAP2)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.t.nk, pch = 21, size = 1, stroke = 0.1, color = "black", fill = "#85378A") +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.highlighting.mature.lymphocytes
ggsave(paste0(current.plots.path, "/UMAP_RNA_T_NK_cell_highlight.tiff"), plot = p.highlighting.mature.lymphocytes, device = "tiff", dpi = 300, width = 5, height = 4, unit = "in")

## Highlighting just the HSCs
md.hsc <- md %>% filter(cluster_celltype == "HSC")
md.na <- md
p.hsc <- ggplot(md.na, aes(x = UMAP1, y = UMAP2)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.hsc, pch = 21, size = 1, stroke = 0.1, color = "black", fill = "#FF7B0B") +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.hsc
ggsave("figures/current_figure_drafts/UMAP_RNA_HSC_highlight.pdf", plot = p.hsc, device = "pdf", width = 5, height = 4)
ggsave("figures/current_figure_drafts/UMAP_RNA_HSC_highlight.tiff", plot = p.hsc, device = "tiff", width = 5, height = 4)

###################################### VEXAS and control data plots #######################################################

## VEXAS vs. control
rna.obj.control <- readRDS("/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pseudotime/vexas_pseudotime_cells.Rds")
md.control <- rna.obj.control[[]]
coords.control <- Embeddings(rna.obj.control[["UMAP"]])
md.control <- cbind(md.control, coords.control) %>% dplyr::rename(`UMAP1` = umap_1) %>% dplyr::rename(`UMAP2` = umap_2)
p.vexas.vs.control <- ggplot(md.control[sample(1:nrow(md.control)), ], aes(x = UMAP1, y = UMAP2, fill = Object)) +
  geom_point(pch = 21, size = 1, stroke = 0.1, color = "black") + ## size controls width of point, stroke controls width of border, color is border color
  scale_fill_manual(values = c(CONTROL = "grey50", VEXAS = "#2B858C")) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.vexas.vs.control
ggsave(paste0(current.plots.path, "/UMAP_RNA_vexas_vs_control_teal_grey.pdf"), plot = p.vexas.vs.control, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Cell type labels, no labels, controls included
p.celltypes.control <- ggplot(md.control, aes(x = UMAP1, y = UMAP2, fill = AllCellType, label = AllCellType)) +
  geom_point(pch = 21, size = 1, stroke = 0.1, color = "black") + ## size controls width of point, stroke controls width of border, color is border color
  scale_fill_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.celltypes.control
ggsave(paste0(current.plots.path, "/UMAP_RNA_celltypes_no_labels_controls_included.pdf"), plot = p.celltypes.control, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")


###################################### Making plots - genotyping UMAPS ################################################

## Print the genotyping fraction (overall)
rna.obj@meta.data %>% 
  count(Genotype_Approx_Match_No_Gene_Thresh) %>% 
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped)

## Print the genotyping fraction (by sample)
rna.obj@meta.data %>% 
  group_by(Donor, Genotype_Approx_Match_No_Gene_Thresh) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped)

## Print the genotyping fraction (by sample, stat summaries)
rna.obj@meta.data %>% 
  group_by(Donor, Genotype_Approx_Match_No_Gene_Thresh) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  ungroup() %>% 
  summarize(mean_frac_genotyped = mean(frac_genotyped), sd = sd(frac_genotyped))

## Genotyping UMAP (individual points)
mut.count <- md %>% filter(Genotype == "MUT") %>% nrow()
wt.count <- md %>% filter(Genotype == "WT") %>% nrow()
na.or.het.count <- md %>% filter(Genotype == "NA") %>% nrow()
mut.title <- paste0("MUT (n = ", format(mut.count, big.mark = ","), ")")
wt.title <- paste0("WT (n = ", format(wt.count, big.mark = ","), ")")
na.title <- paste0("NA (n = ", format(na.or.het.count, big.mark = ","), ")")
md.wt <- md[md$Genotype %in% c("WT"), ]
md.mut <- md[md$Genotype %in% c("MUT"), ]
md.na <- md[md$Genotype %in% c("NA", "Control"), ]

p.genotyping <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, fill = Genotype)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.mut, shape = 21, size = 1, fill = vexas.genotyping.palette[["MUT"]], stroke = 0.1, color = "black") +
  geom_point(data = md.wt, shape = 21, size = 1, fill = vexas.genotyping.palette[["WT"]], stroke = 0.1, color = "black") +
  # geom_point(data = md.genotyped, shape = 21, size = 1) +
  # scale_color_manual(values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none") + coord_fixed()
p.genotyping
ggsave("figures/current_figure_drafts/UMAP_RNA_genotyping.tiff", plot = p.genotyping, device = "tiff", dpi = 300, width = 5, height = 4, unit = "in")


## Genotyping UMAP, split by patient
per.patient.plots <- lapply(levels(rna.obj$Donor), function(x) {
  print(x)
  md.mut.donor <- md.mut %>% filter(Donor == x)
  md.wt.donor <- md.wt %>% filter(Donor == x)
  fraction.gentotyped <- md %>% filter(Donor == x) %>% filter(Genotype %in% c("MUT", "WT")) %>% nrow() /  md %>% filter(Donor == x) %>% nrow()
  percent.genotyped <- round(fraction.gentotyped, digits = 4) *100
  p.donor.genotyped <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, color = Genotype)) +
    geom_point(shape = 19, size = 1.5, color = "black") + 
    geom_point(shape = 19, size = 1, color = "grey80")+
    geom_point(data = md.mut.donor, shape = 19, size = 0.5, color = vexas.genotyping.palette[["MUT"]]  ) +
    geom_point(data = md.wt.donor, shape = 19, size = 0.5, color = vexas.genotyping.palette[["WT"]]) +
    scale_color_manual(values = unlist(vexas.genotyping.palette)) +
    theme_classic() + coord_fixed() +
    theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
    guides(x = axis, y = axis) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    ggtitle(x, subtitle = paste0("Cells genotyped = ", percent.genotyped, "%")) +
    theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + coord_fixed() 
  return(p.donor.genotyped)
}) 
p.combined <- wrap_plots(per.patient.plots, nrow = 2)
ggsave("figures/current_figure_drafts/RNA_UMAP_genotyping_split.tiff", plot = p.combined, device = "tiff", dpi = 300, width = 12, height = 5.5)

## Genotyping UMAP (heatmaps)
umap.lims.x <- c(min(md$UMAP1), max(md$UMAP1)) 
umap.lims.y <- c(min(md$UMAP2), max(md$UMAP2)) 
p.genotyping.wt <- ggplot(md[which(md$Genotype =="WT"),], aes(x=UMAP1, y=UMAP2))+
  stat_density2d(aes(alpha=..level.., fill=..level..), bins=20, geom="polygon") +
  scale_fill_gradient2(low = "white", high = vexas.genotyping.palette[[1]]) +
  scale_alpha(range = c(0.00, 1), guide = FALSE) +
  guides(x = axis, y = axis) +
  coord_fixed(xlim = umap.lims.x, ylim = umap.lims.y) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none")
p.genotyping.wt
ggsave("figures/current_figure_drafts/UMAP_RNA_genotyping_heatmap_wt.pdf", plot = p.genotyping.wt, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")
p.genotyping.mut <- ggplot(md[which(md$Genotype =="MUT"),], aes(x=UMAP1, y=UMAP2))+
  stat_density2d(aes(alpha=..level.., fill=..level..), bins=20, geom="polygon") +
  scale_fill_gradient2(low = "white", high = vexas.genotyping.palette[[2]]) +
  scale_alpha(range = c(0.00, 1), guide = FALSE) +
  guides(x = axis, y = axis) +
  coord_fixed(xlim = umap.lims.x, ylim = umap.lims.y) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none")
p.genotyping.mut
ggsave("figures/current_figure_drafts/UMAP_RNA_genotyping_heatmap_mut.pdf", plot = p.genotyping.mut, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")
p.genotyping.na <- ggplot(md, aes(x=UMAP1, y=UMAP2))+
  geom_point(color = "grey80", size = 1.5) +
  coord_fixed(xlim = umap.lims.x, ylim = umap.lims.y) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))
p.genotyping.na
ggsave(paste0(current.plots.path, "/UMAP_RNA_points_heatmap_na.pdf"), plot = p.genotyping.na, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")



###################################### Plots - MUT cell frequencies ################################################

# median.pseudotime.values <- md %>% 
#   group_by(CellType) %>% 
#   summarize(median_pseudotime = median(monocle3_pseudotime))

overall_mut_fraction <- md %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count(Genotype) %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(overall_geno_frac = MUT / (MUT + WT)) %>% 
  pull(overall_geno_frac)

## Summarize the MUT cell frequency per donor
mut.fraction.per.donor <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  mutate(CellType = gsub("CD4 |CD8 ", "", CellType))  %>%
  filter(n() > 100) %>% ## Only including cell types with > this number genotyped cells
  group_by(CellType, Donor) %>% 
  filter(n() > 10) %>%  ## Only including patient data points with > 10 genotyped cells
  group_by(CellType, Genotype, Donor) %>%
  summarize(n = n()) %>% 
  group_by(CellType, Donor) %>% 
  # mutate(mut_fraction = n / sum(n), overall_mut_fraction = overall_count_genotype/overall_count, normalized_mut_fraction = mut_fraction / overall_mut_fraction) %>% View
  mutate(mut_fraction = n / sum(n), normalized_mut_fraction = mut_fraction / overall_mut_fraction) %>%
  filter(Genotype == "MUT") %>% 
  group_by(CellType) %>% 
  filter(n_distinct(Donor) >= 3)

summarized.freqs <- mut.fraction.per.donor %>% 
  summarize(mean_val = mean(normalized_mut_fraction), sd_val = sd(normalized_mut_fraction), n = n(), se_val = sd_val / sqrt(n), median_val = median(normalized_mut_fraction)) %>% 
  merge(., median.pseudotime.values, group.by = CellType, all.x = T)

## Make the bar plot
p.mutant_cell_fraction.bar.plot <- summarized.freqs %>% 
  ggplot(aes(x = reorder(CellType, -mean_val), y = mean_val, fill = CellType)) +
  geom_bar(stat = "identity", color = "black", size = 0.25) +
  geom_errorbar(aes(x = CellType, ymin = mean_val-se_val, ymax = mean_val+se_val), width=0.4, size = 0.25) +
  geom_point(data = mut.fraction.per.donor, mapping = aes(x = CellType, y = normalized_mut_fraction, shape = Donor), position = position_jitter(width = 0.25), size = 1) +
  scale_shape_manual(values = seq(1, 10, 1)) +
  scale_fill_manual(values = cell.type.palette, guide = "none") +
  # scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept = 1, linetype = "longdash", color = "black") +
  coord_fixed(ylim = c(0, 1.5), ratio = 7) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5)) +
  labs(y = "Normalized mutant-\ncell frequency", x = element_blank())
p.mutant_cell_fraction.bar.plot
ggsave(paste0(current.plots.path, "/normalized_MUT_cell_freq_bar_plot.pdf"), plot = p.mutant_cell_fraction.bar.plot, device = "pdf", dpi = 300, width = 6, height = 4, unit = "in")

## Save as a csv
summarized.freqs %>% 
  write_csv("data/metadata/mut_cell_fractions_RNA.csv")

###################################### Separating T and NK cells ################################################

## RNA

## Based on specific genes
rna.obj.t.nk <- subset(rna.obj, CellType %in% c("CD4 T", "CD8 T", "NK"))
m.t.nk <- rna.obj.t.nk@assays$RNA@data[c("CD3E", "KLRF1"), ] %>% as.matrix() %>% t()
rna.obj.t.nk <- AddMetaData(rna.obj.t.nk, m.t.nk)
making_violin_plot_for_RNA <-  function(x) {
  quantiles = c(0.05, 0.95)
  quantiles.current.mut <- rna.obj.t.nk@meta.data %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "MUT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  quantiles.current.wt <- rna.obj.t.nk@meta.data %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "WT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  p.boxplot <- rna.obj.t.nk@meta.data %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    ggplot(aes(x = Genotype, y = !! sym(x), fill = Genotype)) +
    geom_violin() +
    geom_point(size = 0.1, alpha = 0.5, position = position_jitter(width = 0.25)) +
    # coord_cartesian(ylim = c(min(quantiles.current.mut[[1]], quantiles.current.wt[[1]]), max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]))) +
    scale_y_continuous(expand = c(0.3, 0)) +
    theme_classic() +
    ylab("RNA expression") +
    theme(axis.title.x = element_blank()) +
    # stat_compare_means(label = "p.format", label.y =max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]) * 1.1, size = 3) +
    stat_compare_means(label = "p.format", size = 3) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#2A9098", "#F5B935")) +
    ggtitle(x)
  return(p.boxplot)
}
p.cd3E <- making_violin_plot_for_RNA("CD3E")
p.klrf1 <- making_violin_plot_for_RNA("KLRF1")
combined <- p.cd3E | p.klrf1
combined

## Based on azimuth labels
p.azimuth.labels <- md %>% 
  filter(CellType %in% c("CD4 T", "CD8 T", "NK")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  mutate(Category = if_else(CellType == "NK", "NK", "CD4/CD8 T")) %>% 
  group_by(Donor, Category) %>% 
  filter(n() > 10) %>% ## Only including patient data points with > 10 genotyped cells
  group_by(Donor, Category, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  summarize(mut_cell_fraction = MUT / (MUT + WT)) %>% 
  ggplot(aes(x = Category, y = mut_cell_fraction, fill = Category)) +
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FFCB88", "#FF7B0B")) +
  stat_compare_means(label = "p.format", size = 3)  +
  theme(legend.position = "none") + 
  ylab("Mutant cell fraction") +
  ggtitle("Azimuth labels")
p.azimuth.labels
ggsave("figures/current_figure_drafts/RNA_T_NK_separation_azimuth_labels.pdf", p.azimuth.labels, device = "pdf", dpi = 300, width = 2, height = 3, unit = "in")

## Based on azimuth labels
p.azimuth.labels <- md %>% 
  filter(CellType %in% c("CD4 T", "CD8 T", "NK")) %>% 
  filter(predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "NK")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  mutate(Category = if_else(predicted.celltype.l1 == "NK", "NK", "CD4/CD8 T")) %>% 
  group_by(Donor, Category) %>% 
  filter(n() > 10) %>% ## Only including patient data points with > 10 genotyped cells
  group_by(Donor, Category, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  summarize(mut_cell_fraction = MUT / (MUT + WT)) %>% 
  ggplot(aes(x = Category, y = mut_cell_fraction, fill = Category)) +
  geom_boxplot() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#FFCB88", "#FF7B0B")) +
  stat_compare_means(label = "p.format", size = 3)  +
  theme(legend.position = "none") + 
  ylab("Mutant cell fraction") +
  # facet_grid(cols = vars(Donor)) +
  ggtitle("Azimuth labels")
p.azimuth.labels
ggsave("figures/current_figure_drafts/RNA_T_NK_separation_azimuth_labels_actual_labels_only.pdf", p.azimuth.labels, device = "pdf", dpi = 300, width = 2, height = 3, unit = "in")

## ADT

## Isolate the cells with ADT information in the T/NK cell cluster
adt.counts.per.bc <- colSums(rna.obj@assays$ADT_DSB@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.adt <- subset(rna.obj, cells = bcs.keep)
rna.obj.adt <- ScaleData(rna.obj.adt, assay = "ADT_DSB")
rna.obj.adt.tcell.nk <- subset(rna.obj.adt, CellType %in% c("CD4 T", "CD8 T", "NK"))
adt.matrix <- rna.obj.adt@assays$ADT_DSB@data %>% t()
md.adt <- merge(md, adt.matrix, all.x = FALSE, by = "row.names")
md.no.adt <- md[!(rownames(md) %in% md.adt), ]
quantiles = c(0.01, 0.99)

## Make a helper function to highlight ADT expression on the UMAP
plot_adt_on_umap <- function(x) {
  quantiles.current <- md.adt %>% dplyr::rename(current_tag = x) %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  p.adt.highlight <- md.no.adt %>% 
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    geom_point(color = "grey90", size = 0.5) +
    geom_point(data = md.adt, aes(color = !! sym(x)), size = 0.5) +
    scale_color_viridis(option = "magma", limits = c(quantiles.current[[1]], quantiles.current[[2]]), oob = scales::squish, ) +
    guides(x = axis, y = axis) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
  return(p.adt.highlight)
}

## Plot CD3, CD56 UMAPs
plist <- lapply(c("Hu.CD3-UCHT1", "Hu.CD56"), plot_adt_on_umap)
p.combined <- wrap_plots(plist)
p.combined
ggsave(paste0(current.plots.path, "/RNA_CD3_CD56_on_umap.pdf"), plot = p.combined, device = "pdf", dpi = 300, width = 10, height = 6, unit = "in")

## Make boxplots of CD3, CD56 expression
geno.df <- rna.obj.adt.tcell.nk@meta.data %>% select(Genotype, Sample) %>% filter(Genotype %in% c("MUT", "WT"))
m <- rna.obj.adt.tcell.nk@assays$ADT_DSB@data %>% t() %>% as.data.frame()
m.genotyped <- merge(m, geno.df, by = "row.names", all = F)

making_boxplot_for_adt <-  function(x) {
  quantiles = c(0.05, 0.95)
  quantiles.current.mut <- m.genotyped %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "MUT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  quantiles.current.wt <- m.genotyped %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "WT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  p.boxplot <- m.genotyped %>% 
    ggplot(aes(x = Genotype, y = !! sym(x), fill = Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = c(min(quantiles.current.mut[[1]], quantiles.current.wt[[1]]), max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]))) +
    scale_y_continuous(expand = c(0.3, 0)) +
    theme_classic() +
    ylab("ADT expression") +
    theme(axis.title.x = element_blank()) +
    stat_compare_means(label = "p.format", label.y =max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]) * 1.1, size = 3) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#2A9098", "#F5B935")) +
    ggtitle(gsub("Hu.", "", x) %>% gsub("-.*", "", .))
  return(p.boxplot)
}
p.cd3 <- making_boxplot_for_adt("Hu.CD3-UCHT1")
p.cd56 <- making_boxplot_for_adt("Hu.CD56")
p.combined <- p.cd3 | p.cd56
ggsave(paste0(current.plots.path, "/RNA_ADT_CD3_CD56_boxplot_95th_quantile.pdf"), plot = p.combined, device = "pdf", dpi = 300, width = 3, height = 3, unit = "in")


## Make a scatter plot of CD3, CD56 expression
geno.df <- rna.obj.adt.tcell.nk@meta.data %>% select(Genotype, Sample)
rna.obj.adt.tcell.nk <- ScaleData(rna.obj.adt.tcell.nk, assay = "ADT_DSB")
m <- rna.obj.adt.tcell.nk@assays$ADT_DSB@scale.data %>% t() %>% as.data.frame()
m.genotyped <- merge(m, geno.df, by = "row.names", all = F)
min.cd3 <- m.genotyped %>% filter(Genotype %in% c("MUT", "WT")) %>% pull(`Hu.CD3-UCHT1`)  %>% min()
adt.plot <- m.genotyped %>% 
  filter(Genotype %in% c("MUT", "WT")) %>%
  ggplot(aes(x = `Hu.CD3-UCHT1`, y = `Hu.CD56`, fill = Genotype)) +
  geom_point(pch=21, size = 2) +
  coord_fixed(xlim = c(min.cd3, 2)) +
  # facet_grid(cols = vars(Sample), scales = "free") +
  scale_fill_manual(values = c(MUT = vexas.genotyping.palette[[2]], WT = "grey")) + ggtitle("ADT expression")
adt.plot
ggsave(paste0(current.plots.path, "/ADT_nk_t_cell_separation.pdf"), plot = adt.plot, device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")


##################################### Plots - PERK UPR gene expression ##############################################

## Plotting by expression
rna.obj.genotyped <- subset(rna.obj, Genotype %in% c("MUT", "WT"))
Idents(rna.obj.genotyped) <- "CellType"

celltypes.plot <- c("HSC", "LMPP", "EMP", "MkP")
genes.plot <- list(
  PERK = c("DDIT3", "DDIT4"),
  IRE1 = c("CALR", "XBP1"),
  Inflammation = c("CAT", "HLA-DRA")
  # ATF6 = c("ATF6", "ATF6B", "CREBZF")
  # interferon = c("IRF1", "CXCL3", "IL1B", "HLA-DRA"),
  # chaperones = c("HSPD1", "CCT3", "CCT7")
)
gene.list <- Reduce(union, genes.plot)
gene.df <- data.frame(gene = gene.list, category = NA)
for (category in names(genes.plot)) {
  gene.df$category[gene.df$gene %in% genes.plot[[category]]] <- category
}
m <- rna.obj.genotyped@assays$RNA@data[gene.list, ] %>% as.matrix() %>% t()
m <- merge(rna.obj.genotyped@meta.data %>% dplyr::select(Genotype, CellType, Donor), m, by = "row.names")
p.gene.violin.grid <- m %>% 
  filter(CellType %in% celltypes.plot) %>% 
  pivot_longer(cols = gene.list, names_to = "gene", values_to = "expression") %>% 
  merge(., gene.df, by = "gene") %>% 
  ggplot(aes(x = gene, y = expression, fill = Genotype)) +
  facet_grid(rows = vars(CellType), cols = vars(category), space = "free", scales = "free") +
  # geom_boxplot() +
  geom_violin() +
  geom_point(size = 0.05, position=position_jitterdodge(dodge.width=0.9), alpha = 1) +
  scale_fill_manual(values = vexas.genotyping.palette) +
  stat_compare_means(label = "p.signif", hide.ns = T) +
  coord_cartesian(ylim=c(0, 6)) +
  scale_y_continuous(expand = c(0.2, 0)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(element_blank()) +
  ylab("Normalized expression")
p.gene.violin.grid
ggsave(paste0(current.plots.path, "/RNA_important_gene_violin_grid.pdf"), plot = p.gene.violin.grid, dpi = 300, device = "pdf", width = 5, height = 4)


##################################### Plots - Dotplot of important genes ################################################

cell.type.list <- c("HSC", "LMPP", "EMP", "MkP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT/", x, "_no_mito_no_ribo_renormalized.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})
de.results.combined <- do.call(rbind, de.results)
gene.list <- c("XBP1", "CALR", "HLA-DRA", "CAT", "ATF4", "DDIT4", "DDIT3", "HLA-DRB1", "PCDH9", "TPT1", "EEF1B2")
rna.obj.genotyped.hsc <- subset(rna.obj, cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))

gene.dotplot <- de.results.combined %>% 
  group_by(feature) %>% 
  filter(feature %in% gene.list) %>% 
  mutate(mean_val = mean(avg_log2FC)) %>% 
  ggplot(aes(x = factor(cluster, levels = cell.type.list), y = reorder(feature, mean_val), size = abs(avg_log2FC), color = -log10(fdr)*sign(avg_log2FC))) +
  scale_color_gradient2(low = vexas.genotyping.palette[["WT"]], high = vexas.genotyping.palette[["MUT"]], mid = "white", limits = c(-3, 3),  oob = scales::squish) +
  geom_point(stat = "identity", ) +
  # facet_grid(rows = vars(category), scales = "free") +
  theme_classic() +
  # theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("") + ylab("")
gene.dotplot
ggsave(paste0(current.plots.path, "/RNA_important_genes_dotplot.pdf"), plot = gene.dotplot, dpi = 300, device = "pdf", width = 4.25, height = 4)


###################################### Plots - Violins of important genes ################################################

## Helper violin plot split by donor
rna.obj.genotyped.hsc <- subset(rna.obj, cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))
VlnPlot(rna.obj.genotyped.hsc, features = c("HLA-DRB1"), 
        group.by = "Donor", 
        split.by = "Genotype", 
        assay = "RNA",
        slot = "data", 
        cols = vexas.genotyping.palette) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_y_continuous(expand = c(0.1, 0)) +
  scale_fill_manual(values = vexas.genotyping.palette) + xlab("") +
  stat_summary(fun.y = median, geom='point', colour = "black", shape = 95)

## Violin plots
rna.obj.genotyped <- subset(rna.obj, Genotype %in% c("MUT", "WT"))
Idents(rna.obj.genotyped) <- "CellType"
rna.obj.genotyped$CellType <- factor(rna.obj.genotyped$CellType, levels = c("HSC", "LMPP", "EMP", "MkP"))
plist <- lapply(c("XBP1", "CALR", "HLA-DRA", "CAT"), function(x) {
  VlnPlot(rna.obj.genotyped, features = c(x), 
          group.by = "CellType", 
          split.by = "Genotype", 
          assay = "RNA",
          slot = "data",
          idents = c("HSC", "LMPP", "EMP", "MkP"), 
          pt.size = 0.1,
          alpha = 0.5,
          cols = vexas.genotyping.palette) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_y_continuous(expand = c(0.1, 0)) +
    scale_fill_manual(values = vexas.genotyping.palette) + xlab("") +
    stat_summary(fun.y = median, geom='point', colour = "black", shape = 95)
}) 
p.important.genes <- wrap_plots(plist, ncol = 2) + plot_layout(guides = "collect")
p.important.genes
ggsave(paste0(current.plots.path, "/RNA_important_gene_violins.pdf"), plot = p.important.genes, dpi = 300, device = "pdf", width = 10, height = 7)

## Violin plots - just HSCs
rna.obj.genotyped <- subset(rna.obj, Genotype %in% c("MUT", "WT"))
Idents(rna.obj.genotyped) <- "CellType"
rna.obj.genotyped$CellType <- factor(rna.obj.genotyped$CellType, levels = c("HSC"))
plist <- lapply(c("DDIT3", "DDIT4", "IL1B"), function(x) {
  DefaultAssay(rna.obj.genotyped) <- "RNA"
  md <- rna.obj.genotyped@meta.data
  gene.expression <- FetchData(object = rna.obj.genotyped, vars = c(x), slot = "data", assay = "RNA")
  md <- merge(md, gene.expression, by = "row.names") %>% column_to_rownames("Row.names")
  md %>% 
    dplyr::rename(gene = x) %>% 
    filter(CellType %in% c("HSC")) %>% 
    mutate(CellType = factor(CellType, levels = c("HSC"))) %>% 
    group_by(Genotype) %>% 
    # sample_n(size = 500, replace = T) %>% 
    ggplot(., aes_string(x = "Genotype", y = "gene", fill = "Genotype")) +
    geom_violin() +
    # geom_boxplot() +
    geom_jitter(width = 0.1, size = 0.1) +
    scale_fill_manual(values = vexas.genotyping.palette) +
    facet_grid(rows = vars(Donor)) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_y_continuous(expand = c(0.1, 0)) +
    theme(legend.position = "none") +
    ggtitle(x) +
    xlab("") +
    ylab("Expression Level")
}) 
p.important.genes <- wrap_plots(plist, ncol = 3) + plot_layout(guides = "collect")
p.important.genes
ggsave(paste0(current.plots.path, "/RNA_important_gene_violins_DDIT4_DDIT3_IL1B.pdf"), plot = p.important.genes, dpi = 300, device = "pdf", width = 7, height = 3)


## Violin plots, with medians
rna.obj.genotyped <- subset(rna.obj, Genotype %in% c("MUT", "WT"))
Idents(rna.obj.genotyped) <- "CellType"
plist <- lapply(c("HSPA8", "CAT", "DDIT4", "HLA-DRA"), function(x) {
  DefaultAssay(rna.obj.genotyped) <- "RNA"
  md <- rna.obj.genotyped@meta.data
  gene.expression <- FetchData(object = rna.obj.genotyped, vars = c(x), slot = "data", assay = "RNA")
  md <- merge(md, gene.expression, by = "row.names") %>% column_to_rownames("Row.names")
  md %>% 
    dplyr::rename(gene = x) %>% 
    filter(CellType %in% c("HSC", "LMPP", "EMP", "Prog Mk")) %>% 
    mutate(CellType = factor(CellType, levels = c("HSC", "LMPP", "EMP", "Prog Mk"))) %>% 
    ggplot(., aes_string(x = "Genotype", y = "gene", fill = "Genotype")) +
    geom_violin(draw_quantiles = 0.5) +
    geom_jitter(width = 0.1, size = 0.1) +
    scale_fill_manual(values = c(WT = "lightsteelblue1", MUT = "lightsteelblue3")) +
    facet_grid(cols = vars(CellType)) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_y_continuous(expand = c(0.1, 0)) +
    theme(legend.position = "none") +
    ggtitle(x) +
    xlab("") +
    ylab("Expression Level")
}) 
p.important.genes <- wrap_plots(plist, ncol = 2) + plot_layout(guides = "collect")
ggsave(paste0(current.plots.path, "/RNA_important_gene_violins_medians.pdf"), plot = p.important.genes, dpi = 300, device = "pdf", width = 10, height = 7)


############################################## XBP1s - sample-aware permutation test (unfiltered, HSPC) ###################################

## Save a table of the spliced/unspliced XBP1 UMI counts per donor, cell type, and genotype
rna.obj@meta.data %>% filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  group_by(CellType, Donor, Genotype) %>% summarize(cell_count = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  filter(Genotype != "NA") %>% 
  write_csv("data/xbp1_splicing/xbp1_read_counts_kallisto_output.csv")

donors.keep <- rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1)) %>% 
  group_by(Donor) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  filter(frac_genotyped > 0) %>% 
  pull(Donor)

md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the mean between the two ratios
calculate_difference_of_means <- function(df) {
  df %>% 
    group_by(Genotype) %>% 
    summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
    mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
    select(Genotype, splicing_ratio) %>% 
    pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
    mutate(mean_diff = MUT - WT) %>% 
    pull(mean_diff)
}

## Shuffle the genotype labels and calculate within patients
md.xbp1.shuf <- md.xbp1
mean_diff_results <- sapply(1:10000, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_means(.)
})

## Calculate the actual mean difference
observed.mean.diff <- calculate_difference_of_means(md.xbp1)

## Calculate the fraction of simulated results that are above the observed value
p.value <- sum(mean_diff_results > observed.mean.diff) / length(mean_diff_results)
print(p.value)

## Plot histogram of all iterations
as.data.frame(mean_diff_results) %>% 
  ggplot(aes(x = mean_diff_results, y = ..count..)) +
  geom_histogram() +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = "Sample-aware permutation test, n = 10,000") +
  annotate("text", x = -0.02, y = 1000, label = paste0("Observed\ndifference = ", round(observed.mean.diff, digits = 4), "\np < 0.0001"), hjust = 0, size = 3) +
  # annotate("text", x = 0, y = 800, label = paste0("p = ", p.value)) +
  # annotate("text", x = 0, y = 800, label = "p < 0.0001") +
  geom_vline(xintercept = observed.mean.diff, linetype = "longdash", color = "darkred") + 
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_10k_iterations_kallisto_output.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

############################################## XBP1s - sample-aware permutation test (unfiltered, HSC) ###################################

donors.keep <- rna.obj@meta.data %>% 
  filter(cluster_celltype == "HSC") %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n > 10) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor)
  
md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the mean between the two ratios
calculate_difference_of_means <- function(df) {
  df %>% 
    group_by(Genotype) %>% 
    summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
    mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
    select(Genotype, splicing_ratio) %>% 
    pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
    mutate(mean_diff = MUT - WT) %>% 
    pull(mean_diff)
}

## Shuffle the genotype labels and calculate within patients
md.xbp1.shuf <- md.xbp1
n = 1000
mean_diff_results <- sapply(1:n, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    group_by(Donor, Genotype) %>% 
    sample_n(size = 25) %>% 
    calculate_difference_of_means(.)
})

## Calculate the actual mean difference
observed.mean.diff <- calculate_difference_of_means(md.xbp1)

## Calculate the fraction of simulated results that are above the observed value
p.value <- sum(mean_diff_results > observed.mean.diff) / length(mean_diff_results)
print(p.value)

## Plot histogram of all iterations
as.data.frame(mean_diff_results) %>% 
  ggplot(aes(x = mean_diff_results, y = ..count..)) +
  geom_histogram() +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("Sample-aware permutation test, n = ", n)) +
  annotate("text", x = -0.02, y = 100, label = paste0("Observed\ndifference = ", round(observed.mean.diff, digits = 4)), hjust = 0, size = 3) +
  annotate("text", x = 0, y = 80, label = paste0("p = ", p.value)) +
  # annotate("text", x = 0, y = 800, label = "p < 0.0001") +
  geom_vline(xintercept = observed.mean.diff, linetype = "longdash", color = "darkred") + 
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
# ggsave("figures/current_figure_drafts/sample_aware_permutation_test_10k_iterations_kallisto_output.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)


############################################### XBP1s - sample-aware permutation test boxplots (unfiltered) ##############################################################

## Make a bar plot with XBP1s values
donors.keep <- rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1)) %>% 
  group_by(Donor) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  filter(frac_genotyped > 0) %>% 
  pull(Donor)

## Pseudobulk by donor
p.xbp1.bars <- md %>% 
  filter(Donor %in% donors.keep) %>% 
  # filter(cluster_celltype == "HSC") %>%
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>%
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  summarize(mean_splicing_ratio = mean(splicing_ratio),
            sd = sd(splicing_ratio),
            n = n(),
            se = sd / sqrt(n)) %>% 
  ggplot(aes(x = Genotype, y = mean_splicing_ratio, fill = Genotype)) +
  geom_col(color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_splicing_ratio-se, ymax = mean_splicing_ratio+se), width=0.4, size = 0.25) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = "HSPCs") +
  # annotate("text", x = 1.2, y = 0.4, label=paste0("p = ", p.value), size = 3)
  annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_column_plot_kallisto_output.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)


################################################ Plots - Monocyte volcano, GSEA enrichment plots ###################################

cd14.mono.results <- read_csv("data/differential_expression/VEXAS_MUT_vs_WT/CD14 Mono_min_10_cells_20240505.csv") %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = "Mono")
cd14.mono.results <- read_csv("data/differential_expression/VEXAS_MUT_vs_WT/CD14 Mono_min_20_cells_varfeatures_20240505.csv") %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = "Mono")

## Plot CD14 monocyte plot
genes.highlight = list(
  MUT = cd14.mono.results %>% slice_min(n = 30, order_by = pval) %>% filter(avg_log2FC > 0) %>% pull(feature),
  WT = cd14.mono.results %>% slice_min(n = 30, order_by = pval) %>% filter(avg_log2FC < 0) %>% pull(feature)
  # pathway = pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB
)
p.mono.volcano <- plot_volcano(cd14.mono.results, stat_column = "pval", stat_line = 0.05, effect_line = 0.2, max.overlaps = 30,
                               only_genes = genes.highlight, only_genes_colors = vexas.genotyping.palette) + 
  coord_cartesian(xlim = c(-1.4, 1.4), ylim = c(0, 5)) +
  ggtitle("CD14 Monocytes, MUT vs. WT")
p.mono.volcano
ggsave("figures/current_figure_drafts/RNA_volcano_monocytes_exp_10_cells_all_features_20240305.pdf", plot = p.mono.volcano, device = "pdf", dpi = 300, width = 7, height = 5)

## Plot CD14 monocyte plot - no labels
p.mono.volcano <- plot_volcano(cd14.mono.results, genotyping_column = NA, stat_column = "pval", stat_line = 0.05, effect_line = 0.25, label_genes = F) + 
  coord_cartesian(xlim = c(-1.2, 1.2)) +
  ggtitle("CD14 Monocytes, MUT vs. WT")
p.mono.volcano
ggsave("figures/current_figure_drafts/RNA_volcano_monocytes_no_labels.pdf", plot = p.mono.volcano, device = "pdf", dpi = 300, width = 7, height = 5)

fgsea.out <- run_fgsea_for_list(cd14.mono.results)

## Try re-ranking
de_results.ranked <- cd14.mono.results %>%
  filter(!is.na(pval)) %>% 
  # mutate(rank = -log10(pval)*sign(avg_log2FC)) %>%
  mutate(rank = avg_log2FC) %>% 
  # mutate(rank = -log10(pval)) %>%
  arrange(-rank)
current.ranks.list <- de_results.ranked$rank
names(current.ranks.list) <- de_results.ranked$feature

fgsea.out <- run_fgsea_for_list(ranks.list = current.ranks.list)

fgsea.out %>% 
  write_csv("figures/current_figure_drafts/fgsea_results_for_figure_CD14_mono_10_cells_all_features.csv")

## Plot enrichment - troubleshooting
m_t2g = msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
plotEnrichment(pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB, current.ranks.list)

## Plot enrichment
p1 <- plot_enrichment(ranks.list = current.ranks.list, pathway = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 1)
p2 <- plot_enrichment(cd14.mono.results, ranks.list = current.ranks.list, pathway.name = "HALLMARK_G2M_CHECKPOINT", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 1)
p3 <- plot_enrichment(cd14.mono.results, ranks.list = current.ranks.list, pathway.name = "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 1)
p.combined <- p1 | p2
p.combined
ggsave("figures/current_figure_drafts/RNA_enrichment_monocytes_exp_10_cells_all_features_20240305.pdf", plot = p.combined, device = "pdf", dpi = 300, width = 8, height = 3.5)

## Plot gsea bar plot
p.mono.gsea <- fgsea.out %>% 
  # filter(padj < 0.25) %>% 
  # mutate(status = sign(NES)) %>% 
  # group_by(status) %>% 
  slice_min(n = 10, order_by = pval) %>% 
  mutate(pathway = gsub("_", " ", pathway) %>% str_to_title(.) %>% str_wrap(., width = 40)) %>% 
  dplyr::select(pathway, NES, padj) %>% 
  # column_to_rownames("pathway")
  mutate(`-log10(padj)*sign(NES)` = -log10(padj) * sign(NES)) %>% 
  ggplot(aes(x = reorder(pathway, NES), y = NES, fill = `-log10(padj)*sign(NES)`)) +
  geom_col() +
  xlab(element_blank()) +
  theme_classic() +
  scale_fill_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], limits=c(-3, 3), oob = scales::squish) +
  coord_flip()
p.mono.gsea
ggsave("figures/current_figure_drafts/RNA_gsea_bar_plot_NES_plotted_monocytes_exp_10_cells_all_features_20240305.pdf", plot = p.mono.gsea, device = "pdf", dpi = 300, width = 6, height = 2)

## Plot gsea bar plot by pval
p.mono.gsea <- fgsea.out %>% 
  # filter(padj < 0.25) %>% 
  mutate(status = if_else(sign(NES) > 0, "MUT", "WT")) %>%
  # group_by(status) %>% 
  slice_min(n = 10, order_by = pval) %>% 
  mutate(pathway = gsub("_", " ", pathway) %>% str_to_title(.) %>% str_wrap(., width = 40)) %>% 
  dplyr::select(pathway, NES, padj, status) %>% 
  # column_to_rownames("pathway")
  mutate(`-log10(padj)*sign(NES)` = -log10(padj) * sign(NES)) %>% 
  ggplot(aes(x = reorder(pathway, -log10(padj)*sign(NES)), y = `-log10(padj)*sign(NES)`, fill = status)) +
  scale_fill_manual(values = vexas.genotyping.palette) +
  geom_col() +
  xlab(element_blank()) +
  geom_hline(yintercept = -log10(0.1), linetype = "longdash") +
  geom_hline(yintercept = log10(0.1), linetype = "longdash") +
  theme_classic() +
  # scale_fill_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], limits=c(-1, 1), oob = scales::squish) +
  coord_flip()
p.mono.gsea
ggsave("figures/current_figure_drafts/RNA_gsea_bar_plot_pval_plotted_monocytes_exp_10_cells_all_features_20240305.pdf", plot = p.mono.gsea, device = "pdf", dpi = 300, width = 6, height = 2)


############################################### Plots - HSC volcano, GSEA enrichment plots ##########################################

de.results <- read_csv("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>%
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(cluster = "HSC") %>% 
  mutate(fdr = if_else(fdr < 1e-30, 1e-30, fdr))
  

## Annotated volcanoes
fgsea.out <- run_fgsea_for_list(de.results)

## Try re-ranking
de_results.ranked <- de.results %>%
  filter(!is.na(pval)) %>% 
  # mutate(rank = -log10(pval)*sign(avg_log2FC)) %>%
  mutate(rank = avg_log2FC) %>% 
  # mutate(rank = -log10(pval)) %>%
  arrange(-rank)
current.ranks.list <- de_results.ranked$rank
names(current.ranks.list) <- de_results.ranked$feature

fgsea.out <- run_fgsea_for_list(de.results, ranks.list = current.ranks.list)
fgsea.out2 <- run_fgsea_for_list(de.results, genes.tested = de.results$feature)

## Annotating genes based on Sars's annotations
volcano.gene.list <- list(
  # `all` = c(de.results %>% filter(abs(avg_log2FC) > 0.4 | -log10(fdr) > 10) %>% pull(feature))
  # `all` = c("EEF1B2", "EEF1A1"),
  Translation = c("EEF1B2", "EEF1A1", "EEF1G", "EEF2", "EIF3L", "EIF4B", "EIF3F"),
  Inflammation = c("IFITM1", "HLA-DPB1", "HLA-DRA", "HLA-A", "IFITM2", "ISG20", "IRF1", "CD47", "IL1B", "CXCL3", "CXCL2", "CXCL8"),
  UPR = c("CD99", "DDIT4", "CALCRL", "XBP1", "CALR", "ATF4", "DDIT3"),
  # `UPR` = c("DDIT4", "ATF4", "DDIT3", "CALCRL", "CALR", "XBP1", "CD99", "EEF1A1", "EIF3L", "EIF3F", "EEF2"),
  # `Inflammation` = c("IFITM1", "IFITM2", "ISG20", "IRF1", "CD47", "IL1B", "HLA-DRA", "HLA-DPB1", "HLA-A"),
  `Myeloid` = c("CD38", "MYC", "CAT", "CDK6")
)

## Annotating genes based on overlap with gene modules
fgsea.out <- run_fgsea_for_list(de.results)
volcano.gene.list <- list(
  `UPR` = intersect(fgsea.out %>% filter(pathway == "ATF4_ONLY_HAN") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
  `Inflammation` = union(
    intersect(fgsea.out %>% filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
    intersect(fgsea.out %>% filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature))
    ),
  `Translation` = intersect(fgsea.out %>% filter(pathway == "REACTOME_TRANSLATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
  `Myeloid` = c("CD99", intersect(fgsea.out %>% filter(pathway == "GOBP_LEUKOCYTE_DIFFERENTIATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)))
)

## Pick colors + plot volcano
volcano.gene.list.colors <- c(
  all = "darkgreen",
  `Myeloid` = "blue4",
  `UPR` = "darkgreen", 
  `Inflammation` = "brown4",
  `Translation` = "orange"
)

## Plot volcano with labels
p.volcano <- plot_volcano(de.results, metadata.df = rna.obj@meta.data %>% filter(CellType == "HSC") %>% filter(Genotype %in% c("MUT", "WT")) %>% filter(!(Donor %in% c("PT13", "PT14", "PT25"))), 
                          effect_cutoff = 0.2, effect_line = 0.2, stat_line = 0.05, stat_column = "fdr", title = paste0("HSCs"), 
                          genotyping_column = "Genotype", 
                          only_genes = volcano.gene.list, only_genes_colors = volcano.gene.list.colors) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)") +
  coord_cartesian(xlim = c(-1.5, 1.5))
  # coord_cartesian(ylim=c(0, 30), xlim = c(-1.5, 1.5))
p.volcano
ggsave(paste0(current.plots.path, "/RNA_HSC_diff_exp_20240103.pdf"), plot = p.volcano, dpi = 300, device = "pdf", width = 6.5, height = 5)

## Plot volcano without labels
p.volcano <- plot_volcano(de.results, metadata.df = rna.obj@meta.data %>% filter(CellType == "HSC") %>% filter(Genotype %in% c("MUT", "WT")) %>% filter(!(Donor %in% c("PT13", "PT14", "PT25"))), 
                          effect_cutoff = 0.2, effect_line = 0.2, stat_line = 0.05, stat_column = "fdr", title = paste0("HSCs"), 
                          genotyping_column = "Genotype", 
                          label_genes = F) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)") +
  coord_cartesian(ylim=c(0, 30), xlim = c(-1.5, 1.5))
p.volcano
ggsave(paste0(current.plots.path, "/RNA_HSC_diff_exp_no_labels_20240304.pdf"), plot = p.volcano, dpi = 300, device = "pdf", width = 6.5, height = 5)

## Plot volcano highlighting HLAs
p.volcano <- plot_volcano(de.results, effect_cutoff = 0.05, effect_line = 0.2, stat_line = 0.05, stat_column = "fdr", title = paste0("HSCs"), 
                          only_genes = c(HLA = fgsea.out %>% filter(pathway == "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION") %>% pull(leadingEdge)),
                          only_genes_colors = c("HLA" = "darkred")) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)") +
  coord_cartesian(ylim=c(0, 30), xlim = c(-1.5, 1.5))
p.volcano


## Enrichment plots
p1 <- plot_enrichment(ranks.list = current.ranks.list, pathway = "ATF4 targets (Han et al.)", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 0.5)
p2 <- plot_enrichment(ranks.list = current.ranks.list, pathway = "HALLMARK_INTERFERON_ALPHA_RESPONSE", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 0.5)
p3 <- plot_enrichment(ranks.list = current.ranks.list, pathway = "REACTOME_TRANSLATION", pval_show = "padj") + theme_classic()  + geom_line(color = "#BA242A", size = 0.5)
p4 <- plot_enrichment(ranks.list = current.ranks.list, pathway = "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 0.5)
p.combined <- (p2 | p1 ) / (p4 + p3)
p.combined
ggsave("figures/current_figure_drafts/HSC_RNA_volcano_gsea_plots_linewidth_05.pdf", width = 6.5, height = 5)

plot_enrichment(ranks.list = current.ranks.list, pathway = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 0.5)

###################################### Plots - ADT differential expression ################################################

## NOTE: Only patients PT02, PT12, PT15 were retained (others excluded for low cell counts)
adt.hsc.diff.exp <- read_csv("data/differential_expression/VEXAS_ADT/HSCs_MUT_vs_WT_ADT_20231205.csv")

volcano.gene.list <- list(
  `MUT` = c("CD62L", "CD99", "CD38-HIT2", "HLA-DR", "Hu.CD13"),
  `WT` = c("CD49b", "CD49f", "CD71", "CD62p", "CD107a")
)

## Pick colors + plot volcano
volcano.gene.list <- list(
  Highlight = intersect(c("Hu.CD49b", "Hu.CD62L", "Hu.CD99", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD18", 
                "Hu.CD54", "Hu.HLA.DR", "Hu.CD244", "Hu.CD62P", "Hu.CD119", "Hu.CD11a", 
                "Hu.CD123", "Hu.CD26", "Hu.CD33", "Hu.CD13", "Hu.CD95", "Hu.CD267", "Hu.CD137", 
                "Hu.CLEC12A", "Hu.CD49d", "Hu.CD69", "Hu.CD161", "Hu.CD152", "Hu.CD158e1", 
                "Hu.CD272", "Hu.KLRG1", "Hu.CD25", "Hu.CD29", "Hu.HLA.ABC", "HuMs.CD49f", "Hu.CD107a"), adt.hsc.diff.exp %>% filter(abs(delta) > 1 & fdr < 0.1) %>% pull(feature))
)
volcano.gene.list <- list(
  Highlight = adt.hsc.diff.exp %>% filter(abs(delta) > 1 & fdr < 0.1) %>% pull(feature)
)
volcano.gene.list.colors <- list(Highlight = "darkgreen")

# ## Use this for CLR
# adt.hsc.diff.exp %>%
#   mutate(avg_log2FC = log2(fc)) %>%
#   mutate(feature = gsub("Hu[A-Za-z]*\\.", "", feature)) %>%
#   plot_volcano(unt.results, stat_column = "fdr", effect_column = "avg_log2FC", effect_line = 0.2, title = "HSCs", subtitle_addition = "Differential ADT expression") +
#   coord_cartesian(ylim = c(0, 25)) + theme(legend.position = "none")

## Use this for DSB
p.adt <- plot_volcano(adt.hsc.diff.exp, stat_column = "fdr", effect_column = "delta", effect_line = 0.5, stat_line = 0.05, 
                      title = "HSCs", only_genes = volcano.gene.list, only_genes_colors = volcano.gene.list.colors, max.overlaps = 35) + 
  coord_cartesian(ylim = c(0, 45), xlim = c(-7, 7)) + theme(legend.position = "none") +
  xlab("Delta z-score for ADT expression") +
  ylab("-log10(FDR)")
p.adt
ggsave("figures/current_figure_drafts/RNA_ADT_diff_exp.pdf", plot = p.adt, dpi = 300, device = "pdf", width = 8, height = 7)

## Plot volcano - no labels
p.adt <- plot_volcano(adt.hsc.diff.exp, stat_column = "fdr", effect_column = "delta", effect_line = 0.5, stat_line = 0.05, 
                      title = "HSCs", label_genes = F, max.overlaps = 25) + 
  coord_cartesian(ylim = c(0, 45), xlim = c(-7, 7)) + theme(legend.position = "none") +
  xlab("Delta z-score for ADT expression") +
  ylab("-log10(FDR)")
p.adt
ggsave("figures/current_figure_drafts/ADT_diff_exp_no_labels.pdf", plot = p.adt, dpi = 300, device = "pdf", width = 4, height = 3)


###################################### Plots - Individual ADT boxplots ################################################

## Isolate samples with ADT data
HSC.subset <- subset(rna.obj, subset = Donor %in% c("PT02", "PT12", "PT13", "PT14", "PT15") & cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))
m <- HSC.subset@assays$ADT_DSB@data %>% t()
HSC.subset <- AddMetaData(HSC.subset, m)
quantiles = c(0.01, 0.99)

making_boxplot_for_adt <-  function(x) {
  quantiles.current.mut <- HSC.subset@meta.data %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "MUT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  quantiles.current.wt <- HSC.subset@meta.data %>% dplyr::rename(current_tag = x) %>% filter(Genotype == "WT") %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  if (x == "Hu.CD33") {
    quantiles.current.mut = c(0, 6)
    quantiles.current.wt = c(0, 6)
  }
  if (x == "Hu.CLEC12A") {
    quantiles.current.mut = c(-2, 10)
    quantiles.current.wt = c(-2, 10)
  }
  p.boxplot <- HSC.subset@meta.data %>% 
    ggplot(aes(x = Genotype, y = !! sym(x), fill = Genotype)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = c(min(quantiles.current.mut[[1]], quantiles.current.wt[[1]]), max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]))) +
    scale_y_continuous(expand = c(0.3, 0)) +
    theme_classic() +
    ylab("ADT expression") +
    theme(axis.title.x = element_blank()) +
    stat_compare_means(label = "p.format", label.y =max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]) * 1.1, size = 3) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(WT = "#9FC0C6", MUT = "#00878D")) +
    ggtitle(gsub("Hu.", "", x) %>% gsub("-.*", "", .))
  return(p.boxplot)
}

## Making plots for many tags of interest
plist <- lapply(c("Hu.CD49b", "Hu.CD62L", "Hu.CD99", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD18", 
                  "Hu.CD54", "Hu.HLA.DR", "Hu.CD244", "Hu.CD62P", "Hu.CD119", "Hu.CD11a", 
                  "Hu.CD123", "Hu.CD26", "Hu.CD33", "Hu.CD13", "Hu.CD95", "Hu.CD267", "Hu.CD137", 
                  "Hu.CLEC12A", "Hu.CD49d", "Hu.CD69", "Hu.CD161", "Hu.CD152", "Hu.CD158e1", 
                  "Hu.CD272", "Hu.KLRG1", "Hu.CD25", "Hu.CD29", "Hu.HLA.ABC"), making_boxplot_for_adt)
p.combined <- wrap_plots(plist)
p.combined
ggsave(paste0(current.plots.path, "/RNA_HSC_ADT_boxplots.pdf"), plot = p.combined, dpi = 300, device = "pdf", width = 10, height = 10)

## Make boxplots for select tags
plist <- lapply(c("Hu.CD49b", "Hu.CD99", "Hu.CD38-HIT2", "Hu.HLA.DR", "Hu.CD13"), making_boxplot_for_adt)
p.boxplots <- wrap_plots(plist, nrow = 1) + plot_annotation(title = "HSC cluster", subtitle = "Surface ADT expression")
ggsave(paste0(current.plots.path, "/RNA_HSC_ADT_boxplots.pdf"), plot = p.boxplots, dpi = 300, device = "pdf", width = 6.7, height = 3)


####################### Plots - UBA1 Genotyping QC metrics #########################

## Pull non-normalized UBA1 expression
uba1.expression <- data.frame(UBA1_UMI_count = rna.obj@assays$RNA@counts["UBA1", ])
md.merged <- merge(md, uba1.expression, by = "row.names")

## Exclude any celltype-donor combinations with < 20 cells from analysis
md.merged <- md.merged %>% 
  group_by(CellType, Donor) %>% 
  filter(n() > 20)

# ## Order celltypes by UBA1 expression
# uba1.exp.celltype.order <- md.merged %>% 
#   group_by(CellType, Donor) %>% 
#   summarize(mean_uba1_umi_count = mean(UBA1_UMI_count)) %>% 
#   group_by(CellType) %>% 
#   summarize(mean_val = median(mean_uba1_umi_count)) %>% 
#   arrange(-mean_val) %>% 
#   pull(CellType)

## Order celltypes by genotyped fraction
uba1.exp.celltype.order <- md.merged %>% 
  group_by(CellType, Donor) %>% 
  summarize(fraction_genotyped = sum(Genotype != "NA") / n()) %>% 
  group_by(CellType) %>% 
  summarize(mean_frac_genotyped = mean(fraction_genotyped)) %>% 
  arrange(-mean_frac_genotyped) %>% 
  pull(CellType)

## Percent of cells genotyped
p.genotyping.rate <- md.merged %>% 
  group_by(CellType, Donor) %>% 
  mutate(CellType = factor(CellType, levels = uba1.exp.celltype.order)) %>% 
  summarize(fraction_genotyped = sum(Genotype != "NA") / n(), mut_count = sum(Genotype == "MUT"), wt_count = sum(Genotype == "WT"), total = n()) %>% 
  group_by(CellType) %>% 
  summarize(mean_val = mean(fraction_genotyped), sd_val = sd(fraction_genotyped), n = n(), se_val = sd_val / sqrt(n)) %>% 
  ggplot(aes(x = CellType, y = mean_val, fill = CellType)) + 
  geom_col(color = "black") +
  geom_errorbar(aes(x = CellType, ymin = mean_val-se_val, ymax = mean_val+se_val), width=0.4, size = 0.25) +
  scale_fill_manual(values = cell.type.palette) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent of\ncells genotyped") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "None")
p.genotyping.rate

## Plot the UBA1 expression per cluster (summarized by patient)
p.uba1.exp <- md.merged %>% 
  group_by(CellType, Donor) %>% 
  summarize(mean_uba1_umi_count = mean(UBA1_UMI_count)) %>% 
  group_by(CellType) %>% 
  mutate(CellType = factor(CellType, levels = uba1.exp.celltype.order)) %>% 
  ggplot(aes(x = CellType, y = mean_uba1_umi_count, fill = CellType)) + 
  # geom_violin(scale = "width") +
  # geom_point(position = "jitter") +
  geom_boxplot() +
  ylab("Average UBA1\nUMIs per cell") +
  scale_fill_manual(values = cell.type.palette) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None")

# ## Plot the UBA1 expression per cluster, no summary
# p.uba1.exp <- md.merged %>% 
#   ggplot(aes(x = CellType, y = UBA1_UMI_count, fill = CellType)) + 
#   geom_violin() +
#   # geom_point(position = "jitter") +
#   ylab("Average UBA1 UMI count per cell") +
#   scale_fill_manual(values = cell.type.palette) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None")

p.combined.uba1 <- p.genotyping.rate / p.uba1.exp + plot_layout(heights = c(2, 1)) + plot_annotation(title = "UBA1 expression vs. genotyping rate")
p.combined.uba1
ggsave(paste0(current.plots.path, "/UBA1_genotyping_per_cluster_metrics.pdf"), plot = p.combined.uba1, dpi = 300, device = "pdf", width = 5, height = 5)


########################################### Plots - QC metrics ############################################

## nCount_RNA
rna.obj@meta.data %>% 
  # filter(orig.ident %in% c("SG_VX_BM15a_GEX", "SG_VX_BM15b_GEX", "UPN14_CD34", "UPN15_CD34", "UPN16_CD34", "UPN17_CD34")) %>% 
  ggplot(aes(y = nCount_RNA, x = orig.ident, fill = orig.ident)) +
  geom_violin() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + coord_cartesian(ylim = c(0, 30000))
  
## nFeature_RNA
rna.obj@meta.data %>% 
  # filter(orig.ident %in% c("SG_VX_BM15a_GEX", "SG_VX_BM15b_GEX", "UPN14_CD34", "UPN15_CD34", "UPN16_CD34", "UPN17_CD34")) %>% 
  ggplot(aes(y = nFeature_RNA, x = orig.ident, fill = orig.ident)) +
  geom_violin() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

################################### Plotting cb-sniffer vs. GoT genotyping rates ##############################

## There may be differences in the percents used here and those used in the PPT made previously because
## the PPT was using all cellranger-identified cells and here we are using only cells in our final QC-filtered
## seurat object

## Mapping
NIH_mapping_list <- c(
  UPN14_BMMC = "GSM5858683",
  UPN15_BMMC = "GSM5858684",
  UPN16_BMMC = "GSM5858685",
  UPN17_BMMC = "GSM5858686",
  UPN14_CD34 = "GSM5858670",
  UPN15_CD34 = "GSM5858671",
  UPN16_CD34 = "GSM5858672",
  UPN17_CD34 = "GSM5858673"
)

## Load a list of every barcode with a cb_sniffer output
barcode.lists <- lapply(unique(rna.obj$orig.ident), function(x) {
  cb.sniffer.file <- x
  if (x %in% names(NIH_mapping_list)) {
    cb.sniffer.file <- NIH_mapping_list[[x]]
  }
  cb.sniff.barcodes <- read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/output/", cb.sniffer.file, "/", cb.sniffer.file, "_counts_CB.tsv")) %>% pull(barcode) %>% paste0(x, "_", .)
  return(cb.sniff.barcodes)
})
barcode.list <- do.call(c, barcode.lists)

## Add cb.sniffer output to metadata
md.cb.sniffer <- md %>% 
  rownames_to_column("barcode") %>% 
  mutate(has_cbsniffer_output = barcode %in% barcode.list) %>% 
  mutate(has_GoT_output = Genotype %in% c("MUT", "WT"))

## Look at genotyping differences
md.cb.sniffer %>% count(has_cbsniffer_output, has_GoT_output)
md.cb.sniffer.summarized.sample <- md.cb.sniffer %>% 
  group_by(orig.ident) %>% 
  summarize(fraction_genotyped_cbsniffer = sum(has_cbsniffer_output) / n(), fraction_genotyped_GoT = sum(has_GoT_output) / n())
md.cb.sniffer.summarized.sample %>% arrange(-fraction_genotyped_cbsniffer)

## Pivot to a longer format for plotting
md.cb.sniffer.summarized.sample.formatted <- md.cb.sniffer.summarized.sample %>% 
  pivot_longer(cols = c(fraction_genotyped_cbsniffer, fraction_genotyped_GoT), names_to = c("genotyping_method"), values_to = "genotyping_fraction") %>% 
  mutate(genotyping_method = gsub("fraction_genotyped_", "", genotyping_method)) %>% 
  mutate(`Sample Type` = if_else(grepl("CD34$|BM15|plus|pos", orig.ident), "CD34+", "CD34-/unenriched"))

## Plot as bar plot
p.genotyping.comparison <- md.cb.sniffer.summarized.sample.formatted %>% 
  group_by(genotyping_method) %>% 
  summarize(mean_val = mean(genotyping_fraction), sd_val = sd(genotyping_fraction), n = n(), se_val = sd_val / sqrt(n)) %>% 
  ggplot(aes(x = genotyping_method, y = mean_val)) +
  geom_col(color = "black", fill = "grey80") +
  # scale_fill_manual(values = c("#FDD49D", "#F58E3D")) +
  geom_errorbar(aes(x = genotyping_method, ymin = mean_val-se_val, ymax = mean_val+se_val), width=0.4, size = 0.25, color = "black") +
  geom_point(data = md.cb.sniffer.summarized.sample.formatted, mapping = aes(x = genotyping_method, y = genotyping_fraction,  color = `Sample Type`), alpha = 0.5) +
  # guides(color = FALSE, size = FALSE)
  geom_line(data = md.cb.sniffer.summarized.sample.formatted, mapping = aes(x = genotyping_method, y = genotyping_fraction, group = orig.ident, color = `Sample Type`), linetype = "longdash", alpha = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("black", "#F58E3D")) +
  xlab("Genotyping method") + ylab("Percent of cells genotyped") +
  ggtitle("Genotyping rate", subtitle = "GEX vs. targeted UBA1 GoT library")
p.genotyping.comparison
ggsave("figures/current_figure_drafts/RNA_GEX_GoT_genotyping_comparison.pdf", device = "pdf", dpi = 300, width = 3.5, height = 4.5)

######################### Looking at cell cycle ##################################

md %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  ggplot(aes(x = S.Score, y = G2M.Score, color = Genotype)) +
  geom_point(alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "longdash") +
  geom_vline(xintercept = 0, linetype = "longdash") +
  # coord_cartesian(ylim = c(-0.2, 0.2)) +
  facet_grid(cols = vars(Genotype), rows = vars(cluster_celltype), scales = "free") +
  scale_color_manual(values = vexas.genotyping.palette)

plist <- lapply(c("HSC", "EMP", "LMPP", "MkP", "GMP"), function(x) {
  donors.keep <- md %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    filter(cluster_celltype == x) %>% 
    count(Donor, Genotype) %>% filter(n >= 5) %>% ungroup() %>% count(Donor) %>% filter(n > 1) %>% pull(Donor)
  
  md %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    filter(cluster_celltype == x) %>% 
    filter(Donor %in% donors.keep) %>% 
    group_by(Donor, Genotype) %>% 
    # filter(n() > 10) %>% 
    summarize(count = n(), S_or_G2M_count = sum(Phase %in% c("S", "G2M")), g1_count =sum(Phase == "G1"), S_G2M_fraction = sum(Phase %in% c("S", "G2M")) / n()) %>% 
    ggplot(aes(x = Genotype, y = S_G2M_fraction)) +
    geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Donor)) +
    geom_line(aes(group = Donor, color = Donor)) +
    stat_compare_means() +
    ggtitle(x) +
    scale_y_continuous(expand = c(0.2, 0))
})

wrap_plots(plist)

## Looking at correlation with apoptosis
md

######################### Making DNMT3A overlap plots ##################################

ont.data <- read_csv("data/ont_genotyping/ONT_Genotyping.csv")
ont.data %>% 
  count(DNMT3A_genotype, UMI_WT_DNMT3A, UMI_MUT_DNMT3A)
p.ont <- ont.data %>% 
  filter(!is.na(UBA1_genotype) & !is.na(DNMT3A_genotype)) %>% 
  group_by(UBA1_genotype, DNMT3A_genotype) %>% 
  count() %>% 
  mutate(`Genotype status` = paste0(UBA1_genotype, DNMT3A_genotype)) %>% 
  mutate(`Genotype status` = plyr::mapvalues(`Genotype status`, from = c("MUTMUT", "MUTWT", "WTWT"), to = c("UBA1 MUT, DNMT3A MUT", "UBA1 MUT, DNMTA3A WT", "UBA1 WT, DNMT3A WT"))) %>% 
  ggplot(aes(y = n, axis1 = UBA1_genotype, axis2 = DNMT3A_genotype)) +
  geom_alluvium(aes(fill = `Genotype status`), width = 1/6) +
  geom_stratum(width = 1/3, alpha = 1) +
  scale_fill_manual(values = c("grey10", "grey50", "grey70")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("UBA1_genotype", "DNMT3A_genotype"), expand = c(.05, .05)) +
  ylab("Cell count") +
  ggtitle("PT02, DNMT3A vs. UBA1 Genotype")
p.ont
ggsave("figures/current_figure_drafts/RNA_ONT_DNMT3A_genotyping.pdf", device = "pdf", dpi = 300, width = 6, height = 4)


######################### Save the number of genes expressed in at least 10 cells for each cell type (for ATAC) ##################

for (celltype in c("HSC", "EMP", "LMPP", "MkP", "GMP")) {
  se <- subset(rna.obj, CellType == celltype)
  gene.counts <- rowSums(se@assays$RNA@counts > 0)
  gene.counts <- gene.counts[gene.counts > 20]
  data.frame(gene = names(gene.counts)) %>% 
    write_csv(paste0("data/expressed_gene_lists/expressed_20_cells_", celltype, ".csv"))
}

######################## Resave the seurat object after re-running RunUMAP (for control integration) #######################

DefaultAssay(rna.obj) <-"integrated"
rna.obj <- RunUMAP(rna.obj, reduction = "pca", dims = 1:50, seed.use = 1, return.model = TRUE, reduction.name = "umap_v2")
rna.obj <- FindNeighbors(rna.obj, reduction = "pca", dims = 1:50)
rna.obj <- FindClusters(rna.obj, resolution = 2)
p1 <- DimPlot(rna.obj, group.by = "seurat_clusters", reduction = "umap", label = T, label.box = F) + NoLegend()
p2 <- DimPlot(rna.obj, group.by = "seurat_clusters", reduction = "umap_v2", label = T, label.box = F) + NoLegend()
p3 <- DimPlot(rna.obj, group.by = "cluster_celltype", reduction = "umap_v2", label = T, label.box = F) + NoLegend()
p2 | p3
saveRDS(rna.obj, "~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_UMAP_rerun.rds")
