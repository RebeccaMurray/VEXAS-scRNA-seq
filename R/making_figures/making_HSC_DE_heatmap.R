library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
theme_set(theme_classic())

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
source("~/rscripts/general_helpers/VEXAS_color_palettes.R")

## All available: c("HSC", "EMP", "LMPP", "MkP", "Early Eryth", "GMP 1", "GMP 2", "Late Eryth")
cell.type.list <- c("HSC", "LMPP", "EMP", "MkP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/10p_genotyped_cells_no_RP_MT_renormalized_3000_var_features/", x, ".csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

########################################## Make a heatmap with the top DE genes in HSCs #####################################

## Pull the top N unique genes for MUT/WT
num.markers.to.pull = 50
mut.markers <- de.results[["HSC"]] %>%
  mutate(stat = avg_log2FC * -log10(pval)) %>% 
  filter(avg_log2FC > 0.2) %>% 
  slice_max(n = num.markers.to.pull, order_by = stat) %>%
  pull(feature)
wt.markers <- de.results[["HSC"]] %>%
  mutate(stat = avg_log2FC * -log10(pval)) %>% 
  filter(avg_log2FC < -0.2) %>% 
  slice_min(n = num.markers.to.pull, order_by = stat) %>%
  arrange(-stat) %>% 
  pull(feature)
top.genes <- c(mut.markers, wt.markers)

## Scale the RNA assay for just these genes
se <- subset(rna.obj, CellType == "HSC")
se <- ScaleData(se, features = top.genes, assay = "RNA")

## Downsample to 50 BCs per cluster
downsampled.bc.list <- se@meta.data %>%
  filter(CellType == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>%
  rownames_to_column("Barcode") %>%
  group_by(Genotype) %>%
  arrange(Genotype) %>% 
  slice_sample(n = 300, replace = F) %>%
  pull(Barcode)

## Pull metadata for annotation row
annotation.df <- se@meta.data[downsampled.bc.list, ] %>%
  dplyr::select(Genotype)

## Heatmap
dev.off()
pheatmap::pheatmap(se@assays$RNA@scale.data[top.genes, downsampled.bc.list],
                   annotation_col = annotation.df,
                   # scale = "row",
                   annotation_colors = list(`Genotype` = vexas.genotyping.palette),
                   show_colnames = F, cluster_rows = F, cluster_cols = F, 
                   show_rownames = F,
                   color=colorRampPalette(c("navy", "white", "red"))(50),  breaks = seq(-1.5, 1.5, length.out = 50),
                   main = "HSC cluster\nTop 50 differentially expressed genes")
dev.off()


####################################### Extra - exploring expression ############################

VlnPlot(se, features = "IGLL1", group.by = "Genotype", assay = "RNA", slot = "scale.data")
se <- subset(rna.obj, CellType == "HSC" & Genotype != "NA")

data.frame(exp = se@assays$RNA@data["CXCL8", ] > 0) %>% 
  count(exp)
