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

######################################### Load all the ironthrone outputs ################################

## Mapping dictionary
sample.dict <- list(
  NIH_1 = "UPN14_BMMC",
  NIH_2 = "UPN15_BMMC",
  NIH_3 = "UPN16_BMMC",
  NIH_4 = "UPN17_BMMC",
  NIH_5 = "UPN14_CD34",
  NIH_6 = "UPN15_CD34",
  NIH_7 = "UPN16_CD34",
  NIH_8 = "UPN17_CD34",
  VX_BM2_CD34_minus = "VEXAS_c_BM2_CD34-SR",
  VX_BM2_CD34_plus = "VEXAS_c_BM2_CD34_plus_SR",
  VX_BM1_CD34_minus = "VEXAS_BM1_CD34-_SR",
  VX_BM1_CD34_plus = "VEXAS_BM1_CD34_plus_SR",
  VX_BM1_CD34_plus = "VEXAS_BM1_CD34_plus_SR",
  VX_BM12_A = "VX_BM12_A_GEX",
  VX_BM13_B = "VX_BM13_B_GEX",
  VX_BM14_Fneg = 'VX_BM14_Fneg_GEX',
  VX_BM14_Fpos = "VX_BM14_Fpos_GEX",
  VX_BM15a = "SG_VX_BM15a_GEX",
  VX_BM15b = "SG_VX_BM15b_GEX"
)

genotyping.dirs <- list.dirs("/gpfs/commons/home/rmurray/VEXAS/GoT", recursive = F)
genotyping.files <- lapply(genotyping.dirs, function(x) {
  geno.file <- list.files(x, full.names = T, recursive = F) %>% grep("myGoT.summTable.concat.umi_collapsed.txt", ., value = T)
  df <- read_table(geno.file)
  sample.name <- gsub("/gpfs/commons/home/rmurray/VEXAS/GoT/", "", geno.file) %>% gsub("/myGoT.summTable.concat.umi_collapsed.txt", "", .) %>% 
    gsub("_GoT", "", .)
  df$orig.ident <- sample.name
  df$orig.ident <- plyr::mapvalues(df$orig.ident, from = names(sample.dict), to = unlist(sample.dict))
  df$Barcode <- paste0(df$orig.ident, "_", df$BC, "-1")
  return(df)
})

genotyping.df <- do.call(rbind, genotyping.files)
setdiff(unique(rna.obj$orig.ident), names(table(genotyping.df$orig.ident)))
setdiff(names(table(genotyping.df$orig.ident)), unique(rna.obj$orig.ident))


#################################### Add a "naive" genotyping ###########################

genotyping.df <- genotyping.df %>% 
  mutate(Genotype_simple = if_else(MUT.calls > 0, "MUT", if_else(WT.calls > 0, "WT", NA_character_))) %>% 
  select(Barcode, Genotype_simple) %>% 
  column_to_rownames("Barcode")

rna.obj <- AddMetaData(rna.obj, genotyping.df)

rna.obj@meta.data %>% 
  filter(cluster_celltype == "HSC") %>% 
  count(Genotype)

rna.obj@meta.data %>% 
  filter(cluster_celltype == "HSC") %>% 
  count(Genotype_simple)


################################## Plot normalized mut cell freq with naive genotyping #####################

md <- rna.obj@meta.data

median.pseudotime.values <- md %>% 
  group_by(CellType) %>% 
  summarize(median_pseudotime = median(monocle3_pseudotime))

overall_mut_fraction <- md %>% 
  filter(Genotype_simple %in% c("MUT", "WT")) %>% 
  count(Genotype_simple) %>% 
  pivot_wider(names_from = Genotype_simple, values_from = n) %>% 
  mutate(overall_geno_frac = MUT / (MUT + WT)) %>% 
  pull(overall_geno_frac)

## Summarize the MUT cell frequency per donor
mut.fraction.per.donor <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype_simple %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  mutate(CellType = gsub("GMP [0-9]", "GMP", CellType)) %>% 
  mutate(CellType = gsub("CD4 |CD8 ", "", CellType))  %>%
  filter(n() > 100) %>% ## Only including cell types with > 20 Genotype_simpled cells
  group_by(CellType, Donor) %>% 
  filter(n() > 10) %>%  ## Only including patient data points with > 10 Genotype_simpled cells
  group_by(CellType, Genotype_simple, Donor) %>%
  summarize(n = n()) %>% 
  group_by(CellType, Donor) %>% 
  # mutate(mut_fraction = n / sum(n), overall_mut_fraction = overall_count_Genotype_simple/overall_count, normalized_mut_fraction = mut_fraction / overall_mut_fraction) %>% View
  mutate(mut_fraction = n / sum(n), normalized_mut_fraction = mut_fraction / overall_mut_fraction) %>% 
  filter(Genotype_simple == "MUT") %>% 
  group_by(CellType) %>% 
  filter(n_distinct(Donor) >= 3)

summarized.freqs <- mut.fraction.per.donor %>% 
  summarize(mean_val = mean(normalized_mut_fraction), sd_val = sd(normalized_mut_fraction), n = n(), se_val = sd_val / sqrt(n), median_val = median(normalized_mut_fraction)) %>% 
  merge(., median.pseudotime.values, group.by = CellType, all.x = T)

## Make the bar plot
p.mutant_cell_fraction.bar.plot.naive <- summarized.freqs %>% 
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
p.mutant_cell_fraction.bar.plot.naive
ggsave(paste0(current.plots.path, "/normalized_MUT_cell_freq_bar_plot_naive_genotyping.pdf"), plot = p.mutant_cell_fraction.bar.plot.naive, device = "pdf", dpi = 300, width = 6, height = 4, unit = "in")


########################################## Making UMAP - naive genotyping ###############################

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

## Genotyping UMAP (individual points)
mut.count <- md %>% filter(Genotype_simple == "MUT") %>% nrow()
wt.count <- md %>% filter(Genotype_simple == "WT") %>% nrow()
na.or.het.count <- md %>% filter(Genotype_simple == "NA") %>% nrow()
mut.title <- paste0("MUT (n = ", format(mut.count, big.mark = ","), ")")
wt.title <- paste0("WT (n = ", format(wt.count, big.mark = ","), ")")
na.title <- paste0("NA (n = ", format(na.or.het.count, big.mark = ","), ")")
md.wt <- md[md$Genotype_simple %in% c("WT"), ]
md.mut <- md[md$Genotype_simple %in% c("MUT"), ]
md.na <- md[md$Genotype_simple %in% c("NA", "Control"), ]

p.genotyping <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, fill = Genotype_simple)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.mut, shape = 21, size = 1, fill = vexas.genotyping.palette[["MUT"]], stroke = 0.1, color = "black") +
  geom_point(data = md.wt, shape = 21, size = 1, fill = vexas.genotyping.palette[["WT"]], stroke = 0.1, color = "black") +
  # geom_point(data = md.Genotype_simpled, shape = 21, size = 1) +
  # scale_color_manual(values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none") + coord_fixed()
p.genotyping
ggsave("figures/current_figure_drafts/UMAP_RNA_naive_genotyping.tiff", plot = p.genotyping, device = "tiff", dpi = 300, width = 5, height = 4, unit = "in")


