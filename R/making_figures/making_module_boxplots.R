library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(ggpubr)
library(msigdbr)
library(clusterProfiler)
theme_set(theme_classic())

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")


######################################## Specific module boxplots #####################################

current.cell.type = "HSC"

rna.obj.subset <- subset(rna.obj, cluster_celltype == current.cell.type & Genotype %in% c("MUT", "WT"))

## Compile modules we're interested in
m_t2g <- tibble()
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy"))
pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea

## Add GoT paper modules
got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv") %>% as.list()
got.modules <- lapply(got.modules, function(x) return(x[!is.na(x)]))
pathways <- c(pathways, got.modules)

## Add Han et al modules
han.modules <- list(
  `ATF4_ONLY_HAN` = read_csv("data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
  `CHOP_ONLY_HAN` = read_csv("data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
  `CHOP_AND_ATF4_HAN` = read_csv("data/gene_module_lists/Han_et_al_references/CHOP_and_ATF4.txt") %>% filter(Feature != "N/A") %>% pull(Feature)
)
pathways <- c(pathways, han.modules)

module.list <- list(
  # `GOBP Myeloid cell differentiation` = pathways$GOBP_MYELOID_CELL_DIFFERENTIATION,
  `Hallmark Inflammation` = pathways$HALLMARK_INFLAMMATORY_RESPONSE,
  `Hallmark IFNa Response` = pathways$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  `Hallmark IFNg Response` = pathways$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  `Hallmark Apoptosis` = pathways$HALLMARK_APOPTOSIS,
  `Reactome Apoptosis` = pathways$REACTOME_APOPTOSIS,
  `Hallmark G2M checkpoint` = pathways$HALLMARK_G2M_CHECKPOINT,
  `Anti-apoptosis` = pathways$`Anti-apoptosis module`,
  `Reactome translation` = pathways$REACTOME_TRANSLATION,
  `GOBP ATF6` = pathways$GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
  `GOBP IRE1` = pathways$GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
  # `HAllmark UPR` = pathways$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
  `ATF4 target genes` = pathways$ATF4_ONLY_HAN,
  `CHOP target genes` = pathways$CHOP_ONLY_HAN,
  # # `XBP1 target module` = pathways$`XBP1 target module`,
  # # `XBP1_01` = pathways$XBP1_01
  # # `Antigen presentation` = pathways$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION
  # `KEGG_REGULATION_OF_AUTOPHAGY` = pathways$KEGG_REGULATION_OF_AUTOPHAGY,
  # `REACTOME_AUTOPHAGY` = pathways$REACTOME_AUTOPHAGY,
  # `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = pathways$REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY,
  # HALLMARK_TNFA_SIGNALING_VIA_NFKB = pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
  # HALLMARK_P53_PATHWAY = pathways$HALLMARK_P53_PATHWAY,
  # KEGG_P53_SIGNALING_PATHWAY = pathways$KEGG_P53_SIGNALING_PATHWAY,
  # P53_02 = pathways$P53_02,
  # P53_DECAMER_Q2 = pathways$P53_DECAMER_Q2,
  `PERK module` = pathways$`PERK module`,
  `GOBP PERK` = pathways$GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
  `GOBP PERK (2)` = pathways$GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
  `GOBP PERK (3)` = pathways$GOBP_NEGATIVE_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE
  # KEGG_LYSOSOME = pathways$KEGG_LYSOSOME,
  # KEGG_PEROXISOME = pathways$KEGG_PEROXISOME,
  # KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY = pathways$KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY,
  # REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES = pathways$REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES,
  # ISHIKAWA_STING_SIGNALING = pathways$ISHIKAWA_STING_SIGNALING
  # 
)

## Add as modules to our RNA object
DefaultAssay(rna.obj.subset) <- "RNA"
for (module in names(module.list)) {
  print(module)
  print(module.list[[module]])
  rna.obj.subset <- AddModuleScore(rna.obj.subset, features = list(module.list[[module]]), name = module, assay = "RNA", search = TRUE)
}

## Make boxplots for the modules
plot_module_boxplot <- function(module, color_pal) {
  module.name <- paste0(module, "1")
  min.quantile <- rna.obj.subset@meta.data %>%
    filter(Genotype %in% c("MUT", "WT")) %>%
    dplyr::rename(module_plotting = module.name) %>%
    group_by(Genotype) %>% 
    summarize(quantile = quantile(module_plotting, c(0.01))) %>% 
    slice_min(quantile) %>% 
    pull(quantile)
  max.quantile <- rna.obj.subset@meta.data %>%
    filter(Genotype %in% c("MUT", "WT")) %>%
    dplyr::rename(module_plotting = module.name) %>%
    group_by(Genotype) %>% 
    summarize(quantile = quantile(module_plotting, c(0.99))) %>% 
    slice_max(quantile) %>% 
    pull(quantile)
  
  p1 <- rna.obj.subset@meta.data %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    dplyr::rename(module_plotting = module.name) %>%
    ggplot(aes(x = Genotype, y = module_plotting, fill = Genotype)) +
    # geom_boxplot(outlier.shape = NA) +
    geom_boxplot() +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(gsub("_", " ", module) %>% str_wrap(., width = 15) %>% gsub("\\(", "\n\\(", .)) + 
    # ggtitle(str_wrap(module, width = 15) %>% gsub("\\(", "\n\\(", .)) + 
    # ggtitle(str_wrap(module, width = 17) %>% gsub("\\(", "\n\\(", .)) +
    # scale_fill_manual(values = c(WT = "lightsteelblue1", MUT = "lightsteelblue3")) +
    scale_fill_manual(values = color_pal) +
    # scale_y_continuous(expand = c(0.2, 0), limits = c(min.quantile, max.quantile)) +
    scale_y_continuous(expand = c(0.2, 0)) +
    stat_compare_means(label = "p.format", label.y = Inf, vjust = 1.5, size = 3) +
    ylab("Module score") +
    xlab("")
  return(p1)
}

## Plot boxplots for figure
plist <- lapply(c("GOBP ATF6", "GOBP IRE1", "ATF4 target genes", "CHOP target genes", "Reactome translation"), function(x) {plot_module_boxplot(x, color_pal = c("#FFCB88", "#FF7B0B"))})
p.module.plots <- wrap_plots(plist, nrow = 1) + plot_annotation(title = current.cell.type)
p.module.plots
ggsave(paste0("figures/current_figure_drafts/boxplots_fig3_", current.cell.type, ".pdf"), plot = p.module.plots, device = "pdf", dpi = 300, width = 9, height = 4, unit = "in")

## Additional boxplots
plist <- lapply(c("GOBP Myeloid cell differentiation", "Hallmark Inflammation", "Hallmark IFNa Response"), function(x) {plot_module_boxplot(x, color_pal = c("#FFCB88", "#FF7B0B"))})
p.module.plots <- wrap_plots(plist, nrow = 1)
p.module.plots
ggsave("figures/current_figure_drafts/RNA_module_score_boxplots_inflammation.pdf", plot = p.module.plots, device = "pdf", dpi = 300, width = 5, height = 3, unit = "in")

plist <- lapply(c("Reactome translation", "GOBP ATF6", "GOBP IRE1"), function(x) {plot_module_boxplot(x, color_pal =c(WT = "#FFB4C0", MUT = "#CD2135"))})
p.module.plots <- wrap_plots(plist, nrow = 1)
p.module.plots
ggsave("figures/current_figure_drafts/RNA_module_score_boxplots_atf6_ire1.pdf", plot = p.module.plots, device = "pdf", dpi = 300, width = 5, height = 3, unit = "in")

plist <- lapply(c("GOBP ATF6", "GOBP PERK", "GOBP IRE1"), function(x) {plot_module_boxplot(x, color_pal =c(WT = "#FFB4C0", MUT = "#CD2135"))})
p.module.plots <- wrap_plots(plist, nrow = 1)
p.module.plots
ggsave("figures/current_figure_drafts/RNA_module_score_boxplots_atf6_ire1_perk.pdf", plot = p.module.plots, device = "pdf", dpi = 300, width = 5, height = 3, unit = "in")

plist <- lapply(c("ATF4 target genes", "CHOP target genes", "Reactome translation"), function(x) {plot_module_boxplot(x, color_pal =c(WT = "#E7D4E9", MUT = "#85378A"))})
p.module.plots <- wrap_plots(plist, nrow = 1)
p.module.plots
ggsave("figures/current_figure_drafts/RNA_module_score_boxplots_atf4.pdf", plot = p.module.plots, device = "pdf", dpi = 300, width = 5, height = 3, unit = "in")

plist <- lapply(c("GOBP PERK", "GOBP PERK (2)", "GOBP PERK (3)", "PERK module"), function(x) {plot_module_boxplot(x, color_pal =c(WT = "#E7D4E9", MUT = "#85378A"))})
p.module.plots <- wrap_plots(plist, nrow = 1)
p.module.plots
ggsave("figures/current_figure_drafts/RNA_module_score_boxplots_perk_modules.pdf", plot = p.module.plots, device = "pdf", dpi = 300, width = 7, height = 3, unit = "in")





########################################## Plotting module correlations - cell cycle ###################################

rna.obj.subset@meta.data %>% 
  mutate(cell_cycle_status = if_else(Phase %in% c("S", "G2M"), "S/G2M", "G1")) %>%
  ggplot(aes(x = cell_cycle_status, y = `ATF4 target genes1`, color = Genotype)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.5) +
  facet_grid(cols = vars(Genotype)) +
  stat_compare_means(comparisons = list(c("G1", "S/G2M")), ) +
  scale_color_manual(values = vexas.genotyping.palette) +
  ggtitle("HSCs")

plot_module <- function(module.current) {
  print(module.current)
  rna.obj.subset@meta.data %>% 
    dplyr::rename(current_module = module.current) %>% 
    mutate(cell_cycle_status = if_else(Phase %in% c("S", "G2M"), "S/G2M", "G1")) %>%
    ggplot(aes_string(x = "cell_cycle_status", y = "current_module", color = "Genotype")) +
    geom_boxplot(outlier.shape = NA) +
    ylab(module.current) +
    geom_point(position = position_jitter(width = 0.25), alpha = 0.5) +
    facet_grid(cols = vars(Genotype)) +
    stat_compare_means(comparisons = list(c("G1", "S/G2M")), ) +
    scale_color_manual(values = vexas.genotyping.palette) +
    ggtitle("HSCs") +
    theme(legend.position = "none")
}

plot_module("Hallmark Apoptosis1")
plot_module("Anti-apoptosis1")
plot_module("Reactome Apoptosis1")
plot_module("GOBP IRE11")
plot_module("GOBP ATF61")
plot_module("HAllmark UPR1")
plot_module("Reactome translation1")
plot_module("Hallmark IFNa Response1")

plot_module("KEGG_REGULATION_OF_AUTOPHAGY1") /
plot_module("REACTOME_AUTOPHAGY1") /
plot_module("REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY1")

plot_module("percent.ribo")

plot_module("HALLMARK_TNFA_SIGNALING_VIA_NFKB1") 
plot_module("HALLMARK_P53_PATHWAY1") /
plot_module("KEGG_P53_SIGNALING_PATHWAY1")
plot_module("P53_021") /
plot_module("P53_DECAMER_Q21") 
plot_module("PERK module1")
plot_module("GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE1")
plot_module("KEGG_LYSOSOME1") |
plot_module("KEGG_PEROXISOME1")
plot_module("KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY1") |
plot_module("REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES1") |
plot_module("ISHIKAWA_STING_SIGNALING1")

########################################## Plotting module correlations ###################################

rna.obj.subset@meta.data %>% 
  ggplot(aes(x = `Hallmark IFNa Response1`, y = `GOBP IRE11`, color = Genotype)) +
  geom_point() +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette)
  # ylim(c(0.75, 1.25)) +
  # xlim(c(-0.05, 0.1))

rna.obj.subset@meta.data %>% 
  ggplot(aes(x = `S.Score`, y = `G2M.Score`, color = `Hallmark IFNa Response1`)) +
  geom_point(alpha = 0.2) +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  scale_color_continuous2(low = "blue" high = "red")
  # scale_color_manual(values = vexas.genotyping.palette)
# ylim(c(0.75, 1.25)) +
# xlim(c(-0.05, 0.1))

rna.obj.subset@meta.data %>% 
  ggplot(aes(x = `Hallmark IFNa Response1`, y = `CHOP target genes\n(Han et al.)1`, color = Genotype)) +
  geom_point() +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette)


# ylim(c(0.75, 1.25)) +
# xlim(c(-0.05, 0.1))
rna.obj.subset@meta.data %>% 
  ggplot(aes(x = `Hallmark IFN Alpha1`, y = `GOBP IRE11`, color = Genotype)) +
  geom_point() +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  facet_grid(cols = vars(Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette)


## Plot percent ribo vs. atf4 content overlap
rna.obj.subset@meta.data %>% 
  summarize(quantile = quantile(percent.ribo, c(0.01, 0.99)))
ribo.plot <- rna.obj.subset@meta.data %>% 
  # filter(percent.ribo < 45 & percent.ribo > 24) %>% 
  sample() %>% 
  ggplot(aes(y = percent.ribo, x = `ATF4 target genes [Han et al.]1`)) +
  geom_point(alpha = 0.5, aes(color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  stat_cor()+
  xlim(c(-0.05, 0.1)) +
  ylim(c(7, 60)) +
  # facet_grid(cols = vars(Genotype)) +
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  ylab("Percent ribosomal reads") +
  xlab("ATF4 module score")
ribo.plot
ggsave("figures/current_figure_drafts/RNA_ribo_vs_atf4.pdf", plot = ribo.plot, device = "pdf", dpi = 300, width = 4.5, height = 3.5, unit = "in")

p.ribo.violin <- rna.obj.subset@meta.data %>% 
  ggplot(aes(x = Genotype, y = percent.ribo, fill = Genotype)) +
  geom_violin() +
  # geom_boxplot(width = 0.3) +
  # geom_point(position = position_jitter(width = 0.3), size = 0.05) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.2, 0)) +
  scale_fill_manual(values = vexas.genotyping.palette) +
  stat_compare_means(label = "..p.format..", label.y = 80)  +
  xlab(element_blank()) +
  ylab("Percent ribo reads")
p.ribo.violin
ggsave("figures/current_figure_drafts/RNA_ribo_violins.pdf", plot = p.ribo.violin, device = "pdf", dpi = 300, width = 2.5, height = 3, unit = "in")

## Plot ribosomal content overlap
ribo.plot <- rna.obj.subset@meta.data %>% 
  # filter(percent.ribo < 50 & percent.ribo > 20) %>% 
  sample() %>% 
  ggplot(aes(y = `Hallmark interferon alpha response1`, x = `ATF4 target genes [Han et al.]1`)) +
  geom_point(alpha = 0.5, aes(color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  stat_cor()+
  # facet_grid(cols = vars(Genotype)) +
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  ylab("Percent ribosomal reads") +
  xlab("ATF4 module score")
ribo.plot

## Plot ribosomal content overlap with cd99

ribo.plot.2 <- merge(rna.obj.subset@meta.data, data.frame(CD99 = rna.obj.subset@assays$ADT_DSB@data["Hu.CD99", ]), by = "row.names") %>% 
  ggplot(aes(y = percent.ribo, x = CD99)) +
  geom_point(alpha = 0.5, aes(color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  ylab("Percent ribosomal reads") +
  xlab("CD99 ADT expression")
ggsave("figures/current_figure_drafts/RNA_ribo_vs_cd99.pdf", plot = ribo.plot.2, device = "pdf", dpi = 300, width = 4.5, height = 3.5, unit = "in")

merge(rna.obj.subset@meta.data, data.frame(CD99 = rna.obj.subset@assays$ADT_DSB@data["Hu.CD99", ]), by = "row.names") %>% 
  ggplot(aes(y = `ATF4 target genes [Han et al.]1`, x = CD99)) +
  geom_point(alpha = 0.5, aes(color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  stat_cor()+
  geom_smooth(method = 'lm',se=F, linetype = "solid", color = "black") +
  ylab("Percent ribosomal reads") +
  xlab("CD99 ADT expression")
