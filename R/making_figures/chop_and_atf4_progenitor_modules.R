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

# ## Adding pseudotime
# vexas_pseudotime = readRDS('/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pseudotime/vexas_pseudotime_cells.Rds')
# rna.obj <- AddMetaData(rna.obj, vexas_pseudotime@meta.data %>% select(monocle3_pseudotime))
# rm(vexas_pseudotime)

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "cDC", "CD14 Mono", "B", "pDC", "CD4 T", "CD8 T", "NK", "Plasma")
rna.obj$CellType <- factor(rna.obj$cluster_celltype, levels = celltype.order)
table(rna.obj$CellType, useNA = "always")

## Formatting the sample and donor columns
meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_manually_updated_2024-03-29.csv")
rna.obj$Sample <- plyr::mapvalues(rna.obj$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$updated_name)
rna.obj$Sample <- factor(rna.obj$Sample, levels = meta.data.mapping$updated_name %>% sort())
rna.obj@meta.data %>% count(orig.ident, Sample)
rna.obj$Donor <- gsub("_BMMC|_CD34.*", "", rna.obj$Sample) %>% gsub("_[A|B]", "", .)
rna.obj$Donor <- factor(rna.obj$Donor, levels = unique(rna.obj$Donor) %>% sort())
rna.obj@meta.data %>% count(orig.ident, Sample, Donor)

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

################################ Adding module scores ###########################

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

## Specify modules to use
module.list = list(
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Inflammation",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Inflammation",
  `REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM` = "Inflammation",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` =  "Differentiation",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` =  "Differentiation",
  `REACTOME_APOPTOSIS` = "Survival",
  `REACTOME_TRANSLATION` = "Survival",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Survival",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "Survival",
  `GOBP_PROTEIN_FOLDING` = "Survival",
  `GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "Survival",
  `CHOP_ONLY_HAN` = "Survival",
  `ATF4_ONLY_HAN` = "Survival"
)

rna.obj.subset <- subset(rna.obj, cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "CLP"))

## Add as modules to our RNA object
DefaultAssay(rna.obj.subset) <- "RNA"
for (module in names(module.list)) {
  print(module)
  rna.obj.subset <- AddModuleScore(rna.obj.subset, features = list(pathways[[module]]), name = module, assay = "RNA", search = TRUE)
}


########################################### Plot across progenitors ###########################################

## plot across mut progenitors
md <- rna.obj.subset@meta.data
md$CellType <- factor(md$CellType, levels = c("HSC", "EMP", "LMPP", "MkP", "CLP"))

md %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count(Genotype, CellType) %>% 
  pivot_wider(names_from = Genotype, values_from = n)

p2 <- md %>% 
  filter(Genotype == "MUT") %>% 
  ggplot(aes(x = CellType, y = ATF4_ONLY_HAN1)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("HSC", "MkP"), c("HSC", "CLP")))
p1 <- md %>% 
  filter(Genotype == "MUT") %>% 
  ggplot(aes(x = CellType, y = CHOP_ONLY_HAN1)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("HSC", "MkP"), c("HSC", "CLP")))
p3 <- md %>% 
  filter(Genotype == "MUT") %>% 
  ggplot(aes(x = CellType, y = REACTOME_APOPTOSIS1)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("HSC", "MkP"), c("HSC", "CLP")))
p2 / p1 / p3 + plot_annotation(title = "MUT cells only")


md %>% 
  filter(Genotype == "MUT") %>% 
  group_by(CellType, Donor) %>% 
  filter(n() > 5) %>% 
  ggplot(aes(x = CellType, y = REACTOME_APOPTOSIS1)) +
  geom_boxplot() +
  facet_wrap(~ Donor)
