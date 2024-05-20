library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
theme_set(theme_classic())

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
source("~/rscripts/general_helpers/VEXAS_color_palettes.R")

## All available: c("HSC", "EMP", "LMPP", "MkP", "Early Eryth", "GMP 1", "GMP 2", "Late Eryth")
cell.type.list <- c("HSC","EMP", "LMPP", "MkP", "GMP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/", x, "_5cells_genotyped_20231206.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

################################## Making combined dotplot ###########################

# fgsea.res.all <- lapply(names(de.results), function(x) { 
#   print("Running for:")
#   print(x)
#   fgsea.res <- run_fgsea_for_list(de.results[[x]])
#   fgsea.res[["cluster"]] <- x
#   return(fgsea.res)
# }
# )
# 
# fgsea.res.combined <- do.call(rbind, fgsea.res.all)
# fgsea.res.combined %>% 
#   write_csv("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/fgsea_results.csv")

fgsea.res.combined <- read_csv("figures/current_figure_drafts/fgsea_used_for_heatmap_avglog2FC_ranking_2024-02-27.csv")

## Isolate and name pathways for a more curated dot plot
## Make a GSEA heatmap - NES, targets we pick 
fgsea.pathways.specific = list(
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "IFNa",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "IFNg",
  `REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM` = "Cytokine signaling",
  `HALLMARK_P53_PATHWAY` = "Hallmark P53 pathway",
  `REACTOME_MTOR_SIGNALLING` = "Reactome MTOR signaling",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` = "Leukocyte differentiation",
  `GOBP_MYELOID_CELL_DIFFERENTIATION` = "Myeloid cell differentiation",
  `GOBP_MONONUCLEAR_CELL_DIFFERENTIATION` = "Differentiation",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "Antigen presentation",
  `KEGG_OXIDATIVE_PHOSPHORYLATION` = "KEGG oxidative phosphorylation",
  `HALLMARK_GLYCOLYSIS` = "Hallmark glycolysis",
  `REACTOME_APOPTOSIS` = "Apoptosis",
  `GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES` = "GOBP ROS",
  `REACTOME_TRANSLATION` = "Translation",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Autophagy",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` =   "GOBP UPR - IRE1",
  `GOBP_PROTEIN_FOLDING` =   "GOBP protein folding",
  `ATF6 module` = "ATF6 activation",
  `CHOP targets (Han et al.)` = "CHOP targets (Han et al.)",
  `ATF4 targets (Han et al.)` = "ATF4 targets (Han et al.)"
  # `Shared ATF4/CHOP targets (Han et al.)` = "Shared ATF4/CHOP targets (Han et al.)"
)
names(fgsea.pathways.specific) <- str_wrap(names(fgsea.pathways.specific), width = 30)

fgsea.pathways.df <- data.frame(pathway = names(fgsea.pathways.specific), pathway.name = unlist(fgsea.pathways.specific))

## Specify blocks
fgsea.blocks = list(
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Inflammation",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Inflammation",
  `REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM` = "Inflammation",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` =  "Differentiation",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` =  "Differentiation",
  # `GOBP_MONONUCLEAR_CELL_DIFFERENTIATION` = "Differentiation",
  `REACTOME_APOPTOSIS` = "Survival",
  `REACTOME_TRANSLATION` = "Survival",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "Survival",
  `GOBP_PROTEIN_FOLDING` = "Survival",
  `ATF6 module` = "Survival",
  `CHOP targets (Han et al.)` = "Survival",
  `ATF4 targets (Han et al.)` = "Survival",
  # `Shared ATF4/CHOP targets (Han et al.)` = "2",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Survival"
)
names(fgsea.blocks) <- str_wrap(names(fgsea.blocks), width = 30)
fgsea.blocks.df <- data.frame(pathway = names(fgsea.blocks), block = unlist(fgsea.blocks))

## Plot the fgsea dotplot with specific pathways - all pathways
p.res <- fgsea.res.combined %>% 
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  merge(., fgsea.blocks.df, by.x = "pathway", by.y = "pathway", all.x = F, all.y = T) %>% 
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  mutate(block = factor(block, levels = c("Inflammation", "Differentiation", "Survival"))) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  # mutate(NES = if_else(padj < 0.1, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = cell.type.list), y = pathway.name, size = -log10(padj), color = NES)) +
  scale_colour_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], na.value = "grey80") +
  geom_point(stat = "identity") +
  scale_size_continuous(range = c(3,8)) +
  # geom_label(aes(label = sig), size = 3) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(block), scales = "free", space = "free_y") +
  xlab("") + ylab("")
p.res
ggsave("figures/current_figure_drafts/RNA_gsea_dot_plot_20240227.pdf", plot = p.res, device = "pdf", dpi = 300, width = 4.5, height = 4.5, unit = "in")

## Plot the fgsea dotplot with specific pathways - stars for significance
p.res <- fgsea.res.combined %>% 
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  merge(., fgsea.blocks.df, by.x = "pathway", by.y = "pathway", all.x = F, all.y = T) %>% 
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  mutate(block = factor(block, levels = c("Inflammation", "Differentiation", "Survival"))) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  mutate(sig = if_else(padj < 0.1, "yes", NA_character_)) %>%
  ggplot(aes(x = factor(cluster, levels = cell.type.list), y = pathway.name, size = -log10(padj), color = NES)) +
  scale_colour_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], na.value = "grey80") +
  geom_point(stat = "identity") +
  scale_size_continuous(range = c(3,8)) +
  geom_point(data = . %>%filter(sig == "yes"), mapping = aes(x = factor(cluster, levels = cell.type.list)), color = "black", shape = 8, size = 2) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(block), scales = "free", space = "free_y") +
  xlab("") + ylab("")
p.res
ggsave("figures/current_figure_drafts/RNA_gsea_dot_plot_stars_significance_0_1_20240227.pdf", plot = p.res, device = "pdf", dpi = 300, width = 4.5, height = 4.5, unit = "in")

## Plot the fgsea dotplot with specific pathways - grey for non-significant pathways
p.res <- fgsea.res.combined %>% 
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  merge(., fgsea.blocks.df, by.x = "pathway", by.y = "pathway", all.x = F, all.y = T) %>% 
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  mutate(block = factor(block, levels = c("Inflammation", "Differentiation", "Survival"))) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  mutate(NES = if_else(padj < 0.1, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = cell.type.list), y = pathway.name, size = -log10(padj), color = NES)) +
  scale_colour_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], na.value = "grey80") +
  geom_point(stat = "identity") +
  scale_size_continuous(range = c(3,8)) +
  # geom_point(data = . %>%filter(sig == "yes"), mapping = aes(x = factor(cluster, levels = cell.type.list)), color = "black", shape = 8, size = 2) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(rows = vars(block), scales = "free", space = "free_y") +
  xlab("") + ylab("")
p.res
ggsave("figures/current_figure_drafts/RNA_gsea_dot_plot_significance_grey_20240227.pdf", plot = p.res, device = "pdf", dpi = 300, width = 4.5, height = 4.5, unit = "in")

######################################### Extra ##############################

## Plot the fgsea dotplot - all results
fgsea.res.combined %>% 
  mutate(pathway.name = gsub("_", " ", pathway)) %>% 
  mutate(pathway.name = str_to_title(pathway.name)) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>%
  # mutate(cluster = factor(cluster, levels = unique(cluster))) %>% 
  mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSC", "LMPP", "EMP", "Prog Mk", "GMP")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_color_gradient2(low = "blue", high = "red", mid = "grey80", na.value = "grey") +
  geom_point(stat = "identity") +
  # facet_grid(rows = vars(category), scales = "free") +
  theme_classic() +
  # theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() +
  xlab("") + ylab("")