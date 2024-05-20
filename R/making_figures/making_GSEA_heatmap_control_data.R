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


# list.files("/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/", pattern = ".csv")
# 
# ## Loading the fgsea results - ran by Theo
# cell.type.list <- c("HSC","EMP", "LMPP", "MkP")
# names(cell.type.list) <- cell.type.list
# fgsea.res.all <- lapply(cell.type.list, function (x) {
#   read_table(paste0("/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/fgseaRes_VEXAS_CONTROL_", x, ".csv")) %>% mutate(cluster = x)
# })

## Loading the DE results - ran by RMM
cell.type.list <- c("HSC","EMP", "LMPP", "MkP", "GMP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_vs_control/vexas.vs.control.", x, ".DGE.csv")) %>% dplyr::rename("feature" = "...1")
})


########################################## Make a heatmap with hallmark GSEA sets #####################################

fgsea.res.all <- lapply(names(de.results), function(x) { 
  print("Running for:")
  print(x)
  
  ## Try re-ranking
  de_results.ranked <- de.results[[x]] %>%
    filter(!is.na(p_val)) %>% 
    # mutate(rank = -log10(pval)*sign(avg_log2FC)) %>%
    mutate(rank = avg_log2FC) %>% 
    # mutate(rank = -log10(pval)) %>%
    arrange(-rank)
  current.ranks.list <- de_results.ranked$rank
  names(current.ranks.list) <- de_results.ranked$feature
  
  fgsea.res <- run_fgsea_for_list(de.results, ranks.list = current.ranks.list, nPermSimple = 100000)
  # fgsea.res <- run_fgsea_for_list(de.results[[x]], pval_col = "p_val")
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
}
)

fgsea.res.combined <- do.call(rbind, fgsea.res.all)
fgsea.res.combined %>% write_csv("figures/current_figure_drafts/fgsea_used_for_heatmap_controls_avglog2FC_ranking_2024-02-26.csv")

## Make a GSEA heatmap - all pathways, filtered by pval
m.gsea <- fgsea.res.combined %>% 
  group_by(pathway) %>% 
  filter(sum(padj < 0.2) > 1) %>% arrange(pathway) %>% 
  filter() %>% 
  mutate(plot_value = -log10(padj) * sign(NES)) %>% 
  dplyr::select(pathway, plot_value, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "plot_value") %>% 
  column_to_rownames("pathway")
pheatmap::pheatmap(m.gsea, 
                   cluster_rows = T,
                   color = colorRampPalette(c("navy", "white", "red"))(20), breaks = seq(-3, 3, length.out = 20))

## Order the pathways
fgsea.pathways = list(
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Inflammation",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Inflammation",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "Inflammation",
  # `CIITA_TARGET_GENES` = "Inflammation",
  `HALLMARK_INFLAMMATORY_RESPONSE` = "Inflammation",
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` = "Inflammation",
  `REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM` = "Inflammation",
  `GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS` = "Survival",
  `ATF4 targets (Han et al.)` = "Survival",
  `CHOP targets (Han et al.)` = "Survival",
  `REACTOME_TRANSLATION` = "Survival",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "Survival",
  `GOBP_PROTEIN_FOLDING` = "Survival",
  `ATF6 module` = "Survival",
  `S phase module` = "Differentiation",
  `G2M phase module` = "Differentiation",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` = "Differentiation",
  `HALLMARK_P53_PATHWAY` = "Survival",
  `KEGG_PROTEASOME` = "Survival",
  `KEGG_RNA_DEGRADATION` = "Survival",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Survival",
  `HALLMARK_APOPTOSIS` = "Survival",
  # `Anti-apoptosis module` = "Survival",
  `HALLMARK_GLYCOLYSIS` = "Other",
  `HALLMARK_MTORC1_SIGNALING` = "Other",
  `HALLMARK_MYC_TARGETS_V1` = "Other",
  `KEGG_SPLICEOSOME` = "Other"
)

fgsea.pathway.list.formatted <- gsub("_", " ", names(fgsea.pathways)) %>% str_to_title(.) %>% 
  gsub("Gobp ", "", .) %>% 
  gsub("Kegg ", "", .) %>% 
  gsub("Reactome", "", .) %>% 
  gsub("Hallmark", "", .) %>% 
  gsub("Atf4", "ATF4", .) %>% 
  gsub("Ciita", "CIITA", .) %>% 
  gsub("Et Al.", "et al.", .) %>% 
  gsub("Mtorc1", "MTORC1", .) %>% 
  gsub("G2m", "G2M", .) %>% 
  gsub("Atf6", "ATF6", .) %>% 
  gsub("Myc", "MYC", .) %>% 
  gsub("Chop", "CHOP", .) %>% 
  gsub("Rna", "RNA", .) %>% 
  gsub("Of Ire1", "of IRE1", .) %>% 
  gsub("V1", "v1", .) %>% 
  gsub("Endoplasmic Reticulum", "ER", .) %>% 
  gsub("Tnfa", "TNFa", .) %>% 
  gsub("Nfkb", "NFkB", .) %>% 
  gsub("Via", "via", .) %>% 
  gsub("Unfolded Protein Response", "UPR", .) %>% 
  gsub("^ ", "", .)
# str_wrap(., width = 25) 

## Make a GSEA heatmap - picking pathways, color is pval
m.gsea <- fgsea.res.combined %>% 
  dplyr::filter(pathway %in% names(fgsea.pathways)) %>% 
  mutate(plot_value = -log10(padj) * sign(NES)) %>% 
  dplyr::select(pathway, plot_value, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "plot_value") %>% 
  column_to_rownames("pathway")
rownames(m.gsea) <- plyr::mapvalues(rownames(m.gsea), from = names(fgsea.pathways), to = fgsea.pathway.list.formatted)
max.val = round(mean(abs(as.matrix(m.gsea))) + 1*sd(abs(as.matrix(m.gsea))))
min.val <- -1*max.val
dev.off()
pdf("figures/current_figure_drafts/RNA_GSEA_heatmap_colors_fold_change_controls_20240224.pdf",  width = 3.8, height = 5)
pheatmap::pheatmap(m.gsea[fgsea.pathway.list.formatted, ],
                   gaps_row = c(6, 13),
                   cluster_rows = F,
                   cluster_cols = F,
                   # color = colorRampPalette(c("navy", "white", "red"))(20), 
                   color = colorRampPalette(c("#2A9098", "white", "#F5B935"))(31),
                   # color = colorRampPalette(c("#FF7B0B", "white", "red"))(20), 
                   breaks = seq(min.val, max.val, length.out = 31))
dev.off()



