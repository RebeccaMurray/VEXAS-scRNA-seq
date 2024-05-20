library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(clusterProfiler)

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")


sample.list <- c("HSC", "LMPP", "EMP", "Prog Mk")
names(sample.list) <- sample.list
de.results <- lapply(sample.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells/", x, "_no_mito_no_ribo.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})
de.results.renormalized <- lapply(sample.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT/", x, "_no_mito_no_ribo_renormalized.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

## HSCs
p1 <- plot_volcano(de.results$HSC, metadata.df = NA, effect_line = 0.4, stat_column = "fdr", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 15), xlim = c(-1, 1)) + ggtitle("Mito/Ribo excluded", subtitle = "HSCs")

p2 <- plot_volcano(de.results.renormalized$HSC, metadata.df = NA, effect_line = 0.4, stat_column = "fdr", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 15), xlim = c(-1, 1)) + ggtitle("Mito/Ribo excluded and\ncounts renormalized", subtitle = "HSCs")

p1 | p2

## Prog Mks
p1 <- plot_volcano(de.results$`Prog Mk`, metadata.df = NA, effect_line = 0.4, stat_column = "pval", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 10), xlim = c(-1.5, 1.5)) + ggtitle("Mito/Ribo excluded", subtitle = "Prog Mks")

p2 <- plot_volcano(de.results.renormalized$`Prog Mk`, metadata.df = NA, effect_line = 0.4, stat_column = "pval", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 10), xlim = c(-1.5, 1.5)) + ggtitle("Mito/Ribo excluded and\ncounts renormalized", subtitle = "Prog Mks")

p1 | p2


output <- run_fgsea_for_list(de.results$HSC)
output <- run_fgsea_for_list(de.results.renormalized$HSC)

output <- run_fgsea_for_list(de.results$`Prog Mk`, msigdb.list = "HALLMARK")
output <- run_fgsea_for_list(de.results.renormalized$`Prog Mk`, msigdb.list = "HALLMARK")

output <- run_fgsea_for_list(de.results$`LMPP`, msigdb.list = "HALLMARK")
output <- run_fgsea_for_list(de.results.renormalized$`LMPP`, msigdb.list = "HALLMARK")



###################################### Making gsea dotplot - regular ##########################

fgsea.res.all <- lapply(names(de.results), function(x) { 
  fgsea.res <- run_fgsea_for_list(de.results[[x]])
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
}
)

fgsea.res.combined <- do.call(rbind, fgsea.res.all)

## Isolate and name pathways for a more curated dot plot
fgsea.pathways.specific = list(
  `HALLMARK_INFLAMMATORY_RESPONSE` = "Hallmark inflammatory response",
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Hallmark interferon alpha response",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Hallmark interferon gamma response",
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` = "Hallmark TNFa signaling via NFkB",
  `HALLMARK_MTORC1_SIGNALING` = "Hallmark mTOR signaling",
  `HALLMARK_IL6_JAK_STAT3_SIGNALING` = "Hallmark IL6 JAK STAT3 signaling",
  # `HALLMARK_TGF_BETA_SIGNALING` = "Hallmark TGF-beta signaling",
  `HALLMARK_GLYCOLYSIS` = "Hallmark glycolysis",
  `KEGG_OXIDATIVE_PHOSPHORYLATION` = "KEGG oxidative phosphorylation",
  `HALLMARK_G2M_CHECKPOINT` = "Hallmark G2M checkpoint",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "KEGG antigen processing and presentation",
  `REACTOME_TRANSLATION` = "Reactome translation",
  `HALLMARK_MYC_TARGETS_V2` = "Hallmark MYC targets v2",
  `HALLMARK_APOPTOSIS` = "Hallmark apoptosis",
  `HALLMARK_P53_PATHWAY` = "Hallmark P53 pathway",
  # `GOBP_LEUKOCYTE_DIFFERENTIATION` = "GOBP leukocyte differentiation",
  # `GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION` = "GOBP regulation of myeloid differentiation",
  `GOBP_MYELOID_CELL_DIFFERENTIATION` = "GOBP myeloid cell differentiation",
  `REACTOME_MTOR_SIGNALLING` = "Reactome mTOR signaling",
  `GOCC_VESICLE_MEMBRANE` = "GOCC vesicle membrane",
  `GOCC_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE` = "GOCC intrinsic component of ER membrane",
  `GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS` = "GOBP response to ER stress",
  `CHOP_AND_ATF4_HAN` = "CHOP and ATF4 target genes [Han et al.]",
  `CHOP_ONLY_HAN` = "CHOP target genes [Han et al.]",
  `ATF4_ONLY_HAN` = "ATF4 target genes [Han et al.]",
  # `PERK module` = "PERK module",
  # `ATF6 module` = "ATF6 module",
  # `XBP1 target module` = "XBP1 module",
  # `KEGG_RNA_DEGRADATION` = "KEGG RNA degradation",
  # `REACTOME_AUTOPHAGY` = "Reactome autophagy",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Reactome chaperone-mediated autophagy"
)

## Format as a data frame
fgsea.pathways.df <- data.frame(pathway = names(fgsea.pathways.specific), pathway.name = unlist(fgsea.pathways.specific))

## Plot the fgsea dotplot with specific pathways
fgsea.res.combined %>% 
  filter(cluster != "GMP") %>%
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  mutate(pathway.name = str_wrap(pathway.name, width = 50)) %>% 
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSC", "LMPP", "EMP", "Prog Mk")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_colour_gradient(low = "blue", high = "red2", na.value = "grey80") +
  geom_point(stat = "identity") +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("") + ylab("")

###################################### Making gsea dotplot - renormalized ##########################

fgsea.res.all <- lapply(names(de.results.renormalized), function(x) { 
  fgsea.res <- run_fgsea_for_list(de.results.renormalized[[x]], msigdb.list = "HALLMARK")
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
}
)

fgsea.res.combined <- do.call(rbind, fgsea.res.all)

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

## Isolate and name pathways for a more curated dot plot
fgsea.pathways.specific = list(
  `HALLMARK_INFLAMMATORY_RESPONSE` = "Hallmark inflammatory response",
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Hallmark interferon alpha response",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Hallmark interferon gamma response",
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` = "Hallmark TNFa signaling via NFkB",
  `HALLMARK_MTORC1_SIGNALING` = "Hallmark mTOR signaling",
  `HALLMARK_IL6_JAK_STAT3_SIGNALING` = "Hallmark IL6 JAK STAT3 signaling",
  # `HALLMARK_TGF_BETA_SIGNALING` = "Hallmark TGF-beta signaling",
  `HALLMARK_GLYCOLYSIS` = "Hallmark glycolysis",
  `KEGG_OXIDATIVE_PHOSPHORYLATION` = "KEGG oxidative phosphorylation",
  `HALLMARK_G2M_CHECKPOINT` = "Hallmark G2M checkpoint",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "KEGG antigen processing and presentation",
  `REACTOME_TRANSLATION` = "Reactome translation",
  `HALLMARK_MYC_TARGETS_V2` = "Hallmark MYC targets v2",
  `HALLMARK_APOPTOSIS` = "Hallmark apoptosis",
  `HALLMARK_P53_PATHWAY` = "Hallmark P53 pathway",
  # `GOBP_LEUKOCYTE_DIFFERENTIATION` = "GOBP leukocyte differentiation",
  # `GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION` = "GOBP regulation of myeloid differentiation",
  `GOBP_MYELOID_CELL_DIFFERENTIATION` = "GOBP myeloid cell differentiation",
  `REACTOME_MTOR_SIGNALLING` = "Reactome mTOR signaling",
  `GOCC_VESICLE_MEMBRANE` = "GOCC vesicle membrane",
  `GOCC_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE` = "GOCC intrinsic component of ER membrane",
  `GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS` = "GOBP response to ER stress",
  `CHOP_AND_ATF4_HAN` = "CHOP and ATF4 target genes [Han et al.]",
  `CHOP_ONLY_HAN` = "CHOP target genes [Han et al.]",
  `ATF4_ONLY_HAN` = "ATF4 target genes [Han et al.]",
  # `PERK module` = "PERK module",
  # `ATF6 module` = "ATF6 module",
  # `XBP1 target module` = "XBP1 module",
  # `KEGG_RNA_DEGRADATION` = "KEGG RNA degradation",
  # `REACTOME_AUTOPHAGY` = "Reactome autophagy",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Reactome chaperone-mediated autophagy"
)

## Format as a data frame
fgsea.pathways.df <- data.frame(pathway = names(fgsea.pathways.specific), pathway.name = unlist(fgsea.pathways.specific))

## Plot the fgsea dotplot with specific pathways
fgsea.res.combined %>% 
  filter(cluster != "GMP") %>%
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  mutate(pathway.name = str_wrap(pathway.name, width = 50)) %>% 
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSC", "LMPP", "EMP", "Prog Mk")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_colour_gradient(low = "blue", high = "red2", na.value = "grey80") +
  geom_point(stat = "identity") +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("") + ylab("")
