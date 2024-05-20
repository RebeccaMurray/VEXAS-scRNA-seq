library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(clusterProfiler)

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")


celltype.list <- c("HSC","MkP", "LMPP", "EMP", "CD14 Mono", "B", "CD4 T", "CD8 T")
names(celltype.list) <- celltype.list
de.results <- lapply(celltype.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/", x, ".csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

############################### Plot volcanoes for specific clusters ###################

plot_volcano(de.results$HSC, stat_cutoff = 0.1, effect_cutoff = 0.2, effect_line = 0.2, stat_column = "fdr", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 20), xlim = c(-1, 1))

plot_volcano(de.results$LMPP, effect_cutoff = 0.2, effect_line = 0.2, stat_column = "pval") + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 15), xlim = c(-1, 1))

plot_volcano(de.results$MkP, metadata.df = rna.obj@meta.data %>% filter(cluster_celltype == "MkP") %>% filter(Genotype_Approx_Match_No_Gene_Thresh %in% c("MUT", "WT")), effect_cutoff = 0.2, effect_line = 0.2, stat_column = "pval", title = paste0("MkPs"), 
             genotyping_column = "Genotype_Approx_Match_No_Gene_Thresh") + 
  theme(legend.position = "right") +
  coord_cartesian(ylim=c(0, 10))

plot_volcano(de.results$`CD14 Mono`, stat_cutoff = 0.1, effect_cutoff = 0.2, effect_line = 0.2, stat_column = "pval") + 
  theme(legend.position = "right")

plot_volcano(de.results$`B`, stat_cutoff = 0.1, effect_cutoff = 0.2, effect_line = 0.2, stat_column = "pval") + 
  theme(legend.position = "right")

plot_volcano(de.results$`CD8 T`, stat_cutoff = 0.1, effect_cutoff = 0.2, effect_line = 0.2, stat_column = "pval") + 
  theme(legend.position = "right")



################################ GSEA for specific clusters #####################

fgsea.out <- run_fgsea_for_list(de.results)
fgsea.out.h <- run_fgsea_for_list(de.results[["HSC"]], msigdb.list = "HALLMARK")
fgsea.out <- run_fgsea_for_list(de.results, msigdb.list = "GO:BP", min.size = 0)
fgsea.out <- run_fgsea_for_list(de.results[["HSC"]], msigdb.list = "REACTOME")
fgsea.out.h <- run_fgsea_for_list(de.results[["B"]])
fgsea.out.h <- run_fgsea_for_list(de.results[["CD4 T"]])
fgsea.out.h <- run_fgsea_for_list(de.results[["CD8 T"]])

run_fgsea_for_list(de.results[["HSC"]], msigdb.list = "HALLMARK")
run_fgsea_for_list(de.results[["LMPP"]])

fgsea.out <- run_fgsea_for_list(de.results$`CD14 Mono`, msigdb.list = "REACTOME")

################################ GSEA for UPR #####################

## PERK pathways
GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE
REACTOME_PERK_REGULATES_GENE_EXPRESSION

## XBP1
GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE

## ATF6
ATF6_01

p.1 <- plot_enrichment(de.list = de.results, pathway.name = "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", msigdb.list = "GO:BP")
p.2 <- plot_enrichment(de.list = de.results, pathway.name = "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", msigdb.list = "GO:BP")
p.3 <- plot_enrichment(de.list = de.results, pathway.name = "GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",  msigdb.list = "GO:BP")
plot_enrichment(de.list = de.results, pathway.name = "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",  msigdb.list = "REACTOME")
plot_enrichment(de.list = de.results, pathway.name = "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE_BY_P53_CLASS_MEDIATOR",  msigdb.list = "GO")
p.1 / p.2 / p.3

################################ Making enrichment plots #######################

p.inf <- plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "HALLMARK_INTERFERON_ALPHA_RESPONSE")
p.atf4 <- plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "ATF4_ONLY_HAN", pval_show = "pval")
p.chop <- plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "CHOP_ONLY_HAN", pval_show = "pval")
p.p53 <- plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "HALLMARK_P53_PATHWAY")
p.translation <- plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "REACTOME_TRANSLATION")
plot_enrichment(de.list = de.results[["HSC"]], pathway.name = "ATF4_ONLY_HAN", msigdb.list = "HALLMARK")

p.inf / p.atf4 / p.p53 / p.translation
p.p53
p.translation

p.inf.mkp <- plot_enrichment(de.list = de.results[["MkP"]], pathway.name = "HALLMARK_INTERFERON_ALPHA_RESPONSE")

############################################# Saving genes that correspond to pathways #######################################

## Identify genes upregulated in MUT cells that are part of ROS pathway
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
ros.pathway <- pathways[["GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES"]]
de.results[["HSC"]] %>% filter(pval < 0.05 & avg_log2FC > 0) %>% filter(feature %in% ros.pathway) %>% arrange(pval) %>% 
  write_csv("data/differential_expression/VEXAS_MUT_vs_WT/backups/backups_saved_06_14_23/gene_sets/HSCs_up_MUT_GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES.csv")
ros.pathway <- pathways[["GOBP_RESPONSE_TO_OXIDATIVE_STRESS"]]
de.results[["HSC"]] %>% filter(pval < 0.05 & avg_log2FC > 0) %>% filter(feature %in% ros.pathway) %>% arrange(pval) %>%
  write_csv("data/differential_expression/VEXAS_MUT_vs_WT/backups/backups_saved_06_14_23/gene_sets/HSCs_up_MUT_GOBP_RESPONSE_TO_OXIDATIVE_STRESS.csv")
chaperone.pathway <- pathways[["GOBP_PROTEIN_FOLDING"]]
de.results[["HSC"]] %>% filter(pval < 0.05 & avg_log2FC > 0) %>% filter(feature %in% chaperone.pathway) %>% arrange(pval) %>%
  write_csv("data/differential_expression/VEXAS_MUT_vs_WT/backups/backups_saved_06_14_23/gene_sets/HSCs_up_MUT_GOBP_PROTEIN_FOLDING.csv")
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name)
cytokine.pathway <- pathways$BIOCARTA_CYTOKINE_PATHWAY
de.results[["HSC"]] %>% filter(feature %in% cytokine.pathway) %>% arrange(pval)
de.results[["LMPP"]] %>% filter(feature %in% cytokine.pathway) %>% arrange(pval)
de.results[["EMP"]] %>% filter(feature %in% cytokine.pathway) %>% arrange(pval)
de.results[["Prog Mk"]] %>% filter(feature %in% cytokine.pathway) %>% arrange(pval)