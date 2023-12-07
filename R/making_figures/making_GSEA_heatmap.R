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
cell.type.list <- c("HSC","EMP", "LMPP", "MkP", "GMP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/", x, "_5cells_genotyped_20231206.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

########################################## Make a heatmap with hallmark GSEA sets #####################################

fgsea.res.all <- lapply(names(de.results), function(x) { 
  print("Running for:")
  print(x)
  fgsea.res <- run_fgsea_for_list(de.results[[x]], genes.tested = rna.obj@assays$RNA@data %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., value = T, invert = T))
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
}
)

fgsea.res.combined <- do.call(rbind, fgsea.res.all)
fgsea.res.combined %>% write_csv("figures/current_figure_drafts/fgsea_used_for_heatmap.csv")

## order
fgsea.pathways = list(
  `CIITA_TARGET_GENES` = "Inflammation",
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Inflammation",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Inflammation",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "Inflammation",
  `REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM` = "Inflammation",
  `G2M phase module` = "Differentiation",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` = "Differentiation",
  `GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS` = "Differentiation",
  `ATF4 targets (Han et al.)` = "Survival",
  `GOBP_PROTEIN_FOLDING` = "Survival",
  `GOBP_PROTEIN_REFOLDING` = "Survival",
  `HALLMARK_P53_PATHWAY` = "Survival",
  `KEGG_PROTEASOME` = "Survival",
  `KEGG_RNA_DEGRADATION` = "Survival",
  `REACTOME_AUTOPHAGY` = "Survival",
  `REACTOME_APOPTOSIS` = "Survival",
  `REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY` = "Survival",
  `ATF6 module` = "Survival",
  `CHOP targets (Han et al.)` = "Survival",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "Survival",
  `REACTOME_TRANSLATION` = "Survival",
  `HALLMARK_GLYCOLYSIS` = "Other",
  `HALLMARK_MTORC1_SIGNALING` = "Other",
  `HALLMARK_MYC_TARGETS_V1` = "Other",
  `KEGG_SPLICEOSOME` = "Other",
  `S phase module` = "Other"
)

fgsea.pathway.list.formatted <- str_to_title(names(fgsea.pathways)) %>% str_wrap(., width = 25) %>% gsub("_", " ", .)

## Make a GSEA heatmap - pval
m.gsea <- fgsea.res.combined %>% 
  group_by(pathway) %>% 
  filter(sum(padj < 0.2) > 1) %>% arrange(pathway) %>% 
  filter() %>% 
  mutate(plot_value = -log10(padj) * sign(NES)) %>% 
  select(pathway, plot_value, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "plot_value") %>% 
  column_to_rownames("pathway")
pheatmap::pheatmap(m.gsea, 
                   cluster_rows = T,
                   color = colorRampPalette(c("navy", "white", "red"))(20), breaks = seq(-3, 3, length.out = 20))

## Make a GSEA heatmap - picking pathways, color is pval
m.gsea <- fgsea.res.combined %>% 
  filter(pathway %in% names(fgsea.pathways)) %>% 
  mutate(plot_value = -log10(padj) * sign(NES)) %>% 
  select(pathway, plot_value, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "plot_value") %>% 
  column_to_rownames("pathway")
rownames(m.gsea) <- plyr::mapvalues(rownames(m.gsea), from = names(fgsea.pathways), to = fgsea.pathway.list.formatted)
pdf("figures/current_figure_drafts/RNA_GSEA_heatmap_colors_v1.pdf", width = 5.2, height = 5)
pheatmap::pheatmap(m.gsea[fgsea.pathway.list.formatted, ],
                   gaps_row = c(5, 8, 20),
                   cluster_rows = F,
                   cluster_cols = F,
                   # color = colorRampPalette(c("navy", "white", "red"))(20), 
                   color = colorRampPalette(c("#2A9098", "white", "#F5B935"))(20),
                   # color = colorRampPalette(c("#FF7B0B", "white", "red"))(20), 
                   breaks = seq(-3, 3, length.out = 20))
dev.off()

## Make a GSEA heatmap - NES, targets we pick 
fgsea.pathways.specific = list(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "REACTOME_MTOR_SIGNALLING",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_LEUKOCYTE_DIFFERENTIATION",
  "GOBP_MYELOID_CELL_DIFFERENTIATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "REACTOME_TRANSLATION",
  "GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "CHOP targets (Han et al.)",
  "ATF4 targets (Han et al.)",
  "Shared ATF4/CHOP targets (Han et al.)",
  "REACTOME_AUTOPHAGY",
)
## Format as a data frame
fgsea.pathways.df <- data.frame(pathway = names(fgsea.pathways.specific), pathway.name = unlist(fgsea.pathways.specific))

m.gsea <- fgsea.res.combined %>% 
  group_by(pathway) %>% 
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  # mutate(NES = ifelse(padj > 0.25, NA, NES)) %>% 
  dplyr::select(pathway.name, NES, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "NES") %>% 
  column_to_rownames("pathway.name")
pheatmap::pheatmap(m.gsea[fgsea.pathways.df$pathway.name, ], 
                   cluster_rows = F,
                   gaps_row = c(4, 8, 11, 15),
                   color = colorRampPalette(c("navy", "white", "red"))(20), breaks = seq(-2, 2, length.out = 20))


################################### Extra #################################

## Make a GSEA heatmap - NES, targets we pick 
fgsea.pathways.specific = list(
  `HALLMARK_INTERFERON_ALPHA_RESPONSE` = "Hallmark interferon alpha response",
  `HALLMARK_INTERFERON_GAMMA_RESPONSE` = "Hallmark interferon gamma response",
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` = "Hallmark TNFa signaling via NFkB",
  `HALLMARK_TGF_BETA_SIGNALING` = "Hallmark TGF-beta signaling",
  `HALLMARK_P53_PATHWAY` = "Hallmark P53 pathway",
  `REACTOME_MTOR_SIGNALLING` = "Reactome MTOR signaling",
  `KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION` = "KEGG antigen processing and presentation",
  `GOBP_LEUKOCYTE_DIFFERENTIATION` = "GOBP leukocyte differentiation",
  `GOBP_MYELOID_CELL_DIFFERENTIATION` = "GOBP myeloid cell differentiation",
  `KEGG_OXIDATIVE_PHOSPHORYLATION` = "KEGG oxidative phosphorylation",
  `HALLMARK_GLYCOLYSIS` = "Hallmark glycolysis",
  `GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES` = "GOBP ROS",
  `HALLMARK_G2M_CHECKPOINT` = "Hallmark G2M checkpoint",
  `HALLMARK_MYC_TARGETS_V1` = "Hallmark MYC targets v1",
  `HALLMARK_MYC_TARGETS_V2` = "Hallmark MYC targets v2",
  `REACTOME_TRANSLATION` = "Reactome translation",
  `GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` =   "GOBP UPR - IRE1",
  # `GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` = "GOBP UPR - PERK",
  # `GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE` =   "GOBP UPR - ATF6",
  `CHOP targets (Han et al.)` = "CHOP targets (Han et al.)",
  `ATF4 targets (Han et al.)` = "ATF4 targets (Han et al.)",
  `Shared ATF4/CHOP targets (Han et al.)` = "Shared ATF4/CHOP targets (Han et al.)",
  `REACTOME_AUTOPHAGY` = "Reactome autophagy"
)
## Format as a data frame
fgsea.pathways.df <- data.frame(pathway = names(fgsea.pathways.specific), pathway.name = unlist(fgsea.pathways.specific))

m.gsea <- fgsea.res.combined %>% 
  group_by(pathway) %>% 
  filter(pathway %in% names(fgsea.pathways.specific)) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>%
  mutate(pathway.name = factor(pathway.name, levels = rev(unlist(fgsea.pathways.specific)))) %>% 
  # mutate(NES = ifelse(padj > 0.25, NA, NES)) %>% 
  dplyr::select(pathway.name, NES, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "NES") %>% 
  column_to_rownames("pathway.name")
pheatmap::pheatmap(m.gsea[fgsea.pathways.df$pathway.name, ], 
                   cluster_rows = F,
                   gaps_row = c(4, 8, 11, 15),
                   color = colorRampPalette(c("navy", "white", "red"))(20), breaks = seq(-2, 2, length.out = 20))

## Make a GSEA heatmap - NES
m.gsea <- fgsea.res.combined %>% 
  group_by(pathway) %>% 
  filter(sum(padj < 0.2) > 1) %>% arrange(pathway) %>% 
  dplyr::select(pathway, NES, cluster) %>% 
  pivot_wider(names_from = "cluster", values_from = "NES") %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway) %>% gsub("_", " ", .) %>% str_to_sentence(.)) %>% 
  column_to_rownames("pathway")
pheatmap::pheatmap(m.gsea, 
                   cluster_rows = T,
                   color = colorRampPalette(c("navy", "white", "red"))(20), breaks = seq(-2, 2, length.out = 20))

################# Looking for overlapping genes ####################

de.results.combined <- do.call(rbind, de.results)
de.results.combined %>% 
  filter(fdr < 0.05) %>% 
  group_by(feature) %>% filter(n() == 3)


gene.dotplot <- de.results.combined %>% 
  filter(fdr < 0.2) %>% 
  group_by(feature) %>% filter(n() == 4) %>% 
  mutate(mean_val = mean(avg_log2FC)) %>% 
  ggplot(aes(x = factor(cluster, levels = rev(c("HSC", "LMPP", "EMP", "MkP"))), y = reorder(feature, mean_val), size = -log10(fdr), color = avg_log2FC)) +
  scale_color_gradient2(low = vexas.genotyping.palette[["WT"]], high = vexas.genotyping.palette[["MUT"]], mid = "grey", na.value = "grey") +
  geom_point(stat = "identity") +
  # facet_grid(rows = vars(category), scales = "free") +
  theme_classic() +
  # theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_flip() +
  xlab("") + ylab("")
gene.dotplot
ggsave(paste0(current.plots.path, "/RNA_genes_overlapping_de.csv"), plot = gene.dotplot, dpi = 300, device = "pdf", width = 5, height = 3)


################## Hypergeometric test for HSCs #####################

gene.list.up <- de.results[["HSC"]] %>% filter(pval < 0.05 & avg_log2FC > 0.2) %>% pull(feature)
gene.list.down <- de.results[["HSC"]] %>% filter(pval < 0.05 & avg_log2FC < 0.2) %>% pull(feature)
universe <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., invert = T, value = T)

ego <- enrichGO(gene = gene.list.up,
  universe = universe,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  minGSSize = 30,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)
ego@result %>% arrange(pvalue) %>% head(20)

## Plot the results
p.res <- ego@result %>% 
  slice_min(pvalue, n = 10) %>% 
  mutate(stat_plot = -log10(p.adjust)) %>% 
  # mutate(Description = str_wrap(Description, width = 30)) %>% 
  ggplot(aes(x = reorder(Description, stat_plot), y = stat_plot)) +
  geom_col(fill = vexas.genotyping.palette[["MUT"]], color = "black") +
  coord_flip() +
  ylab("-log10(p.adjust)") +
  xlab(element_blank()) +
  ggtitle("GO:CC terms")
p.res
ggsave(paste0(current.plots.path, "/HSCs_GOCC_bar_plot_hypergeometric.pdf"), plot = p.res, device = "pdf", dpi = 300, width = 4, height = 2.5, unit = "in")



################# Bar plot with top HSC pathways ###############

hsc.fgsea.reactome <- run_fgsea_for_list(de.results$HSC, msigdb.list = "KEGG", min.size = 20)
hsc.fgsea.reactome %>% filter(pval < 0.05 & NES > 0) %>% arrange(pval) %>% head(n = 10)

hsc.fgsea.reactome <- run_fgsea_for_list(de.results$HSC, msigdb.list = "REACTOME", min.size = 20)
hsc.fgsea.reactome %>% filter(pval < 0.05 & NES > 0) %>% arrange(pval) %>% head(n = 10)

hsc.fgsea <- run_fgsea_for_list(de.results$HSC, msigdb.list = "GO:BP", min.size = 20)
hsc.fgsea %>% filter(pval < 0.05 & NES > 0) %>% arrange(pval) %>% head(n = 10)
hsc.fgsea %>% filter(pval < 0.05 & NES > 0) %>% arrange(-NES) %>% head(n = 10)

hsc.fgsea <- run_fgsea_for_list(de.results$HSC, msigdb.list = "GO:CC", min.size = 20)
hsc.fgsea %>% filter(pval < 0.05 & NES > 0) %>% arrange(pval) %>% head(n = 20)

hsc.fgsea.tf <- run_fgsea_for_list(de.results$HSC, msigdb.list = "TF", min.size = 20)
hsc.fgsea.tf %>% filter(pval < 0.05 & NES < 0) %>% arrange(pval) %>% head(n = 20)

hsc.fgsea.all <- run_fgsea_for_list(de.results$HSC, msigdb.list = "ALL", min.size = 20)
hsc.fgsea.all %>% filter(pval < 0.05 & NES < 0) %>% arrange(pval) %>% head(n = 20)
