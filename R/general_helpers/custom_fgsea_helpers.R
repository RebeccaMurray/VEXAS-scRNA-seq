library(fgsea)
library(msigdbr)

run_fgsea_for_list <- function(de.list = data.frame(), ranks.list = list(), pval_col = "pval", genes.tested = list(), msigdb.list = "none", min.size = 0, nPermSimple = 1000) {
  
  if (length(ranks.list) == 0) {
    if (!("feature" %in% colnames(de.list))) {
      de.list$feature <- rownames(de.list)
    }
    
    de_results.ranked <- de.list %>%
      dplyr::rename(pval = pval_col) %>% 
      filter(!is.na(pval)) %>% 
      filter(is.finite(-log10(pval))) %>% 
      mutate(rank = avg_log2FC * -log10(pval)) %>%
      arrange(-rank)
    
    ranks.list <- de_results.ranked$rank
    names(ranks.list) <- de_results.ranked$feature
    
  }
  
  print(head(ranks.list))
  print(tail(ranks.list))
  if (msigdb.list == "KEGG") {
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "GO") {
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C5")
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "GO:BP") {
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "GO:CC") {
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "REACTOME") {
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") 
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "HALLMARK") {
    m_t2g <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "ALL") {
    m_t2g <- msigdbr(species = "Homo sapiens")
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else if (msigdb.list == "TF") {
    m_t2g <- tibble()
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD"))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy"))
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
  } else {
    ## Use a custom set of msigdb gene sets, if not specified
    m_t2g <- tibble()
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% filter(gs_name %in% fgsea.pathways))
    m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% filter(gs_name %in% fgsea.pathways))
    pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    ## Add GoT paper modules
    got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv") %>% as.list()
    got.modules <- lapply(got.modules, function(x) return(x[!is.na(x)]))
    pathways <- c(pathways, got.modules)
    
    ## Add Han et al modules
    han.modules <- list(
      `ATF4 targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
      `CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
      `Shared ATF4/CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_and_ATF4.txt") %>% filter(Feature != "N/A") %>% pull(Feature)
    )
    pathways <- c(pathways, han.modules)
  }
  
  ## Keep only genes tested (if specified)
  if (length(genes.tested) > 0) {
    
    # pathways <- lapply(pathways, function(pathway) grep("^RP|^MT-", pathway, value = T, invert = T))
    pathways <- lapply(pathways, function(pathway) {intersect(genes.tested, pathway)})
  }
  
  fgseaRes <- fgsea(pathways = pathways, stats = ranks.list, minSize = min.size, nPermSimple = nPermSimple)
  print("NES > 0:")
  print(fgseaRes %>% filter(NES > 0 & pval < 0.05) %>% arrange(pval))
  print("NES < 0:")
  print(fgseaRes %>% filter(NES < 0 & pval < 0.05) %>% arrange(pval))
  return(fgseaRes)
}

plot_enrichment <- function(de.list = data.frame(), ranks.list = list(), pathway.name, pval_col = "pval", pval_show = "padj", msigdb.list = "none") {
  if (length(ranks.list) == 0) {
    if (!("feature" %in% colnames(de.list))) {
      de.list$feature <- rownames(de.list)
    }
    
    de_results.ranked <- de.list %>%
      dplyr::rename(pval = pval_col) %>% 
      filter(!is.na(pval)) %>% 
      filter(is.finite(-log10(pval))) %>% 
      mutate(rank = avg_log2FC * -log10(pval)) %>%
      arrange(-rank)
    
    ranks.list <- de_results.ranked$rank
    names(ranks.list) <- de_results.ranked$feature
  }
    
    head(ranks.list)
    tail(ranks.list)
    if (msigdb.list == "KEGG") {
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    } else if (msigdb.list == "GO") {
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C5")
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    } else if (msigdb.list == "REACTOME") {
      m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") 
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    } else if (msigdb.list == "HALLMARK") {
      m_t2g <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    } else if (msigdb.list == "TF") {
      m_t2g <- tibble()
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD"))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy"))
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
    } else {
      ## Use a custom set of msigdb gene sets, if not specified
      m_t2g <- tibble()
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% filter(gs_name %in% fgsea.pathways))
      m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% filter(gs_name %in% fgsea.pathways))
      pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
      ## Add GoT paper modules
      got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv") %>% as.list()
      got.modules <- lapply(got.modules, function(x) return(x[!is.na(x)]))
      pathways <- c(pathways, got.modules)

      ## Add Han et al modules
      han.modules <- list(
        `ATF4 targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
        `CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
        `Shared ATF4/CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_and_ATF4.txt") %>% filter(Feature != "N/A") %>% pull(Feature)
      )
      pathways <- c(pathways, han.modules)
    }
  
  fgseaRes <- fgsea(pathways = pathways, stats = ranks.list)
    
  pval <- fgseaRes %>% filter(pathway == pathway.name) %>% pull(pval_show) %>% formatC(., format = "e", digits = 2)
  ES <- fgseaRes %>% filter(pathway == pathway.name) %>% pull(ES)
  y.pos <- if_else(ES > 0, ES - 0.25, ES + 0.25)
  
  p1 <- plotEnrichment(pathways[[pathway.name]], ranks.list) + annotate("text", label = paste0(pval_show, " = ", pval), x = length(ranks.list) * 0.8, y = y.pos) + 
    ggtitle(gsub("_", " ", pathway.name) %>% str_to_title(.) %>% str_wrap(., width = 30)) + theme_classic() + xlab("Rank") + ylab("Enrichment score")
  return(p1)
  
}

# hallmark.pathways <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% 
fgsea.pathways = list(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_APOPTOSIS",
  "KEGG_BASE_EXCISION_REPAIR",
  "KEGG_CALCIUM_SIGNALING_PATHWAY",
  "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  "KEGG_CELL_CYCLE",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_CITRATE_CYCLE_TCA_CYCLE",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_DNA_REPLICATION",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_LYSOSOME",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_MTOR_SIGNALING_PATHWAY",
  "KEGG_NOTCH_SIGNALING_PATHWAY",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_P53_SIGNALING_PATHWAY",
  "KEGG_PEROXISOME",
  "KEGG_PROTEASOME",
  "KEGG_REGULATION_OF_AUTOPHAGY",
  "KEGG_RNA_DEGRADATION",
  "KEGG_SPLICEOSOME",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
  "KEGG_WNT_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_LIPID_BIOSYNTHETIC_PROCESS",
  "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING",
  "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION",
  "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION",
  "GOBP_REGULATION_OF_LEUKOCYTE_DIFFERENTIATION",
  "GOBP_MYELOID_CELL_DIFFERENTIATION",
  "GOBP_LEUKOCYTE_DIFFERENTIATION",
  "GOBP_PROTEIN_FOLDING",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_PROTEIN_REFOLDING",
  "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
  "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS",
  "REACTOME_TRANSLATION",
  "REACTOME_PERK_REGULATES_GENE_EXPRESSION",
  "REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT",
  "REACTOME_AUTOPHAGY",
  "REACTOME_MTOR_SIGNALLING",
  "REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",
  "REACTOME_APOPTOSIS",
  "P53_02",
  "P53_DECAMER_Q2",
  "ATF4_Q2",
  "CHOP_01",
  "CIITA_TARGET_GENES",
  "GATA1_01",
  "GATA1_02",
  "GATA1_03",
  "GATA1_04",
  "GATA1_05",
  "MYC_Q2",
  "MYCMAX_01",
  "MYCMAX_02",
  "MYCMAX_03",
  "MYCMAX_B",
  "NMYC_01"
)

# fgsea.pathways = list(
#   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   "HALLMARK_TGF_BETA_SIGNALING",
#   "HALLMARK_FATTY_ACID_METABOLISM",
#   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#   "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
#   "HALLMARK_HYPOXIA",
#   "HALLMARK_GLYCOLYSIS",
#   "HALLMARK_P53_PATHWAY",
#   "HALLMARK_G2M_CHECKPOINT",
#   "HALLMARK_MYC_TARGETS_V1",
#   "HALLMARK_MYC_TARGETS_V2",
#   "HALLMARK_APOPTOSIS",
#   "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
#   "HALLMARK_INFLAMMATORY_RESPONSE",
#   "HALLMARK_MTORC1_SIGNALING",
#   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
#   "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
#   "KEGG_APOPTOSIS",
#   "KEGG_BASE_EXCISION_REPAIR",
#   "KEGG_CALCIUM_SIGNALING_PATHWAY",
#   "KEGG_CELL_ADHESION_MOLECULES_CAMS",
#   "KEGG_CELL_CYCLE",
#   "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
#   "KEGG_CITRATE_CYCLE_TCA_CYCLE",
#   "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
#   "KEGG_DNA_REPLICATION",
#   "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
#   "KEGG_INSULIN_SIGNALING_PATHWAY",
#   "KEGG_JAK_STAT_SIGNALING_PATHWAY",
#   "KEGG_LYSOSOME",
#   "KEGG_MAPK_SIGNALING_PATHWAY",
#   "KEGG_MTOR_SIGNALING_PATHWAY",
#   "KEGG_NOTCH_SIGNALING_PATHWAY",
#   "KEGG_OXIDATIVE_PHOSPHORYLATION",
#   "KEGG_P53_SIGNALING_PATHWAY",
#   "KEGG_PEROXISOME",
#   "KEGG_PROTEASOME",
#   "KEGG_REGULATION_OF_AUTOPHAGY",
#   "KEGG_RNA_DEGRADATION",
#   "KEGG_SPLICEOSOME",
#   "KEGG_TGF_BETA_SIGNALING_PATHWAY",
#   "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
#   "KEGG_WNT_SIGNALING_PATHWAY",
#   "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#   "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
#   "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
#   "GOBP_LIPID_BIOSYNTHETIC_PROCESS",
#   "GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
#   "GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR",
#   "GOBP_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE_BY_P53_CLASS_MEDIATOR",
#   "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING",
#   "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
#   "GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
#   "GOBP_REGULATION_OF_MYELOID_LEUKOCYTE_DIFFERENTIATION",
#   "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION",
#   "GOBP_REGULATION_OF_LEUKOCYTE_DIFFERENTIATION",
#   "GOBP_MYELOID_CELL_DIFFERENTIATION",
#   "GOBP_LEUKOCYTE_DIFFERENTIATION",
#   "GOBP_PROTEIN_FOLDING",
#   "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#   "GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
#   "GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
#   "GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE",
#   "GOBP_PROTEIN_REFOLDING",
#   "GOCC_VESICLE_MEMBRANE",
#   "GOCC_ER_TO_GOLGI_TRANSPORT_VESICLE_MEMBRANE",
#   "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
#   "GOCC_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE",
#   "GOCC_ENDOCYTIC_VESICLE",
#   "GOCC_PHAGOCYTIC_VESICLE",
#   "GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS",
#   "REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS",
#   "REACTOME_TRANSLATION",
#   "REACTOME_PERK_REGULATES_GENE_EXPRESSION",
#   "REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#   "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
#   "REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT",
#   "REACTOME_AUTOPHAGY",
#   "REACTOME_MTOR_SIGNALLING",
#   "REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",
#   "P53_02",
#   "P53_DECAMER_Q2",
#   "ATF4_Q2",
#   "CHOP_01",
#   "CIITA_TARGET_GENES",
#   "GATA1_01",
#   "GATA1_02",
#   "GATA1_03",
#   "GATA1_04",
#   "GATA1_05",
#   "MYC_Q2",
#   "MYCMAX_01",
#   "MYCMAX_02",
#   "MYCMAX_03",
#   "MYCMAX_B",
#   "NMYC_01"  
# )
