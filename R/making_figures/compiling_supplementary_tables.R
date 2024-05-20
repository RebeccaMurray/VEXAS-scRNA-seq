library(tidyverse)

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

##################################### Compiling GSEA modules used ###################################

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

pathways.merged <- lapply(pathways, function(x) {
  paste0(x, collapse = " ")
})

pathways.df <- as.data.frame(pathways.merged) %>% t() %>% as.data.frame() %>% rownames_to_column("pathway")
colnames(pathways.df) <- c("pathway", "gene list")

write_csv(pathways.df, file = 'data/gene_module_list.csv')
# 
# pathways %>%
#   map(as_tibble) %>%
#   reduce(bind_rows)
# pathways.list <- lapply(pathways, as.list)
# 
# purrr::map_df(pathways, tibble::as_tibble)
# 
# pathways.test <- pathways[1:5]
# 
# do.call(cbind, pathways.test) %>% View
# 
# do.call(rbind.data.frame, pathways)
# pathways %>% 
#   enframe() %>% 
#   unnest_longer(value)
#   
