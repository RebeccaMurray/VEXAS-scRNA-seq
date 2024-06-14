library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)

run_enrichr_for_list <- function(gene.list) {
  m_t2g <- tibble()
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% filter(gs_name %in% fgsea.pathways))
  m_t2g <- m_t2g %>% dplyr::select(gs_name, human_entrez_gene) %>% dplyr::rename(TermID = gs_name) %>% dplyr::rename(geneID = human_entrez_gene)
  
  ## Add GoT paper modules
  got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv")
  additional.modules <- got.modules %>% pivot_longer(cols = colnames(got.modules), names_to = "TermID", values_to = "geneID") %>% 
    dplyr::filter(!is.na(geneID))
  
  ## Add Han et al modules
  atf4.df <- data.frame(geneID = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% dplyr::rename(geneID = Feature))
  atf4.df$TermID = "ATF4 only (Han et al.)"
  atf4.df <- atf4.df %>% dplyr::select(TermID, geneID)
  chop.df <- data.frame(geneID = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% dplyr::rename(geneID = Feature))
  chop.df$TermID = "CHOP only (Han et al.)"
  chop.df <- chop.df %>% dplyr::select(TermID, geneID)
  
  additional.modules <- rbind(additional.modules, atf4.df)
  additional.modules <- rbind(additional.modules, chop.df)
  
  hs <- org.Hs.eg.db
  mapping <- AnnotationDbi::select(hs, 
                                   keys = additional.modules$geneID,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  additional.modules <- merge(additional.modules, mapping, by.x = "geneID", by.y = "SYMBOL")
  additional.modules <- additional.modules %>% dplyr::select(TermID, ENTREZID) %>% dplyr::rename(geneID = ENTREZID) %>% 
    filter(!is.na(geneID))
  m_t2g <- rbind(m_t2g, additional.modules)
  
  ## Map input gene list to entrezid
  gene.list.entrez.ids <- AnnotationDbi::select(hs, 
                                                keys = gene.list,
                                                columns = c("ENTREZID", "SYMBOL"),
                                                keytype = "SYMBOL")
  
  
  res.custom <- enricher(gene.list.entrez.ids$ENTREZID, gson = NULL, TERM2GENE = m_t2g, pvalueCutoff = 1)
  res.custom@result %>% dplyr::filter(pvalue < 0.01) %>% print()
  return(res.custom@result)
}