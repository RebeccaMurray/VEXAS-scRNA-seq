library(Seurat)
library(Signac)
library(msigdbr)

## Extract the raw counts matrix, remove mito/ribo genes
non.mito.ribo.genes <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^[MT-|RP]", ., value = T, invert = T)
m <- rna.obj@assays$RNA@counts[non.mito.ribo.genes, ]

## Replace the RNA assay
DefaultAssay(rna.obj) <- "ADT_DSB" ## Set another random assay as default
rna.obj[["RNA"]] <- NULL
rna.obj[["RNA"]] <- CreateAssayObject(counts = m)
DefaultAssay(rna.obj) <- "RNA"
rna.obj <- NormalizeData(rna.obj, assay = "RNA", normalization.method = "LogNormalize")

rna.obj.HSC <- subset(rna.obj, cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))


m_t2g <- tibble()
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
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
  `Hallmark interferon alpha response` = pathways$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  `Hallmark interferon gamma response` = pathways$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  `Hallmark tnfa` = pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
  `Reactome translation` = pathways$REACTOME_TRANSLATION,
  `ATF4 target genes [Han et al.]` = pathways$ATF4_ONLY_HAN,
  `CHOP target genes [Han et al.]` = pathways$CHOP_ONLY_HAN,
  `Hallmark UPR` = pathways$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
  `Hallmark P53 pathway` = pathways$HALLMARK_P53_PATHWAY,
  KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY = pathways$KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY,
  GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS = pathways$GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS
)

## Add as modules to our RNA object
DefaultAssay(rna.obj.HSC) <- "RNA"
for (module in names(module.list)) {
  print(module)
  print(module.list[[module]])
  rna.obj.HSC <- AddModuleScore(rna.obj.HSC, features = list(module.list[[module]]), name = module, assay = "RNA", search = TRUE)
}

## Make boxplots for the modules
plist <- lapply(names(module.list), function(module) {
  module.name <- paste0(module, "1")
  p1 <- rna.obj.HSC@meta.data %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    dplyr::rename(module_plotting = module.name) %>% 
    ggplot(aes(x = Genotype, y = module_plotting, fill = Genotype)) +
    geom_violin() +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(str_wrap(module, width = 15)) +
    # scale_fill_manual(values = c(WT = "lightsteelblue1", MUT = "lightsteelblue3")) +
    scale_fill_manual(values = vexas.genotyping.palette) +
    scale_y_continuous(expand = c(0.2, 0)) +
    stat_compare_means(label = "p.format", label.y = Inf, vjust = 1.5) +
    ylab("Module score") +
    xlab("")
  return(p1)
})


rna.obj <- AddModuleScore(rna.obj, features = list(pathways$`ATF6 module`), name = "ATF6_module", assay = "RNA", search = TRUE)
FeaturePlot(rna.obj, features = "ATF6_module1")
VlnPlot(rna.obj, features = "ATF6_module1", group.by = "cluster_celltype")

wrap_plots(plist)
Idents(rna.obj.HSC) <- "Genotype"
VlnPlot(rna.obj.HSC, features = "IL1B", split.by = "Genotype", group.by = "Donor", idents = c("MUT", "WT"))

VlnPlot(rna.obj.HSC, features = "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS1", split.by = "Genotype", group.by = "Donor", idents = c("MUT", "WT"))

high.ribo.hscs <- rna.obj.HSC@meta.data %>% filter(percent.ribo > 50) %>% rownames()

## Correlation?
md <- rna.obj.HSC@meta.data
md %>% 
  ggplot(aes(x = `ATF4 target genes [Han et al.]1`, y = `Reactome translation1`, color = Genotype)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = vexas.genotyping.palette) +
  facet_grid(cols = vars(Genotype))

module.names <- paste0(names(module.list), "1")
md <- rna.obj.HSC@meta.data %>% filter(Genotype == "MUT")
m <- md[, module.names]
cor.m <- cor(m)

translation.low.cells <- rna.obj.HSC@meta.data %>% filter(`Reactome translation1` < 0.25) %>% rownames()
DimPlot(rna.obj, cells.highlight = translation.low.cells)

############## running fgsea for all msigdb modules #################

de.list <- de.results$HSC

  if (!("feature" %in% colnames(de.list))) {
    de.list$feature <- rownames(de.list)
  }
  
  de_results.ranked <- de.list %>%
    dplyr::rename(pval = pval) %>% 
    filter(!is.na(pval)) %>% 
    filter(is.finite(-log10(pval))) %>% 
    mutate(rank = avg_log2FC * -log10(pval)) %>%
    arrange(-rank)
  
  ranks.list <- de_results.ranked$rank
  names(ranks.list) <- de_results.ranked$feature
  
head(ranks.list)
tail(ranks.list)  
  
## Compile modules we're interested in
m_t2g <- tibble()
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5"))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
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

fgseaRes <- fgsea(pathways = pathways, stats = ranks.list, minSize = 5)
