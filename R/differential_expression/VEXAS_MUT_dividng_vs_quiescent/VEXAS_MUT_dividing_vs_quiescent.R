# Sys.setenv(RETICULATE_PYTHON="/gpfs/commons/home/rmurray/.virtualenvs/r-reticulate-gotcha/bin/python")

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
# library(Gotcha)
theme_set(theme_bw())

source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")

## Extract the raw counts matrix, remove mito/ribo genes
non.mito.ribo.genes <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., value = T, invert = T)
print("Excluding genes:")
print(setdiff(rna.obj@assays$RNA@counts %>% rownames(), non.mito.ribo.genes))
m <- rna.obj@assays$RNA@counts[non.mito.ribo.genes, ]

## Extract the raw counts matrix, remove cell cycle genes
cell.cycle.genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
print("Excluding genes:")
print(setdiff(rna.obj@assays$RNA@counts %>% rownames(), cell.cycle.genes))
m <- rna.obj@assays$RNA@counts[setdiff(rna.obj@assays$RNA@counts %>% rownames(), cell.cycle.genes) ]


## Replace the RNA assay
DefaultAssay(rna.obj) <- "ADT_DSB" ## Set another random assay as default
rna.obj[["RNA"]] <- NULL
rna.obj[["RNA"]] <- CreateAssayObject(counts = m)
DefaultAssay(rna.obj) <- "RNA"
rna.obj <- NormalizeData(rna.obj, assay = "RNA", normalization.method = "LogNormalize")

# Helper function to filter for minimum expression
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    # object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@assays$RNA@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@assays$RNA@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@assays$RNA@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      # object@assays$RNA@data <- object@assays$RNA@data[genes.use, ] # original way of subsetting
      
      ## modified way of subsetting, removing from entire obj
      counts <- GetAssayData(object, assay = "RNA")
      counts <- counts[(which(rownames(counts) %in% genes.use)),]
      object <- subset(object, features = rownames(counts))
      return(object)
    } else {
      return(object)
    }
  }



se <- subset(rna.obj, cluster_celltype == "HSC" & Genotype_Approx_Match_No_Gene_Thresh %in% c("MUT"))

se$status <- "G1"
se$status[se$Phase %in% c("S", "G2M")] <- "S/G2M"

se@meta.data %>% 
  count(status, Phase)
  
## Only keep donors with at least 5 cells of each genotype
donors.keep <- se@meta.data %>% count(Donor, status) %>% filter(n >= 5) %>% ungroup() %>% count(Donor) %>% filter(n > 1) %>% pull(Donor)
if (length(donors.keep) < 3) {
    print(paste0("Fewer than 3 donors for cluster ", celltype))
}
se <- subset(se, Donor %in% donors.keep)
print(se@meta.data %>% count(status))
print(se@meta.data %>% count(Donor, status))

se <- FilterGenes(object = se, min.value = 0, min.cells = floor(nrow(se@meta.data)*0.05))
  

  
m <- se@assays$RNA@data
diffLMM.results = DiffLMMPseudocount(metadata = se@meta.data,
                                       provided.matrix = m,
                                       sample.column = "Donor",
                                       cluster.column = "cluster_celltype",
                                       selected.clusters = "HSC",
                                       treatment.column = "status",
                                       treatment.levels = c("G1","S/G2M"),
                                       ncores = 18)  
  
print("Done!")
write_delim(diffLMM.results, file = "data/differential_expression/VEXAS_MUT_dividing_vs_quiescent/HSC_5cells_genotyped_20231222.csv", delim = ",")

###################################### Look at results #######################################

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

de.results <- read_csv("data/differential_expression/VEXAS_MUT_dividing_vs_quiescent/HSC_5cells_genotyped_20231222.csv") %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = "HSC")

## Remove cell cycle genes
de.results.filtered <- de.results %>% 
  filter(!(feature %in% cc.genes$s.genes)) %>% 
  filter(!(feature %in% cc.genes$g2m.genes)) 

plot_volcano(de.results.filtered, stat_cutoff = 0.05, effect_line = 0.1, stat_column = "fdr", title = paste0("HSCs")) + 
  theme(legend.position = "right")

fgsea.out <- run_fgsea_for_list(de.results.filtered, msigdb.list = "REACTOME")
