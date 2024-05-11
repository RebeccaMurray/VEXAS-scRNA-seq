# Sys.setenv(RETICULATE_PYTHON="/gpfs/commons/home/rmurray/.virtualenvs/r-reticulate-gotcha/bin/python")

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(parallel)
# library(Gotcha)
theme_set(theme_bw())
source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
source("~/rscripts/general_helpers/custom_hypergeometric_helpers.R")
source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")

## Make sure umap is correct
DimPlot(rna.obj, group.by = "cluster_celltype")
dev.off()

print(table(rna.obj$Genotype_Approx_Match_No_Gene_Thresh, rna.obj$Donor))

######################## Set up object and run the DE ######################

## Extract the raw counts matrix, remove mito/ribo genes
non.mito.ribo.genes <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., value = T, invert = T)
print("Excluding genes:")
print(setdiff(rna.obj@assays$RNA@counts %>% rownames(), non.mito.ribo.genes))
m <- rna.obj@assays$RNA@counts[non.mito.ribo.genes, ]

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

celltype = "CD14 Mono"
  
print("Running for cluster...")
print(celltype)
  
se <- subset(rna.obj, cluster_celltype == celltype & Genotype_Approx_Match_No_Gene_Thresh %in% c("MUT", "WT"))

## Exclude any cells with T and B cell marker genes from the monocyte DE
print("Excluding cells with non-zero B and T cell marker gene expression...")
cd3.expression.level <- se@assays$RNA@data[c("CD3D", "CD3E", "CD19"), ] %>% colSums()
cd3.negative <- cd3.expression.level[cd3.expression.level == 0] %>% names()
se <- subset(se, cells = cd3.negative)
  
## Only keep donors with at least 5 cells of each genotype
donors.keep <- se@meta.data %>% count(Donor, Genotype) %>% filter(n >= 5) %>% ungroup() %>% count(Donor) %>% filter(n > 1) %>% pull(Donor)
se <- subset(se, Donor %in% donors.keep)

print(table(se$Donor, se$Genotype))
print(table(se$Genotype))
rna.obj@meta.data %>%
  filter(cluster_celltype == celltype) %>% 
  group_by(Genotype_Approx_Match_No_Gene_Thresh) %>%
  count() %>%
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>%
  mutate(total = sum(MUT, WT, `NA`)) %>%
  mutate(frac_genotyped = sum(MUT, WT) / total) %>%
  arrange(-frac_genotyped) %>% print

## Filter to only keep genes expressed in at least X% of all genotyped cells for that cell type
DefaultAssay(se) <- "RNA"
se <- FilterGenes(object = se, min.value = 0, min.cells = 10)

  
# Option 3: Just keep all genes
m <- se@assays$RNA@data

diffLMM.results = DiffLMMPseudocount(metadata = se@meta.data,
                                       provided.matrix = m,
                                       sample.column = "Donor",
                                       cluster.column = "cluster_celltype",
                                       selected.clusters = celltype,
                                       treatment.column = "Genotype_Approx_Match_No_Gene_Thresh",
                                       treatment.levels = c("WT","MUT"),
                                       ncores = detectCores())
diffLMM.results <- diffLMM.results %>% mutate(avg_log2FC = log2(fc))


write_delim(diffLMM.results, file = paste0("data/differential_expression/VEXAS_MUT_vs_WT/", celltype, "_20240511.csv"), delim = ",")

