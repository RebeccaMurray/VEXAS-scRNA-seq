# Sys.setenv(RETICULATE_PYTHON="/gpfs/commons/home/rmurray/.virtualenvs/r-reticulate-gotcha/bin/python")

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(parallel)
# library(Gotcha)
theme_set(theme_bw())

source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_ncount_regressed_NormalizeData_soupx_mito_10_sd_PC_30_var_genes_3000_20231001.rds")
DefaultAssay(rna.obj) <- "RNA"

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

mclapply(c("HSC", "EMP", "LMPP", "Prog Mk"), function(sample.type) {
  
  print("Running for cluster...")
  print(sample.type)
  
  se <- subset(rna.obj, cluster_celltype == sample.type & Genotype %in% c("MUT", "WT"))
  
  ## Filter to only keep genes expressed in at least 5% of all genotyped cells for that cell type
  se <- FilterGenes(object = se, min.value = 0, min.cells = floor(nrow(se@meta.data)*0.05))
  
  # Remove mito/ribo genes
  non.mito.ribo.genes <- rownames(se@assays$RNA@data) %>%
    grep("^RP", ., value = T, invert = T) %>%
    grep("^MT-", ., value = T, invert = T)
  
  counts <- GetAssayData(se, assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% non.mito.ribo.genes)),]
  se <- subset(se, features = rownames(counts))
  
  m <- se@assays$RNA@data
  diffLMM.results = DiffLMMPseudocount(metadata = se@meta.data,
                                       provided.matrix = m,
                                       sample.column = "Donor",
                                       cluster.column = "cluster_celltype",
                                       selected.clusters = sample.type,
                                       treatment.column = "Genotype",
                                       treatment.levels = c("WT","MUT"),
                                       ncores = 18)  
  
  print("Saving output...")
  write_delim(diffLMM.results, file = paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells/", sample.type, "_no_mito_no_ribo.csv"), delim = ",")
  print("Done!")
}, mc.cores = detectCores())

