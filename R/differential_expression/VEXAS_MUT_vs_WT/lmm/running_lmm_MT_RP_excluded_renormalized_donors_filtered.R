# Sys.setenv(RETICULATE_PYTHON="/gpfs/commons/home/rmurray/.virtualenvs/r-reticulate-gotcha/bin/python")

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
# library(Gotcha)
theme_set(theme_bw())

source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")

## Make sure umap is correct
DimPlot(rna.obj, group.by = "cluster_celltype")
dev.off()

print(table(rna.obj$Genotype_Approx_Match_No_Gene_Thresh, rna.obj$Donor))

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

mclapply(c("HSC", "EMP", "LMPP", "MkP", "GMP", "Early Eryth", "CD14 Mono"), function(celltype) {
  
  print("Running for cluster...")
  print(celltype)
  
  se <- subset(rna.obj, cluster_celltype == celltype & Genotype_Approx_Match_No_Gene_Thresh %in% c("MUT", "WT"))
  
  if (celltype == "CD14 Mono") {
    ## Exclude any cells with T and B cell marker genes from the monocyte DE
    print("Excluding cells with non-zero B and T cell marker gene expression...")
    cd3.expression.level <- se@assays$RNA@data[c("CD3D", "CD3E", "CD19"), ] %>% colSums()
    cd3.negative <- cd3.expression.level[cd3.expression.level == 0] %>% names()
    se <- subset(se, cells = cd3.negative)
  }
  
  ## Only keep donors with at least 5 cells of each genotype
  donors.keep <- se@meta.data %>% count(Donor, Genotype) %>% filter(n >= 5) %>% ungroup() %>% count(Donor) %>% filter(n > 1) %>% pull(Donor)
  if (length(donors.keep) < 3) {
    print(paste0("Fewer than 3 donors for cluster ", celltype))
  }
  se <- subset(se, Donor %in% donors.keep)
  print(se@meta.data %>% count(Genotype))
  print(se@meta.data %>% count(Donor, Genotype))
  
  # ## Only keep donors with at least 20 cells genotyped
  # donors.keep <- intersect(
  #   se@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "MUT") %>% count(Donor) %>% filter(n >= 10) %>% pull(Donor),
  #   se@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "WT") %>% count(Donor) %>% filter(n >= 10) %>% pull(Donor)
  # )
  # se <- subset(se, Donor %in% donors.keep)
  # print(table(se$Donor))
  # 
  ## Filter to only keep genes expressed in at least 10% of all genotyped cells for that cell type
  se <- FilterGenes(object = se, min.value = 0, min.cells = floor(nrow(se@meta.data)*0.05))
  
  # ## Only keep top 3000 variable features
  # DefaultAssay(se) <- "RNA"
  # se <- FindVariableFeatures(se, assay = "RNA", nfeatures = 3000)
  # m <- se@assays$RNA@data[VariableFeatures(se, assay = "RNA"), ]
  
  m <- se@assays$RNA@data
  diffLMM.results = DiffLMMPseudocount(metadata = se@meta.data,
                                       provided.matrix = m,
                                       sample.column = "Donor",
                                       cluster.column = "cluster_celltype",
                                       selected.clusters = celltype,
                                       treatment.column = "Genotype_Approx_Match_No_Gene_Thresh",
                                       treatment.levels = c("WT","MUT"),
                                       ncores = 18)  
  
  print("Saving output...")
  write_delim(diffLMM.results, file = paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/", celltype, "_5cells_genotyped_20231206.csv"), delim = ",")
  print("Done!")
}, mc.cores = detectCores())

