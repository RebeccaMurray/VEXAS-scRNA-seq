## Script to pre-process each sample by removing doublets.

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(DoubletFinder)
library(parallel)
set.seed(1234)

theme_set(theme_classic())

##################################### Load all samples ###########################################
cellranger.output.dir <- "/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/soupx_output/"
samples <- list.dirs(cellranger.output.dir, recursive = F, full.names = F) %>% grep("Healthydonor",., invert = T, value = T)
sample.list <- mclapply(samples, function(x) {
  print("Loading sample:")
  print(x)
  se <- CreateSeuratObject(counts = Read10X(data.dir = paste0(cellranger.output.dir, x, "/", x, "/")), project = x, min.cells = 3, min.features = 200)
  se <- RenameCells(se, add.cell.id = x)
  se$percent.mito <- PercentageFeatureSet(se, pattern = "^MT-")
  se$percent.ribo <- PercentageFeatureSet(se, pattern = "^RPL|^RPS")
  return(se)
}, mc.cores = detectCores())
names(sample.list) <- samples


##################################### Run some preliminary QC filtering ###########################################

## Note: We are only performing a "basic" QC filtering here to not remove cells with higher numbers of 
sample.list <- mclapply(sample.list, mc.cores = detectCores(), FUN = function(sample) {
  subset(sample, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent.mito < 20)
})

print("Done with filtering")
print("Printing samples...")
for(sample in sample.list) {
  print(sample)
}


########################################## Helper function to add doublet scores #############################

## Function to add doublet scores to seurat object
add_doublet_scores <- function(se.name) {
  
  print("Running for:")
  print(se.name)
  se <- sample.list[[se.name]]
  DefaultAssay(se) <- "RNA"
  
  ## Do a quick clustering
  se <- NormalizeData(se)
  se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
  se <- ScaleData(se)
  se <- RunPCA(se)
  se <- RunUMAP(se, dims = 1:20)
  
  ## Find the max PK value
  sweep.res.list <- paramSweep_v3(se, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_results <- find.pK(sweep.stats)
  
  ## make the PK plot
  p.max.pk <- ggplot(bcmvn_results, aes(pK, BCmetric, group = 1)) + 
    geom_point() +
    geom_line() +
    ggtitle(se.name) + 
    theme_classic()
  print(p.max.pk)
  ggsave(paste0("figures/doublet_removal/", se.name, ".png"), device = "png", width = 10, height = 4)
  
  ## Select the max from the pK plot programmatically
  pK <- bcmvn_results %>% 
    filter(BCmetric == max(BCmetric)) %>% 
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  print(paste0("pK value selected: ", pK))
  
  ######################## Fixed multiplet rate calulation based on 10x manual #####################
  
  ## Calculate the homotypic doublet rate based on the number of recovered cells
  ## Values from 10x v3 manual: https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
  # cell.count <- floor(nrow(se@meta.data)/1000)*1000
  # multiplets <- list(
  #   `0` = 0.004,
  #   `1000` = 0.008,
  #   `2000` = 0.016,
  #   `3000` = 0.023,
  #   `4000` = 0.031,
  #   `5000` = 0.039,
  #   `6000` = 0.046,
  #   `7000` = 0.054,
  #   `8000` = 0.061,
  #   `9000` = 0.069,
  #   `10000` = 0.076,
  #   `11000` = 0.076,
  #   `12000` = 0.076,
  #   `13000` = 0.076,
  #   `14000` = 0.076,
  #   `15000` = 0.076,
  #   `16000` = 0.076,
  #   `17000` = 0.076,
  #   `18000` = 0.076
  # )
  # multiplet.estimate <- multiplets[[as.character(cell.count)]]
  
  ########################## Mutiplet estimation by rate #########################
  
  ## Multiplet estimate by rate
  ## 0.8% for every 1,000 cells targeted 0 https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing-
  multiplet.estimate <- nrow(se@meta.data) * 8e-6
  
  ########################## Model homotypic doublets #########################
  
  print("Cell count:")
  print(nrow(se@meta.data))
  print("Multiplet estimate used:")
  print(multiplet.estimate)
  
  ## Calculate number of expected homotypic doublets (which the tool can't detect as well)
  homotypic.prop <- modelHomotypic(se@meta.data$seurat_clusters)
  nExp_poi <- round(multiplet.estimate*nrow(se@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  print("Running DoubletFinder...")
  
  ## Add doublet classifications to seurat object and return
  se <- doubletFinder_v3(se, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  print("Done!")
  df <- se@meta.data
  pANN.col <- colnames(df) %>% grep("pANN", ., value = T)
  DF.classifications.col <- colnames(df) %>% grep("DF.classifications_", ., value = T)
  df <- df %>% 
    dplyr::rename(pANN = pANN.col) %>% 
    dplyr::rename(DF.classification = DF.classifications.col) %>% 
    select(pANN, DF.classification)
  return(df)
}

########################################## Running for all samples #############################

doublets.df.list <- mclapply(names(sample.list), add_doublet_scores, mc.cores = detectCores())
doublets.df <- do.call(rbind, doublets.df.list)

## Save as a csv
doublets.df %>% 
  rownames_to_column("BC") %>% 
  write_csv("data/qc_filtering/doubletfinder_output_20231103_v2.csv")
