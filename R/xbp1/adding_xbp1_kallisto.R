library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(viridis)
library(pheatmap)
theme_set(theme_classic())


rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")

## Clear previous XBP1s fields
rna.obj$Spliced_XBP1 <- NULL
rna.obj$Unspliced_XBP1 <- NULL
rna.obj$log10_splicing_ratio <- NULL
rna.obj$splicing_fraction <- NULL
rna.obj$splicing_ratio <- NULL


############################################ Load XBP1s data - kallisto outputs #################################

xbp1.kallisto.output.dictionary <- list(
  `UPN14_BMMC` = "NIH_1",
  `UPN15_BMMC` = "NIH_2",
  `UPN16_BMMC` = "NIH_3",
  `UPN17_BMMC` = "NIH_4",
  `UPN14_CD34` = "NIH_5",
  `UPN15_CD34` = "NIH_6",
  `UPN16_CD34` = "NIH_7",
  `UPN17_CD34` = "NIH_8",
  `VEXAS_BM1_CD34_plus_SR` = "VEXAS_BM1_CD34_plus",
  `VEXAS_BM1_CD34-_SR` = "VEXAS_BM1_CD34_negative",
  `VEXAS_c_BM2_CD34_plus_SR` = "VEXAS_BM2_CD34_plus",
  `VEXAS_c_BM2_CD34-SR` = "VEXAS_BM2_CD34_negative",
  `VX_BM14_Fpos_GEX` = "VX_BM14_Fpos",
  `VX_BM14_Fneg_GEX` = "VX_BM14_Fneg",
  `VX_BM12_A_GEX` = "VX_BM12_A",
  `VX_BM13_B_GEX` = "VX_BM13_B",
  `SG_VX_BM15a_GEX` = "VX_BM15a",
  `SG_VX_BM15b_GEX` = "VX_BM15b"
)

pull.xbp1.counts <- function(sample.name) {
  print(sample.name)
  sample.name.modified <- sample.name
  if (sample.name %in% names(xbp1.kallisto.output.dictionary)) {
    sample.name.modified <- xbp1.kallisto.output.dictionary[[sample.name]]
  }
  
  base.kallisto.path <- paste0("~/VEXAS/XBP1_splicing/kallisto/", sample.name.modified, "/counts_unfiltered/")
  adt_data <- Matrix::readMM(paste0(base.kallisto.path, "cells_x_features.mtx"))
  adt_bcs <- read.csv(paste0(base.kallisto.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(sample.name, "_", V1, "-1")) %>% pull(V1) 
  adt_tag_names <- read.csv(paste0(base.kallisto.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
  rownames(adt_data) <- adt_bcs
  colnames(adt_data) <- adt_tag_names
  return(adt_data)
}

xbp1.counts <- lapply(rna.obj$orig.ident %>% unique(), pull.xbp1.counts)
xbp1.counts.df <- do.call(rbind, xbp1.counts)

rna.obj <- AddMetaData(rna.obj, xbp1.counts.df)

## Add splicing ratio, fraction
rna.obj$splicing_fraction <- rna.obj$Spliced_XBP1 / (rna.obj$Spliced_XBP1 + rna.obj$Unspliced_XBP1)
rna.obj$splicing_ratio <- rna.obj$Spliced_XBP1 / rna.obj$Unspliced_XBP1
rna.obj$log10_splicing_ratio <- log10(rna.obj$splicing_ratio)

