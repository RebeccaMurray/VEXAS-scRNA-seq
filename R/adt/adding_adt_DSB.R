library(Seurat)
library(tidyverse)
library(dsb)

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_NormalizeData_soupx_mito_20_PC_30_20230926.rds")

DefaultAssay(rna.obj) <- "RNA"
# rna.obj@assays$ADT_DSB <- NULL

######################## Load ADT data for each sample, run DSB, then run PCA on the ADT assay #########################

# Load adt data for both real and empty droplets

print("Beginning ADT")

kallisto.output.dictionary <- list(
  `VEXAS_c_BM2_CD34_plus_SR` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_plus/counts_unfiltered/",
  `VEXAS_c_BM2_CD34-SR` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_minus/counts_unfiltered/",
  `VX_BM14_Fpos_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fpos/counts_unfiltered/",
  `VX_BM14_Fneg_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fneg/counts_unfiltered/",
  `VX_BM12_A_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM12_A/counts_unfiltered/",
  `VX_BM13_B_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM13_B/counts_unfiltered/",
  `SG_VX_BM15a_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15a/counts_unfiltered/",
  `SG_VX_BM15b_GEX` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15b/counts_unfiltered/"
)

pull.dsb.counts <- function(sample.name) {
  
  print("Running DSB for:")
  print(sample.name)
  
  base.kallisto.path <- kallisto.output.dictionary[[sample.name]]
  adt_data <- Matrix::readMM(paste0(base.kallisto.path, "cells_x_features.mtx"))
  adt_bcs <- read.csv(paste0(base.kallisto.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(sample.name, "_", V1, "-1")) %>% pull(V1) 
  adt_tag_names <- read.csv(paste0(base.kallisto.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
  rownames(adt_data) <- adt_bcs
  colnames(adt_data) <- adt_tag_names
  
  # Load the barcodes from the cellranger output
  real.bcs.cellranger <- read_csv(paste0("/gpfs/commons/home/rmurray/VEXAS/", sample.name, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"), col_names = c("BC")) %>%
    mutate(BC = paste0(sample.name, "_", BC)) %>% 
    pull(BC)
  all.bcs.cellranger <- read_csv(paste0("/gpfs/commons/home/rmurray/VEXAS/", sample.name, "/outs/raw_feature_bc_matrix/barcodes.tsv.gz"), col_names = c("BC")) %>%
    mutate(BC = paste0(sample.name, "_", BC)) %>% 
    pull(BC)
  
  # Define empty droplet barcodes based on cellranger output
  empty.droplet.bcs <- setdiff(all.bcs.cellranger, real.bcs.cellranger)
  
  # Define real barcodes based on the QC filtered seurat object
  real.bcs <- rna.obj@meta.data %>% filter(orig.ident == sample.name) %>% rownames()
  
  # Filter both lists for only barcodes present in the ADT output
  empty.droplet.bcs <- empty.droplet.bcs[(empty.droplet.bcs %in% adt_bcs)]
  real.bcs <- real.bcs[(real.bcs %in% adt_bcs)]
  
  # # Filter empty droplet barcodes based on previously identified background drop distributions
  # empty.droplet.bcs.within.boundaries <- read_csv(paste0("data/qc_filtering/dsb_background/", sample.name, ".csv")) %>%
  #    mutate(barcode = paste0(prefix.dictionary[[sample.name]], barcode, suffix.dictionary[[sample.name]])) %>%
  #    pull(barcode)
  # empty.droplet.bcs <- empty.droplet.bcs[empty.droplet.bcs %in% empty.droplet.bcs.within.boundaries]
  
  # Separate real and empty cells
  adt.cells = adt_data[real.bcs,] %>% as.matrix()
  adt.empty = adt_data[empty.droplet.bcs,] %>% as.matrix()
  Prot = t(adt.cells)
  empty = t(adt.empty)
  print(dim(Prot))
  print(dim(empty))
  
  isotypes <- adt_tag_names[grepl("CTRL|Isotype", adt_tag_names)]
  print(isotypes)
  
  # Run DSB
  normProt = DSBNormalizeProtein(
    cell_protein_matrix = Prot,
    empty_drop_matrix = empty,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = isotypes
  )
  
  return(normProt)
}


## Load ADT values and run DSB for each sample
dsb.counts.list <- lapply(names(kallisto.output.dictionary), pull.dsb.counts)
dsb.counts <- do.call(cbind, dsb.counts.list)
dsb.counts <- t(dsb.counts)

## Add NA counts for remaining samples
# If there are any atac barcodes not present in adt, add them to the adt output
rna.obj.bcs <- rna.obj@meta.data %>% rownames()
bc_missing_from_adt <- rna.obj.bcs[!(rna.obj.bcs %in% rownames(dsb.counts))]
print("Number of barcodes missing from dsb.counts data for:")
print(length(bc_missing_from_adt))
missing_matrix <- matrix(NA_real_, ncol=ncol(dsb.counts), nrow=length(bc_missing_from_adt), dimnames=list(bc_missing_from_adt, colnames(dsb.counts)))
dsb.counts <- rbind(dsb.counts, missing_matrix)

## Add as a new assay, normalize, and scale
dsb.counts <- dsb.counts[rna.obj@meta.data %>% rownames, ]
rna.obj[["ADT_DSB"]] <- CreateAssayObject(data = t(dsb.counts)) ## Important for DSB to add as the "data" slot
DefaultAssay(rna.obj) <- "ADT_DSB"
rna.obj <- ScaleData(rna.obj)

## Correct CD90 for BM2
BM2.barcodes <- rna.obj@meta.data %>% filter(orig.ident == "VEXAS_BM2_CD34_plus" | orig.ident == "VEXAS_BM2_CD34_negative") %>% rownames()
rna.obj@assays$ADT_DSB@data["Hu.CD90", BM2.barcodes] <- NA_real_
rna.obj@assays$ADT_DSB@scale.data["Hu.CD90", BM2.barcodes] <- NA_real_

## Done!

saveRDS(rna.obj, "~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_NormalizeData_soupx_mito_20_PC_30_20230926.rds")
