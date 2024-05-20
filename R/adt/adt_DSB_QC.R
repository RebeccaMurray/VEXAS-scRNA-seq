library(Seurat)
library(tidyverse)
library(patchwork)

kallisto.output.dictionary <- list(
  `VEXAS_BM2_CD34_plus` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_plus/counts_unfiltered/",
  `VEXAS_BM2_CD34_negative` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VEXAS_c_BM2_CD34_minus/counts_unfiltered/",
  `VX_BM14_Fpos` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fpos/counts_unfiltered/",
  `VX_BM14_Fneg` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM14_Fneg/counts_unfiltered/",
  `VX_BM12_A` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM12_A/counts_unfiltered/",
  `VX_BM13_B` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM13_B/counts_unfiltered/",
  `VX_BM15a` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15a/counts_unfiltered/",
  `VX_BM15b` = "/gpfs/commons/home/rmurray/VEXAS_CITE/kallisto_v2/VX_BM15b/counts_unfiltered/"
)

cellranger.dictionary <- list(
  `VEXAS_BM2_CD34_plus` = "VEXAS_c_BM2_CD34_plus_SR",
  `VEXAS_BM2_CD34_negative` = "VEXAS_c_BM2_CD34-SR",
  `VX_BM14_Fpos` = "VX_BM14_Fpos_GEX",
  `VX_BM14_Fneg` = "VX_BM14_Fneg_GEX",
  `VX_BM12_A` = "VX_BM12_A_GEX",
  `VX_BM13_B` = "VX_BM13_B_GEX",
  `VX_BM15a` = "SG_VX_BM15a_GEX",
  `VX_BM15b` = "SG_VX_BM15b_GEX"
)


dsb.background.params <- list(
  `VEXAS_BM2_CD34_plus` = list(prot.min = 1.5, prot.max = 3, ngene.min = 1, ngene.max = 2),
  `VEXAS_BM2_CD34_negative` = list(prot.min = 1.5, prot.max = 3, ngene.min = 1.5, ngene.max = 2.5),
  `VX_BM14_Fpos` = list(prot.min = 1.5, prot.max = 3, ngene.min = 1.5, ngene.max = 2.25),
  `VX_BM14_Fneg` = list(prot.min = 1.5, prot.max = 3.5, ngene.min = 1.5, ngene.max = 2.5),
  `VX_BM12_A` = list(prot.min = 1, prot.max = 2.5, ngene.min = 1, ngene.max = 2),
  `VX_BM13_B` = list(prot.min = 1, prot.max = 2.5, ngene.min = 1, ngene.max = 2),
  `VX_BM15a` = list(prot.min = 1, prot.max = 2.5, ngene.min = 1.5, ngene.max = 2.5),
  `VX_BM15b` = list(prot.min = 1, prot.max = 2.5, ngene.min = 1, ngene.max = 2)
)

plot_background_for_sample <- function(sample.name) {
  print("Running for:")
  print(sample.name)
  ## adt data, all cells
  base.kallisto.path <- kallisto.output.dictionary[[sample.name]]
  adt_data <- Matrix::readMM(paste0(base.kallisto.path, "cells_x_features.mtx"))
  adt_bcs <- read.csv(paste0(base.kallisto.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% pull(V1) %>% paste0(., "-1")
  adt_tag_names <- read.csv(paste0(base.kallisto.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
  rownames(adt_data) <- adt_bcs
  colnames(adt_data) <- adt_tag_names
  adt_data <- t(adt_data %>% as.matrix())
  
  # Load the RNA data, with and without background
  cellranger.raw <- Read10X(paste0("/gpfs/commons/home/rmurray/VEXAS/", cellranger.dictionary[[sample.name]], "/outs/raw_feature_bc_matrix/"))
  cellranger.filtered <- Read10X(paste0("/gpfs/commons/home/rmurray/VEXAS/", cellranger.dictionary[[sample.name]], "/outs/filtered_feature_bc_matrix/"))
  
  ## Filter for only RNA data (cells + empty droplets) with ADT data
  cellranger.raw <- cellranger.raw[, adt_bcs]
  
  ## Separate the real cells and the empty droplets
  stained_cells = colnames(cellranger.filtered)
  background = setdiff(colnames(cellranger.raw), stained_cells)
  
  mtgene = grep(pattern = "^MT-", rownames(cellranger.raw), value = TRUE) 
  md = data.frame(
    rna.size = log10(Matrix::colSums(cellranger.raw)), 
    prot.size = log10(Matrix::colSums(adt_data)), 
    n.gene = Matrix::colSums(cellranger.raw > 0), 
    mt.prop = Matrix::colSums(cellranger.raw[mtgene, ]) / Matrix::colSums(cellranger.raw)
  )
  md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  md = md[md$rna.size > 0 & md$prot.size > 0, ]
  
  p1 <- ggplot(md, aes(y = log10(n.gene), x = prot.size)) +
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") +
    facet_wrap(~drop.class) +
    ggtitle(sample.name)
  
  ## Load thresholds for sample
  thresholds <- dsb.background.params[[sample.name]]
  
  bc.background <- md %>% filter(drop.class == "background") %>%
    mutate(n.gene = log10(n.gene)) %>%
    # filter(n.gene > thresholds$ngene.min & n.gene < thresholds$ngene.max) %>%
    # filter(prot.size > thresholds$prot.min & prot.size < thresholds$prot.max) %>%
    filter(n.gene < thresholds$ngene.max) %>%
    filter(prot.size < thresholds$prot.max) %>%
    rownames_to_column("barcode") %>%
    select(barcode)
  
  bc.background %>% 
    write_csv(., file = paste0("data/qc_filtering/dsb_background/", sample.name, ".csv"))
  
  p1 <- p1 + 
    # geom_hline(yintercept = thresholds$ngene.min, color = "red") +
    geom_hline(yintercept = thresholds$ngene.max, color = "red") +
    # geom_vline(xintercept = thresholds$prot.min, color = "red") +
    geom_vline(xintercept = thresholds$prot.max, color = "red")
  
  return(p1)
}

sample.list <- names(cellranger.dictionary)
names(sample.list) <- sample.list
plist <- lapply(sample.list, plot_background_for_sample)

wrap_plots(plist)

