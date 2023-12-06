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

se <- subset(rna.obj, subset = Donor %in% c("PT02", "PT12", "PT15") & cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))
DefaultAssay(se) <- "ADT_DSB"
# se <- subset(se, Donor %in% c("VEXAS_BM1", "VEXAS_BM2", "VEXAS_BM12", "VEXAS_BM15")) ## Remove donors with < 5 cells for either genotype
m <- se@assays$ADT_DSB@data
se$cluster_const <- "HSC"
  
# For use with Gotcha helper function included in package
unt.results = DiffLMMPseudocount(metadata = se@meta.data,
                                 provided.matrix = m,
                                 sample.column = "Donor",
                                 cluster.column = "cluster_const",
                                 selected.clusters = "HSC",
                                 treatment.column = "Genotype",
                                 treatment.levels = c("WT","MUT"),
                                 pseudocount.value = 1,
                                 ncores = 18)  

print("Saving output...")
write_delim(unt.results, file = paste0("data/differential_expression/VEXAS_ADT/HSCs_MUT_vs_WT_ADT_20231205.csv"), delim = ",")

## Plot
unt.results %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  plot_volcano(unt.results, stat_column = "fdr", effect_column = "delta", effect_line = 0.2) + coord_cartesian(ylim = c(0, 25))
  
