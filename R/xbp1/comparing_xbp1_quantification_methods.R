library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(viridis)
library(pheatmap)
theme_set(theme_classic())

source("~/rscripts/global_vexas_configs/color_palettes.R")
source("~/rscripts/general_helpers/adding_module_scores.R")


xbp1.kallisto.output.dictionary <- list(
  `BM1_CD34_plus` = "VEXAS_BM1_CD34_plus/counts_unfiltered/",
  `BM1_CD34_negative` = "VEXAS_BM1_CD34_minus/counts_unfiltered/",
  `BM2_CD34_plus` = "VEXAS_c_BM2_CD34_plus/counts_unfiltered/",
  `BM2_CD34_negative` = "VEXAS_c_BM2_CD34_minus/counts_unfiltered/",
  `VX_BM14_Fpos` = "VX_BM14_Fpos/counts_unfiltered/",
  `VX_BM14_Fneg` = "VX_BM14_Fneg/counts_unfiltered/",
  `VX_BM12_A` = "VX_BM12_A/counts_unfiltered/",
  `VX_BM13_B` = "VX_BM13_B/counts_unfiltered/",
  `VX_BM15a` = "VX_BM15a/counts_unfiltered/",
  `VX_BM15b` = "VX_BM15b/counts_unfiltered/"
)

## Pull suffixes and prefixes for sample
suffix.df <- samples.vexas.and.control@meta.data %>% rownames_to_column("Barcode") %>% group_by(orig.ident) %>% summarize(suffix = gsub("landau_.*-1", "", Barcode) %>% unique())
suffix.dictionary <- suffix.df$suffix
names(suffix.dictionary) <- suffix.df$orig.ident


prefix.df <- samples.vexas.and.control@meta.data %>% rownames_to_column("Barcode") %>% group_by(orig.ident) %>% summarize(prefix = gsub("[ACGT]*-1_[0-9]", "", Barcode) %>% unique())
prefix.dictionary <- prefix.df$prefix
names(prefix.dictionary) <- prefix.df$orig.ident

sample.name = "VX_BM15a"

## Original run
base.kallisto.path <- "~/VEXAS/XBP1_splicing/kallisto/VX_BM15a/counts_unfiltered/"
xbp1.data.original <- Matrix::readMM(paste0(base.kallisto.path, "cells_x_features.mtx"))
xbp1.barcodes <- read.csv(paste0(base.kallisto.path, "cells_x_features.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(prefix.dictionary[[sample.name]], V1, "-1", suffix.dictionary[[sample.name]])) %>% pull(V1) 
xbp1.tag.names <- read.csv(paste0(base.kallisto.path, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
rownames(xbp1.data.original) <- xbp1.barcodes
colnames(xbp1.data.original) <- xbp1.tag.names


## "Standard" run
base.kallisto.path <- "~/VEXAS/XBP1_splicing/kallisto_standard/VX_BM15a/counts_unfiltered/"
xbp1.data.standard <- Matrix::readMM(paste0(base.kallisto.path, "cells_x_genes.mtx"))
xbp1.barcodes <- read.csv(paste0(base.kallisto.path, "cells_x_genes.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(prefix.dictionary[[sample.name]], V1, "-1", suffix.dictionary[[sample.name]])) %>% pull(V1) 
xbp1.tag.names <- read.csv(paste0(base.kallisto.path, "cells_x_genes.genes.txt"),  header = FALSE) %>% pull(V1)
rownames(xbp1.data.standard) <- xbp1.barcodes
colnames(xbp1.data.standard) <- xbp1.tag.names


## RNA velocity run
base.kallisto.path <- "~/VEXAS/XBP1_splicing/kallisto_rna_velocity_method/VX_BM15a/counts_unfiltered/"
xbp1.data.unspliced <- Matrix::readMM(paste0(base.kallisto.path, "unspliced.mtx"))
xbp1.barcodes <- read.csv(paste0(base.kallisto.path, "unspliced.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(prefix.dictionary[[sample.name]], V1, "-1", suffix.dictionary[[sample.name]])) %>% pull(V1) 
xbp1.tag.names <- read.csv(paste0(base.kallisto.path, "unspliced.genes.txt"),  header = FALSE) %>% pull(V1)
rownames(xbp1.data.unspliced) <- xbp1.barcodes
colnames(xbp1.data.unspliced) <- xbp1.tag.names
xbp1.data.spliced <- Matrix::readMM(paste0(base.kallisto.path, "spliced.mtx"))
xbp1.barcodes <- read.csv(paste0(base.kallisto.path, "spliced.barcodes.txt"), header = FALSE) %>% mutate(V1 = paste0(prefix.dictionary[[sample.name]], V1, "-1", suffix.dictionary[[sample.name]])) %>% pull(V1) 
xbp1.tag.names <- read.csv(paste0(base.kallisto.path, "spliced.genes.txt"),  header = FALSE) %>% pull(V1)
rownames(xbp1.data.spliced) <- xbp1.barcodes
colnames(xbp1.data.spliced) <- xbp1.tag.names

xbp1.data.rna.velocity <- data.frame(Unspliced = xbp1.data.unspliced[, "ENSG00000100219.16"], xbp1.data.spliced[, "ENSG00000100219.16"])

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- colnames(xbp1.data.spliced)
df$id <- NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "description"),values=genes,mart= mart)

## Compare

xbp1.data.original[1:5, ]
xbp1.data.standard[1:5, ]
hist(log10(xbp1.data.original[, "Spliced_XBP1"]))
hist(log10(xbp1.data.standard[, "Spliced_XBP1"]))
hist(log10(xbp1.data.rna.velocity[, "Unspliced"]))
