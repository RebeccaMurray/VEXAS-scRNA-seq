library(tidyverse)
library(Seurat)
library(patchwork)
library(parallel)


###################################### Loading the CD34+ NIH samples from scratch ##############################

## Load the NIH samples
sample.list.dict <- list(
  # `GSM5858665` = "UPN6_CD34",
  # `GSM5858666` = "UPN11_CD34",
  # `GSM5858667` = "UPN1_CD34",
  # `GSM5858668` = "UPN10_CD34",
  # `GSM5858669` = "UPN13_CD34",
  `GSM5858670` = "UPN14_CD34",
  `GSM5858671` = "UPN15_CD34",
  `GSM5858672` = "UPN16_CD34",
  `GSM5858673` = "UPN17_CD34"
)

se.list <- lapply(names(sample.list.dict), function(x) {
  cellranger.data <- Read10X(data.dir = paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/", x, "/outs/filtered_feature_bc_matrix/"))
  se <- CreateSeuratObject(counts = cellranger.data, project = sample.list.dict[[x]], min.cells = 0, min.features = 0)
  se[["percent.mt"]] <- PercentageFeatureSet(se, pattern = "^MT-")
  
  ## Add the cb_sniffer results
  cb.sniffer.output <- read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/", x, "_counts_CB.tsv"))
  mut.barcodes <- cb.sniffer.output %>% filter(alt_count > 0) %>% pull(barcode)
  wt.barcodes <- cb.sniffer.output %>% filter(alt_count == 0 & ref_count > 0) %>% pull(barcode)
  cb.sniffer.output <- cb.sniffer.output %>% 
    mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
    select(barcode, genotype_status_cbsniff) %>% 
    column_to_rownames("barcode")
  
  if (nrow(cb.sniffer.output) > 0) { ## one sample has 0 reads in GEX covering UBA1
    se <- AddMetaData(se, cb.sniffer.output) 
  }
  
  return(se)
})
names(se.list) <- unlist(sample.list.dict)


## Print the genotyping counts for each sample
lapply(se.list, function(x) {
  print(x@project.name)
  if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
    x@meta.data %>% 
      count(genotype_status_cbsniff) %>% print()
    print("genotyping fraction:")
    print(x@meta.data %>% filter(!is.na(genotype_status_cbsniff)) %>% nrow() /  x@meta.data %>% nrow())
  } else {
    print("No reads, continuing...")
  }
})

# Visualize QC metrics as a violin plot
p1 <- VlnPlot(se.list[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(se.list[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(se.list[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(se.list[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

p1 | p2
p3 | p4

## Look at UBA1 expression
lapply(se.list, function(x) {
  print(x@project.name)
  uba1.expression.df <- data.frame(UBA1.counts = x@assays$RNA@counts["UBA1", ])
  x <- AddMetaData(x, uba1.expression.df, col.name = "UBA1_UMI_count")
  print("Fraction of cells with non-zero UBA1 expression:")
  print(x@meta.data %>% filter(UBA1_UMI_count > 0) %>% nrow())
  print(x@meta.data %>% filter(UBA1_UMI_count > 0) %>% nrow() / x@meta.data %>% nrow())
  
  # if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
  #   x@meta.data %>% 
  #     mutate(nonzero_UBA1_expression = UBA1_UMI_count > 0) %>% 
  #     mutate(has_cbsniffer_genotype = genotype_status_cbsniff %in% c("MUT", "WT")) %>% 
  #     count(nonzero_UBA1_expression, has_cbsniffer_genotype) %>% 
  #     print()
  #   
  #   print("") 
  # }
  
})


###################################### Loading the CD34+ NIH samples from provided count matrices ##############################

## Load reprocessed count matrices
load("~/rscripts/VEXAS_RNA_seurat/data/NIH_data/CD34.rawData.RData")
cd34.rawData <- rawData
cd34.se <- CreateSeuratObject(counts = cd34.rawData, project = "NIH_CD34")
# load("~/rscripts/VEXAS_RNA_seurat/data/NIH_data/BM.rawData.RData")
# BM.rawData <- rawData
# rm(rawData)
# BM.se <- CreateSeuratObject(counts = BM.rawData, project = "NIH_BM")

# metadata.df.1 <- read_csv("data/NIH_data/metaDatatSNECellType_ALiceManual.csv")
metadata.df.2 <- read_csv("data/NIH_data/metaDatatSNECellType_ALiceManual (1).csv")

cd34.se <- AddMetaData(cd34.se, metadata.df.2)
cd34.se[["percent.mt"]] <- PercentageFeatureSet(cd34.se, pattern = "^MT-")
# BM.se <- AddMetaData(BM.se, metadata.df.1)

# cd34.se <- subset(cd34.se, group == "PT")
table(cd34.se$group)
table(cd34.se$subject)

## Look for barcode formatting issues? - the "-1" at the end is sometimes different numbers
length(rownames(md))
rownames(md) %>% length()
rownames(md) %>% unique() %>% length()
md$orig.ident %>% gsub("^[A|C|G|T]{16}", "", .) %>% unique()
md$barcode_unformatted %>% gsub("^[A|C|G|T]{16}", "", .) %>% unique()
md %>% 
  mutate(suffix = gsub("^[A|C|G|T]{16}", "", orig.ident)) %>% 
  count(subject, suffix)


## Identifying samples
md <- cd34.se@meta.data %>% 
  mutate(barcode_unformatted = gsub("-[0-9]*$", "-1", orig.ident))

## Printing for each sample loaded from fastqs, how many barcode overlap it has
## With each sample identifier in the count matrix seurat object
lapply(se.list, function(x) {
  print(x@project.name)
  barcode.list <- x@meta.data %>% rownames
  md$has_current_sample_bc <- md$barcode_unformatted %in% barcode.list
  md %>% 
    filter(has_current_sample_bc == T) %>% 
    count(subject, has_current_sample_bc) %>% 
    arrange(-n) %>% 
    head(n=5) %>% 
    print()
})

## UPN14_CD34 = PT6
## UPN15_CD34 = PT7
## UPN16_CD34 = PT8
## UPN17_CD34 = PT9


###################################### Adding QC filtering #######################################################

## Filtering instructions from their manuscript
se.list.filtered <- lapply(se.list, function(sample) {
  
  sample.filtered <- subset(sample,
                            subset = nFeature_RNA > 500 &
                              nFeature_RNA < 3000 &
                              percent.mt > 0 &
                              percent.mt < 10)
  return(sample.filtered)
})

lapply(c(se.list.filtered$UPN14_CD34, se.list.filtered$UPN15_CD34, se.list.filtered$UPN16_CD34, se.list.filtered$UPN17_CD34), function(x) {
  print(x@project.name)
  if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
    x@meta.data %>% 
      count(genotype_status_cbsniff) %>% print()
    print("genotyping fraction:")
    print(x@meta.data %>% filter(!is.na(genotype_status_cbsniff)) %>% nrow() /  x@meta.data %>% nrow())
  } else {
    print("No reads, continuing...")
  }
})



## Filtering used in my rna integration script
bm15.list.filtered <- lapply(c(vx.BM15a, vx.BM15b), function(sample) {
  
  med_nGene <- median(sample$nFeature_RNA)
  sd_nGene <- mad(sample$nFeature_RNA)
  nGene.low <- max(med_nGene-3*sd_nGene, 250)
  nGene.high <- med_nGene+3*sd_nGene
  
  med_nUMI <- median(sample$nCount_RNA)
  sd_nUMI <- mad(sample$nCount_RNA)
  nUMI.low <- max(med_nUMI-3*sd_nUMI, 500)
  nUMI.high <- med_nUMI+3*sd_nUMI
  
  mito.low <- 0
  mito.high <- 10
  
  print("\n")
  print(sample@project.name)
  print("nGene.low")
  print(nGene.low)
  print("nGene high")
  print(nGene.high)
  print("nUMI low")
  print(nUMI.low)
  print("nUMI high")
  print(nUMI.high)
  print("\n")
  
  sample.filtered <- subset(sample,
                            subset = nFeature_RNA>nGene.low &
                              nFeature_RNA<nGene.high &
                              nCount_RNA > nUMI.low &
                              nCount_RNA < nUMI.high &
                              percent.mt >mito.low &
                              percent.mt < mito.high)
  return(sample.filtered)
})


## Print the genotyping counts for each sample
lapply(c(se.list$UPN14_CD34, se.list$UPN15_CD34, se.list$UPN16_CD34, se.list$UPN17_CD34), function(x) {
  print(x@project.name)
  if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
    x@meta.data %>% 
      count(genotype_status_cbsniff) %>% print()
    print("genotyping fraction:")
    print(x@meta.data %>% filter(!is.na(genotype_status_cbsniff)) %>% nrow() /  x@meta.data %>% nrow())
  } else {
    print("No reads, continuing...")
  }
})

lapply(bm15.list.filtered, function(x) {
  print(x@project.name)
  if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
    x@meta.data %>% 
      count(genotype_status_cbsniff) %>% print()
    print("genotyping fraction:")
    print(x@meta.data %>% filter(!is.na(genotype_status_cbsniff)) %>% nrow() /  x@meta.data %>% nrow())
  } else {
    print("No reads, continuing...")
  }
})

###################################### Comparing with Landau Lab fresh sample (BM15) ##############################


## Comparing with BM15
load("~/rscripts/VEXAS_v2_CD34_pos/data/gex/180/VX_BM15a_MUT_0_WT_0/VX_BM15a_gex_UMI_filtered_UMI_threshold_180_MUT_0_WT_0.Rdata")
vx.BM15a <- gex_UMI_filtered
load("~/rscripts/VEXAS_v2_CD34_pos/data/gex/100/VX_BM15b_MUT_0_WT_0/VX_BM15b_gex_UMI_filtered_UMI_threshold_100_MUT_0_WT_0.Rdata")
vx.BM15b <- gex_UMI_filtered

x = "SG_VX_BM15a_GEX"
cb.sniffer.output <- read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/", x, "_counts_CB.tsv"))
mut.barcodes <- cb.sniffer.output %>% filter(alt_count > 0) %>% pull(barcode)
wt.barcodes <- cb.sniffer.output %>% filter(alt_count == 0 & ref_count > 0) %>% pull(barcode)
cb.sniffer.output <- cb.sniffer.output %>% 
  mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
  select(barcode, genotype_status_cbsniff) %>% 
  column_to_rownames("barcode")
vx.BM15a <- AddMetaData(vx.BM15a, cb.sniffer.output)
vx.BM15a[["percent.mt"]] <- PercentageFeatureSet(vx.BM15a, pattern = "^MT-")

x = "SG_VX_BM15b_GEX"
cb.sniffer.output <- read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/", x, "_counts_CB.tsv"))
mut.barcodes <- cb.sniffer.output %>% filter(alt_count > 0) %>% pull(barcode)
wt.barcodes <- cb.sniffer.output %>% filter(alt_count == 0 & ref_count > 0) %>% pull(barcode)
cb.sniffer.output <- cb.sniffer.output %>% 
  mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
  select(barcode, genotype_status_cbsniff) %>% 
  column_to_rownames("barcode")
vx.BM15b <- AddMetaData(vx.BM15b, cb.sniffer.output)
vx.BM15b[["percent.mt"]] <- PercentageFeatureSet(vx.BM15b, pattern = "^MT-")

lapply(c(vx.BM15a, vx.BM15b), function(x) {
  print(x@project.name)
  if ("genotype_status_cbsniff" %in% colnames(x@meta.data)) {
    x@meta.data %>% 
      count(genotype_status_cbsniff) %>% print()
    print("genotyping fraction:")
    print(x@meta.data %>% filter(!is.na(genotype_status_cbsniff)) %>% nrow() /  x@meta.data %>% nrow())
  } else {
    print("No reads, continuing...")
  }
  
})

## Look at UBA1 expression
lapply(c(vx.BM15a, vx.BM15b), function(x) {
  print(x@project.name)
  uba1.expression.df <- data.frame(UBA1.counts = x@assays$RNA@counts["UBA1", ])
  x <- AddMetaData(x, uba1.expression.df, col.name = "UBA1_UMI_count")
  print("Fraction of cells with non-zero UBA1 expression:")
  print(x@meta.data %>% filter(UBA1_UMI_count > 0) %>% nrow())
  print(x@meta.data %>% filter(UBA1_UMI_count > 0) %>% nrow() / x@meta.data %>% nrow())
  print("\n")
  # print(x@meta.data %>% count(UBA1_UMI_count))
})
