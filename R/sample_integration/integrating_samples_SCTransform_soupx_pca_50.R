library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(Azimuth)
library(gridExtra)
library(parallel)
library(DoubletFinder)
library(dsb)
set.seed(1)

## Load all the samples
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

print("Beginning filtering...")

## Run QC filtering
sample.list <- mclapply(sample.list, mc.cores = detectCores(), FUN = function(sample) {

  med_nGene <- median(sample$nFeature_RNA)
  sd_nGene <- mad(sample$nFeature_RNA)
  nGene.low <- max(med_nGene-3*sd_nGene, 250)
  nGene.high <- med_nGene+3*sd_nGene

  med_nUMI <- median(sample$nCount_RNA)
  sd_nUMI <- mad(sample$nCount_RNA)
  nUMI.low <- max(med_nUMI-3*sd_nUMI, 1000)
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
                              percent.mito>mito.low &
                              percent.mito<mito.high)
  return(sample.filtered)
})


print("Done with filtering")
print("Printing samples...")
for(sample in sample.list) {
  print(sample)
}


########################################## Removing previously-identified doublets #############################


doublets.df <- read_csv("data/qc_filtering/doubletfinder_output_20231103.csv") %>%
  column_to_rownames("BC")

print("Removing previously-defined doublets")
sample.list <- lapply(sample.list, function(x) {
  x <- AddMetaData(x, doublets.df)
  x <- subset(x, DF.classification == "Singlet")
  return(x)
})

print("Done with removing doublets")
print("Printing samples...")
for(sample in sample.list) {
  print(sample)
}

###################################### Merge the CD34+ and - compartment for each VEXAS sample #########################################


## Merge the CD34- and CD34+ compartments of each sample
sample.list.merged <- list()
sample.list.merged[["BM1"]] <- merge(sample.list$`VEXAS_BM1_CD34_plus_SR`, sample.list$`VEXAS_BM1_CD34-_SR`, project = "BM1")
sample.list.merged[["BM2"]] <- merge(sample.list$`VEXAS_c_BM2_CD34_plus_SR`, sample.list$`VEXAS_c_BM2_CD34-SR`, project = "BM2")
sample.list.merged[["BM14"]] <- merge(sample.list$VX_BM14_Fpos_GEX, sample.list$VX_BM14_Fneg_GEX, project = "BM14")
sample.list.merged[["B15a"]] <- sample.list$SG_VX_BM15a_GEX
sample.list.merged[["B15b"]] <- sample.list$SG_VX_BM15b_GEX ## Keep BM15 separate, since they do require merging based on looking at HSPC cluster
sample.list.merged[["B12"]] <- sample.list$VX_BM12_A_GEX
sample.list.merged[["B13"]] <- sample.list$VX_BM13_B_GEX
sample.list.merged[["UPN14"]] <- merge(sample.list$UPN14_CD34, sample.list$UPN14_BMMC)
sample.list.merged[["UPN15"]] <- merge(sample.list$UPN15_CD34, sample.list$UPN15_BMMC)
sample.list.merged[["UPN16"]] <- merge(sample.list$UPN16_CD34, sample.list$UPN16_BMMC)
sample.list.merged[["UPN17"]] <- merge(sample.list$UPN17_CD34, sample.list$UPN17_BMMC)

rm(sample.list)
sample.list <- sample.list.merged
rm(sample.list.merged)

print("Samples have been merged...")
for(sample in sample.list) {
  print(sample)
}

###################################### Adding cell cycle genes ###################################################

print("Adding cell cycle genes...")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sample.list <- mclapply(sample.list, mc.cores = detectCores(), FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
})


############################################# Normalize, merge, integrate ###################################

## Print all the info for all samples to check in
print("Printing samples prior to normalizing...")
for(sample in sample.list) {
  print(sample)
}

## Normalize and integrate the data using SCTransform
print("Normalizing w/ regression of percent.mito, S.Score, G2M.Score")
sample.list <- mclapply(X = sample.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, vst.flavor = "v2", vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"), method = "glmGamPoi", verbose = TRUE, seed.use = 1)
  DefaultAssay(x) <- "SCT"
  return(x)
}, mc.cores = detectCores())

print("Selecting integration features...")

integration.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = integration.features)
sample.list <- mclapply(X = sample.list, mc.cores = detectCores(), FUN = RunPCA, features = integration.features, seed.use = 1) ## Necessary to use rPCA

integration.anchors <- FindIntegrationAnchors(object.list = sample.list, 
                                              k.anchor = 2, 
                                              k.filter = 5, 
                                              k.score = 30, 
                                              reduction = "rpca",
                                              normalization.method = "SCT", 
                                              anchor.features = integration.features,
                                              dims = 1:50)
samples.vexas.integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT", dims = 1:50)

print("Integration complete! Starting downstream analysis...")

## Continue with downstream analysis
DefaultAssay(samples.vexas.integrated) <- "integrated"
samples.vexas.integrated <- RunPCA(samples.vexas.integrated, verbose = FALSE, seed.use = 1, npcs = 70)
ElbowPlot(samples.vexas.integrated, ndims = 70)
samples.vexas.integrated <- RunUMAP(samples.vexas.integrated, reduction = "pca", dims = 1:50, seed.use = 1, return.model = TRUE)
samples.vexas.integrated <- FindNeighbors(samples.vexas.integrated, reduction = "pca", dims = 1:50)
samples.vexas.integrated <- FindClusters(samples.vexas.integrated, resolution = 1)


## Make sure that the RNA assay is normalized and scaled as well
DefaultAssay(samples.vexas.integrated) <- "RNA"
samples.vexas.integrated <- NormalizeData(samples.vexas.integrated, assay = "RNA")


############################################ Add Azimuth cell type annotations #####################################

## Adding Azimuth predictions programmatically
print("Starting Azimuth...")
DefaultAssay(samples.vexas.integrated) <- "SCT"
samples.vexas.integrated <- RunAzimuth(query = samples.vexas.integrated, reference = "bonemarrowref")
samples.vexas.integrated@reductions$ref.umap <- NULL ## Remove the azimuth UMAP, it keeps defaulting to this while plotting downstream

p1 <- DimPlot(samples.vexas.integrated, group.by = "predicted.celltype.l2", label = T, label.box = T) + NoLegend()
ggsave("current_umap.pdf", plot = p1, device = "pdf")

## Adding bone marrow mapping (previously run)
mapping.data <- read_csv('data/BoneMarrowMap_annotations/querydata_projected_labeled_ccdiff_20230924.csv') %>%
  column_to_rownames("Cell")
samples.vexas.integrated <- AddMetaData(samples.vexas.integrated, mapping.data)
samples.vexas.integrated$predicted_CellType[is.na(samples.vexas.integrated$predicted_CellType)] <- "NA"


############################################ Add UMI-filtered genotyping counts #####################################

## Mapping dictionary
sample.dict <- list(
  NIH_1 = "UPN14_BMMC",
  NIH_2 = "UPN15_BMMC",
  NIH_3 = "UPN16_BMMC",
  NIH_4 = "UPN17_BMMC",
  NIH_5 = "UPN14_CD34",
  NIH_6 = "UPN15_CD34",
  NIH_7 = "UPN16_CD34",
  NIH_8 = "UPN17_CD34",
  VEXAS_BM2_CD34_minus = "VEXAS_c_BM2_CD34-SR",
  VEXAS_BM2_CD34_plus = "VEXAS_c_BM2_CD34_plus_SR",
  VEXAS_BM1_CD34_negative = "VEXAS_BM1_CD34-_SR",
  VEXAS_BM1_CD34_plus = "VEXAS_BM1_CD34_plus_SR",
  VEXAS_BM1_CD34_plus = "VEXAS_BM1_CD34_plus_SR"
)

## Read all the genotyping dfs
genotyping.files <- list.files("/gpfs/commons/home/rmurray/rscripts/running_GoT_UMI_filtering_VEXAS/data/gex/", recursive = T, full.names = T)
genotyping.dfs <- lapply(genotyping.files, function(x) {
  df <- read_csv(x) %>% 
    select(Barcode, orig.ident, Genotype_Approx_Match_No_Gene_Thresh)
  df$orig.ident <- plyr::mapvalues(df$orig.ident, from = names(sample.dict), to = unlist(sample.dict))
  df$Barcode <- paste0(df$orig.ident, "_", df$Barcode)
  return(df)
})
genotyping.df <- do.call(rbind, genotyping.dfs) %>% column_to_rownames("Barcode")

## Get an initial count of the genotyping fraction per sample
genotyping.df %>% 
  mutate(is.genotyped = !is.na(Genotype_Approx_Match_No_Gene_Thresh)) %>% 
  group_by(orig.ident) %>% 
  summarize(frac_genotyped = sum(is.genotyped) / n()) %>% 
  arrange(-frac_genotyped)

samples.vexas.integrated <- AddMetaData(samples.vexas.integrated, genotyping.df)

## Get the genotyping fraction per sample
samples.vexas.integrated@meta.data %>% 
  group_by(orig.ident) %>% 
  count(Genotype_Approx_Match_No_Gene_Thresh) %>% 
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped)


############################ Nicely format some fields so the object is figure-ready #######################

sample.mapping <- list(
  `VEXAS_BM1_CD34-_SR` = "VEXAS_BM1_CD34_minus",
  `VEXAS_BM1_CD34_plus_SR` = "VEXAS_BM1_CD34_plus",
  `VEXAS_c_BM2_CD34-SR` = "VEXAS_BM2_CD34_minus",
  `VEXAS_c_BM2_CD34_plus_SR` = "VEXAS_BM2_CD34_plus",
  `VX_BM12_A_GEX` = "VEXAS_BM12",
  `VX_BM13_B_GEX` = "VEXAS_BM13",
  `VX_BM14_Fpos_GEX` = "VEXAS_BM14_plus",
  `VX_BM14_Fneg_GEX` = "VEXAS_BM14_minus",
  `SG_VX_BM15a_GEX` = "VEXAS_BM15_A",
  `SG_VX_BM15b_GEX` = "VEXAS_BM15_B"
)

samples.vexas.integrated$Sample <- plyr::mapvalues(samples.vexas.integrated$orig.ident, from = names(sample.mapping), to = unlist(sample.mapping))
samples.vexas.integrated$Source <- if_else(grepl("UPN", samples.vexas.integrated$orig.ident), "NIH", "Landau")
samples.vexas.integrated$Sample <- gsub("UPN1", "PT2", samples.vexas.integrated$Sample) ## Give NIH samples the annotation PT2...
samples.vexas.integrated$Sample <- gsub("VEXAS_BM", "PT", samples.vexas.integrated$Sample) ## Give Landau samples the annotation PT1...

samples.vexas.integrated$Donor <- gsub("(_CD34_minus|_CD34_plus|_plus|_minus|_A|_B|-T.*|_CD34|_BMMC)$", "", samples.vexas.integrated$Sample)
samples.vexas.integrated$Donor <- plyr::mapvalues(samples.vexas.integrated$Donor, from = c("PT1", "PT2"), to = c("PT01", "PT02"))
samples.vexas.integrated$Donor <- factor(samples.vexas.integrated$Donor, levels = c("PT01", "PT02", "PT12", "PT13", "PT14", "PT15", "PT24", "PT25", "PT26", "PT27"))

samples.vexas.integrated$Genotype <- samples.vexas.integrated$Genotype_Approx_Match_No_Gene_Thresh
samples.vexas.integrated$Genotype[is.na(samples.vexas.integrated$Genotype)] <- "NA"
samples.vexas.integrated$Genotype <- factor(samples.vexas.integrated$Genotype, levels = c("WT", "MUT", "NA"))

samples.vexas.integrated$MDS <- "No"
samples.vexas.integrated$MDS[samples.vexas.integrated$Donor %in% c("PT02", "PT13", "PT14", "PT24", "PT27")] <- "Yes"


######################## Load ADT data for each sample and run DSB #########################

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
  real.bcs <- samples.vexas.integrated@meta.data %>% filter(orig.ident == sample.name) %>% rownames()
  
  # Filter both lists for only barcodes present in the ADT output
  empty.droplet.bcs <- empty.droplet.bcs[(empty.droplet.bcs %in% adt_bcs)]
  real.bcs <- real.bcs[(real.bcs %in% adt_bcs)]
  
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
samples.vexas.integrated.bcs <- samples.vexas.integrated@meta.data %>% rownames()
bc_missing_from_adt <- samples.vexas.integrated.bcs[!(samples.vexas.integrated.bcs %in% rownames(dsb.counts))]
print("Number of barcodes missing from dsb.counts data for:")
print(length(bc_missing_from_adt))
missing_matrix <- matrix(NA_real_, ncol=ncol(dsb.counts), nrow=length(bc_missing_from_adt), dimnames=list(bc_missing_from_adt, colnames(dsb.counts)))
dsb.counts <- rbind(dsb.counts, missing_matrix)

## Add as a new assay, normalize, and scale
dsb.counts <- dsb.counts[samples.vexas.integrated@meta.data %>% rownames, ]
samples.vexas.integrated[["ADT_DSB"]] <- CreateAssayObject(data = t(dsb.counts)) ## Important for DSB to add as the "data" slot
DefaultAssay(samples.vexas.integrated) <- "ADT_DSB"
samples.vexas.integrated <- ScaleData(samples.vexas.integrated)

## Correct CD90 for BM2
BM2.barcodes <- samples.vexas.integrated@meta.data %>% filter(orig.ident == "VEXAS_c_BM2_CD34_plus_SR" | orig.ident == "VEXAS_c_BM2_CD34-SR") %>% rownames()
samples.vexas.integrated@assays$ADT_DSB@data["Hu.CD90", BM2.barcodes] <- NA_real_
samples.vexas.integrated@assays$ADT_DSB@scale.data["Hu.CD90", BM2.barcodes] <- NA_real_


############################################ Load XBP1s data #################################

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

xbp1.counts <- lapply(samples.vexas.integrated$orig.ident %>% unique(), pull.xbp1.counts)
xbp1.counts.df <- do.call(rbind, xbp1.counts)

samples.vexas.integrated <- AddMetaData(samples.vexas.integrated, xbp1.counts.df)

## Add splicing ratio, fraction
samples.vexas.integrated$splicing_fraction <- samples.vexas.integrated$Spliced_XBP1 / (samples.vexas.integrated$Spliced_XBP1 + samples.vexas.integrated$Unspliced_XBP1)
samples.vexas.integrated$splicing_ratio <- samples.vexas.integrated$Spliced_XBP1 / samples.vexas.integrated$Unspliced_XBP1
samples.vexas.integrated$log10_splicing_ratio <- log10(samples.vexas.integrated$splicing_ratio)

######################################## Save ##############################################

print("Done! Saving...")
saveRDS(samples.vexas.integrated, file = "data/seurat_objects/vexas_final_RNA_data.rds")
print("Object is saved")


