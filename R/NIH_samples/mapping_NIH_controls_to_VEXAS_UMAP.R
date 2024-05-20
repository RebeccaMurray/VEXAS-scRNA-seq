library(Seurat)
# library(sctransform)
library(tidyverse)
library("ggplot2")
library("viridis")
library("tibble")
library("gridExtra")
library("stringr")
# library(readr)
# library(ggpmisc)
library(ggrepel)
library(data.table)
library(patchwork)
# library(sleuthALR)
library(ggpubr)
# library(rstatix)
# library("scater")
library(stats)
library(ggh4x)
library(parallel)
# library(EnhancedVolcano)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

print("Starting integration script...")

vexas <- readRDS("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_only_sscore_g2mscore_mito_regressed_doublets_filtered_NormalizeData_cell_types.rds")
DefaultAssay(vexas)<-"RNA"
vexas=FindVariableFeatures(vexas, nfeatures = 3000)

vexas$Object='Vexas'

# CREATE NIH_CONTROLS OBJECTS
sra.ID.df <- read_csv("data/NIH_data/SraRunTable.txt") %>%
  filter(AvgSpotLen != 310) %>% ## Remove the TCR/BCR libraries
  filter(grepl("Healthydonor", `subject_id/status`)) %>% ## Keep only healthy donors
  select(`GEO_Accession (exp)`, `subject_id/status`, `tissue/cell_type`) %>% distinct() %>%
  mutate(suffix = `tissue/cell_type`) %>%
  mutate(suffix = gsub("FACS sorted Linjeage\\-CD34\\+ bone marrow", "_CD34", suffix)) %>%
  mutate(suffix = gsub("Human Bone Marrow", "_BMMC", suffix)) %>%
  mutate(sample_name = paste0(`subject_id/status`, suffix)) %>%
  select(`GEO_Accession (exp)`, sample_name)

print("Selected control IDs...")

sra.IDs <- sra.ID.df$sample_name
names(sra.IDs) <- sra.ID.df$`GEO_Accession (exp)`

nih.cellranger.output.dir <- "/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/"
nih.sample.list <- mclapply(names(sra.IDs), mc.cores = detectCores(), function(x) {
  print("Loading NIH sample:")
  print(x)
  nih.sample <- CreateSeuratObject(counts = Read10X(
    data.dir = paste0(nih.cellranger.output.dir, x, "/outs/filtered_feature_bc_matrix/")),
    project = sra.IDs[[x]], min.cells = 3, min.features = 200
  )
  return(nih.sample)
})
names(nih.sample.list) <- sra.IDs

print("Created controls seruat objects")

nih.sample.list <- mclapply(nih.sample.list, mc.cores = detectCores(), FUN = function(x) {
  
  print(x)
  
  # CREATE NIH_CONTROLS OBJECTS
  id = as.vector(unique(x$orig.ident))
  
  # QC NIH_CONTROLS OBJECTS
  x$percent.mito <- PercentageFeatureSet(x, pattern = "^MT-")
    med_nGene <- median(x$nFeature_RNA)
    sd_nGene <- mad(x$nFeature_RNA)
    nGene.low <- max(med_nGene-3*sd_nGene, 250)
    nGene.high <- med_nGene+3*sd_nGene

    med_nUMI <- median(x$nCount_RNA)
    sd_nUMI <- mad(x$nCount_RNA)
    nUMI.low <- max(med_nUMI-3*sd_nUMI, 500)
    nUMI.high <- med_nUMI+3*sd_nUMI

    mito.low <- 0
    mito.high <- 10

    print("\n")
    print(x@project.name)
    print("nGene.low")
    print(nGene.low)
    print("nGene high")
    print(nGene.high)
    print("nUMI low")
    print(nUMI.low)
    print("nUMI high")
    print(nUMI.high)
    print("\n")

    x.filtered <- subset(x, subset = nFeature_RNA>nGene.low &
                                nFeature_RNA<nGene.high &
                                nCount_RNA > nUMI.low &
                                nCount_RNA < nUMI.high &
                                percent.mito>mito.low &
                                percent.mito<mito.high)
  VlnPlot(x, features = c('nCount_RNA','nFeature_RNA','percent.mito'), pt.size = 0.01)
  
  
  # PRE-PROCESS NIH_CONTROLS
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  #x <- SCTransform(x, vst.flavor = "v2", vars.to.regress = c("percent.mito", "CC.Difference"))
  x <- FindVariableFeatures(x, nfeatures=2000)
  # x <- ScaleData(x, vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"))
  # x <- RunPCA(x, npcs = 30)
  
  
  gc()
  
  # PROJECTION
  anchors <- FindTransferAnchors(reference = vexas, reference.assay = 'integrated', query = x,
                                 dims = 1:30, reference.reduction = "pca", normalization.method = 'LogNormalize', query.assay = 'RNA')
  
  predictions <- TransferData(reference = vexas, query = x, anchorset = anchors, refdata = vexas$cluster_celltype,
                              dims = 1:30)
  
  x <- MapQuery(anchorset = anchors, reference = vexas, query = x,
                refdata = list(cluster_celltype = "cluster_celltype"), reference.reduction = "pca", reduction.model = "umap")
  
  
  gc()
  # MERGE DATASETS
  
  x$Object='CONTROL'
  
  merged = merge(vexas, x)
  
  merged@reductions[["umap"]]=vexas@reductions[["umap"]]
  
  merged@reductions[["umap"]]@cell.embeddings = rbind(vexas@reductions[["umap"]]@cell.embeddings,x@reductions[["ref.umap"]]@cell.embeddings)
  
  vexas=merged
  
  rm(merged)
  
  vexas$Allcluster_celltype = ifelse(is.na(vexas$cluster_celltype),vexas$predicted.cluster_celltype, vexas$cluster_celltype)
  
  # saveRDS(x, file=paste0(id,'_processed.Rds'))
  
  gc()
  
  return(x)
})


gc()

merged = merge(vexas, nih.sample.list)

table(merged$Object)  

umap_df=rbind(vexas@reductions[["umap"]]@cell.embeddings,HD1@reductions[["ref.umap"]]@cell.embeddings,HD2@reductions[["ref.umap"]]@cell.embeddings,HD3@reductions[["ref.umap"]]@cell.embeddings,HD4@reductions[["ref.umap"]]@cell.embeddings)

merged[["umap"]] <- SeuratObject::CreateDimReducObject(umap_df, assay="RNA")

merged$Allcluster_celltype = ifelse(is.na(merged$cluster_celltype),merged$predicted.cluster_celltype, merged$cluster_celltype)

saveRDS(merged, file='data/seurat_objects/vexas.ctrl.integrated.Rds')

p1=DimPlot(merged, group.by = 'cluster_celltype',label = T, label.box = T)+NoLegend()+ggtitle('VEXAS + NIH-CONTROLS')
p2=DimPlot(merged, group.by = 'Object',label = T, label.box = T)+NoLegend()+ggtitle('VEXAS + NIH-CONTROLS')
p3=DimPlot(merged, group.by = 'orig.ident',label = T, label.box = T)+NoLegend()+ggtitle('VEXAS + NIH-CONTROLS')

p1+p2
ggsave('figures/vexas.ctrl.integrated.Object.pdf',width=15,height=10)

p1+p3
ggsave('figures/vexas.ctrl.integrated.Ident.pdf',width=15,height=10)


gc()


