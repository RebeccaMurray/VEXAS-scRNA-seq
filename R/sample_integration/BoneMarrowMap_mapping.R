library(Seurat)
library(tidyverse)
library(symphony)
library(ggpubr)
library(patchwork)
library(BoneMarrowMap)

# Set directory to store projection reference files
projection_path = 'data/BoneMarrowMap_annotations/'

## Complete this section with something prior to Seurat v5 (I used seurat v4).
## Can use conda ENV signac_seurat_env_06_04_23
# 
# # Download Bone Marrow Reference - 344 Mb
# curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_SymphonyRef.rds', 
#                     destfile = paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# # Download uwot model file - 221 Mb
# curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
#                     destfile = paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot'))
# 
# # Load Symphony reference
# ref <- readRDS(paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# 
# # Set uwot path for UMAP projection
# ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')
# 
# # Visualize Bone Marrow reference
# ReferenceSeuratObj <- create_ReferenceObject(ref)
# DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', label=TRUE, label.size = 4)
# 
# # Save the seurat object
# saveRDS(ReferenceSeuratObj, paste0(projection_path, "ReferenceSeuratObj.rds"))

## Complete this section with Seurat v5
## Can use conda ENV signac_seurat_env_06_04_23_seurat_5

# Load Symphony reference - redo these two steps
ref <- readRDS(paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')

ReferenceSeuratObj <- readRDS(paste0(projection_path, "ReferenceSeuratObj.rds"))

## Load the current VEXAS seurat object
rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")

# batch variable to correct in the query data, set as NULL if no batches in query
batchvar <- 'orig.ident'

# Map query dataset using Symphony 
query <- map_Query(
  exp_query = rna.obj@assays$RNA@counts, 
  metadata_query = rna.obj@meta.data,
  ref_obj = ref,
  vars = batchvar
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref) 
plot_MappingErrorQC(query)

# Predict Hematopoietic Cell Types by KNN classification
query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap', group.by = c('predicted_CellType'), label=TRUE, label.size = 4, label.box = T, repel = T) + NoLegend()

# Predict Pseudotime values by KNN
query <- predict_Pseudotime(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_Pseudotime', 
  final_label = 'predicted_Pseudotime'
)
FeaturePlot(subset(query, mapping_error_QC == 'Pass'), features = c('predicted_Pseudotime'))

print("Query is done! Saving csv...")

# Save the output as a csv
save_ProjectionResults(
  query_obj = query, 
  celltype_label = 'predicted_CellType', 
  celltype_KNNprob_label = 'predicted_CellType_prob', 
  pseudotime_label = 'predicted_Pseudotime', 
  file_name = 'data/BoneMarrowMap_annotations/boneMarrowMap_annotations_20231108.csv')

# # Plot
# celltypes.exclude <- rna.obj@meta.data %>% count(predicted_CellType) %>% filter(n < 25) %>% pull(predicted_CellType)
# rna.obj$predicted_CellType_plot <- rna.obj$predicted_CellType
# rna.obj$predicted_CellType_plot[rna.obj$predicted_CellType_plot %in% celltypes.exclude] <- "NA"
# DimPlot(rna.obj, group.by = "predicted_CellType_plot", label = T, label.box = T, repel = T) + NoLegend()

