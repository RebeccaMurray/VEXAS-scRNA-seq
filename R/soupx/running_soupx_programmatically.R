## Conda ENV: myR4_RMM_3_17_23_soupx

library(Seurat)
library(SoupX)
library(DropletUtils)
library(tidyverse)
library(parallel)

## Load sample names, cellranger out path for landau VEXAS samples
landau.sample.list <- list.dirs("/gpfs/commons/home/rmurray/VEXAS/", recursive = F, full.names = F) %>% grep("^[VEXAS_|VX_]", ., value = T) %>% 
  grep("GOT|XBP1", ., value = T, invert = T) %>% 
  grep("BM3", ., value = T, invert = T)
landau.dir.list <- paste0("/gpfs/commons/home/rmurray/VEXAS/", landau.sample.list)
names(landau.dir.list) <- landau.sample.list

## Load sample names, cellranger out path for NIH VEXAS samples
nih.sample.df <- read_csv("data/NIH_data/SraRunTable.txt") %>%
  filter(AvgSpotLen != 310) %>% ## Remove the TCR/BCR libraries
  # filter(grepl("UPN", `subject_id/status`)) %>% ## Keep only patients for now
  # filter(`subject_id/status` %in% c("UPN14", "UPN15", "UPN16", "UPN17")) %>% ## Keep only 5' samples for now
  filter(`subject_id/status` %in% c("Healthydonor1", "Healthydonor2", "Healthydonor3", "Healthydonor4")) %>% 
  select(`GEO_Accession (exp)`, `subject_id/status`, `tissue/cell_type`) %>% distinct() %>%
  mutate(suffix = `tissue/cell_type`) %>%
  mutate(suffix = gsub("FACS sorted Linjeage\\-CD34\\+ bone marrow", "_CD34", suffix)) %>%
  mutate(suffix = gsub("Human Bone Marrow", "_BMMC", suffix)) %>%
  mutate(sample_name = paste0(`subject_id/status`, suffix)) %>%
  # filter(`GEO_Accession (exp)` != "GSM5858678") %>% ## Temporary while we finish processing this sample
  select(`GEO_Accession (exp)`, sample_name)
nih.dir.list <- paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/", nih.sample.df$`GEO_Accession (exp)`)
names(nih.dir.list) <- nih.sample.df$sample_name

# dir.list <- c(landau.dir.list, nih.dir.list)
dir.list <- c(nih.dir.list)
dir.list

run_soupx_on_cellranger_dir <- function(sample.name) {
  
  sample.path <- dir.list[[sample.name]]
  
  print("Running for...")
  print(sample.name)
  
  print("Path is...")
  print(sample.path)
  
  pdf(paste0("figures/soupx/", sample.name, ".pdf"))
  
  filt.matrix <- Read10X_h5(paste0(sample.path, "/outs/filtered_feature_bc_matrix.h5"))
  raw.matrix <- Read10X_h5(paste0(sample.path, "/outs/raw_feature_bc_matrix.h5"))
  
  ## Perform a quick clustering of just this pool
  srat <- CreateSeuratObject(counts = filt.matrix)
  srat <- SCTransform(srat, verbose = F)
  srat <- RunPCA(srat, verbose = F)
  srat <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat <- FindClusters(srat, verbose = T)
  
  ## Create the soupx object
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)
  
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  cluster.maps <- setNames(meta$seurat_clusters, rownames(meta))
  soup.channel <- setClusters(soup.channel, cluster.maps)
  soup.channel  <- setDR(soup.channel, umap)
  head(meta)
  
  soup.channel  <- autoEstCont(soup.channel)
  
  print(head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20))
  
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  output.dir <- paste0("data/soupx_output/", sample.name)
  dir.create(output.dir)
  DropletUtils:::write10xCounts(paste0(output.dir, "/", sample.name), adj.matrix)
  
  print("done!\n\n")
  
  dev.off()
  
}

mclapply(names(dir.list), run_soupx_on_cellranger_dir, mc.cores = detectCores())
