library(Seurat)
library(tidyverse)

rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes.rds")


############################################ Add UMI-filtered XBP1s counts #####################################

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
genotyping.files <- list.files("/gpfs/commons/home/rmurray/rscripts/running_GoT_UMI_filtering_VEXAS_XBP1/data/gex/", recursive = T, full.names = T) %>% grep("local_min", ., invert = T, value = T)
genotyping.dfs <- lapply(genotyping.files, function(x) {
  df <- read_csv(x) %>% 
    dplyr::select(Barcode, orig.ident, WT.calls.Approx_Match_No_Gene_Thresh, MUT.calls.Approx_Match_No_Gene_Thresh) %>% 
    dplyr::rename(Unspliced_XBP1_umi_filtered = WT.calls.Approx_Match_No_Gene_Thresh) %>% 
    dplyr::rename(Spliced_XBP1_umi_filtered = MUT.calls.Approx_Match_No_Gene_Thresh)
  df$orig.ident <- plyr::mapvalues(df$orig.ident, from = names(sample.dict), to = unlist(sample.dict))
  df$Barcode <- paste0(df$orig.ident, "_", df$Barcode)
  return(df)
})

genotyping.df <- do.call(rbind, genotyping.dfs) %>% column_to_rownames("Barcode")
## Get an initial count of the genotyping fraction per sample
genotyping.df %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1_umi_filtered)) %>% 
  group_by(orig.ident) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  arrange(-frac_genotyped)

rna.obj <- AddMetaData(rna.obj, genotyping.df)

## Get the genotyping fraction per sample
rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1_umi_filtered)) %>% 
  group_by(orig.ident) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  arrange(-frac_genotyped)

####################### Making boxpots ###################
# 
# ## Looking at monocytes
# rna.obj.HSC@meta.data %>% 
#   filter(!is.na(Unspliced_XBP1_umi_filtered)) %>%
#   # filter(seurat_clusters %in% c(3, 34, 37, 28)) %>% 
#   count(Genotype, orig.ident) %>% filter(Genotype != "NA")
# 
# rna.obj.HSC@meta.data %>% 
#   filter(Unspliced_XBP1_umi_filtered > 1) %>% 
#   mutate(spliced_ratio_umi_filtered = log2(Spliced_XBP1_umi_filtered / Unspliced_XBP1_umi_filtered)) %>% 
#   filter(Genotype %in% c("MUT", "WT"))  %>% 
#   ggplot(aes(x = spliced_ratio_umi_filtered, y = ..count..)) + 
#   geom_histogram() +
#   facet_grid(rows = vars(Donor))
#   
#   
# ## Just XBP1s reads
# rna.obj.HSC@meta.data %>% 
#   # filter(Donor == "UPN24") %>% 
#   filter(Unspliced_XBP1_umi_filtered > 0) %>% 
#   # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
#   mutate(spliced_ratio_umi_filtered = log2(Spliced_XBP1_umi_filtered / Unspliced_XBP1_umi_filtered)) %>% 
#   ggplot(aes(x = spliced_ratio_umi_filtered)) +
#   geom_histogram() +
#   facet_grid(cols = vars(Donor))
# 
# rna.obj.HSC@meta.data %>% 
#   filter(Donor != "UPN25") %>%
#   filter(Unspliced_XBP1 > 0) %>%
#   # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
#   # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
#   mutate(spliced_ratio = log2(Spliced_XBP1 / Unspliced_XBP1)) %>%
#   ggplot(aes(x = spliced_ratio)) +
#   geom_histogram() +
#   facet_grid(rows = vars(Donor), scales = "free")
# 
# 
# p1 <- rna.obj.HSC@meta.data %>% 
#     # filter(Donor == "UPN24") %>% 
#     filter(Unspliced_XBP1_umi_filtered > 0) %>% 
#     # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
#     mutate(spliced_ratio_umi_filtered = log2(Spliced_XBP1_umi_filtered / Unspliced_XBP1_umi_filtered)) %>% 
#     mutate(is_high_splicing = spliced_ratio_umi_filtered > -2) %>% 
#     ggplot(aes(x = is_high_splicing, y = `Hallmark IFNa Response1`)) +
#     geom_boxplot() +
#     facet_grid(cols = vars(Donor)) +
#     stat_compare_means()
# 
# p2 <- rna.obj.HSC@meta.data %>% 
#   # filter(Donor == "UPN24") %>% 
#   filter(Unspliced_XBP1 > 0) %>% 
#   # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
#   # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
#   mutate(spliced_ratio = log2(Spliced_XBP1 / Unspliced_XBP1)) %>%
#   mutate(is_high_splicing = spliced_ratio > -2) %>% 
#   ggplot(aes(x = is_high_splicing, y = `Hallmark IFNa Response1`)) +
#   geom_boxplot() +
#   facet_grid(cols = vars(Donor)) +
#   stat_compare_means(label = "p.format")
# 
# p1 | p2
# 
# 
# ## Previous
# test.table <- rna.obj.HSC@meta.data %>% 
#   group_by(Genotype) %>%
#   filter(Genotype != "NA") %>% 
#   filter(Donor == "UPN24") %>% 
#   summarize(
#     Unspliced_XBP1_total = sum(Unspliced_XBP1, na.rm = T),
#     Spliced_XBP1_total = sum(Spliced_XBP1, na.rm = T)
#   ) %>% 
#   column_to_rownames("Genotype")
# print(test.table)
# test.results <-  test.table %>% 
#   fisher.test(.)
# print("p_value:")
# print(test.results$p.value)
# print("odds ratio:")
# print(test.results$estimate)  
# 
# 
# ## UMI filtered
# test.table <- rna.obj.HSC@meta.data %>% 
#   group_by(Genotype) %>%
#   filter(Genotype != "NA") %>% 
#   filter(Donor == "UPN24") %>% 
#   summarize(
#     Unspliced_XBP1_total = sum(Unspliced_XBP1_umi_filtered, na.rm = T),
#     Spliced_XBP1_total = sum(Spliced_XBP1_umi_filtered, na.rm = T)
#   ) %>% 
#   column_to_rownames("Genotype")
# print(test.table)
# test.results <-  test.table %>% 
#   fisher.test(.)
# print("p_value:")
# print(test.results$p.value)
# print("odds ratio:")
# print(test.results$estimate)  
# 
# ################# Genotyping-agnostic analysis
# 
# 
# 
# 
