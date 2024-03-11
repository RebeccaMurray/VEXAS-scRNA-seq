library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
theme_set(theme_classic())
set.seed(1)

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
current.plots.path <- c("figures/current_figure_drafts/")

################################ Perform a sample-aware permutation test ##################################

donors.keep <- rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1_umi_filtered)) %>% 
  group_by(Donor) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  filter(frac_genotyped > 0) %>% 
  pull(Donor)

md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% c("HSC")) %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function to calculate the mean difference
calculate_difference_of_means <- function(df) {
  df %>% 
  group_by(Genotype) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), unspliced_sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  select(Genotype, splicing_ratio) %>% 
  pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
  mutate(mean_diff = MUT - WT) %>% 
  pull(mean_diff)
}

## Shuffle the genotype labels and calculate within patients
md.xbp1.shuf <- md.xbp1
mean_diff_results <- sapply(1:1000, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_means(.)
})
hist(mean_diff_results)

## Calculate the actual mean difference
observed.mean.diff <- calculate_difference_of_means(md.xbp1)

sum(mean_diff_results > observed.mean.diff) / length(mean_diff_results)

