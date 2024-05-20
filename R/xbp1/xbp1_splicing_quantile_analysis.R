library(Seurat)
library(tidyverse)

# hsc.clusters <- c(8, 13, 5, 0, 31, 40)
hsc.clusters <- c(8, 13, 5, 0, 31, 40, 21, 15, 6)
rna.obj.HSC <- subset(rna.obj, seurat_clusters %in% hsc.clusters)

## Splicing fracrion
rna.obj@meta.data %>% 
  filter(seurat_clusters %in% hsc.clusters) %>% 
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 1) %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  ggplot(aes(x = splicing_fraction, y = ..count..)) +
  geom_histogram() +
  facet_grid(rows = vars(Donor), cols = vars(Genotype), scales = "free")

## SPlicing ratio
mean.vals <- rna.obj@meta.data %>% 
  filter(seurat_clusters %in% hsc.clusters) %>% 
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 1) %>% 
  mutate(log2_splicing_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>% 
  group_by(Donor) %>% 
  summarize(median_log2_splicing_ratio = median(log2_splicing_ratio))

rna.obj@meta.data %>% 
  filter(seurat_clusters %in% hsc.clusters) %>% 
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 1) %>% 
  mutate(log2_splicing_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>% 
  group_by(Donor) %>% 
  mutate(median_log2_splicing_ratio = median(log2_splicing_ratio)) %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  ggplot(aes(x = log2_splicing_ratio, y = ..count..)) +
  geom_histogram() +
  facet_grid(rows = vars(Donor), cols = vars(Genotype), scales = "free") +
  ggplot(data = mean.vals) +
  geom_vline(xintercept = median_log2_splicing_ratio)


rna.obj.HSC@meta.data %>% 
  filter(Genotype == "MUT") %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  mutate(log2_splicing_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>% 
  group_by(Donor) %>% 
  mutate(avg_log2_splicing = mean(log2_splicing_ratio)) %>% 
  mutate(is_high_splicing = log2_splicing_ratio > avg_log2_splicing) %>% 
  ggplot(aes(x = is_high_splicing, y = `ATF4 target genes\n(Han et al.)1`)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.HSC@meta.data %>% 
  filter(Genotype == "MUT") %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  mutate(log2_splicing_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>% 
  group_by(Donor) %>% 
  mutate(avg_log2_splicing = mean(log2_splicing_ratio)) %>% 
  mutate(is_high_splicing = log2_splicing_ratio > avg_log2_splicing) %>% 
  ggplot(aes(x = is_high_splicing, y = `GOBP IRE11`)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

## Genotype-agnostic
rna.obj.HSC@meta.data %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  mutate(log2_splicing_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>% 
  group_by(Donor) %>% 
  mutate(avg_log2_splicing = mean(log2_splicing_ratio)) %>% 
  mutate(is_high_splicing = log2_splicing_ratio > avg_log2_splicing) %>% 
  ggplot(aes(x = is_high_splicing, y = `Hallmark IFNa Response1`)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

## Just XBP1s reads
rna.obj@meta.data %>% 
  filter(seurat_clusters %in% hsc.clusters) %>% 
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  ggplot(aes(x = log10(Spliced_XBP1), y = ..count..)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(rows = vars(Donor), cols = vars(Genotype), scales = "free")

## Looking at just HAS xbp1s reads or not
rna.obj@meta.data %>% 
  filter(seurat_clusters %in% hsc.clusters) %>% 
  filter(Genotype == "MUT") %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 1) %>%
  mutate(has_xbp1 = Spliced_XBP1 > 0) %>% 
  count(Donor, has_xbp1)

rna.obj.HSC@meta.data %>% 
  filter(Genotype == "MUT") %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 1) %>%
  mutate(has_xbp1s = Spliced_XBP1 > 0) %>% 
  ggplot(aes(x = has_xbp1s, y = `Hallmark IFNa Response1`)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means()

rna.obj.HSC@meta.data %>% 
  filter(Genotype == "MUT") %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  filter(Unspliced_XBP1 > 1) %>%
  mutate(has_xbp1s = Spliced_XBP1 > 0) %>% 
  ggplot(aes(x = has_xbp1s, y = `Hallmark IFNa Response1`)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means()



## Looking at just MUT cells
rna.obj.HSC@meta.data %>% 
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  mutate(log2_splicing_ratio = log2(Spliced_XBP1/Unspliced_XBP1)) %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  group_by(Donor) %>% 
  summarize(splicing_median = median(log2_splicing_ratio))
  # summarize(xbp1_quantile = quantile(log2_splicing_ratio, c(0.5, 0.8)), q = c(0.2, 0.8))
  mutate(log2_splicing_ratio = log2(Spliced_XBP1/Unspliced_XBP1)) %>% 
  ggplot(aes(x = log2_splicing_ratio, y = ..count..)) +
  geom_histogram() +
  facet_grid(rows = vars(Donor), cols = vars(Genotype), scales = "free")


################## Genotyping-agnostic


df.medians <- rna.obj.HSC@meta.data %>% 
  filter(Donor != "UPN25") %>%
  filter(Unspliced_XBP1 > 0) %>%
  # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
  # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  group_by(Donor) %>% 
  summarize(mean_spliced_ratio = median(spliced_ratio))
  
rna.obj.HSC@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
  # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  ggplot(aes(x = spliced_ratio)) +
  geom_histogram() +
  facet_grid(rows = vars(Donor), scales = "free") +
  geom_vline(aes(xintercept = mean_spliced_ratio, group = Donor), color = "green", linetype = 'longdash') 

rna.obj.HSC@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `GOBP IRE11`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.HSC@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `XBP1 target module1`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.HSC@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `XBP1_011`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")
  
################## mut cells only
  rna.obj.HSC@meta.data %>% 
    filter(Donor != "UPN25") %>%
    filter(Genotype == "MUT") %>% 
    filter(Unspliced_XBP1 > 0) %>%
    count(Donor)
  
  df.medians <- rna.obj.HSC@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Genotype == "MUT") %>% 
    filter(Unspliced_XBP1 > 0) %>% 
    # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
    # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    group_by(Donor) %>% 
    summarize(mean_spliced_ratio = mean(spliced_ratio))
  
  rna.obj.HSC@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Genotype == "MUT") %>% 
    filter(Unspliced_XBP1 > 0) %>%
    merge(., df.medians, by = "Donor") %>% 
    # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
    # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    ggplot(aes(x = spliced_ratio)) +
    geom_histogram() +
    facet_grid(rows = vars(Donor), scales = "free") +
    geom_vline(aes(xintercept = mean_spliced_ratio, group = Donor), color = "green", linetype = 'longdash') 
  
  rna.obj.HSC@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Unspliced_XBP1 > 0) %>% 
    filter(Genotype == "MUT") %>% 
    group_by(Donor) %>% 
    filter(n() > 100) %>% 
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    merge(., df.medians, by = "Donor") %>% 
    mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
    mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
    ggplot(aes(x = splicing_category, y = `Hallmark IFNa Response1`, fill = splicing_category)) +
    geom_boxplot() +
    facet_grid(cols = vars(Donor)) +
    stat_compare_means(label = "p.format")
  
  rna.obj.HSC@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Unspliced_XBP1 > 0) %>% 
    filter(Genotype == "MUT") %>% 
    group_by(Donor) %>% 
    filter(n() > 100) %>% 
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    merge(., df.medians, by = "Donor") %>% 
    mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
    mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
    ggplot(aes(x = splicing_category, y = `XBP1 target module1`, fill = splicing_category)) +
    geom_boxplot() +
    facet_grid(cols = vars(Donor)) +
    stat_compare_means(label = "p.format")
  
  rna.obj.HSC@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Unspliced_XBP1 > 0) %>% 
    filter(Genotype == "MUT") %>% 
    group_by(Donor) %>% 
    filter(n() > 100) %>% 
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    merge(., df.medians, by = "Donor") %>% 
    mutate(splicing_category = if_else(spliced_ratio > mean_spliced_ratio, "High", "Low")) %>% 
    mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
    ggplot(aes(x = splicing_category, y = `XBP1_011`, fill = splicing_category)) +
    geom_boxplot() +
    facet_grid(cols = vars(Donor)) +
    stat_compare_means(label = "p.format")
  
  
  
  