library(Seurat)
library(tidyverse)


rna.obj$spliced_ratio <- log2((rna.obj$Spliced_XBP1_umi_filtered + 0.1)/(rna.obj$Unspliced_XBP1_umi_filtered + 0.1))
rna.obj.sub <- rna.obj.HSC
rna.obj.sub <- subset(rna.obj.HSC, cluster_celltype == "HSC")

################## Genotyping-agnostic - any splicing vs none

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  count(Donor, has_XBP1s)

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  ggplot(aes(x = has_XBP1s, y = `GOBP IRE11`, fill = has_XBP1s)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`Yes` = "lightyellow2", `No` = "grey80"))

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  ggplot(aes(x = has_XBP1s, y = `Hallmark IFNa Response1`, fill = has_XBP1s)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`Yes` = "lightyellow2", `No` = "grey80"))

rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(Spliced_XBP1_umi_filtered > 0, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `Hallmark Inflammation1`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(Spliced_XBP1_umi_filtered > 0, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `ATF4 target genes1`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")




################## MUT cells only - any splicing vs none

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  count(Donor, has_XBP1s)

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  ggplot(aes(x = has_XBP1s, y = `GOBP IRE11`, fill = has_XBP1s)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`Yes` = "lightyellow2", `No` = "grey80"))

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  ggplot(aes(x = has_XBP1s, y = `Hallmark IFNa Response1`, fill = has_XBP1s)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`Yes` = "lightyellow2", `No` = "grey80"))

rna.obj.sub@meta.data %>% 
  # filter(!Donor %in% c("UPN25")) %>%
  filter(Genotype == "MUT") %>% 
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(has_XBP1s = if_else(Spliced_XBP1_umi_filtered > 0, "Yes", "No")) %>% 
  mutate(has_XBP1s = factor(has_XBP1s, levels = c("No", "Yes"))) %>% 
  ggplot(aes(x = has_XBP1s, y = `Hallmark Inflammation1`, fill = has_XBP1s)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`Yes` = "lightyellow2", `No` = "grey80"))


################## Genotyping-agnostic

df.medians <- rna.obj.sub@meta.data %>% 
  filter(Donor != "UPN25") %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  group_by(Donor) %>% 
  summarize(median_spliced_ratio = median(spliced_ratio))
  
rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1_umi_filtered > 0) %>%
  merge(., df.medians, by = "Donor") %>% 
  ggplot(aes(x = spliced_ratio)) +
  geom_histogram() +
  facet_grid(rows = vars(Donor), scales = "free") +
  geom_vline(aes(xintercept = median_spliced_ratio, group = Donor), color = "green", linetype = 'longdash') 

rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > median_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `GOBP IRE11`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > median_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `XBP1 target module1`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")

rna.obj.sub@meta.data %>% 
  filter(!Donor %in% c("UPN25")) %>%
  filter(Unspliced_XBP1 > 0) %>%
  mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
  merge(., df.medians, by = "Donor") %>% 
  mutate(splicing_category = if_else(spliced_ratio > median_spliced_ratio, "High", "Low")) %>% 
  mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
  ggplot(aes(x = splicing_category, y = `XBP1_011`, fill = splicing_category)) +
  geom_boxplot() +
  facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format")
  
################## mut cells only
  rna.obj.sub@meta.data %>% 
    filter(Donor != "UPN25") %>%
    filter(Genotype == "MUT") %>% 
    filter(Unspliced_XBP1 > 0) %>%
    count(Donor)
  
  df.medians <- rna.obj.sub@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Genotype == "MUT") %>% 
    filter(Unspliced_XBP1 > 0) %>% 
    group_by(Donor) %>% 
    summarize(mean_spliced_ratio = median(spliced_ratio))
  
  rna.obj.sub@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Genotype == "MUT") %>% 
    filter(Spliced_XBP1_umi_filtered > 0) %>%
    merge(., df.medians, by = "Donor") %>% 
    # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
    # count(Unspliced_XBP1_umi_filtered, Spliced_XBP1_umi_filtered)
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    ggplot(aes(x = spliced_ratio)) +
    geom_histogram() +
    facet_grid(rows = vars(Donor), scales = "free") +
    geom_vline(aes(xintercept = mean_spliced_ratio, group = Donor), color = "green", linetype = 'longdash') 
  
  rna.obj.sub@meta.data %>% 
    filter(!Donor %in% c("UPN25", "UPN12", "UPN13")) %>%
    filter(Unspliced_XBP1 > 0) %>% 
    filter(Genotype == "MUT") %>% 
    group_by(Donor) %>% 
    filter(n() > 100) %>% 
    mutate(spliced_ratio = log2((Spliced_XBP1 + 0.1)/(Unspliced_XBP1 + 0.1))) %>%
    merge(., df.medians, by = "Donor") %>% 
    mutate(splicing_category = if_else(spliced_ratio > median, "High", "Low")) %>% 
    mutate(splicing_category = factor(splicing_category, levels = c("Low", "High"))) %>% 
    ggplot(aes(x = splicing_category, y = `Hallmark IFNa Response1`, fill = splicing_category)) +
    geom_boxplot() +
    facet_grid(cols = vars(Donor)) +
    stat_compare_means(label = "p.format")
  
  rna.obj.sub@meta.data %>% 
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
  
  rna.obj.sub@meta.data %>% 
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
  
  
  
  