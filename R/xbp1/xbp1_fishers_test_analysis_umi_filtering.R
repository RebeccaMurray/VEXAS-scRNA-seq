library(rstatix)
library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(ggpubr)
theme_set(theme_classic())


celltype.list <- c("HSC", "EMP", "LMPP", "MkP", "GMP")

donor.list <- levels(rna.obj$Donor)
names(donor.list) <- donor.list


rna.obj@meta.data %>% filter(CellType %in% celltype.list) %>% count(Donor, Genotype) %>% filter(Genotype != "NA")
rna.obj@meta.data %>% filter(CellType %in% celltype.list) %>% 
  group_by(Donor, CellType, Genotype, ) %>% summarize(cell_count = n(), unspliced_total = sum(Unspliced_XBP1_umi_filtered, na.rm = T), spliced_total = sum(Spliced_XBP1_umi_filtered, na.rm = T)) %>% 
  filter(Genotype != "NA") %>% 
  filter(cell_count > 10) %>% 
  group_by(Donor, CellType) %>% 
  filter(n() > 1) %>% 
  write_csv("data/xbp1_splicing/xbp1_read_counts_umi_filtered.csv")

run_fishers_test <- function(donor.name, celltype = NA, celltype.list = list()) {
  print(donor.name)
  
  md <- rna.obj@meta.data %>% filter(Donor == donor.name) %>% filter(Genotype %in% c("MUT", "WT"))
  
  if (length(celltype.list) != 0) {
    print("Multiple celltypes")
    md <- md %>% filter(cluster_celltype %in% celltype.list)
  } else {
    md <- md %>% filter(cluster_celltype == celltype)
  }
  
  ## If we have fewer than 10 genotyped cells for that sample for that cell type, return NULL
  if(nrow(md %>% filter(Genotype == "MUT")) < 10 | nrow(md %>% filter(Genotype == "WT")) < 10) {
    return(NULL)
  }
  
  test.table <- md %>% 
    group_by(Genotype) %>% 
    summarize(
      Unspliced_XBP1_total = sum(Unspliced_XBP1_umi_filtered, na.rm = T),
      Spliced_XBP1_total = sum(Spliced_XBP1_umi_filtered, na.rm = T)
    ) %>% 
    column_to_rownames("Genotype")
  print(test.table)
  test.results <-  test.table %>% 
      fisher.test(.)
  print("p_value:")
  print(test.results$p.value)
  print("odds ratio:")
  print(test.results$estimate)
  test.results.df <- data.frame(
    sample = donor.name, 
    p_value = test.results$p.value, 
    odds_ratio = test.results$estimate, 
    conf_interval_low = test.results$conf.int[[1]],
    conf_interval_high = test.results$conf.int[[2]],
    cluster = celltype
  )
  return(test.results.df)
}

## Run fisher's test for progenitors individually
names(celltype.list) <- celltype.list
fishers.test.results <- lapply(celltype.list, function(celltype) {
  print(paste0("Celltype: ", celltype))
  do.call(rbind, lapply(donor.list, function(donor.name) {run_fishers_test(donor.name = donor.name, celltype = celltype)}))
})
fishers.test.results
do.call(rbind, fishers.test.results) %>% 
  write_csv("data/xbp1_splicing/xbp1_fishers_test_results_umi_filtered.csv")


## Run fisher's test for progenitors combined
fishers.test.results.all.progenitors <- do.call(rbind, lapply(donor.list, function(donor.name) {run_fishers_test(donor.name = donor.name, celltype.list = c("HSC", "EMP", "LMPP", "MkP", "GMP"))}))
fishers.test.results.all.progenitors

## Plot results
fishers.test.results.plots <- lapply(names(fishers.test.results), function(celltype) {
  fishers.test.results[[celltype]] %>%
    mutate(label = if_else(p_value < 0.005, "**", if_else(p_value< 0.05, "*", ""))) %>% 
    ggplot(aes(x = reorder(sample, -odds_ratio), y = odds_ratio)) +
    geom_col() +
    geom_errorbar(aes(ymin = conf_interval_low, ymax = conf_interval_high), width = 0.25) +
    geom_text(aes(label = label), vjust = -0.5) +
    theme(axis.title.x = element_blank()) +
    ylab("Odds ratio, \n Fisher's exact test") + 
    scale_y_continuous(expand = c(0.1, 0)) +
    scale_x_discrete(drop = FALSE) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom") +
    ggtitle(celltype)
})
names(fishers.test.results.plots) <- names(fishers.test.results)

## Plot HSCs separately
fishers.test.results.plots$HSC + theme(legend.position = "right")
fishers.test.results.plots$MkP + theme(legend.position = "right")

## Plot other celltypes
wrap_plots(fishers.test.results.plots, nrow = 2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

## Plot combined progenitors analysis
p.fishers.test <- fishers.test.results.all.progenitors %>% 
  # mutate(`MDS status` = if_else(sample %in% mds.sample.list, "MDS", "None")) %>% 
  # mutate(`MDS status` = factor(`MDS status`, levels = c("None", "MDS"))) %>% 
  # mutate(label = if_else(p_value < 0.005, "**", if_else(p_value< 0.05, "*", ""))) %>% 
  mutate(`pval < 0.05` = p_value < 0.05) %>% 
  ggplot(aes(x = reorder(sample, -odds_ratio), y = odds_ratio, fill = `pval < 0.05`)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = conf_interval_low, ymax = conf_interval_high), width = 0.25) +
  # geom_text(aes(label = label), vjust = -0.5) +
  theme(axis.title.x = element_blank()) +
  ylab("Odds ratio, \n Fisher's exact test") + 
  scale_y_continuous(expand = c(0.1, 0)) +
  coord_cartesian(ylim = c(0, 18)) +
  scale_x_discrete(drop = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE` = "darkseagreen", `FALSE` = "grey80")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "right", legend.justification = "top") +
  ggtitle("XBP1s ratio vs. Genotype\nAll HSPCs")
p.fishers.test
ggsave("figures/current_figure_drafts/XBP1_fishers_bar_plot.pdf", plot = p.fishers.test, device = "pdf", dpi = 300, width = 4, height = 3, unit = "in")


## Make a combined bar plot with values for all patients
fishers.test.results.combined <- do.call(rbind, fishers.test.results)
fishers.test.results.combined %>% 
  mutate(cluster = factor(cluster, levels = celltype.list)) %>% 
  ggplot(aes(x = cluster, y = odds_ratio, fill = cluster)) +
  geom_violin() +
  geom_point(position = "jitter") +
  scale_fill_manual(values = cell.type.palette)

