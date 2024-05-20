
############################################# Additional XBP1s explorations ############################################################

## Single cell - together
md %>% 
  # filter(!(Donor %in% c("PT03", "PT04", "PT05", "PT08"))) %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(Unspliced_XBP1 > 0) %>% 
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) + geom_violin() +scale_fill_manual(values = vexas.genotyping.palette) +
  geom_jitter() +
  stat_compare_means(label = "p.format", method = "wilcox.test")


## Single cell - by patient
rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(Unspliced_XBP1 > 0) %>% 
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Genotype, y = log2(splicing_ratio), fill = Genotype)) + geom_violin() +scale_fill_manual(values = vexas.genotyping.palette) +
  geom_jitter() + facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format", method = "wilcox.test")

## Look at PT07
DimPlot(rna.obj, cells.highlight = c("PT07 HSCs" = md %>% filter(CellType == "HSC" & Donor == "PT07") %>% rownames()))
DimPlot(rna.obj, cells.highlight = c("PT02 HSCs" = md %>% filter(CellType == "HSC" & Donor == "PT02") %>% rownames()))

## Number of cells with XBP1s data?
md %>% 
  filter(CellType == "HSC") %>% 
  group_by(Donor) %>% 
  # summarize(has_XBP1_info_frac = sum(!is.na(Unspliced_XBP1))/n()) %>% 
  summarize(has_XBP1_info_frac = sum(!is.na(Unspliced_XBP1) & Unspliced_XBP1 > 0)/n()) %>% 
  arrange(has_XBP1_info_frac)
  
md %>% 
  filter(CellType == "HSC") %>% 
  mutate(Total_XBP1 = Unspliced_XBP1 + Unspliced_XBP1) %>%
  ggplot(aes(x = Spliced_XBP1)) +
  geom_histogram() +
  facet_wrap(~ Donor, scales = "free")


############################################# Recreating GoT paper ####################################################################


rna.obj.subset <- read_rds(paste0("/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pseudotime/vexas_meyloid_pseudotime_cells_leidengraph.Rds")) %>% subset(., Object == "VEXAS")
monocle3.data <- rna.obj.subset@meta.data %>% select(monocle3_pseudotime) %>% transmute(Pseudotime = monocle3_pseudotime)
md.pseudo.all <- merge(md, monocle3.data, by = "row.names", all.x = T) %>% column_to_rownames("Row.names")
md.pseudo.has.values <- md.pseudo.all %>% filter(!is.na(Pseudotime))

p1 <- md.pseudo.has.values %>% 
  # filter(Donor == "PT02") %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP", "CLP", "GMP")) %>% 
  filter(Unspliced_XBP1 > 1) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  mutate(Total_XBP1 = Spliced_XBP1 + Unspliced_XBP1) %>% 
  ggplot(., aes(x = Pseudotime, group = Genotype, y = Total_XBP1, color = Genotype)) +
  geom_smooth() +
  geom_smooth(aes(x = Pseudotime, y = Spliced_XBP1, group = Genotype, color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette)

p2 <- md.pseudo.has.values %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP", "CLP", "GMP")) %>% 
  filter(Unspliced_XBP1 > 1) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Pseudotime, group = CellType, color = CellType, fill = CellType)) +
  geom_density(alpha = 0.5) +
  scale_y_reverse()

p1 / p2 + plot_layout(heights = c(2, 1)) + plot_annotation("Myeloid lineage")

############################################## UMI-filtered values? ###################################################################
rna.obj@meta.data %>% head()

xbp1.celltypes.list <- c("HSC")
rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% 
  filter(!(Donor %in% c("PT03", "PT04", "PT08"))) %>% 
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), unspliced_sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme() +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = 7 donors")) 
  # annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)


############################################## XBP1s - sample-aware permutation test (kallisto, HSC) ###################################


xbp1.celltypes.list <- c("HSC")
min.genotyped.cell.count <- 20
# xbp1.celltypes.list <- c("HSC", "EMP", "LMPP", "MkP")
# min.genotyped.cell.count <- 50

rna.obj@meta.data %>% 
  # filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  pivot_wider(names_from = Genotype, values_from = n)


donors.keep <- rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n > min.genotyped.cell.count) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor) %>% as.character()
print(donors.keep)

## Save a table of the spliced/unspliced XBP1 UMI counts per donor, cell type, and genotype
rna.obj@meta.data %>% filter(CellType %in% xbp1.celltypes.list) %>% 
  filter(Donor %in% donors.keep) %>% 
  group_by(CellType, Donor, Genotype) %>% summarize(cell_count = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  filter(Genotype != "NA") %>% 
  write_csv("data/xbp1_splicing/xbp1_read_counts_kallisto_output_HSC_only_20240502.csv")


md.xbp1 <- rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the difference between the two ratios
calculate_difference_of_ratios <- function(df) {
  df %>% 
    group_by(Genotype) %>% 
    summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
    mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
    select(Genotype, splicing_ratio) %>% 
    pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
    mutate(ratio_diff = MUT - WT) %>% 
    pull(ratio_diff)
}

## Shuffle the genotype labels and calculate within patients
n = 1000
scrambled_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>% ## Randomly downsample to 25 MUT and 25 WT cells per donor
    slice_sample(n = min.genotyped.cell.count, replace = T) %>%
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% ## Shuffle the genotype labels within donors 
    calculate_difference_of_ratios(.)
})

actual_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>% ## Randomly downsample to 25 MUT and 25 WT cells per donor
    slice_sample(n = min.genotyped.cell.count, replace = T) %>%
    calculate_difference_of_ratios(.)
})

## Convert ratios into data frames and plot
scrambled.df <- data.frame(values = scrambled_diff_results)
actual.df <- data.frame(values = actual_diff_results)

## Calculate the actual mean difference
observed.ratio.diff <- calculate_difference_of_ratios(md.xbp1)
observed.ratio.diff.downsampled <- md.xbp1 %>% group_by(Donor, Genotype) %>% slice_sample(n = min.genotyped.cell.count, replace = T) %>% calculate_difference_of_ratios(.)
observed.ratio.diff.downsampled.mean <- mean(actual_diff_results)

print(observed.ratio.diff)
print(observed.ratio.diff.downsampled)
print(observed.ratio.diff.downsampled.mean)

## Calculate the fraction of simulated results that are above the observed value
p.value <- sum(scrambled_diff_results > observed.ratio.diff.downsampled.mean) / length(scrambled_diff_results)
print(p.value)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "dodgerblue") +
  geom_histogram(data = actual.df, alpha = 0.2, fill = "salmon") +
  geom_vline(xintercept = observed.ratio.diff.downsampled.mean, linetype = "longdash", color = "darkred") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("HSCs, sample-aware permutation test, n = ", n)) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_1k_iterations_downsampling_kallisto_output_HSC_only_20240502.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

## Plot bar plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  summarize(mean_splicing_ratio = mean(splicing_ratio),
            sd = sd(splicing_ratio),
            n = n(),
            se = sd / sqrt(n)) %>% 
  ggplot(aes(x = Genotype, y = mean_splicing_ratio, fill = Genotype)) +
  geom_col(color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_splicing_ratio-se, ymax = mean_splicing_ratio+se), width=0.4, size = 0.25) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

## Plot box plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_box_plot_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

## Plot line plot with splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme() +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_line_plot_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)


## Paired t-test
t.test.df <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% select(Donor, Genotype, splicing_ratio) %>% pivot_wider(names_from = Genotype, values_from = splicing_ratio)
t.test(t.test.df$WT, t.test.df$MUT, paired = TRUE, alternative = "two.sided")

## Normalized number of spliced transcripts?
md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(normalized_spliced_sum = spliced_sum / cell_count) %>% 
  arrange(Donor, Genotype) %>% 
  ggplot(aes(x = Genotype, y = normalized_spliced_sum, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  theme() +
  ggtitle("XBP1s", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors"))

md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(normalized_spliced_sum = spliced_sum / cell_count) %>% 
  mutate(normalized_unspliced_sum = unspliced_sum / cell_count) %>% 
  arrange(Donor, Genotype) %>% 
  ggplot(aes(x = Genotype, y = normalized_unspliced_sum, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  theme() +
  ggtitle("XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors"))

t.test.df <- md.xbp1 %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(normalized_spliced_sum = spliced_sum / cell_count) %>% 
  select(Donor, Genotype, normalized_spliced_sum) %>% 
  pivot_wider(names_from = Genotype, values_from = normalized_spliced_sum)
t.test(t.test.df$WT, t.test.df$MUT, paired = TRUE, alternative = "two.sided")

############# other plots

## Single cell - together
rna.obj@meta.data %>% 
  filter(!(Donor %in% c("PT03", "PT04", "PT05", "PT08"))) %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Genotype, y = log2(splicing_ratio), fill = Genotype)) + geom_violin() +scale_fill_manual(values = vexas.genotyping.palette) +
  geom_jitter() +
  stat_compare_means(label = "p.format", method = "wilcox.test")


## Single cell - by patient
rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Genotype, y = log2(splicing_ratio), fill = Genotype)) + geom_violin() +scale_fill_manual(values = vexas.genotyping.palette) +
  geom_jitter() + facet_grid(cols = vars(Donor)) +
  stat_compare_means(label = "p.format", method = "wilcox.test")

rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype == "HSC") %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(Donor, Genotype) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), count = n()) %>% 
  mutate(normalized_spliced_sum = spliced_sum / count, normalized_unspliced_sum = unspliced_sum / count ) %>% 
  mutate(normalized_spliced_ratio = normalized_spliced_sum / normalized_unspliced_sum)
  ggplot(aes(x = Genotype, y = normalized_spliced_ratio, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme() +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
  
  
  
############### IRE1 module score

  ## Extract the raw counts matrix, remove mito/ribo genes
  non.mito.ribo.genes <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., value = T, invert = T)
  print("Excluding genes:")
  print(setdiff(rna.obj@assays$RNA@counts %>% rownames(), non.mito.ribo.genes))
  m <- rna.obj@assays$RNA@counts[non.mito.ribo.genes, ]
  
  ## Replace the RNA assay
  DefaultAssay(rna.obj) <- "ADT_DSB" ## Set another random assay as default
  rna.obj[["RNA"]] <- NULL
  rna.obj[["RNA"]] <- CreateAssayObject(counts = m)
  DefaultAssay(rna.obj) <- "RNA"
  rna.obj <- NormalizeData(rna.obj, assay = "RNA", normalization.method = "LogNormalize")
  
  
pathways = readRDS('/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pathways.Rds')

# `Reactome translation` = pathways$REACTOME_TRANSLATION,
# `GOBP ATF6` = pathways$GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
# `GOBP IRE1` = pathways$GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE,
# # `HAllmark UPR` = pathways$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
# `ATF4 target genes` = pathways$ATF4_ONLY_HAN,
# `CHOP target genes` = pathways$CHOP_ONLY_HAN,
# # # `XBP1 target module` = pathways$`XBP1 target m
# 

reactome_xbp1_module <- c("ACADVL", "ADD1", "ARFGAP1", "ATP6V0D1", "CTDSP2", "CUL7", "CXXC1", "DCTN1", "DDX11", "DDX12P", "DNAJB11", "DNAJB9", "DNAJC3", "EDEM1", "EXTL3", "FKBP14", "GOSR2", "GSK3A", "HDGF", "HYOU1", "KDELR3", "KLHDC3", "LMNA", "MYDGF", "PDIA5", "PDIA6", "PPP2R5B", "PREB", "SEC31A", "SERP1", "SHC1", "SRPRA", "SRPRB", "SSR1", "SULT1A3", "SULT1A4", "SYVN1", "TATDN2", "TLN1", "TPP1", "TSPYL2", "WFS1", "WIPI1", "XBP1", "YIF1A", "ZBTB17")
module.name <- "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE"

rna.obj <- AddModuleScore(rna.obj, assay = "RNA", features = list(pathways[[module.name]]), name = "current_module", search = T)
# rna.obj <- AddModuleScore(rna.obj, assay = "RNA", features = list(reactome_xbp1_module), name = "current_module", search = T)
md <- rna.obj@meta.data
md %>% 
  filter(!(Donor %in% c("PT03", "PT04", "PT08"))) %>% 
  filter(CellType == "HSC" & Genotype %in% c("MUT",  "WT")) %>% 
  ggplot(aes(x = Genotype, y = current_module1, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = vexas.genotyping.palette) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  ggtitle(module.name) + 
  theme(legend.position = "none")

md %>% 
  filter(!(Donor %in% c("PT03", "PT04", "PT05", "PT08"))) %>% 
  filter(CellType == "HSC" & Genotype %in% c("MUT",  "WT")) %>% 
  ggplot(aes(x = Donor, y = current_module1, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = vexas.genotyping.palette) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  ggtitle(module.name) + 
  theme(legend.position = "none")

## Correlate with inflammation?
md %>% 
  filter(!(Donor %in% c("PT03", "PT04", "PT08"))) %>% 
  filter(CellType == "HSC" & Genotype == "MUT") %>% 
  filter(Unspliced_XBP1 > 0) %>% 
  ggplot(aes(x = log2(splicing_ratio), y = current_module1)) +
  geom_point()
  

############################################## XBP1s - sample-aware permutation test (kallisto, HSC, separating limits) ###################################


xbp1.celltypes.list <- c("HSC")

rna.obj@meta.data %>% 
  # filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  pivot_wider(names_from = Genotype, values_from = n)


md.xbp1 <- rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% 
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the difference between the two ratios
calculate_difference_of_ratios <- function(df) {
  df %>% 
    group_by(Genotype) %>% 
    summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
    mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
    select(Genotype, splicing_ratio) %>% 
    pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
    mutate(ratio_diff = MUT - WT) %>% 
    pull(ratio_diff)
}

## Shuffle the genotype labels and calculate within patients
n = 1000
scrambled_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>% ## Randomly downsample to 25 MUT and 25 WT cells per donor
    slice_sample(n = 50, replace = T) %>%
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% ## Shuffle the genotype labels within donors 
    calculate_difference_of_ratios(.)
})

actual_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>% ## Randomly downsample to 25 MUT and 25 WT cells per donor
    slice_sample(n = min.genotyped.cell.count, replace = T) %>%
    calculate_difference_of_ratios(.)
})

## Convert ratios into data frames and plot
scrambled.df <- data.frame(values = scrambled_diff_results)
actual.df <- data.frame(values = actual_diff_results)

## Calculate the actual mean difference
observed.ratio.diff <- calculate_difference_of_ratios(md.xbp1)
observed.ratio.diff.downsampled <- md.xbp1 %>% group_by(Donor, Genotype) %>% slice_sample(n = min.genotyped.cell.count, replace = T) %>% calculate_difference_of_ratios(.)
observed.ratio.diff.downsampled.mean <- mean(actual_diff_results)

print(observed.ratio.diff)
print(observed.ratio.diff.downsampled)
print(observed.ratio.diff.downsampled.mean)

## Calculate the fraction of simulated results that are above the observed value
p.value <- sum(scrambled_diff_results > observed.ratio.diff.downsampled.mean) / length(scrambled_diff_results)
print(p.value)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "dodgerblue") +
  geom_histogram(data = actual.df, alpha = 0.2, fill = "salmon") +
  geom_vline(xintercept = observed.ratio.diff.downsampled.mean, linetype = "longdash", color = "darkred") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("HSCs, sample-aware permutation test, n = ", n)) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_1k_iterations_downsampling_kallisto_output_HSC_only_20240502.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

## Plot bar plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  summarize(mean_splicing_ratio = mean(splicing_ratio),
            sd = sd(splicing_ratio),
            n = n(),
            se = sd / sqrt(n)) %>% 
  ggplot(aes(x = Genotype, y = mean_splicing_ratio, fill = Genotype)) +
  geom_col(color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_splicing_ratio-se, ymax = mean_splicing_ratio+se), width=0.4, size = 0.25) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

## Plot box plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_box_plot_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

## Plot line plot with splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_point(aes(color = Donor)) +
  geom_line(aes(color = Donor, group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme() +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_HSC_only_line_plot_20240502.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)



## Plot bar plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  summarize(mean_splicing_ratio = mean(splicing_ratio),
            sd = sd(splicing_ratio),
            n = n(),
            se = sd / sqrt(n)) %>% 
  ggplot(aes(x = Genotype, y = mean_splicing_ratio, fill = Genotype)) +
  geom_col(color = "black", size = 0.25) +
  geom_errorbar(aes(ymin = mean_splicing_ratio-se, ymax = mean_splicing_ratio+se), width=0.4, size = 0.25) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSPCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.5, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_20240506.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)


## Plot box plot with distribution of splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_box_plot_20240506.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

## Plot line plot with splicing ratios for donors used in the analysis
p.xbp1.bars <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_point() +
  geom_line(aes(group = Donor)) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSPCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
# annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_bar_plot_kallisto_output_line_plot_20240506.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)


