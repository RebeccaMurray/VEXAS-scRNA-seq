
set.seed(1234)

rna.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  count(Unspliced_XBP1)

rna.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Unspliced_XBP1 == 0 & Spliced_XBP1 == 0) %>% 
  count(Donor, Spliced_XBP1)

rna.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Unspliced_XBP1 == 0) %>% 
  count(Donor, Spliced_XBP1)


## Check how many cells we have per donor - HSPCs
rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  pivot_wider(names_from = Genotype, values_from = n)

## Check how many cells we have per donor - HSCs
rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% c("HSC")) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  pivot_wider(names_from = Genotype, values_from = n)


############################################## XBP1s - sample-aware permutation test (kallisto, HSPC) ###################################

xbp1.celltypes.list <- c("HSC", "EMP", "LMPP", "MkP")
min.genotyped.cell.count <- 25

donors.keep <- rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n >= min.genotyped.cell.count) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor) %>% as.character()
print(donors.keep)

## Print the cell type composition per donor
rna.obj@meta.data %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  filter(Donor %in% donors.keep) %>% 
  count(Donor, CellType)  %>% 
  group_by(Donor) %>% 
  mutate(frac = n / sum(n)) %>% 
  ggplot(aes(x = Donor, y = frac, fill = CellType)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  ylab("Cell Type Composition (%)")
ggsave("figures/current_figure_drafts/distribution_of_progenitor_subtypes_XBP1s.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)


## Save a table of the spliced/unspliced XBP1 UMI counts per donor, cell type, and genotype
rna.obj@meta.data %>% filter(CellType %in% xbp1.celltypes.list) %>% 
  filter(Donor %in% donors.keep) %>% 
  group_by(CellType, Donor, Genotype) %>% summarize(cell_count = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  filter(Genotype != "NA") %>% 
  write_csv("data/xbp1_splicing/xbp1_read_counts_HSPC_kallisto_output_20240506.csv")


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
num_cells_to_sample = 50
scrambled_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>%
    slice_sample(n = num_cells_to_sample, replace = T) %>% ## Randomly select X cells, with replacement
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% ## Shuffle the genotype labels within donors 
    calculate_difference_of_ratios(.)
})

actual_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>%
    slice_sample(n = num_cells_to_sample, replace = T) %>% ## Randomly select X cells, with replacement
    calculate_difference_of_ratios(.) ## Calculate difference WITHOUT shuffling genotype labels
})

## Convert the output lists into data frames
scrambled.df <- data.frame(values = scrambled_diff_results)
actual.df <- data.frame(values = actual_diff_results)

## Save the output data frames
saveRDS(scrambled.df, file = "data/xbp1_splicing_scrambled_df_HSPC.rds")
saveRDS(actual.df, file = "data/xbp1_splicing_actual_df_HSPC.rds")
# scrambled.df <- readRDS("data/xbp1_splicing_scrambled_df_HSPC.rds")
# actual.df <- readRDS("data/xbp1_splicing_actual_df_HSPC.rds")

## Calculate the actual mean difference
observed.ratio.diff <- calculate_difference_of_ratios(md.xbp1)
observed.ratio.diff.downsampled.mean <- mean(actual.df$values)

print(observed.ratio.diff)
print(observed.ratio.diff.downsampled.mean)

## Calculate the fraction of scrambled genotype label results that are above the mean of the real genotype label results
p.value <- sum(scrambled.df$values > observed.ratio.diff.downsampled.mean) / length(scrambled.df$values)
print(p.value)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "dodgerblue") +
  geom_histogram(data = actual.df, alpha = 0.2, fill = "salmon") +
  geom_vline(xintercept = observed.ratio.diff.downsampled.mean, linetype = "longdash", color = "darkred") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("HSPCs, sample-aware permutation test, n = ", n)) +
  annotate("text", x = observed.ratio.diff.downsampled.mean, y = 1500, 
           label=paste0("Mean = ", round(observed.ratio.diff.downsampled.mean, digits = 3)), size = 3, color = "darkred", hjust = -0.5) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_1k_iterations_downsampling_HSPCs_kallisto_output_20240506.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

## Plot the boxplot with lines overlaid
p.xbp1.boxplot.lines <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  geom_point() +
  geom_line(aes(group = Donor)) +
  ylab("log2(XBP1s/XBP1u)") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSPCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
p.xbp1.boxplot.lines
ggsave("figures/current_figure_drafts/RNA_xbp1s_boxplots_dots_overlaid_HSPC_20240506.pdf", plot = p.xbp1.boxplot.lines, device = "pdf", dpi = 300, width = 1.5, height = 2.5)



############################################## XBP1s - sample-aware permutation test (kallisto, HSC) ###################################

xbp1.celltypes.list <- c("HSC")
min.genotyped.cell.count <- 25

donors.keep <- rna.obj@meta.data %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n >= min.genotyped.cell.count) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor) %>% as.character()
print(donors.keep)

## Print the number of cells per donor
rna.obj@meta.data %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(!is.na(Unspliced_XBP1)) %>% ## Remove cells without xbp1s info
  filter(cluster_celltype %in% xbp1.celltypes.list) %>% 
  filter(Donor %in% donors.keep) %>% 
  count(Donor)  %>% 
  group_by(Donor) %>% 
  ggplot(aes(x = Donor, y = n, fill = Donor)) +
  geom_col(position = "dodge")

## Save a table of the spliced/unspliced XBP1 UMI counts per donor, cell type, and genotype
rna.obj@meta.data %>% filter(CellType %in% xbp1.celltypes.list) %>% 
  filter(Donor %in% donors.keep) %>% 
  group_by(CellType, Donor, Genotype) %>% summarize(cell_count = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
  filter(Genotype != "NA") %>% 
  write_csv("data/xbp1_splicing/xbp1_read_counts_kallisto_output_HSC_20240506.csv")


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
num_cells_to_sample = 50
scrambled_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>%
    slice_sample(n = num_cells_to_sample, replace = T) %>% ## Randomly select X cells, with replacement
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% ## Shuffle the genotype labels within donors 
    calculate_difference_of_ratios(.)
})

actual_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor, Genotype) %>%
    slice_sample(n = num_cells_to_sample, replace = T) %>% ## Randomly select X cells, with replacement
    calculate_difference_of_ratios(.) ## Calculate difference WITHOUT shuffling genotype labels
})

## Convert output lists into dataframes
scrambled.df <- data.frame(values = scrambled_diff_results)
actual.df <- data.frame(values = actual_diff_results)

# ## Save the output data frames
# saveRDS(scrambled.df, file = "data/xbp1_splicing_scrambled_df_HSC.rds")
# saveRDS(actual.df, file = "data/xbp1_splicing_actual_df_HSC.rds")

## Calculate the actual mean difference
observed.ratio.diff <- calculate_difference_of_ratios(md.xbp1)
observed.ratio.diff.downsampled.mean <- mean(actual_diff_results)

print(observed.ratio.diff)
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
  annotate("text", x = observed.ratio.diff.downsampled.mean, y = 1500, 
           label=paste0("Mean = ", round(observed.ratio.diff.downsampled.mean, digits = 3)), size = 3, color = "darkred", hjust = -0.5) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_1k_iterations_downsampling_HSCs_kallisto_output_20240506.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

## Plot the boxplot with lines overlaid
p.xbp1.boxplot.lines <- md.xbp1 %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T), cell_count = n()) %>% 
  mutate(splicing_ratio = log2(spliced_sum/unspliced_sum)) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_boxplot() +
  scale_fill_manual(values = c(WT = "#FFCB88", MUT = "#FF7B0B")) +
  geom_point() +
  geom_line(aes(group = Donor)) +
  ylab("XBP1s/XBP1u") +
  theme(legend.position = "none") +
  ggtitle("XBP1s/XBP1u", subtitle = paste0("HSCs\nn = ", length(donors.keep), " donors")) +
  annotate("text", x = 1.2, y = 0.07, label=paste0("p = ", p.value), size = 3)
p.xbp1.boxplot.lines
ggsave("figures/current_figure_drafts/RNA_xbp1s_boxplots_dots_overlaid_HSC_20240506.pdf", plot = p.xbp1.boxplot.lines, device = "pdf", dpi = 300, width = 1.5, height = 2.5)
