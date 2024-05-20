############################################## XBP1s - sample-aware permutation test (unfiltered, HSC) ###################################

allowed.cell.types <- c("HSC", "EMP", "LMPP", "MkP", "GMP")
allowed.cell.types <- c("HSC")

donors.keep <- rna.obj@meta.data %>% 
  filter(cluster_celltype %in% allowed.cell.types) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n > 50) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor) %>% as.character()

md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% allowed.cell.types) %>% 
  filter(Genotype %in% c("MUT", "WT"))

md.xbp1 %>%
  group_by(Donor,Genotype) %>% 
  summarize(n = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum)

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the mean between the two ratios
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
md.xbp1.shuf <- md.xbp1
n = 1000
scrambled_diff_ratios <- sapply(1:n, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor, Genotype) %>% 
    slice_sample(n = 25) %>%
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_ratios(.)
})
actual_diff_ratios <- sapply(1:n, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor, Genotype) %>% 
    slice_sample(n = 25, replace = F) %>% 
    group_by(Donor) %>% 
    calculate_difference_of_ratios(.)
})

scrambled.df <- data.frame(values = scrambled_diff_ratios)
actual.df <- data.frame(values = actual_diff_ratios)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "dodgerblue") +
  geom_histogram(data = actual.df, alpha = 0.2, fill = "salmon") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("Sample-aware permutation test, n = ", n)) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")
# ggsave("figures/current_figure_drafts/sample_aware_permutation_test_10k_iterations_kallisto_output.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

## No downsampling - actual values
calculate_difference_of_ratios(md.xbp1)

## Downsample one time - actual values
md.xbp1 %>% 
  group_by(Donor, Genotype) %>% 
  slice_sample(n = 25, replace = T) %>% 
  group_by(Donor) %>% 
  calculate_difference_of_ratios(.)

p.value <- sum(scrambled_diff_ratios > mean(actual_diff_ratios)) / length(scrambled_diff_ratios)
print(p.value)



f########################################## Bar plots ##########################################


## Pseudobulk by donor
p.xbp1.bars <- md %>% 
  filter(Donor %in% donors.keep) %>% 
  # filter(cluster_celltype == "HSC") %>%
  filter(cluster_celltype %in% allowed.cell.types) %>%
  filter(Genotype %in% c("MUT", "WT")) %>% 
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
  ggtitle("XBP1s/XBP1u", subtitle = "HSPCs") +
  # annotate("text", x = 1.2, y = 0.4, label=paste0("p = ", p.value), size = 3)
  annotate("text", x = 1.2, y = 0.4, label="p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_column_plot_kallisto_output.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)




###################################### Assessment together ######################

## Shuffle the genotype labels and calculate within patients
md.xbp1.shuf <- md.xbp1
n = 1000
diff_ratios <- sapply(1:n, function(iter) {
  md.downsample <- md.xbp1.shuf %>% 
    group_by(Donor, Genotype) %>% 
    slice_sample(n = 25)
  scrambled.diff <- md.downsample %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_ratios(.)
  real.diff <- md.downsample  %>% 
    calculate_difference_of_ratios(.)
  # print(scrambled.diff)
  # print(real.diff) 
  return(scrambled.diff > real.diff)
})

actual_diff_ratios <- sapply(1:n, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor, Genotype) %>% 
    slice_sample(n = 25) %>% 
    group_by(Donor) %>% 
    calculate_difference_of_ratios(.)
})


################ Downsampling once ###############################

donors.keep <- rna.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC")) %>% 
  count(Donor, Genotype) %>% 
  filter(Genotype != "NA") %>% 
  filter(n > 25) %>% 
  count(Donor) %>% 
  filter(n > 1) %>% 
  pull(Donor) %>% as.character(.)

md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% c("HSC")) %>% 
  # filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT"))

md.xbp1 %>%
  filter(Donor %in% donors.keep) %>% 
  group_by(Donor,Genotype) %>% 
  summarize(n = n(), spliced_sum = sum(Spliced_XBP1, na.rm = T), unspliced_sum = sum(Unspliced_XBP1, na.rm = T)) %>% 
  mutate(splicing_ratio = spliced_sum/unspliced_sum)

selection.values <- md.xbp1 %>%
  filter(Donor %in% donors.keep) %>% 
  group_by(Donor,Genotype) %>% 
  summarize(n = n()) %>% group_by(Donor) %>% 
  summarize(val_to_select = min(n))
selection.list <- selection.values$val_to_select
names(selection.list) <- selection.values$Donor

md.xbp1.downsampled.list <- lapply(donors.keep, function(donor.current) {
  print(donor.current)
  print(selection.list[[donor.current]])
  md.xbp1 %>% 
    filter(Donor == donor.current) %>% 
    group_by(Genotype) %>% 
    slice_sample(n = selection.list[[donor.current]])
})

md.xbp1.combined <- do.call(rbind, md.xbp1.downsampled.list)

## Helper function: pseudobulks by genotype, calculates the spliced/unspliced ratio within genotype,
## then takes the mean between the two ratios
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
md.xbp1.shuf <- md.xbp1.combined
n = 1000
scrambled_diff_ratios <- sapply(1:n, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor, Genotype) %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_ratios(.)
})

scrambled.df <- data.frame(values = scrambled_diff_ratios)

p.value <- sum(scrambled_diff_ratios > 0.01762931) / length(scrambled_diff_ratios)
print(p.value)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "dodgerblue") +
  geom_histogram(data = actual.df, alpha = 0.2, fill = "salmon") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("Sample-aware permutation test, n = ", n)) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Count")




## dynamically set downsampling number
md.xbp1.shuf <- md.xbp1
n = 100
scrambled_diff_ratios <- sapply(1:n, function(iter) {
  if (iter %% 10 == 0) {
    print(iter)
  }
  
  scrambled.dfs <- lapply(donors.keep, function(donor) {
    print("using donor:")
    print(donor)
    print(selection.list[[donor]])
    md.xbp1.shuf %>% 
      filter(Donor == donor) %>% 
      group_by(Genotype) %>% 
      slice_sample(n = selection.list[[donor]]) %>% 
      mutate(Genotype = sample(Genotype, replace = F))
  })
  # do.call(rbind, scrambled.dfs) %>% 
  #   group_by(Donor, Genotype) %>% 
  #   count() %>% 
  #   print()
  do.call(rbind, scrambled.dfs) %>% 
    calculate_difference_of_ratios(.)
})

actual_diff_ratios <- sapply(1:n, function(iter) {
  if (iter %% 10 == 0) {
    print(iter)
  }
  actual.dfs <- lapply(donors.keep, function(donor) {
    md.xbp1.shuf %>% 
      filter(Donor == donor) %>% 
      group_by(Genotype) %>% 
      mutate(count = n()) %>% 
      slice_sample(n = selection.list[[donor]])
  })
  # do.call(rbind, actual.dfs) %>% 
  #   group_by(Donor, Genotype) %>% 
  #   count() %>% 
  #   print()
  do.call(rbind, actual.dfs) %>% 
    calculate_difference_of_ratios(.)
})

ks.test(scrambled_diff_ratios, actual_diff_ratios)