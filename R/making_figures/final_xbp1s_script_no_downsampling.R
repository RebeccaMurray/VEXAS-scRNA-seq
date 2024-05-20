
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


############################################## XBP1s - sample-aware permutation test  ###################################

xbp1.celltypes.list <- c("HSC")

md.xbp1 <- rna.obj@meta.data %>% 
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
n = 10000
scrambled_diff_results <- sapply(1:n, function(iter) {
  md.xbp1 %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% ## Shuffle the genotype labels within donors 
    calculate_difference_of_ratios(.)
})

## Convert the output lists into data frames
scrambled.df <- data.frame(values = scrambled_diff_results)

## Calculate the actual mean difference
observed.ratio.diff <- calculate_difference_of_ratios(md.xbp1)
print(observed.ratio.diff)

## Calculate the fraction of scrambled genotype label results that are above the mean of the real genotype label results
p.value <- sum(scrambled.df$values > observed.ratio.diff) / length(scrambled.df$values)
print(p.value)

## Plot histogram of all iterations
scrambled.df %>% 
  ggplot(aes(x = values, y = ..count..)) +
  geom_histogram(alpha = 0.2, fill = "grey20") +
  geom_vline(xintercept = observed.ratio.diff, linetype = "longdash", color = "darkred") +
  ggtitle("Pseudobulked splicing ratios, MUT - WT", subtitle = paste0("HSCs, sample-aware permutation test, n = ", n)) +
  annotate("text", x = observed.ratio.diff, y = 200, 
           label=paste0("Mean = ", round(p.value, digits = 3)), size = 3, color = "darkred", hjust = -0.5) +
  xlab("XBP1s/XPB1u[MUT] - XBP1s/XPB1u[WT]") +
  ylab("Frequency")
ggsave("figures/current_figure_drafts/sample_aware_permutation_test_1k_iterations_downsampling_HSPCs_kallisto_output_20240506.pdf", device = "pdf", width = 4.5, height = 3, dpi = 300)

