library(tidyverse)


mapping <- list(
  PT6 = "GSM5858670",
  PT7 = "GSM5858671",
  PT8 = "GSM5858672",
  PT9 = "GSM5858673"
)

nih.calls <- lapply(names(mapping), function(x) {
  read_tsv(paste0("data/NIH_data/cb_sniffer_data_from_Shouguo/cb_sniffer_Cornell_UBA1/", x, "_CD34_barcode_counts_CB.tsv")) %>% 
    mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
    mutate(sample = mapping[[x]])
})
names(nih.calls) <- mapping

landau.calls <- lapply(mapping, function(x) {
  print(x)
  read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/", x, "_counts_CB.tsv")) %>% 
    mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
    mutate(sample = mapping[[x]])
})

names(landau.calls) <- mapping

## GSM5858670
nih.calls$GSM5858670 %>% 
  filter(start == "47199051") %>% 
  count(start, ref_count, alt_count)
landau.calls$GSM5858670 %>% 
  filter(start == "47199051") %>% 
  count(start, ref_count, alt_count)

## GSM5858671
nih.calls$GSM5858671 %>% 
  filter(start == "47199051") %>% 
  count(start, ref_count, alt_count)
landau.calls$GSM5858671 %>% 
  filter(start == "47199051") %>%
  count(start, ref_count, alt_count)

## GSM5858672
nih.calls$GSM5858672 %>% 
  filter(start == "47199052") %>% 
  count(start, ref_count, alt_count)
landau.calls$GSM5858672 %>% 
  filter(start == "47199052") %>%
  count(start, ref_count, alt_count)


## GSM5858673
nih.calls$GSM5858673 %>% 
  filter(start == "47199052") %>%
  count(start, ref_count, alt_count)
landau.calls$GSM5858673 %>% 
  filter(start == "47199052") %>%
  count(start, ref_count, alt_count)

## Summarize the counts per sample
nih.calls <- lapply(names(mapping), function(x) {
  read_tsv(paste0("data/NIH_data/cb_sniffer_data_from_Shouguo/cb_sniffer_Cornell_UBA1/", x, "_CD34_barcode_counts_CB.tsv")) %>% 
    mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
    select(-c(chrm, end, type, gene)) %>% 
    mutate(sample = mapping[[x]]) %>% 
    mutate(source = "NIH")
})
names(nih.calls) <- mapping

landau.calls <- lapply(mapping, function(x) {
  read_tsv(paste0("/gpfs/commons/home/rmurray/VEXAS/cb_sniffer/", x, "_counts_CB.tsv")) %>% 
    mutate(genotype_status_cbsniff = if_else(alt_count > 0, "MUT", if_else(ref_count > 0, "WT", NA_character_))) %>% 
    select(-c(chrm, end, type, gene)) %>% 
    mutate(sample = x) %>% 
    # dplyr::rename(ref_count_landau = ref_count) %>% 
    # dplyr::rename(alt_count_landau = alt_count) %>% 
    mutate(source = "Landau")
})

names(landau.calls) <- mapping


nih.calls.total <- do.call(rbind, nih.calls)
landau.calls.total <- do.call(rbind, landau.calls)

## Remove the other base calls for the non-mutation sites for NIH
nih.calls.total <- nih.calls.total %>% 
  filter((sample == "GSM5858670" & start == "47199051") |
      (sample == "GSM5858671" & start == "47199051") |
      (sample == "GSM5858672" & start == "47199052") |
      (sample == "GSM5858673" & start == "47199052"))

calls.merged <- merge(nih.calls.total, landau.calls.total, by = c("start", "ref", "alt", "barcode", "sample"), suffixes = c(".nih", ".landau"), all = T)


### Just look at the number of counts
nih.summary <- nih.calls.total %>% 
  count(sample, genotype_status_cbsniff) %>% 
  pivot_wider(names_from = genotype_status_cbsniff, values_from = n)

landau.summary <- landau.calls.total %>% 
  count(sample, genotype_status_cbsniff) %>% 
  pivot_wider(names_from = genotype_status_cbsniff, values_from = n)

counts.merged <- merge(nih.summary, landau.summary, by = "sample", suffixes = c("_NIH", "_Landau")) %>% 
  select(sample, MUT_NIH, MUT_Landau, WT_NIH, WT_Landau)
