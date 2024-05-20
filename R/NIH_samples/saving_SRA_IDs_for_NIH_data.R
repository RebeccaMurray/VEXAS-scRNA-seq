library(tidyverse)

sra.data <- read_csv("data/NIH_data/SraRunTable.txt")


###################### By GSM - VEXAS samples ################################

path.to.fastqs <- "/gpfs/commons/groups/landau_lab/VEXAS/fastqs/NIH_data/"

## Keep only the VEXAS samples
sra.data <- sra.data %>% filter(grepl("UPN", `subject_id/status`))

## Pick GSM numbers to process
gsm.numbers <- sra.data %>% 
  filter(source_name %in% c("FACS sorted Linjeage-CD34+ bone marrow", "Human Bone Marrow")) %>% 
  filter(AvgSpotLen != 310) %>% ## Remove the TCR/BCR libraries
  filter(grepl("UPN", `subject_id/status`)) %>% ## for VEXAS samples
  pull(`GEO_Accession (exp)`) %>% unique()

## Save the GSM numbers to a file
data.frame(geo_numbers = gsm.numbers) %>% 
  write_csv(., "/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/GEO_accession_list_CD34_and_BMMC.csv", col_names = F)

## Create a file with the fastqs corresponding to each GSM number
lapply(gsm.numbers, function(x) {
  ## Get the associated SRA numbers
  sra.numbers <- sra.data %>%
    filter(`GEO_Accession (exp)` == x) %>% 
    select(Run)
  
  ## Get the associated UPN
  upn <- sra.data %>%
    filter(`GEO_Accession (exp)` == x) %>% 
    pull(`subject_id/status`) %>% unique()
  
  ## Save a list of the SRAs
  sra.numbers %>% 
    write_csv(., paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/", x, "_SRA_numbers.csv"), 
              col_names = F)
  
  ## Save a list of the fastqs corresponding to each SRA
  sra.numbers %>% 
    mutate(Run = paste0(path.to.fastqs, "/", upn, "/", Run, "/")) %>% 
    write_csv(., paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/", x, "_SRA_numbers_fastq_paths.csv"), 
              col_names = F)
  
})

###################### By GSM - healthy donors ################################

path.to.fastqs <- "/gpfs/commons/groups/landau_lab/VEXAS/fastqs/NIH_data/CONTROLS/"

## Keep only the healthy donors
sra.data <- sra.data %>% filter(grepl("Healthydonor", `subject_id/status`))

## Pick GSM numbers to process
gsm.numbers <- sra.data %>% 
  filter(source_name %in% c("FACS sorted Linjeage-CD34+ bone marrow", "Human Bone Marrow")) %>% 
  filter(AvgSpotLen != 310) %>% ## Remove the TCR/BCR libraries
  pull(`GEO_Accession (exp)`) %>% unique()

## Save the GSM numbers to a file
data.frame(geo_numbers = gsm.numbers) %>% 
  write_csv(., "/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/GEO_accession_list_healthy_donor_CD34_and_BMMC.csv", col_names = F)

## Create a file with the fastqs corresponding to each GSM number
lapply(gsm.numbers, function(x) {
  ## Get the associated SRA numbers
  sra.numbers <- sra.data %>%
    filter(`GEO_Accession (exp)` == x) %>% 
    select(Run)
  
  ## Get the associated UPN
  upn <- sra.data %>%
    filter(`GEO_Accession (exp)` == x) %>% 
    pull(`subject_id/status`) %>% unique() %>% gsub("Healthydonor", "HD", .)
  
  ## Save a list of the SRAs
  sra.numbers %>% 
    write_csv(., paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/", x, "_SRA_numbers.csv"), 
              col_names = F)
  
  ## Save a list of the fastqs corresponding to each SRA
  sra.numbers %>% 
    mutate(Run = paste0(path.to.fastqs, "/", upn, "/", Run, "/fastqs/")) %>% 
    write_csv(., paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/", x, "_SRA_numbers_fastq_paths.csv"), 
              col_names = F)
  
})



###################### By patient ID ################################

## By patient ID
# lapply(unique(sra.data$`subject_id/status`), function(x) {
#   df <- sra.data %>%
#     filter(`subject_id/status` == x) %>% 
#     filter(source_name == "FACS sorted Linjeage-CD34+ bone marrow") %>% 
#     select(Run)
#   df %>% 
#     write_csv(., 
#             paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/NIH_", x, "_CD34_geneexpression.csv"), 
#             col_names = F)
#   df %>% 
#     mutate(Run = paste0(path.to.fastqs, x, "/", Run, "/")) %>% 
#     write_csv(., 
#             paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/NIH_SRA_files/NIH_", x, "_CD34_geneexpression_fastq_paths.csv"), 
#             col_names = F)
#   
# })
