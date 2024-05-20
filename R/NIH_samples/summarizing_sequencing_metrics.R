library(tidyverse)

## Load the NIH samples
dir.list <- list.dirs("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs", recursive = F, full.names = F) %>% 
  grep("GSM585867[0|1|2|3]", ., value = T)
metrics.dfs <- lapply(dir.list, function(x) {
  read_csv(paste0("/gpfs/commons/groups/landau_lab/VEXAS/cellranger_outs/", x, "/outs/metrics_summary.csv")) %>% 
    mutate(sample = x)
})

## Load BM15
metrics.dfs <- c(metrics.dfs, lapply(c("SG_VX_BM15a_GEX", "SG_VX_BM15b_GEX"), function(x) {
  read_csv(paste0("/gpfs/commons/home/rmurray/VEXAS/", x, "/outs/metrics_summary.csv")) %>% 
    mutate(sample = x)
}))

metrics.df <- do.call(rbind, metrics.dfs)

metrics.df %>% select(sample, `Estimated Number of Cells`, `Number of Reads`, `Mean Reads per Cell`, `Median UMI Counts per Cell`)  %>% 
  write_csv("data/NIH_data/cellranger_qc_metrics.csv")
