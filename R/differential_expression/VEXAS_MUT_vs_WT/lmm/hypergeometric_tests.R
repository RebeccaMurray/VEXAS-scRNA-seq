library(clusterProfiler)
library(org.Hs.eg.db)

cell.type.list <- c("HSC", "LMPP", "EMP", "MkP")
names(cell.type.list) <- cell.type.list
de.results <- lapply(cell.type.list, function (x) {
  read_csv(paste0("data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT/", x, "_no_mito_no_ribo_renormalized.csv")) %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = x)
})

ego.output <- lapply(names(de.results), function(x) {
  de.df <- de.results[[x]]
  gene.list.up <- de.df %>% filter(pval < 0.05 & avg_log2FC > 0.2) %>% pull(feature)
  gene.list.down <- de.df %>% filter(pval < 0.05 & avg_log2FC < 0.2) %>% pull(feature)
  universe <- rna.obj@assays$RNA@counts %>% rownames() %>% grep("^MT-|^RPL|^RPS", ., invert = T, value = T)
  
  ego.up <- enrichGO(gene = gene.list.up,
                  universe      = universe,
                  keyType = "SYMBOL",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  minGSSize = 30,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  ego.down <- enrichGO(gene = gene.list.down,
                     universe      = universe,
                     keyType = "SYMBOL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     minGSSize = 30,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  ego.up <- ego.up@result
  ego.up$direction = "Up"
  ego.down <- ego.down@result
  ego.down$direction = "Down"
  return.df <- rbind(ego.up, ego.down)
  return.df$cluster <- x
  return(return.df)

})

ego.output.df <- do.call(rbind, ego.output)

ego.output.df %>% filter(direction == "Up") %>% filter(p.adjust < 0.2) %>% group_by(Description) %>% mutate(cluster_count = n()) %>% arrange(-n)



