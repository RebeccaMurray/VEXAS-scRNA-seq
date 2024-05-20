##################################### Sub-clustering the HSCs ##########################################

## Select the HSCs and recluster
rna.obj.HSC <- subset(rna.obj, CellType == "HSC")
DimPlot(rna.obj.HSC)

## Set up the new HSC-only RNA object
rna.obj.HSC <- RunUMAP(rna.obj.HSC, reduction = "pca", dims = 1:50, seed.use = 1, return.model = TRUE)
rna.obj.HSC <- FindNeighbors(rna.obj.HSC, reduction = "pca", dims = 1:50)
rna.obj.HSC <- FindClusters(rna.obj.HSC, resolution = 1)
DimPlot(rna.obj.HSC)
md.hsc <- rna.obj.HSC[[]]
coords <- Embeddings(rna.obj.HSC[["umap"]])
md.hsc <- cbind(md.hsc, coords) %>% dplyr::rename(`UMAP1` = umap_1) %>% dplyr::rename(`UMAP2` = umap_2)
md.hsc <- md.hsc[sample(1:nrow(md.hsc)), ] ## Shuffle the rows in the dataframe, so plotting order is random

## Genotyping heatmap, just the HSCs
md.hsc.wt <- md.hsc[md.hsc$Genotype %in% c("WT"), ]
md.hsc.mut <- md.hsc[md.hsc$Genotype %in% c("MUT"), ]
md.hsc.na <- md.hsc[md.hsc$Genotype %in% c("NA", "Control"), ]
p.genotyping <- ggplot(md.hsc.na, aes(x = UMAP1, y = UMAP2, color = Genotype)) +
  geom_point(shape = 19, size = 1, color = "grey80") + 
  geom_point(data = md.hsc.mut, shape = 19, size = 1, color = vexas.genotyping.palette[["MUT"]]  ) +
  geom_point(data = md.hsc.wt, shape = 19, size = 1, color = vexas.genotyping.palette[["WT"]]) +
  # geom_point(data = md.hsc.genotyped, shape = 21, size = 1) +
  scale_color_manual(values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.genotyping
ggsave(paste0(current.plots.path, "/UMAP_RNA_genotyping_dots_HSCs.pdf"), plot = p.genotyping, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Isolating ADT data and adding to dataframe
adt.counts.per.bc <- colSums(rna.obj.HSC@assays$ADT_DSB@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.HSC.adt <- subset(rna.obj.HSC, cells = bcs.keep)
rna.obj.HSC.adt <- ScaleData(rna.obj.HSC.adt, assay = "ADT_DSB")
adt.matrix <- rna.obj.HSC.adt@assays$ADT_DSB@data %>% t()
md.hsc.adt <- merge(md.hsc, adt.matrix, all.x = FALSE, by = "row.names")
md.hsc.no.adt <- md.hsc[!(rownames(md.hsc) %in% md.hsc.adt), ]
quantiles = c(0.05, 0.95)

## Make a helper function to highlight ADT expression on the UMAP
plot_adt_on_umap <- function(x) {
  quantiles.current <- md.hsc.adt %>% dplyr::rename(current_tag = x) %>% summarise(x = quantile(current_tag, quantiles), q = quantiles) %>% pull(x)
  md.hsc.adt <- md.hsc.adt[md.hsc.adt$Genotype != "NA", ]
  md.hsc.no.adt <- md.hsc.no.adt[md.hsc.no.adt$Genotype != "NA", ]
  p.adt.highlight <- md.hsc.no.adt %>% 
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    geom_point(color = "grey90", size = 0.5) +
    geom_point(data = md.hsc.adt, aes(color = !! sym(x)), size = 0.5) +
    # scale_color_viridis(limits = c(-2, 2), oob = scales::squish) +
    # scale_color_gradient(limits = c(0, 10), low = "white", high = "darkgreen", oob = scales::squish) +
    scale_color_gradient(limits = c(-quantiles.current[[1]], quantiles.current[[2]]), low = "white", high = "darkgreen", oob = scales::squish) +
    guides(x = axis, y = axis) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed() +
    facet_grid(cols = vars(Genotype), drop = T)
  return(p.adt.highlight)
}
plist <- lapply(c("Hu.CD11a", "Hu.CD26", "Hu.CD49d", "Hu.CD29", "Hu.CD13", "Hu.CD33", "Hu.HLA.DR", "Hu.CD71", "Hu.CD18"), plot_adt_on_umap)
p.combined <- wrap_plots(plist)
ggsave(paste0(current.plots.path, "/UMAP_RNA_HSCs_only_ADT_markers.pdf"), plot = p.combined, device = "pdf", dpi = 300, width = 12, height = 8, unit = "in")

