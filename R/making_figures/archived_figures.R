###################################### Plots - MUT cell frequencies ################################################

median.pseudotime.values <- md %>% 
  group_by(CellType) %>% 
  summarize(median_pseudotime = median(monocle3_pseudotime))

## Summarize the MUT cell frequency per donor
mut.fraction.per.donor <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  mutate(CellType = gsub("GMP [0-9]", "GMP", CellType)) %>% 
  mutate(CellType = gsub("CD4 |CD8 ", "", CellType))  %>%
  filter(n() > 100) %>% ## Only including cell types with > 20 genotyped cells
  group_by(CellType, Donor) %>% 
  filter(n() > 10) %>%  ## Only including patient data points with > 10 genotyped cells
  group_by(CellType, Genotype, Donor) %>%
  summarize(n = n()) %>% 
  group_by(CellType, Donor) %>% 
  mutate(mut_fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT") %>% 
  group_by(CellType) %>% 
  filter(n_distinct(Donor) >= 3)

summarized.freqs <- mut.fraction.per.donor %>% 
  summarize(mean_val = mean(mut_fraction), sd_val = sd(mut_fraction), n = n(), se_val = sd_val / sqrt(n), median_val = median(mut_fraction)) %>% 
  merge(., median.pseudotime.values, group.by = CellType, all.x = T)

## Make the bar plot
p.mutant_cell_fraction.bar.plot <- summarized.freqs %>% 
  ggplot(aes(x = reorder(CellType, -mean_val), y = mean_val, fill = CellType)) +
  geom_bar(stat = "identity", color = "black", size = 0.25) +
  geom_errorbar(aes(x = CellType, ymin = mean_val-se_val, ymax = mean_val+se_val), width=0.4, size = 0.25) +
  geom_point(data = mut.fraction.per.donor, mapping = aes(x = CellType, y = mut_fraction, shape = Donor), position = position_jitter(width = 0.25), size = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_shape_manual(values = seq(1, 10, 1)) +
  scale_fill_manual(values = cell.type.palette, guide = "none") +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Mutant cell fraction", x = element_blank()) + coord_fixed(9)
p.mutant_cell_fraction.bar.plot
ggsave(paste0(current.plots.path, "/MUT_cell_freq_bar_plot.pdf"), plot = p.mutant_cell_fraction.bar.plot, device = "pdf", dpi = 300, width = 6, height = 4, unit = "in")

## Save as a csv
summarized.freqs %>% 
  write_csv("data/metadata/mut_cell_fractions_RNA.csv")


###################################### Plots - Cytokine bar plot ################################################

rna.obj.genotyped <- subset(rna.obj, Genotype %in% c("MUT", "WT"))
Idents(rna.obj.genotyped) <- "CellType"

celltypes.plot <- c("HSC", "LMPP", "EMP", "MkP")
cytokines <- c("IL1B", "IL12A", "CXCL8", "CXCL3","TGFB1", "DDIT3")

m <- rna.obj.genotyped@assays$RNA@data[cytokines, ] %>% as.matrix() %>% t()
m <- merge(rna.obj.genotyped@meta.data %>% dplyr::select(Genotype, CellType, Donor), m, by = "row.names")

## Plotting % expressed logFC between MUT and WT per patient, splitting by genotyping
m %>% 
  filter(CellType %in% celltypes.plot) %>% 
  pivot_longer(cols = cytokines, names_to = "gene", values_to = "expression") %>% 
  merge(., gene.df, by = "gene") %>% 
  group_by(CellType, Genotype, gene, Donor) %>% 
  summarize(count = n(), mean_exp = mean(expression), nonzero_pct = sum(expression > 0) / n()) %>% View
filter(count > 5) %>% 
  dplyr::select(-count) %>% 
  pivot_wider(names_from = Genotype, values_from = c(nonzero_pct, mean_exp)) %>% 
  mutate(nonzero_pct_fc = log2(nonzero_pct_MUT / nonzero_pct_WT)) %>% 
  # filter(CellType == "HSC", gene == "IL1B") %>% 
  ggplot(aes(x = gene, y = nonzero_pct_fc)) +
  geom_boxplot() +
  facet_grid(rows = vars(CellType)) +
  geom_hline(yintercept = 0, linetype = "dashed")

filter(count > 5) %>% 
  ggplot(aes(x = gene, y = nonzero_pct, fill = Genotype)) +
  # coord_cartesian(ylim=c(0, 0.8)) 
  # geom_col(position = "dodge") +
  # geom_boxplot() +
  geom_violin() +
  geom_point(position = "jitter") +
  facet_grid(rows = vars(CellType)) +
  stat_compare_means(label = "p.signif") +
  scale_fill_manual(values = vexas.genotyping.palette) +
  theme(legend.position = "none") +
  ylab("Percent of cells with expression")

## Plotting % expressed, no splitting by genotyping, scaling color by expression
m %>% 
  filter(CellType %in% celltypes.plot) %>% 
  pivot_longer(cols = gene.list, names_to = "gene", values_to = "expression") %>% 
  merge(., gene.df, by = "gene") %>% 
  group_by(CellType, category, Genotype, gene) %>% 
  summarize(mean_exp = mean(expression), nonzero_pct = sum(expression > 0) / n()) %>% 
  ggplot(aes(x = gene, y = nonzero_pct, fill = mean_exp)) +
  coord_cartesian(ylim=c(0, 0.8)) +
  facet_grid(cols = vars(CellType), rows = vars(Genotype), space = "free", scales = "free") +
  geom_col(position = "dodge") +
  # scale_fill_manual(values = vexas.genotyping.palette) +
  theme(legend.position = "none")

###################################### Plots - correlating ADTs ################################################

adt.list <- c("Hu.CD49b", "Hu.CD99", "Hu.CD38-HIT2")
HSC.subset <- subset(rna.obj, subset = Donor %in% c("UPN2", "UPN12", "UPN13", "UPN14", "UPN15") & cluster_celltype == "HSC" & Genotype %in% c("MUT", "WT"))

# Add ADT expression
m <- HSC.subset@assays$ADT_DSB@data[adt.list, ] %>% t()
HSC.subset <- AddMetaData(HSC.subset, m)

# Add gene expression
gene.list <- c("ITGA2", "CD99", "CD38")
m <- HSC.subset@assays$RNA@data[gene.list, ] %>% t()
HSC.subset <- AddMetaData(HSC.subset, m)

HSC.subset@meta.data %>% 
  ggplot(aes(x = Hu.CD99, y = CD99, color = Genotype)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = vexas.genotyping.palette)

HSC.subset@meta.data %>% 
  ggplot(aes(x = `Hu.CD38-HIT2`, y = CD38, color = Genotype)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = vexas.genotyping.palette)


## Plot VEXAS vs. control
VEXAS.count <- rna.obj@meta.data %>% filter(SampleType == "VEXAS") %>% nrow()
control.count <- rna.obj@meta.data %>% filter(SampleType == "Control") %>% nrow()
vexas.title <- paste0("VEXAS (n = ", format(VEXAS.count, big.mark = ","), ")")
control.title <- paste0("Control (n = ", format(control.count, big.mark = ","), ")")
p.vexas.vs.control <- ggplot(md, aes(x = UMAP_1, y = UMAP_2, fill = SampleType)) +
  geom_point(shape = 21, size = 1) + theme_classic() + coord_fixed() +
  scale_fill_manual(labels = c(VEXAS = vexas.title, Control = control.title), values = c("VEXAS" = vexas.color, "Control" = control.color)) +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 14)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))
p.vexas.vs.control


## Mut cell frquency - no crosshairs
p.mutant_cell_fraction.scatter.plot <- mut.freqs.merged %>% 
  ggplot(aes(x = mean_val_RNA, y = mean_val_ATAC)) +
  geom_point(size = 0) +
  geom_abline(linetype = "longdash", color = "grey70") +
  stat_cor(aes(x = mean_val_RNA, y = mean_val_ATAC), method = "pearson") +
  # geom_point(aes(color = CellType), size = 4) +
  geom_point(aes(color = CellType), size = 1) +
  scale_color_manual(aes(color = CellType, label = CellType), values = cell.type.palette) +
  geom_text_repel(aes(color = CellType, label = CellType), point.padding = 10) +
  # geom_crossbar(mapping = aes(x = mean_val_RNA, 
  #                      ymin = mean_val_ATAC - sd_val_ATAC, 
  #                      ymax = mean_val_ATAC + sd_val_ATAC, color = CellType)) + 
  geom_crossbar(mapping = aes(y = mean_val_ATAC, 
                              xmin = (mean_val_RNA - sd_val_RNA), 
                              xmax = (mean_val_RNA + sd_val_RNA), color = CellType)) + 
  theme(legend.position = "none") + 
  ggtitle("Mutant cell frequencies", subtitle = "Genotyping from RNA (GoT) vs. ATAC (GoTChA)") +
  coord_fixed() +
  xlab("Mean MUT cell fraction,\nRNA genotyping") +
  ylab("Mean MUT cell fraction,\nATAC genotyping") +
  xlim(c(0, 1)) + ylim(c(0, 1))


############################################### XBP1s ###############################################

## Option 1: filter for 50 MUT and WT cells
pt.keep <- md %>% 
  group_by(Donor, Genotype) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count() %>% 
  filter(n > 50) %>% 
  group_by(Donor) %>% 
  count() %>% 
  filter(n > 1) %>% 
  pull(Donor)

## Option 2: filter for at least 20 genotyped cells
pt.keep <- md %>% 
  group_by(Donor, Genotype) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count() %>% 
  filter(n > 5) %>% 
  group_by(Donor) %>% 
  count() %>% 
  filter(n > 1) %>% 
  filter(Donor != "PT12") %>% ## Has 0 XBP1s reads detected
  pull(Donor)



## Line plot
p.xbp1s <- md %>% 
  filter(Donor %in% pt.keep) %>%
  mutate(Donor = factor(Donor, levels = pt.keep)) %>% 
  group_by(Donor, Genotype) %>% 
  # filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP", "Early Eryth", "Late Eryth", "GMP")) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  summarize(Spliced_Sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), Unspliced_Sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T)) %>% 
  mutate(spliced_ratio = log2(Spliced_Sum / Unspliced_Sum)) %>% 
  ggplot(aes(x = Genotype, y = spliced_ratio, color = Donor)) +
  geom_point() +
  scale_color_manual(values = ggsci::pal_aaas()(10)) +
  ylab("Log2(XBP1s/XBP1u)") +
  geom_line(aes(group = Donor)) +
  ggtitle("HSPCs")
p.xbp1s
ggsave(paste0(current.plots.path, "/RNA_XBP1s_line_plot_all_donors.pdf"), plot = p.xbp1s, dpi = 300, device = "pdf", width = 3, height = 4)


## Other plot ideas
md %>% filter(Donor %in% pt.keep) %>% 
  group_by(Genotype) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  summarize(Spliced_Sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), Unspliced_Sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T)) %>% 
  mutate(spliced_frac = Spliced_Sum / (Spliced_Sum + Unspliced_Sum))

md %>% 
  # filter(Donor %in% pt.keep) %>% 
  group_by(Donor, Genotype) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  summarize(Spliced_Avg = sum(Spliced_XBP1_umi_filtered, na.rm = T)/n(), Unspliced_Avg = sum(Unspliced_XBP1_umi_filtered, na.rm = T)/n()) %>% 
  mutate(ratio_avg = Spliced_Avg/Unspliced_Avg) %>% 
  ggplot(aes(x = Genotype, y = ratio_avg, color = Donor)) +
  geom_point() +
  geom_line(aes(group = Donor))

md %>% 
  # filter(Donor %in% pt.keep) %>% 
  group_by(Donor, Genotype) %>% 
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  # filter(Unspliced_XBP1_umi_filtered > 0) %>% 
  # mutate(Spliced_Ratio = Spliced_XBP1_umi_filtered / Unspliced_XBP1_umi_filtered) %>% 
  ggplot(aes(x = Donor, y = Unspliced_XBP1_umi_filtered, fill = Genotype)) +
  geom_violin(scale = "width") +
  stat_compare_means(label = "p.format")


## UMAP
rna.obj.splicing <- subset(rna.obj, Unspliced_XBP1_umi_filtered > 0)
rna.obj.splicing$log2_splicing_ratio_umi_filtered <- log2(rna.obj.splicing$Spliced_XBP1_umi_filtered / rna.obj.splicing$Unspliced_XBP1_umi_filtered)
rna.obj.splicing$splicing_ratio_umi_filtered <- (rna.obj.splicing$Spliced_XBP1_umi_filtered / rna.obj.splicing$Unspliced_XBP1_umi_filtered)
FeaturePlot(rna.obj.splicing, features = "splicing_ratio_umi_filtered", max.cutoff = 1, cols = rev(viridis_pal()(2))) +
  theme_classic() + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + ggtitle("XBP1s/XBP1u")
FeaturePlot(rna.obj.splicing, features = "log10_splicing_ratio", min.cutoff = , max.cutoff = 0, cols = rev(viridis_pal()(2)))
FeaturePlot(rna.obj.splicing, features = "splicing_ratio", min.cutoff = 0, max.cutoff = 5)


## HSPCs
p.xbp1s <- rna.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>%
  filter(Unspliced_XBP1 > 0) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(!(Donor %in% c("UPN12", "UPN13", "UPN14", "UPN25"))) %>% 
  ggplot(aes(x = Genotype, y = splicing_ratio, fill = Genotype)) +
  geom_violin() +stat_compare_means(label = "p.format", label.y = 4.5)+ scale_fill_manual(values = vexas.genotyping.palette) +
  # geom_point(position = position_jitterdodge(), size = 0.1) + 
  coord_cartesian(ylim = c(0, 5)) +
  theme(legend.position = "none") +
  ylab("XBP1s/XBP1u") + ggtitle("HSPCs") +
  facet_grid(cols = vars(Donor), scales = "free")
p.xbp1s
ggsave(paste0(current.plots.path, "/XBP1s_volcano_HSC.pdf"), plot = p.xbp1s, dpi = 300, device = "pdf", width = 2.5, height = 4)


################################ ADT highlights ############################

## Highlight CD3, CD56
p.adt.highlight <- FeaturePlot(rna.obj.adt, features = c("Hu.CD3-UCHT1"), slot = "data", min.cutoff = "q10", max.cutoff = "q90", cols = scales::viridis_pal()(2)) +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  ggtitle("Hu.CD3") +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.adt.highlight
ggsave(paste0(current.plots.path, "/UMAP_ADT_cd3_highlight.pdf"), plot = p.adt.highlight, device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")

p.adt.highlight <- FeaturePlot(rna.obj.adt, features = c("Hu.CD56"), slot = "data", min.cutoff = "q10", max.cutoff = "q90", cols = scales::viridis_pal()(2)) +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  ggtitle("Hu.CD56") +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.adt.highlight
ggsave(paste0(current.plots.path, "/UMAP_ADT_cd56_highlight.pdf"), plot = p.adt.highlight, device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")

############################### Genotyping ###################################
## Plot genotyping
mut.count <- md %>% filter(Genotype == "MUT") %>% nrow()
wt.count <- md %>% filter(Genotype == "WT") %>% nrow()
na.or.het.count <- md %>% filter(Genotype == "NA") %>% nrow()
mut.title <- paste0("MUT (n = ", format(mut.count, big.mark = ","), ")")
wt.title <- paste0("WT (n = ", format(wt.count, big.mark = ","), ")")
na.title <- paste0("NA (n = ", format(na.or.het.count, big.mark = ","), ")")
md.wt <- md[md$Genotype %in% c("WT"), ]
md.mut <- md[md$Genotype %in% c("MUT"), ]
md.genotyped <- md[md$Genotype %in% c("MUT", "WT"), ]
md.na <- md[md$Genotype %in% c("NA", "Control"), ]

p.genotyping <- ggplot(md.na, aes(x = umap_1, y = umap_2, fill = Genotype)) +
  geom_point(shape = 19, size = 0.25, color = "grey80") + 
  geom_point(data = md.mut, shape = 21, size = 2, color = "grey") +
  geom_point(data = md.wt, shape = 21, size = 2, color = "grey") +
  # geom_point(data = md.genotyped, shape = 21, size = 1) +
  scale_fill_manual(labels = c(
    `VEXAS MUT` = mut.title,
    `VEXAS WT` = wt.title,
    `VEXAS unknown` = na.title,
    `Control` = na.title
  ), values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 8)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.genotyping
## Second version
p.genotyping <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, color = Genotype)) +
  geom_point(shape = 19, size = 0.25, color = "grey80") + 
  geom_point(data = md.mut, shape = 19, size = 1, color = "red") +
  geom_point(data = md.wt, shape = 19, size = 1, color = "blue") +
  # geom_point(data = md.genotyped, shape = 21, size = 1) +
  # scale_color_manual(labels = c(
  #   `MUT` = mut.title,
  #   `WT` = wt.title,
  #   `NA` = na.title
  # ), values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 8)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.genotyping


## Other attempts at the genotyping heatmap
p.genotyping.mut <- ggplot(md[which(md$Genotype =="MUT"),], aes(x=UMAP1, y=UMAP2))+
  stat_density2d(aes(alpha=..level.., fill=..level..),
                 geom="polygon", breaks = seq(0, 0.01, 0.00001)) +
  scale_fill_gradient2(low = "white", high = vexas.genotyping.palette[[2]]) +
  # scale_alpha(range = c(0, 0.1), guide = FALSE) +
  coord_fixed(xlim = umap.lims.x, ylim = umap.lims.y) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none")
p.genotyping.mut
p.genotyping.mut <- ggplot(md[which(md$Genotype =="MUT"),], aes(x=UMAP1, y=UMAP2))+
  stat_density2d(aes(alpha=log10(..level..), fill=log10(..level..)), bins=5, geom="polygon", size = 1) +
  scale_fill_gradient2(low = "white", high = vexas.genotyping.palette[[2]]) +
  # scale_alpha_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, 0.001))+
  coord_fixed(xlim = umap.lims.x, ylim = umap.lims.y) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "right")
p.genotyping.mut
ggsave(paste0(current.plots.path, "/UMAP_RNA_genotyping_heatmap_mut.pdf"), plot = p.genotyping.mut, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

############################### ADTs ##########################################


## Create an seurat object with only cells with non-zero ADT expression
adt.counts.per.bc <- colSums(rna.obj@assays$ADT_CLR@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.adt <- subset(rna.obj, cells = bcs.keep)

## Dot plot to validate cell types with ADTs - DSB-normalized, by cluster assignments
all.tags <- rna.obj@assays$ADT_DSB %>% rownames()
to.plot.tags <- c(grep("Isotype|CTRL", all.tags, value = T), 
                  "Hu.CD19", "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD14-M5E2",  "Hu.CD16", "Hu.CD56",
                  "Hu.CD7", "Hu.CD41", "Hu.CD71", "Hu.HLA.DR", "Hu.CD36", "Hu.CD99", "Hu.CD90", "Hu.CD49b", "HuMs.CD49f", "Hu.CD38-HIT2", "Hu.CD45RA", "Hu.CD123")
DefaultAssay(rna.obj.adt) <- "ADT_DSB"
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = to.plot.tags, scale = F) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "cluster_celltype", features = to.plot.tags) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Dot plot to validate cell types with ADTs - DSB-normalized, by cell type label
rna.obj.adt.celltype <- subset(rna.obj.adt, predicted.celltype.l2.score > 0.5)
DotPlot(rna.obj.adt.celltype, assay = "ADT_DSB", group.by = "predicted.celltype.l2", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt.celltype, assay = "ADT_DSB", group.by = "predicted.celltype.l2", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




############################### Cell types ###################################
## Plot cell types (no dot outline)
p1 <- DimPlot(rna.obj, group.by = "CellType", label = T, label.box = T, repel = T, cols = cell.type.palette) + NoLegend()
p.celltypes <- p1 + ggtitle("") + coord_fixed() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))
p.celltypes


## Plotting mut/wt cell frequency per cluster - bar plot with dots
mut.cell.fractions.per.sample <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType, Sample) %>% 
  filter(n() > 10) %>% 
  group_by(CellType, Genotype, Sample) %>%
  count() %>%
  group_by(CellType, Sample) %>% 
  mutate(fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT")
p.mut.cell.fractions <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  filter(n() > 20) %>% 
  group_by(CellType, Genotype) %>%
  count() %>%
  group_by(CellType) %>% 
  mutate(fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT") %>% 
  ggplot(aes(x = reorder(CellType, -fraction), y = fraction, fill = CellType)) +
  geom_col() +
  scale_shape_manual(values = c(seq(1:10))) +
  scale_fill_manual(values = cell.type.palette) +
  geom_point(data = mut.cell.fractions.per.sample, aes(x = CellType, y = fraction, shape = Sample)) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(fill = F) +
  labs(y = "Mutant cell fraction", x = element_blank()) + coord_fixed(9)
p.mut.cell.fractions

## Plot celltypes - with cell outlines
p.celltypes <- ggplot(md, aes(x = UMAP1, y = UMAP2, fill = CellType)) +
  geom_point(pch = 21, size = 2, stroke = 0.1, color = "black") + ## size controls width of point, stroke controls width of border, color is border color
  scale_fill_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 8)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.celltypes



############################### ADTs ###################################

############# Checking cell type assignments with ADT ################
all.tags <- rna.obj@assays$ADT_DSB %>% rownames()
to.plot.tags <- c(grep("Isotype|CTRL", all.tags, value = T), 
                  "Hu.CD19", "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD14-M5E2",  "Hu.CD16", "Hu.CD56",
                  "Hu.CD7", "Hu.CD41", "Hu.CD71", "Hu.HLA.DR", "Hu.CD36", "Hu.CD99", "Hu.CD90", "Hu.CD49b", "HuMs.CD49f", "Hu.CD38-HIT2", "Hu.CD45RA", "Hu.CD123")

## CLR normalization
adt.counts.per.bc <- colSums(rna.obj@assays$ADT_CLR@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.adt <- subset(rna.obj, cells = bcs.keep)

DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "cluster_celltype", features = to.plot.tags) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "cluster_celltype", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "cluster_celltype", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
rna.obj.adt <- subset(rna.obj.adt, predicted.celltype.l2.score > 0.5)
DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "predicted.celltype.l2", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_CLR", group.by = "predicted.celltype.l2", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## DSB normalization
adt.counts.per.bc <- colSums(rna.obj@assays$ADT_DSB@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
rna.obj.adt <- subset(rna.obj, cells = bcs.keep)

DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = to.plot.tags) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
rna.obj.adt <- subset(rna.obj.adt, predicted.celltype.l2.score > 0.5)
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "predicted.celltype.l2", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(rna.obj.adt, assay = "ADT_DSB", group.by = "predicted.celltype.l2", features = all.tags[75:138], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

atac.obj.adt.hpscs <- subset(atac.obj.adt, cluster_celltype %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk"))
DotPlot(atac.obj.adt.hpscs, assay = "ADT_DSB", features = to.plot.tags, group.by = "cluster_celltype", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
atac.obj.adt.hpscs <- subset(atac.obj.adt, predicted.l2 %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk") & predicted.l2.score > 0.5)
DotPlot(atac.obj.adt.hpscs, assay = "ADT_DSB", features = to.plot.tags, group.by = "predicted.l2", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



## Plotting ADTs - DSB
rna.obj.adt <- subset(rna.obj, SampleType == "VEXAS" & Donor != "VEXAS_BM1")
DefaultAssay(rna.obj.adt) <- "ADT_DSB"
avg.expression.per.cluster <- list()
for (cluster in c("LMPP", "EMP", "HSC", "Prog Mk")) {
  mut.barcodes <- rna.obj.adt@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "MUT" & cluster_celltype == cluster) %>% rownames()
  wt.barcodes <- rna.obj.adt@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "WT" & cluster_celltype == cluster) %>% rownames()
  avg_expression <- data.frame(MUT = rowMeans(rna.obj.adt@assays$ADT_DSB@data[, mut.barcodes], na.rm = T), WT = rowMeans(rna.obj.adt@assays$ADT_DSB@data[, wt.barcodes], na.rm = T))
  avg_expression$CellType <- cluster
  avg.expression.per.cluster[[cluster]] <- avg_expression
}
plist <- lapply(names(avg.expression.per.cluster), function(cluster.name) {
  df <- avg.expression.per.cluster[[cluster.name]]
  as.data.frame(df) %>% 
    rownames_to_column("protein") %>% 
    ggplot(aes(x = WT, y = MUT, label = protein)) + 
    geom_point() +
    theme_bw() +
    geom_label_repel(max.overlaps = 15) + 
    geom_abline() + ggtitle(cluster.name, subtitle = "DSB normalized ADT expression") + 
    theme_classic() +  
    theme(legend.position = "none")
})

wrap_plots(plist)


## Plotting ADTs to annotate clusters - CLR
markers.list <- c("Hu.CD3-UCHT1", "Hu.CD19", "Hu.CD14-M5E2", "Hu.CD71", "Hu.CD27", "Hu.CD90", "Hu.CD16", "Hu.CD64")
rna.obj.adt <- subset(rna.obj, SampleType == "VEXAS" & Donor != "VEXAS_BM1" & Donor != "VEXAS_BM2")
cluster.list <- levels(rna.obj.adt$cluster_celltype)
names(cluster.list) <- cluster.list
exp.per.cluster <- lapply(cluster.list, function(cluster.name) {
  cluster.bcs <- rna.obj.adt@meta.data %>% filter(cluster_celltype == cluster.name) %>% rownames()
  adt.exp.matrix <-rowMeans(rna.obj.adt@assays$ADT_DSB@data[, cluster.bcs], na.rm = T)
  return(adt.exp.matrix)  
})
exp.per.cluster <- as.data.frame(exp.per.cluster)

pheatmap::pheatmap(exp.per.cluster[markers.list, ], 
                   scale = "row",
                   color=colorRampPalette(c("white", "navy"))(50))
Idents(rna.obj.adt) <- "cluster_celltype"
DefaultAssay(rna.obj.adt) <- "ADT_DSB"
rna.obj.adt.sub <- subset(rna.obj.adt, subset = orig.ident == "VX_BM15a")
VlnPlot(rna.obj.adt.sub, features = "Hu.CD90", assay = "ADT_DSB", group.by = "cluster_celltype")
FeaturePlot(rna.obj.adt, features= "Hu.CD90", min.cutoff = "q10", max.cutoff = "q90")

## Plotting ADTs
rna.obj.adt <- subset(rna.obj, SampleType == "VEXAS" & Donor != "VEXAS_BM1")
DefaultAssay(rna.obj.adt) <- "ADT_CLR"
avg.expression.per.cluster <- list()
for (cluster in c("LMPP", "EMP", "HSPC", "Prog Mk")) {
  mut.barcodes <- rna.obj.adt@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "MUT" & cluster_celltype == cluster) %>% rownames()
  wt.barcodes <- rna.obj.adt@meta.data %>% filter(Genotype_Approx_Match_No_Gene_Thresh == "WT" & cluster_celltype == cluster) %>% rownames()
  avg_expression <- data.frame(MUT = rowMeans(rna.obj.adt@assays$ADT_CLR@data[, mut.barcodes], na.rm = T), WT = rowMeans(rna.obj.adt@assays$ADT_CLR@data[, wt.barcodes], na.rm = T))
  avg_expression$CellType <- cluster
  avg.expression.per.cluster[[cluster]] <- avg_expression
}

plist <- lapply(names(avg.expression.per.cluster), function(cluster.name) {
  df <- avg.expression.per.cluster[[cluster.name]]
  as.data.frame(df) %>% 
    rownames_to_column("protein") %>% 
    ggplot(aes(x = WT, y = MUT, label = protein)) + 
    geom_point() +
    theme_bw() +
    geom_label_repel(max.overlaps = 15) + 
    geom_abline() + ggtitle(cluster.name, subtitle = "CLR normalized ADT expression") + theme_classic() +  theme(legend.position = "none")  + coord_cartesian(xlim = c(0, 3), ylim = c(0, 3))
})

wrap_plots(plist)



## Plot module scores for inflammation, etc
make_module_boxplot <- function(celltype, module.name) {
  
  md.filtered <- md %>% filter(predicted.celltype.l2.score > 0.6 & predicted.celltype.l2 == celltype)
  md.filtered <- md
  
  sample.counts <- md.filtered %>% filter(CellType == celltype) %>% pull(SampleType) %>% table()
  vexas.label <- paste0("VEXAS\n(n = ", format(sample.counts[["VEXAS"]], big.mark = ","), ")")
  control.label <- paste0("Control\n(n = ", format(sample.counts[["Control"]], big.mark = ","), ")")
  
  md.filtered %>% 
    dplyr::rename(module_to_test = module.name) %>% 
    filter(CellType == celltype) %>%
    ggplot(aes(x = SampleType, y = log2(module_to_test), fill = SampleType)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(x = "", y = "log2(module score)") +
    # scale_fill_manual(values = vexas.vs.control.palettes[[celltype]]) +
    ggpubr::stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0) +
    ggtitle(module.name) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c(VEXAS = vexas.label, Control = control.label)) +
    expand_limits(y = 0.5)
}

make_module_boxplot_for_celltype <- function(celltype.value) {
  plist <- list()
  plist[[1]] <- make_module_boxplot(celltype.value, "InterferonAlphaHallmark1") + ggtitle("Interferon alpha module")
  plist[[2]] <- make_module_boxplot(celltype.value, "InterferonGammaHallmark1") + ggtitle("Interferon gamma module")
  plist[[3]] <- make_module_boxplot(celltype.value, "InflammatoryHallmark1") + ggtitle("Inflammatory module")
  plist[[4]] <- make_module_boxplot(celltype.value, "XBP1_target_module1") + ggtitle("XBP1 module")
  return(plist)
}

make_got_paper_boxplot_for_celltype <- function(celltype.value) {
  plist <- list()
  plist[[1]] <- make_module_boxplot(celltype.value, "PERK_module1") + ggtitle("PERK module")
  plist[[2]] <- make_module_boxplot(celltype.value, "ATF6_module1") + ggtitle("AFT6 module")
  plist[[3]] <- make_module_boxplot(celltype.value, "XBP1_target_module1") + ggtitle("XBP1 module")
  plist[[4]] <- make_module_boxplot(celltype.value, "Anti.apoptosis_module1") + ggtitle("Anti-apoptosis module")
  plist[[5]] <- make_module_boxplot(celltype.value, "OxPhosHallmark1") + ggtitle("OxPhos module")
  return(plist)
}

boxplot.list <- list("HSC", "EMP", "GMP", "LMPP", "Early Eryth", "CD14 Mono")
names(boxplot.list) <- boxplot.list

hallmark.module.boxplots <- lapply(boxplot.list, function(x) {
  wrap_plots(make_module_boxplot_for_celltype(x), nrow = 1) + plot_annotation(x)
})

got.paper.module.boxplots <- lapply(boxplot.list, function(x) {
  wrap_plots(make_got_paper_boxplot_for_celltype(x), nrow = 1) + plot_annotation(x)
})


wrap_plots(hallmark.module.boxplots, ncol = 1)
wrap_plots(hsc.boxplots, nrow = 1) / wrap_plots(emp.boxplots, nrow = 1) / wrap_plots(gmp.boxplots, nrow = 1) / wrap_plots(lmpp.boxplots, nrow = 1)
# combined.boxplots <- hsc.boxplots / emp.boxplots



## Plot CD90
cd49b.subset <- subset(rna.obj, subset = SampleType == "VEXAS" & Donor != "VEXAS_BM1" & cluster_celltype %in% c("HSPC", "EMP", "LMPP", "GMP", "Prog Mk") & Genotype %in% c("MUT", "WT"))
m <- cd49b.subset@assays$ADT_CLR@data["Hu.CD90", ] %>% data.frame(CD49b = .)
cd49b.subset <- AddMetaData(cd49b.subset, m)
cd49b.subset@meta.data %>% 
  ggplot(aes(x = cluster_celltype, y = CD49b, fill = cluster_celltype)) +
  geom_boxplot() +
  stat_compare_means(label.y = 4) +
  theme(legend.position = "none")


## Plot NK cell genotypes - as stacked bar plots
md %>% 
  filter(predicted.celltype.l1 %in% c("NK", "CD4 T", "CD8 T")) %>% 
  mutate(Category = if_else(predicted.celltype.l1 %in% c("CD4 T", "CD8 T"), "CD4/CD8 T", "NK")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(orig.ident) %>% 
  filter(n() > 5) %>% 
  ungroup() %>% 
  count(Genotype, Category) %>%
  group_by(Category) %>% 
  mutate(MUT_fraction = n / sum(n)) %>% 
  # filter(Genotype == "MUT") %>%
  ggplot(aes(Category, y = MUT_fraction, fill = Genotype)) +
  geom_col() + ggtitle("") +
  scale_fill_manual(values = vexas.genotyping.palette) +
  # stat_compare_means(label.y = 1.2) +
  xlab("Azimuth label") + ylab("MUT cell fraction (%)") +
  coord_cartesian(ylim=c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

## Plot NK cell genotypes - as box plots, based on azimuth label
md %>% 
  filter(CellType %in% c("CD4/CD8 T cells", "NK cells")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count(predicted.celltype.l1, orig.ident, Genotype)
md %>% 
  filter(predicted.celltype.l1 %in% c("NK", "CD4 T", "CD8 T")) %>% 
  mutate(Category = if_else(predicted.celltype.l1 %in% c("CD4 T", "CD8 T"), "CD4/CD8 T", "NK")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(orig.ident) %>% 
  filter(n() > 5) %>% 
  ungroup() %>% 
  count(Genotype, orig.ident, Category) %>%
  group_by(orig.ident, Category) %>% 
  mutate(MUT_fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT") %>%
  ggplot(aes(Category, y = MUT_fraction, fill = Category)) +
  geom_boxplot() + ggtitle("") +
  scale_fill_manual(values = c(`CD4/CD8 T` = "#DE9ED6", NK = "#CE6DBD")) + 
  stat_compare_means(label.y = 1.2) +
  xlab("Azimuth label") + ylab("MUT cell fraction") +
  coord_cartesian(ylim=c(0, 1.2)) +
  theme(legend.position = "none")

## Plot NK cell genotypes - as box plots, based on cluster assignment
md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  filter(n() > 20) %>% ## Only including cell types with > 20 genotyped cells
  group_by(CellType, Sample) %>% 
  filter(n() > 5) %>% 
  group_by(CellType, Genotype, Sample) %>%
  # filter(n() > 10) %>%  ## Only including patient data points with > 10 genotyped cells
  count() %>% 
  group_by(CellType, Sample) %>% 
  mutate(mut_fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT") %>% 
  filter(CellType %in% c("CD4/CD8 T cells", "NK cells")) %>% 
  group_by(CellType) %>% 
  ggplot(aes(x = CellType, y = mut_fraction, fill = CellType))+
  geom_boxplot() +
  stat_compare_means()


############################### VEXAS vs Control ###################################
## VEXAS vs. control - plotting lfc of a few key genes - using averaging function
gene.list <- c("EEF1A1", "UBA1", "TNF", "LTA", "IL1B", "CXCL8", "ISG15", "IFNG", "IL6")
cluster.list <- c("HSPC", "LMPP", "EMP", "GMP", "Prog Mk")
Idents(rna.obj) <- "cluster_celltype"

DefaultAssay(rna.obj) <- "RNA"
vex.rna.obj <- subset(rna.obj, SampleType == "VEXAS")
as.data.frame() %>% 
  mutate(SampleType = "VEXAS") %>% 
  rownames_to_column("feature")
ctrl.rna.obj <- subset(rna.obj, SampleType == "Control" & Donor != "N-02")
ctrl.avg.exp <- AverageExpression(ctrl.rna.obj, features = gene.list, group.by = "cluster_celltype", assay = "RNA")[[1]] %>% 
  as.data.frame() %>% 
  mutate(SampleType = "Control") %>% 
  rownames_to_column("feature")
avg.exp <- rbind(vexas.avg.exp, ctrl.avg.exp)
avg.exp %>%
  filter(feature != "EEF1A1") %>% 
  pivot_longer(cols = levels(rna.obj$cluster_celltype), names_to = "cluster") %>%
  filter(cluster %in% cluster.list) %>% 
  group_by(feature, cluster) %>% 
  mutate(cluster = factor(cluster, levels = c("HSPC", "LMPP", "EMP", "GMP", "Prog Mk"))) %>% 
  summarize(fc = (value[SampleType == "VEXAS"]) / (value[SampleType == "Control"])) %>% 
  vexas.avg.exp <- AverageExpression(vex.rna.obj, features = gene.list, group.by = "cluster_celltype", assay = "RNA")[[1]] %>% 
  ggplot(aes(x = feature, y = log10(fc), fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cell.type.palette) +
  facet_grid(rows = vars(cluster)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


################################## Monocytes ###############################


p.mono.gsea <- fgsea.out %>% 
  filter(padj < 0.25) %>% 
  mutate(status = sign(NES)) %>% 
  group_by(status) %>% 
  slice_min(n = 5, order_by = pval) %>% 
  mutate(pathway = gsub("_", " ", pathway) %>% str_to_title(.) %>% str_wrap(., width = 40)) %>% 
  dplyr::select(pathway, NES, padj) %>% 
  # column_to_rownames("pathway")
  mutate(`-log10(padj)*sign(NES)` = -log10(padj) * sign(NES)) %>% 
  ggplot(aes(x = reorder(pathway, NES), y = NES, fill = `-log10(padj)*sign(NES)`)) +
  geom_col() +
  xlab(element_blank()) +
  scale_fill_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], limits=c(-1, 1), oob = scales::squish) +
  coord_flip()
p.mono.gsea
ggsave("figures/current_figure_drafts/RNA_gsea_monocytes.pdf", plot = p.mono.gsea, device = "pdf", dpi = 300, width = 6, height = 2)

p.mono.dotplot <- fgsea.out %>% 
  filter(padj < 0.25) %>% 
  mutate(status = sign(NES)) %>% 
  group_by(status) %>% 
  slice_min(n = 5, order_by = pval) %>% 
  mutate(pathway = gsub("_", " ", pathway) %>% str_to_title(.) %>% str_wrap(., width = 15)) %>% 
  dplyr::select(pathway, NES, padj) %>% 
  # column_to_rownames("pathway")
  mutate(`-log10(padj)*sign(NES)` = -log10(padj) * sign(NES)) %>% 
  ggplot(aes(y = reorder(pathway, NES), x = "", fill = `-log10(padj)*sign(NES)`)) +
  geom_point(aes(size = abs(NES)), shape = 21) +
  ylab(element_blank()) +
  scale_size_continuous(range = c(3,8)) +
  xlab(element_blank()) +
  scale_fill_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], limits=c(-1.5, 1.5), oob = scales::squish) 
p.mono.dotplot
ggsave("figures/current_figure_drafts/RNA_dotplot_monocytes_custom.pdf", plot = p.mono.dotplot, device = "pdf", dpi = 300, width = 4, height = 3.5)

## Plot as heatmap
fgsea.out %>% 
  filter(padj < 0.1) %>% 
  mutate(pathway = gsub("_", " ", pathway) %>% str_to_title(.) %>% str_wrap(., width = 18)) %>% 
  dplyr::select(pathway, NES) %>%
  # mutate(`-log10(padj)*sign(NES)` = -log10(padj) * sign(NES)) %>% 
  # column_to_rownames("pathway")
  ggplot(aes(x = reorder(pathway, -`NES`), y = "const", fill = NES)) +
  geom_tile() +
  scale_fill_gradient2(low = vexas.genotyping.palette[[1]], mid = "white", high = vexas.genotyping.palette[[2]], limits = c(0, 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank()) + ylab("") 


############################################## XBP1s - sample-aware permutation test (umi-filtered) ###################################

donors.keep <- rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1_umi_filtered)) %>% 
  group_by(Donor) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  filter(frac_genotyped > 0) %>% 
  pull(Donor)

md.xbp1 <- rna.obj@meta.data %>% 
  filter(Donor %in% donors.keep) %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  filter(Genotype %in% c("MUT", "WT"))

## Helper function to calculate the mean difference
calculate_difference_of_means <- function(df) {
  df %>% 
    group_by(Genotype) %>% 
    summarize(spliced_sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), unspliced_sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T)) %>% 
    mutate(splicing_ratio = spliced_sum/unspliced_sum) %>% 
    select(Genotype, splicing_ratio) %>% 
    pivot_wider(values_from = splicing_ratio, names_from = Genotype) %>% 
    mutate(mean_diff = MUT - WT) %>% 
    pull(mean_diff)
}

## Shuffle the genotype labels and calculate within patients
md.xbp1.shuf <- md.xbp1
mean_diff_results <- sapply(1:10000, function(iter) {
  md.xbp1.shuf %>% 
    group_by(Donor) %>% 
    mutate(Genotype = sample(Genotype, replace = F)) %>% 
    calculate_difference_of_means(.)
})
hist(mean_diff_results)

## Calculate the actual mean difference
observed.mean.diff <- calculate_difference_of_means(md.xbp1)

sum(mean_diff_results > observed.mean.diff) / length(mean_diff_results)

############################################### Plots - XBP1s (umi-filtered) ##############################################################

## Make a bar plot with XBP1s values
donors.keep <- rna.obj@meta.data %>% 
  mutate(has_xbp1s_info = !is.na(Unspliced_XBP1_umi_filtered)) %>% 
  group_by(Donor) %>% 
  summarize(frac_genotyped = sum(has_xbp1s_info) / n()) %>% 
  filter(frac_genotyped > 0) %>% 
  pull(Donor)


## Pseudobulk by donor
p.xbp1.bars <- md %>% 
  filter(Donor %in% donors.keep) %>% 
  # filter(cluster_celltype == "HSC") %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  group_by(Genotype, Donor) %>% 
  summarize(spliced_sum = sum(Spliced_XBP1_umi_filtered, na.rm = T), unspliced_sum = sum(Unspliced_XBP1_umi_filtered, na.rm = T), cell_count = n()) %>% 
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
  annotate("text", x = 1.2, y = 0.4, label="  p < 0.0001", size = 3)
p.xbp1.bars
ggsave("figures/current_figure_drafts/RNA_xbp1s_column_plot.pdf", plot = p.xbp1.bars, device = "pdf", dpi = 300, width = 1.5, height = 2.5)

