
library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
theme_set(theme_classic())
set.seed(1)

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
current.plots.path <- c("figures/current_figure_drafts/")


############################################ Making figures - VEXAS + control data #########################################

rna.obj.control <- readRDS("/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pseudotime/vexas_pseudotime_cells.Rds")

## Normalize the RNA assay
rna.obj.control <- NormalizeData(rna.obj.control, assay = "RNA")
rna.obj.control@assays$RNA@counts[1:5, 1:5]
rna.obj.control@assays$RNA@data[1:5, 1:5]
DefaultAssay(rna.obj.control) <- "RNA"

# ## Format the "Sample" field
# # rna.obj.control@meta.data %>% count(orig.ident) %>% arrange(orig.ident) %>% write_csv("data/sample_ids.csv")
# meta.data.mapping <- read_csv("data/sample_ids.csv")
# rna.obj.control$Sample <- plyr::mapvalues(rna.obj.control$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$formatted)
# rna.obj.control@meta.data %>% count(orig.ident, Sample)
# 
# ## Update the "Donor" column for the healthy controls
# rna.obj.control$Donor[is.na(rna.obj.control$Donor)] <- rna.obj.control$orig.ident[is.na(rna.obj.control$Donor)] 
# healthy.sample.list <- unique(rna.obj.control$orig.ident) %>% grep("HD", ., value = T)
# rna.obj.control$Donor <- plyr::mapvalues(rna.obj.control$Donor, from = healthy.sample.list, to = gsub("_.*", "", healthy.sample.list)) %>% gsub("HD", "HD0", .)
# rna.obj.control@meta.data %>% count(Sample, Donor, orig.ident)

## Format the "Sample" field
meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_manually_updated_2024-03-29.csv")
rna.obj.control$Sample <- plyr::mapvalues(rna.obj.control$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$updated_name)
rna.obj.control@meta.data %>% count(orig.ident, Sample)

# ## Update the "Donor" column for the healthy controls
rna.obj.control$Donor <- gsub("_BMMC|_CD34.*", "", rna.obj.control$Sample) %>% gsub("_[A|B]", "", .)
rna.obj.control@meta.data %>% count(orig.ident, Sample, Donor)

## Saving metadata
rna.obj.control@meta.data %>% 
  group_by(Donor, Genotype_Approx_Match_No_Gene_Thresh) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype_Approx_Match_No_Gene_Thresh, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`, na.rm = T)) %>% 
  mutate(frac_genotyped = sum(MUT, WT, na.rm = T) / total) %>% 
  arrange(-frac_genotyped) %>% 
  write_csv("data/metadata/cell_counts_donor_genotype_controls_included.csv")

## nCount_RNA
sample.order <- rna.obj.control@meta.data %>% group_by(Sample) %>% summarize(mean_nCount_RNA = mean(nCount_RNA)) %>% arrange(-mean_nCount_RNA) %>% pull(Sample)
rna.obj.control@meta.data %>%
  mutate(Sample = factor(Sample, levels = sample.order)) %>%
  ggplot(aes(x = Sample, y = log10(nCount_RNA), fill = Object)) +
  geom_boxplot() +
  xlab(element_blank()) +
  scale_fill_manual(values = c(CONTROL = "#F5B935", VEXAS = "#2B858C")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
ggsave("figures/current_figure_drafts/qc_plots_nCount_RNA_vexas_and_control_20240328.pdf", device = "pdf", dpi = 300, width = 7, height = 4)

## nFeature_RNA
sample.order <- rna.obj.control@meta.data %>% group_by(Sample) %>% summarize(mean_nFeature_RNA = mean(nFeature_RNA)) %>% arrange(-mean_nFeature_RNA) %>% pull(Sample)
rna.obj.control@meta.data %>%
  mutate(Sample = factor(Sample, levels = sample.order)) %>%
  ggplot(aes(x = Sample, y = log10(nFeature_RNA), fill = Object)) +
  geom_boxplot() +
  scale_fill_manual(values = c(CONTROL = "#F5B935", VEXAS = "#2B858C")) +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
ggsave("figures/current_figure_drafts/qc_plots_nFeature_RNA_vexas_and_control_20240328.pdf", device = "pdf", dpi = 300, width = 7, height = 4)

## percent.mito
sample.order <- rna.obj.control@meta.data %>% group_by(Sample) %>% summarize(mean_pct_mito = mean(percent.mito)) %>% arrange(-mean_pct_mito) %>% pull(Sample)
rna.obj.control@meta.data %>%
  mutate(Sample = factor(Sample, levels = sample.order)) %>%
  ggplot(aes(x = Sample, y = percent.mito, fill = Object)) +
  geom_violin(scale = "width") +
  scale_fill_manual(values = c(CONTROL = "#F5B935", VEXAS = "#2B858C")) +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
ggsave("figures/current_figure_drafts/qc_plots_percent_mito_RNA_vexas_and_control_20240328.pdf", device = "pdf", dpi = 300, width = 7, height = 4)

## Sample of origin - split into separate panes (all samples)
p.sample.origin.split <- DimPlot(rna.obj.control, split.by = "Donor", group.by = "Donor", ncol = 7, pt.size = 0.1, raster = F) + NoLegend() +
  coord_fixed() +
  ggtitle(element_blank()) + theme(plot.title = element_text(size=2))
p.sample.origin.split
ggsave("figures/current_figure_drafts/qc_plots_rna_UMAP_orig_ident_split_vexas_and_controls_20240328.tiff", plot = p.sample.origin.split, device = "tiff", dpi = 300, width = 11, height = 4, unit = "in")

# ## Sample of origin - split into separate panes (controls only)
# p.sample.origin.split <- DimPlot(subset(rna.obj.control, Object == "CONTROL"), split.by = "Donor", group.by = "Donor", ncol = 9, pt.size = 0.05) + NoLegend() +
#   coord_fixed() +
#   ggtitle(element_blank()) + theme(plot.title = element_text(size=2))
# p.sample.origin.split
# ggsave("figures/current_figure_drafts/qc_plots_rna_UMAP_orig_ident_split_vexas_and_controls.pdf", plot = p.sample.origin.split, device = "pdf", dpi = 300, width = 8, height = 4, unit = "in")

## Bar plot of cell type proportions per donor
cell.types.include <- c("HSC", "EMP", "LMPP", "MkP", "GMP")
donor.order <- rna.obj.control@meta.data %>% 
  filter(AllCellType %in% cell.types.include) %>% 
  group_by(Donor, AllCellType) %>% 
  summarize(cell_type_count = n()) %>% 
  group_by(Donor) %>% 
  mutate(fraction = cell_type_count / sum(cell_type_count)) %>% 
  filter(AllCellType == "HSC") %>% 
  arrange(-fraction) %>% 
  pull(Donor)
p.celltype.fractions.per.donor <- rna.obj.control@meta.data %>% 
  mutate(Donor = factor(Donor, levels = donor.order)) %>% 
  filter(AllCellType %in% cell.types.include) %>% 
  group_by(Donor, AllCellType) %>% 
  summarize(cell_type_count = n()) %>% 
  group_by(Donor) %>% 
  mutate(fraction = cell_type_count / sum(cell_type_count)) %>% 
  ggplot(aes(x = Donor, y = fraction, fill = AllCellType)) +
  geom_col() +
  scale_fill_manual(values = cell.type.palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(element_blank()) +
  ylab("Percent of total cells") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "right") + guides(fill=guide_legend(ncol=3,byrow=TRUE))
p.celltype.fractions.per.donor
ggsave(paste0(current.plots.path, "/celltype_proportions_per_donor_20240328.pdf"), plot = p.celltype.fractions.per.donor, device = "pdf", dpi = 300, width = 6, height = 2.5, unit = "in")




####################################### Module boxplots #####################################

# BOXPLOTS Module Score
pathways = readRDS('/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pathways.Rds')

rna.obj.control$Object=factor(rna.obj.control$Object,levels = c('CONTROL','VEXAS'))

HALLMARK_APOPTOSIS = pathways[['HALLMARK_APOPTOSIS']]
HALLMARK_INFLAMMATORY_RESPONSE = pathways[['HALLMARK_INFLAMMATORY_RESPONSE']]
GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION=pathways[['GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION']]
REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR	=pathways[['REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR']]
REACTOME_TRANSLATION	=pathways[['REACTOME_TRANSLATION']]

rna.obj.control=AddModuleScore(rna.obj.control, assay = "RNA", features = list(HALLMARK_APOPTOSIS,HALLMARK_INFLAMMATORY_RESPONSE,REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR	,GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION, REACTOME_TRANSLATION), name=c('HALLMARK_APOPTOSIS','HALLMARK_INFLAMMATORY_RESPONSE','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION', 'REACTOME_TRANSLATION'),  seed=1234)

rna.obj.control_hsc=rna.obj.control[,rna.obj.control$AllCellType%in%c('HSC')]
rna.obj.control_emp=rna.obj.control[,rna.obj.control$AllCellType%in%c('EMP')]
rna.obj.control_lmpp=rna.obj.control[,rna.obj.control$AllCellType%in%c('LMPP')]
rna.obj.control_gmp=rna.obj.control[,rna.obj.control$AllCellType%in%c('GMP')]

met_hsc=rna.obj.control_hsc@meta.data[,c('AllCellType','Object','HALLMARK_APOPTOSIS1','HALLMARK_INFLAMMATORY_RESPONSE2','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR3','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION4', 'REACTOME_TRANSLATION5')]
met_emp=rna.obj.control_emp@meta.data[,c('AllCellType','Object','HALLMARK_APOPTOSIS1','HALLMARK_INFLAMMATORY_RESPONSE2','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR3','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION4', 'REACTOME_TRANSLATION5')]
met_lmpp=rna.obj.control_lmpp@meta.data[,c('AllCellType','Object','HALLMARK_APOPTOSIS1','HALLMARK_INFLAMMATORY_RESPONSE2','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR3','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION4', 'REACTOME_TRANSLATION5')]
met_gmp=rna.obj.control_gmp@meta.data[,c('AllCellType','Object','HALLMARK_APOPTOSIS1','HALLMARK_INFLAMMATORY_RESPONSE2','REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR3','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION4', 'REACTOME_TRANSLATION5')]

met=do.call("rbind", list(met_hsc, met_emp,met_lmpp,met_gmp))
met$AllCellType=factor(met$AllCellType, levels = c('HSC','EMP','LMPP','GMP'))

for (i in c('HSC','EMP','LMPP','GMP')){
  met_df=met[met$AllCellType==i,]
  
  ggplot(met_df, aes(x=Object, y=HALLMARK_APOPTOSIS1, group=Object,fill=Object)) + 
    geom_boxplot()+theme_classic()+ggtitle(paste0('HALLMARK_APOPTOSIs\n ',i))+ylab('Module Score')+xlab('')+
    scale_fill_manual(values=c("#E7D4E9", "#85378A"))+NoGrid()+
    NoLegend()+theme(plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(label = "p.format")
  ggsave(paste0('figures/current_figure_drafts/figure1_boxplot_apoptosis_',i,'.pdf'),dpi = 300, width = 2, height = 3, unit = "in")
  
}

for (i in c('HSC','EMP','LMPP','GMP')){
  met_df=met[met$AllCellType==i,]
  
  ggplot(met_df, aes(x=Object, y=HALLMARK_INFLAMMATORY_RESPONSE2, group=Object,fill=Object)) + 
    geom_boxplot()+theme_classic()+ggtitle(paste0('HALLMARK_INFLAMMATORY_RESPONSE\n ',i))+scale_y_continuous(expand = c(0.2, 0)) +
    ylab('Module Score')+xlab('')+scale_fill_manual(values=c("#FFCB88", "#FF7B0B"))+NoGrid()+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(label = "p.format")
  ggsave(paste0('figures/current_figure_drafts/figure1_boxplot_inflammatory_',i,'.pdf'),dpi = 300, width =2, height = 3, unit = "in")
  
}


for (i in c('HSC','EMP','LMPP','GMP')){
  met_df=met[met$AllCellType==i,]
  
  ggplot(met_df, aes(x=Object, y=REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR3, group=Object,fill=Object)) + 
    geom_boxplot()+theme_classic()+ggtitle(paste0('REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR\n',i))+
    ylab('Module Score')+xlab('')+scale_fill_manual(values=c("#FFB4C0", "#CD2135"))+
    NoGrid()+NoLegend()+theme(plot.title = element_text(hjust = 0.5)) + stat_compare_means(label = "p.format")
  ggsave(paste0('figures/current_figure_drafts/figure1_boxplot_UPR_',i,'.pdf'),dpi = 300, width = 2, height = 3, unit = "in")
  
}

for (i in c('HSC','EMP','LMPP','GMP')){
  met_df=met[met$AllCellType==i,]
  
  ggplot(met_df, aes(x=Object, y=GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION4, group=Object,fill=Object)) + 
    geom_boxplot()+theme_classic()+ggtitle(paste0('GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION\n',i))+ylab('Module Score')+
    xlab('')+scale_fill_manual(values=c("#A0CBE1", "#419CCC"))+NoGrid()+NoLegend()+theme(plot.title = element_text(hjust = 0.5)) + stat_compare_means(label = "p.format")
  ggsave(paste0('figures/current_figure_drafts/figure1_boxplot_MYELOID_',i,'.pdf'),dpi = 300, width =2, height = 3, unit = "in")
  
}

for (i in c('HSC','EMP','LMPP','GMP')){
  met_df=met[met$AllCellType==i,]
  
  ggplot(met_df, aes(x=Object, y=REACTOME_TRANSLATION5, group=Object,fill=Object)) + 
    geom_boxplot()+theme_classic()+ggtitle(paste0('REACTOME_TRANSLATION\n',i))+ylab('Module Score')+
    xlab('')+scale_fill_manual(values=c("#A0CBE1", "#419CCC"))+NoGrid()+NoLegend()+theme(plot.title = element_text(hjust = 0.5)) + stat_compare_means(label = "p.format")
  ggsave(paste0('figures/current_figure_drafts/figure1_boxplot_TRANSLATION_',i,'.pdf'),dpi = 300, width =2, height = 3, unit = "in")
  
}


## Making the myeloid bias density plot
GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION=pathways[['GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION']]
HOEBEKE_LYMPHOID_STEM_CELL_DN=pathways[['HOEBEKE_LYMPHOID_STEM_CELL_DN']]
HOEBEKE_LYMPHOID_STEM_CELL_UP=pathways[['HOEBEKE_LYMPHOID_STEM_CELL_UP']]

rna.obj.control=AddModuleScore(rna.obj.control, 
                   features = list(GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION,HOEBEKE_LYMPHOID_STEM_CELL_DN,HOEBEKE_LYMPHOID_STEM_CELL_UP), 
                   name=c('GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION','HOEBEKE_LYMPHOID_STEM_CELL_DN','HOEBEKE_LYMPHOID_STEM_CELL_UP'),  
                   assay = "RNA",
                   seed=1234)

hscs=rna.obj.control[,colnames(rna.obj.control)[rna.obj.control$AllCellType %in% c('HSC')]]

df=hscs@meta.data[,c('Object','GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION1','HOEBEKE_LYMPHOID_STEM_CELL_DN2','HOEBEKE_LYMPHOID_STEM_CELL_UP3')]
df$donors=df$Object
df$Object=ifelse(df$donors %in% c('CONTROL'),'Healthy','VEXAS')
df$Object=factor(df$Object, levels = c('Healthy','VEXAS'))
df1=data.frame('donors'=df$Object, 'diff.score'=df[,'GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION1'])
df=data.frame('donors'=df$Object, 'diff.score'=df[,'GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION1']-df[,'HOEBEKE_LYMPHOID_STEM_CELL_UP3'])


x_limits=c(min(df1$diff.score),max(df1$diff.score))
ggplot(data = df1, aes(x = diff.score, group = donors, fill = donors)) +
  geom_density(adjust = 1.5, alpha = 0.8)+
  scale_fill_manual(values = c(Healthy = "#F5B935", VEXAS = "#2B858C")) +
  theme_classic(base_size = 15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(x_limits)+xlab('')+ylab('Density')+ 
  theme(text = element_text(face = "bold")) +
  xlab("Myeloid Differentiation Score\n \n    (GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION)")
ggsave('figures/current_figure_drafts/figure1_density_scores_COMMON_MYE.pdf',width = 8, height=5)

ggplot() +
  theme_classic(base_size = 15) +
  theme(text = element_text(face = "bold")) + ylab('')+
  geom_boxplot(data = df1, aes(x = diff.score, y = donors, group = donors, fill = donors))+ scale_fill_manual(values = c(Healthy = "#F5B935", VEXAS = "#2B858C")) +
  xlab("Myeloid Differentiation Score\n \n    (GOBP_COMMON_MYELOID_PROGENITOR_CELL_PROLIFERATION)")
ggsave('figures/current_figure_drafts/figure1_density_scores_COMMON_MYE_boxplot_only.pdf',width = 8, height=2.5)

## Alternate colors: c('Healthy' = '#2B858C', 'VEXAS' = '#FF7B0B')

###################################### Run the differential expression testing #############################

DefaultAssay(rna.obj.control)='RNA'


# Remove the genes we are NOT testing
blacklist <- read_csv("/gpfs/commons/home/tbotella/VEXAS/RNA/MDS_splicing_project/blacklist_chry_chrx_mtgene_ribosome.csv")$blacklist
rna.obj.control=rna.obj.control[!rownames(rna.obj.control)%in%blacklist,]
rna.obj.control=subset(rna.obj.control, subset=AllCellType %in% c('HSC','LMPP','EMP','MkP', "GMP"))

## Renormalize after removing genes
rna.obj.control <- NormalizeData(rna.obj.control, normalization.method = 'LogNormalize', assay = "RNA")

celltypes=unique(rna.obj.control$AllCellType)

# Helper function to filter for minimum expression
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    # object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@assays$RNA@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@assays$RNA@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@assays$RNA@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      # object@assays$RNA@data <- object@assays$RNA@data[genes.use, ] # original way of subsetting
      
      ## modified way of subsetting, removing from entire obj
      counts <- GetAssayData(object, assay = "RNA")
      counts <- counts[(which(rownames(counts) %in% genes.use)),]
      object <- subset(object, features = rownames(counts))
      return(object)
    } else {
      return(object)
    }
  }

## Run the DE
## We considered all genes detected in at least 10% of the AML cells with a log2 fold change larger than 0.25 or smaller than −0.25 and a Bonferroni-adjusted P < 0.05 as differentially expressed.
#mclapply(X = celltypes,  FUN = function(celltype){
for (celltype in celltypes){
  print(celltype)
  
  prog.rna.obj.control = subset(rna.obj.control, subset = AllCellType == celltype)
  
  Idents(prog.rna.obj.control)=as.factor(prog.rna.obj.control$Donor)
  prog.rna.obj.control <- subset(prog.rna.obj.control, downsample=200)
  table(prog.rna.obj.control$Donor)
  
  ## Remove genes not detected in 10% of VEXAS cells
  prog.rna.obj.control.vexas <- subset(prog.rna.obj.control, Object == "VEXAS")
  prog.rna.obj.control.vexas <- FilterGenes(object = prog.rna.obj.control.vexas, min.value = 0, min.cells = floor(nrow(prog.rna.obj.control.vexas@meta.data)*0.1))
  features.test <- rownames(prog.rna.obj.control.vexas@assays$RNA@data)
  
  vexas.vs.control <- FindMarkers(prog.rna.obj.control, assay = "RNA", ident.1 = "VEXAS", ident.2 = "CONTROL", group.by = 'Object', test.use='wilcox', features = features.test, logfc.threshold = 0)
  vexas.vs.control=vexas.vs.control[order(vexas.vs.control$avg_log2FC, decreasing = T),]
  write.csv(vexas.vs.control,file = paste0('data/differential_expression/VEXAS_vs_control/vexas.vs.control.',celltype,'.DGE.csv'))
}  
#}, mc.cores = detectCores())


############################### Make the DE volcano plot #########################################

de.results <- read_csv("data/differential_expression/VEXAS_vs_control/vexas.vs.control.HSC.DGE.csv") %>% mutate(cluster = "HSC") %>% dplyr::rename(feature = "...1")

# pathways = readRDS('/gpfs/commons/home/tbotella/VEXAS/PAPER/Dec_2023/pathways.Rds')
# fgsea.out <- run_fgsea_for_list(de.results, pval_col = "p_val")
# volcano.gene.list <- list(
#   `UPR` = intersect(fgsea.out %>% filter(pathway == "ATF4 targets (Han et al.)") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#   `Inflammation` = union(
#     intersect(fgsea.out %>% filter(pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#     intersect(fgsea.out %>% filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature))
#     ),
#   `Apoptosis` = intersect(fgsea.out %>% filter(pathway == "HALLMARK_APOPTOSIS") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#   `Myeloid` = intersect(fgsea.out %>% filter(pathway == "GOBP_MYELOID_CELL_DIFFERENTIATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#   `Translation` = intersect(fgsea.out %>% filter(pathway == "REACTOME_TRANSLATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature))
# )
# 
# ## Pick colors + plot volcano
# volcano.gene.list.colors <- c(
#   Apoptosis = "dodgerblue",
#   `Inflammation` = "brown4",
#   `Myeloid` = "blue4",
#   `UPR` = "darkgreen", 
#   Translation = "orange"
# )

volcano.gene.list <- list(
  `highlight` = de.results %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>% pull(feature)
)

## Pick colors + plot volcano
volcano.gene.list.colors <- c(
  highlight = vexas.genotyping.palette[[2]]
)


## Plot volcano with labels
p.volcano <- plot_volcano(de.results, effect_line = 0.5, stat_line = 0.05, stat_column = "p_val_adj", title = paste0("HSCs"), max.overlaps = 25, only_genes = volcano.gene.list, only_genes_colors = volcano.gene.list.colors) + 
  theme(legend.position = "right") +
  ylab("-log10(adjusted p-value)") +
  xlab("Average log2(FC)")
ggsave(paste0(current.plots.path, "/RNA_control_comparison_HSC_diff_exp_gene_categories_20231218.pdf"), plot = p.volcano, dpi = 300, device = "pdf", width = 6.5, height = 5)
ggsave(paste0(current.plots.path, "/RNA_control_comparison_HSC_diff_exp_gene_categories_zoomed_20231218_zoomed.pdf"), plot = p.volcano, dpi = 300, device = "pdf", width = 10, height = 8)

p.volcano <- plot_volcano(de.results, effect_line = 0.5, stat_line = 0.05, stat_column = "p_val_adj", title = paste0("HSCs"), label_genes = F, only_genes = volcano.gene.list, only_genes_colors = volcano.gene.list.colors) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)")
ggsave(paste0(current.plots.path, "/RNA_control_comparison_HSC_diff_exp_20231209_no_labels.pdf"), plot = p.volcano, dpi = 300, device = "pdf", width = 6.5, height = 5)



############################ Genotype aware analysis with controls ###########################


DefaultAssay(rna.obj.control)='RNA'

## Update the "Donor" column for the healthy controls
rna.obj.control$Sample[is.na(rna.obj.control$Sample)] <- rna.obj.control$orig.ident[is.na(rna.obj.control$Sample)] 
rna.obj.control$Donor[is.na(rna.obj.control$Donor)] <- rna.obj.control$orig.ident[is.na(rna.obj.control$Donor)] 
healthy.sample.list <- unique(rna.obj.control$orig.ident) %>% grep("HD", ., value = T)
rna.obj.control$Donor <- plyr::mapvalues(rna.obj.control$Donor, from = healthy.sample.list, to = gsub("_.*", "", healthy.sample.list)) 
rna.obj.control@meta.data %>% count(Sample, Donor, orig.ident)


# Remove the genes we are NOT testing
blacklist <- read_csv("/gpfs/commons/home/tbotella/VEXAS/RNA/MDS_splicing_project/blacklist_chry_chrx_mtgene_ribosome.csv")$blacklist
rna.obj.control=rna.obj.control[!rownames(rna.obj.control)%in%blacklist,]
rna.obj.control=subset(rna.obj.control, subset=AllCellType %in% c('HSC','LMPP','EMP','MkP', "GMP"))

## Renormalize after removing genes
rna.obj.control <- NormalizeData(rna.obj.control, normalization.method = 'LogNormalize', assay = "RNA")


# Helper function to filter for minimum expression
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    # object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@assays$RNA@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@assays$RNA@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@assays$RNA@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      # object@assays$RNA@data <- object@assays$RNA@data[genes.use, ] # original way of subsetting
      
      ## modified way of subsetting, removing from entire obj
      counts <- GetAssayData(object, assay = "RNA")
      counts <- counts[(which(rownames(counts) %in% genes.use)),]
      object <- subset(object, features = rownames(counts))
      return(object)
    } else {
      return(object)
    }
  }


######################## WT vs control #######################

## Run the DE
prog.rna.obj.control = subset(rna.obj.control, subset = AllCellType == "HSC")
prog.rna.obj.control@meta.data %>% count(Donor, Genotype, Object) %>% 
  filter(Object == "CONTROL" | Genotype == "WT")
prog.rna.obj.control <- subset(prog.rna.obj.control, Object == "CONTROL" | Genotype == "WT")
prog.rna.obj.control$Status = paste0(prog.rna.obj.control$Object, "_", prog.rna.obj.control$Genotype)
prog.rna.obj.control@meta.data %>% 
  count(Donor, Status)

bc.keep <- prog.rna.obj.control@meta.data %>% 
  rownames_to_column("BC") %>% 
  group_by(Donor) %>% 
  filter(n() > 40) %>%
  slice_sample(n = 40, replace = F) %>% 
  pull(BC)

se <- subset(prog.rna.obj.control, cells = bc.keep)
se@meta.data %>% 
  count(Donor, Status, Genotype)

## Remove genes not detected in 10% of VEXAS cells
se.vexas <- subset(se, Object == "VEXAS")
se.vexas <- FilterGenes(object = se.vexas, min.value = 0, min.cells = floor(nrow(se.vexas@meta.data)*0.1))
features.test <- rownames(se.vexas@assays$RNA@data)
  
vexas.WT.vs.control <- FindMarkers(se, assay = "RNA", ident.1 = "VEXAS_WT", ident.2 = "CONTROL_NA", group.by = 'Status', test.use='wilcox', features = features.test, logfc.threshold = 0)
vexas.WT.vs.control=vexas.WT.vs.control[order(vexas.WT.vs.control$avg_log2FC, decreasing = T),]

vexas.WT.vs.control %>% 
  write_csv("data/differential_expression/VEXAS_WT_vs_control/WT_vs_control_40_cells_per_donor_DE.csv")

## Plot volcano with labels
p.volcano <- plot_volcano(vexas.WT.vs.control, effect_line = 0.5, stat_line = 0.05, stat_column = "p_val_adj", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)") +
  scale_color_manual(labels = c("VEXAS WT", "controls"), values = c("dodgerblue3", "black"))
p.volcano

fgsea.res <- run_fgsea_for_list(vexas.WT.vs.control, pval_col = "p_val")
fgsea.res %>% 
  write_csv("data/differential_expression/VEXAS_WT_vs_control/WT_vs_control_40_cells_per_donor_fgsea.csv")

######################## Looking at cGAS and STING #######################

rna.obj.control$category <- rna.obj.control$Genotype
rna.obj.control$category[rna.obj.control$Object == "CONTROL"] <- "Control"
table(rna.obj.control$category)
rna.obj.control$category <- factor(rna.obj.control$category, levels = c("Control", "WT", "MUT"))

## Look at HSCs
rna.HSC <- subset(rna.obj.control, AllCellType == "HSC" & category %in% c("Control", "MUT", "WT"))
VlnPlot(rna.HSC, features = "CGAS", group.by = "category", assay = "RNA", cols = c(Control = "grey80", MUT = "#BA242A", WT = "dodgerblue"), y.max = 3) + stat_compare_means(comparisons = list(c("WT", "MUT"))) + NoLegend() |
  VlnPlot(rna.HSC, features = "TMEM173", group.by = "category", assay = "RNA",  cols = c(Control = "grey80", MUT = "#BA242A", WT = "dodgerblue"), y.max = 3)  + stat_compare_means(comparisons = list(c("WT", "MUT"))) + NoLegend()
gene.list <- rna.HSC


## Look at monocytes
rna.mono <- subset(rna.obj.control, AllCellType == "CD14 Mono" & category %in% c("Control", "MUT", "WT"))
VlnPlot(rna.mono, features = "CGAS", group.by = "category", assay = "RNA", cols = c(Control = "grey80", MUT = "#BA242A", WT = "dodgerblue"), y.max = 3.25) + stat_compare_means(comparisons = list(c("WT", "MUT"))) + NoLegend() |
  VlnPlot(rna.mono, features = "TMEM173", group.by = "category", assay = "RNA",  cols = c(Control = "grey80", MUT = "#BA242A", WT = "dodgerblue"), y.max = 3.25) + stat_compare_means(comparisons = list(c("WT", "MUT"))) + NoLegend()

######################## Looking at MHC II signature #######################

signature.genes <- c("IL10", "LAG3", "CD274", "HAVCR2", "CD226", "MAF")
rna.obj.control <- AddModuleScore(rna.obj.control, features = list(signature.genes), name = c("MHCII"), assay = "RNA")

rna.cd4t <- subset(rna.obj.control, AllCellType == "CD4 T")
bc.keep <- rna.cd4t@meta.data %>% rownames_to_column("BC") %>% group_by(Object) %>% sample_n(size = 2500, replace = T) %>% pull(BC)
rna.cd4t <- subset(rna.cd4t, cells = bc.keep)
VlnPlot(rna.cd4t, features = "MHCII1", group.by = "Object", cols = c(Control = "grey", VEXAS = "tan")) + stat_compare_means() + NoLegend()
VlnPlot(rna.cd4t, features = c("IL10", "CD274"), group.by = "Object", cols = c(Control = "grey", VEXAS = "tan")) + NoLegend()


######################## MUT vs control #######################

## Run the DE
prog.rna.obj.control = subset(rna.obj.control, subset = AllCellType == "HSC")
prog.rna.obj.control@meta.data %>% count(Donor, Genotype, Object) %>% 
  filter(Object == "CONTROL" | Genotype == "MUT")
prog.rna.obj.control <- subset(prog.rna.obj.control, Object == "CONTROL" | Genotype == "MUT")
prog.rna.obj.control$Status = paste0(prog.rna.obj.control$Object, "_", prog.rna.obj.control$Genotype)
prog.rna.obj.control@meta.data %>% 
  count(Donor, Status)

bc.keep <- prog.rna.obj.control@meta.data %>% 
  rownames_to_column("BC") %>% 
  group_by(Donor) %>% 
  filter(n() > 200) %>%
  slice_sample(n = 200, replace = F) %>% 
  pull(BC)

se <- subset(prog.rna.obj.control, cells = bc.keep)
se@meta.data %>% 
  count(Donor, Status, Genotype)

## Remove genes not detected in 10% of VEXAS cells
se.vexas <- subset(se, Object == "VEXAS")
se.vexas <- FilterGenes(object = se.vexas, min.value = 0, min.cells = floor(nrow(se.vexas@meta.data)*0.1))
features.test <- rownames(se.vexas@assays$RNA@data)

vexas.MUT.vs.control <- FindMarkers(se, assay = "RNA", ident.1 = "VEXAS_MUT", ident.2 = "CONTROL_NA", group.by = 'Status', test.use='wilcox', features = features.test, logfc.threshold = 0)
vexas.MUT.vs.control=vexas.MUT.vs.control[order(vexas.MUT.vs.control$avg_log2FC, decreasing = T),]

vexas.MUT.vs.control %>% 
  write_csv("data/differential_expression/VEXAS_MUT_vs_control/MUT_vs_control_200_cells_per_donor_DE.csv")

## Plot volcano with labels
p.volcano <- plot_volcano(vexas.MUT.vs.control, effect_line = 0.5, stat_line = 0.05, stat_column = "p_val_adj", title = paste0("HSCs")) + 
  theme(legend.position = "right") +
  ylab("-log10(FDR)") +
  xlab("Average log2(FC)") +
  scale_color_manual(labels = c("VEXAS MUT", "controls"), values = c("darkred", "black"))
p.volcano

fgsea.res <- run_fgsea_for_list(vexas.MUT.vs.control, pval_col = "p_val")
fgsea.res %>% 
  write_csv("data/differential_expression/VEXAS_MUT_vs_control/MUT_vs_control_40_cells_per_donor_fgsea.csv")


##################################### Double-checking sex of donor ##########################

VlnPlot(rna.obj, features = c("UBA1"), assay = "RNA", group.by = "orig.ident") + NoLegend()
VlnPlot(rna.obj, features = c("ZFY"), assay = "RNA", group.by = "orig.ident") + NoLegend()
VlnPlot(rna.obj.control, features = c("ZFY"), assay = "RNA", group.by = "orig.ident") + NoLegend()
VlnPlot(rna.obj.control, features = c("UBA1"), assay = "RNA", group.by = "orig.ident") + NoLegend()
VlnPlot(rna.obj.control, features = c("BTK"), assay = "RNA", group.by = "orig.ident") + NoLegend()
