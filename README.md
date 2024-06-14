# VEXAS-scRNA-seq
scRNA-seq analysis scripts to accompany the Landau Lab's recent preprint, "Single-cell genotype-phenotype mapping identifies therapeutic vulnerabilities in VEXAS syndrome".
Read the preprint here: https://doi.org/10.1101/2024.05.19.594376
Scripts from here and from https://github.com/RebeccaMurray/VEXAS-scATAC-seq can be used to reproduce the figures from the manuscript.

# Installation and setup
Source code can be downloaded using the following commands:
```
git clone https://github.com/RebeccaMurray/VEXAS-scRNA-seq.git
git clone https://github.com/RebeccaMurray/VEXAS-scATAC-seq.git
```
Raw sequencing data can be downloaded from the EGA (deposit info to come). Sequencing data can be processed with [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) and [Cell Ranger ATAC](https://www.10xgenomics.com/support/single-cell-atac) via the "count" commands. (For scRNA-seq data, the 5' mode was used). 
The cellranger outputs for both data sets can then be used as input to the source code above.

Downloading the source code should take seconds on a standard computer. The analysis scripts vary in time to run, with the most computationally intesive being the sample preprocessing and integration scripts [integrating_samples_SCTransform_soupx_pca_50.R](https://github.com/RebeccaMurray/VEXAS-scRNA-seq/blob/main/R/sample_integration/integrating_samples_SCTransform_soupx_pca_50.R) for the scRNA-seq data and [integrating_samples.R](https://github.com/RebeccaMurray/VEXAS-scATAC-seq/blob/main/R/sample_integration/integrating_samples.R) for the scATAC-seq paper. Both scripts took several hours using 150G of memory and 20 cores.

## Genotyping data
Genotyping data can be processed using the previously published [Genotyping of Transcriptomes](https://github.com/dan-landau/IronThrone-GoT) and [Genotyping of Targeted Loci with Chromatin Accessibility](https://github.com/landau-lab/Gotcha) pipelines.


# Software Dependencies
AnnotationDbi	1.60.0
AnnotationFilter	1.22.0
Azimuth	0.4.6.9004
base	4.2.2
Biobase	2.58.0
BiocGenerics	0.44.0
Biostrings	2.66.0
bonemarrowref.SeuratData	1.0.0
BSgenome	1.66.3
BSgenome.Hsapiens.UCSC.hg38	1.4.5
clusterProfiler	4.6.0
DoubletFinder	2.0.3
dsb	1.0.3
EnsDb.Hsapiens.v86	2.99.0
ensembldb	2.22.0
fgsea	1.24.0
future	1.29.0
GenomeInfoDb	1.34.9
GenomicFeatures	1.50.2
GenomicRanges	1.50.0
ggh4x	0.2.3
ggpubr	0.5.0
ggrepel	0.9.2
ggsci	2.9
gridExtra	2.3
IRanges	2.32.0
JASPAR2020	0.99.10
Matrix	1.5.3
msigdbr	7.5.1
org.Hs.eg.db	3.16.0
pals	1.7
patchwork	1.1.2
pheatmap	1.0.12
rtracklayer	1.58.0
S4Vectors	0.36.0
Seurat	4.9.9.9049
SeuratData	0.2.2
SeuratObject	4.9.9.9085
shinyBS	0.61.1
Signac	1.9.0.9000
sp	1.5.1
TFBSTools	1.36.0
tidyverse	1.3.2
tximport	1.26.0
viridis	0.6.2
viridisLite	0.4.1
XVector	0.38.0
SoupX	1.6.2, run with R 4.0.5
kallisto	0.48.0
bustools	0.41.0
kb_python	0.27.2
Amulet	1.1
ArchR	1.0.2, run with R 4.0.5
chromVAR	1.20.0
Salmon alevin	1.10.1
bowtie2	2.5.2
samtools	1.9
MarkDuplicates	2.27.4
macs2 callpeak	2.2.7.1
annotatr	1.24.0
cb_sniffer	None
Monocle3	v1.3.1

All software was run on a Linux x86_64 machine. No non-standard hardware is required.
Software has only been tested with the version numbers above.
