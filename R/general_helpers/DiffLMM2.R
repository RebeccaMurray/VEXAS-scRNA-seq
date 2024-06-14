library(parallel)
library(dplyr)
# library(rslurm)
library(readr)
library(lme4)
library(stats)

DiffLMM2 = function (metadata, provided.matrix, subset = NULL, sample.column,
                     cluster.column, selected.clusters, treatment.column, treatment.levels,
                     reg = 1e-3,
                     ncores = 18)
{
  if (class(metadata) != "data.frame") {
    stop("Specified metadata is not data.frame class")
  }
  if (!sample.column %in% colnames(metadata)) {
    stop("Specified sample column not found in @cellColData slot")
  }
  if (!cluster.column %in% colnames(metadata)) {
    stop("Specified pseudotime column not found in @cellColData slot")
  }
  if (!(treatment.column %in% colnames(metadata))) {
    stop("Specified treatment column not found in @cellColData slot.")
  }
  if (!sum(treatment.levels %in% unique(metadata[, treatment.column])) ==
      2) {
    stop("Specified treatment levels not found in specified treatment column in @cellColData slot")
  }
  if (sum(selected.clusters %in% unique(metadata[, cluster.column])) !=
      length(selected.clusters)) {
    message(paste0(setdiff(selected.clusters, unique(metadata[,
                                                              cluster.column])), " not found in specified cluster column. Proceeding with: ",
                   intersect(selected.clusters, unique(metadata[, cluster.column]))))
    if (!sum(selected.clusters %in% unique(metadata[, cluster.column])) ==
        0)
      stop("None of the selected clusters found in specified cluster.column")
  }
  m <- provided.matrix
  if (!is.null(subset)) {
    if (nrow(m[rownames(m) %in% subset, ]) == 0) {
      stop("None of the specified features in subset is found in matrix")
    }
    else {
      m = m[rownames(m) %in% subset, ]
    }
  }
  message("------- Running differential analysis based on linear mixture models (LMM) -------")
  out = mclapply(rownames(m), function(x) {
    dev = as.vector(t(m[x, ]))
    names(dev) = colnames(m)
    temp = metadata
    temp = temp[, c(sample.column, cluster.column, treatment.column)]
    temp[, x] = dev[rownames(temp)]
    temp = temp[temp[, treatment.column] %in% treatment.levels,
    ]
    temp = temp[temp[, cluster.column] %in% selected.clusters,
    ]
    temp = temp[complete.cases(temp), ]
    temp = data.frame(motif = temp[, x], treatment = temp[,
                                                          treatment.column], sample = temp[, sample.column])
    lmm1 = suppressMessages(lmer(data = temp, formula = motif ~
                                   treatment + (1 | sample), REML = F))
    lmm2 = suppressMessages(lmer(data = temp, formula = motif ~
                                   (1 | sample), REML = F))
    out = anova(lmm1, lmm2)
    pval = out$`Pr(>Chisq)`[2]
    delta = mean(temp$motif[temp$treatment == treatment.levels[2]],
                 na.rm = T) - mean(temp$motif[temp$treatment == treatment.levels[1]],
                                   na.rm = T)
    fc = (mean(temp$motif[temp$treatment == treatment.levels[2]],
               na.rm = T)+reg)/(mean(temp$motif[temp$treatment == treatment.levels[1]],
                                     na.rm = T)+reg)
    final = c(pval = pval, delta = delta, fc = fc)
    return(final)
  }, mc.cores = ncores)
  names(out) = rownames(m)
  out = as.data.frame(do.call(rbind, out))
  out$fdr = p.adjust(out$pval, method = "fdr")
  out$feature = rownames(out)
  message("------- Done! -------")
  return(out)
}

## Example to run
#
# # Clean up the matrix
# mtx = integrated@assays$GeneScores@data[,!is.na(integrated@assays$GeneScores@data[1,])]
# mtx = mtx[rowSums(mtx) != 0,]
# mtx = mtx[rownames(mtx) %in% expressed_genes,]
# geneAccess_untreated = list()
# selected.clusters = list(HSC_mapped = "HSC")
# message("Running DiffLMM...")
# for(i in names(selected.clusters)){
#   message("Untreated: ",i)
#   cond = Cells(integrated)[integrated$predicted.l2 %in% selected.clusters[[i]] & integrated$Treatment == "Untreated"]
#   # cond = cond[!is.na(cond)]
#   meta = integrated@meta.data[cond,]
#   geneAccess_untreated[[i]] = DiffLMM2(metadata = meta,
#                                        provided.matrix = mtx,
#                                        sample.column = "Sample",
#                                        cluster.column = "predicted.l2",
#                                        selected.clusters = selected.clusters[[i]],
#                                        treatment.column = "JAK2_Genotype",
#                                        treatment.levels = c("WT","MUT"),
#                                        ncores = 18,
#                                        subset = NULL,
#                                        reg = 1e-3)
# }