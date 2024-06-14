library(parallel)
library(dplyr)
# library(rslurm)
library(readr)
library(lme4)
library(stats)

DiffLMMPseudocount = function(metadata,
                   provided.matrix,
                   subset = NULL,
                   sample.column,
                   cluster.column,
                   selected.clusters,
                   treatment.column,
                   treatment.levels,
                   pseudocount.value = 1,
                   ncores = 1){
  
  if(class(metadata) != "data.frame"){
    stop("Specified metadata is not data.frame class")
  }
  
  if(!sample.column %in% colnames(metadata)){
    stop("Specified sample column not found in @cellColData slot")
  }
  
  if(!cluster.column %in% colnames(metadata)){
    stop("Specified pseudotime column not found in @cellColData slot")
  }
  
  if(!(treatment.column %in% colnames(metadata))){
    stop("Specified treatment column not found in @cellColData slot.")
  }
  
  if(!sum(treatment.levels %in% unique(metadata[,treatment.column])) == 2){
    stop("Specified treatment levels not found in specified treatment column in @cellColData slot")
  }
  
  if(sum(selected.clusters %in% unique(metadata[,cluster.column])) != length(selected.clusters)){
    message(paste0(setdiff(selected.clusters, unique(metadata[,cluster.column])), " not found in specified cluster column. Proceeding with: ",intersect(selected.clusters, unique(metadata[,cluster.column]))))
    if(!sum(selected.clusters %in% unique(metadata[,cluster.column])) == 0)
      stop("None of the selected clusters found in specified cluster.column")
  }
  
  
  m <- provided.matrix
  
  
  if(!is.null(subset)){
    if(nrow(m[rownames(m) %in% subset,]) == 0){
      stop("None of the specified features in subset is found in matrix")
    }else{
      m = m[rownames(m) %in% subset,]
    }
  }
  
  message("------- Running differential analysis based on linear mixture models (LMM) -------")
  out = mclapply(rownames(m), function(x){
    dev = as.vector(t(m[x,]))
    names(dev) = colnames(m)
    temp = metadata
    temp = temp[,c(sample.column, cluster.column, treatment.column)]
    temp[,x] = dev[rownames(temp)]
    temp = temp[temp[,treatment.column] %in% treatment.levels,]
    temp = temp[temp[,cluster.column] %in% selected.clusters,]
    temp = temp[complete.cases(temp),]
    
    temp = data.frame(motif = temp[,x], treatment = temp[,treatment.column], sample = temp[,sample.column])
    
    lmm1 = suppressMessages(lmer(data = temp, formula = motif ~ treatment + (1|sample), REML = F))
    lmm2 = suppressMessages(lmer(data = temp, formula = motif ~ (1|sample), REML = F))
    out = anova(lmm1,lmm2)
    
    pval = out$`Pr(>Chisq)`[2]
    delta = mean(temp$motif[temp$treatment == treatment.levels[2]], na.rm = T) - mean(temp$motif[temp$treatment == treatment.levels[1]], na.rm = T)
    
    ######## Changed by RM to add expm1, pseudocount ##############
    fc = (mean(expm1(temp$motif[temp$treatment == treatment.levels[2]]), na.rm = T) + pseudocount.value) / (mean(expm1(temp$motif[temp$treatment == treatment.levels[1]]), na.rm = T) + pseudocount.value)
    ###############################################################
    
    final = c(pval = pval, delta = delta, fc = fc)
    return(final)
  }, mc.cores = ncores)
  
  names(out) = rownames(m)
  out = as.data.frame(do.call(rbind,out))
  out$fdr = p.adjust(out$pval, method = "fdr")
  out$feature = rownames(out)
  message("------- Done! -------")
  
  return(out)
}