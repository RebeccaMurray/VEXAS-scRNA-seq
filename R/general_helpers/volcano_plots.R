library(ggrepel)
library(tidyverse)

plot_volcano <- function(de, 
                                metadata.df = NA,
                                stat_cutoff = NA, 
                                effect_line = 0.5, 
                                effect_cutoff = NA, 
                                stat_line = 0.25, extra_genes = c(), 
                                only_genes = list(),
                                only_genes_colors = list(),
                                title = "Volcano",
                                subtitle_addition = "",
                                genotyping_column = NA,
                                stat_column = "fdr",
                                label_genes = T,
                                max.overlaps = 15,
                                effect_column = "avg_log2FC"
                                ) {
  
  if (!("feature" %in% colnames(de))) {
    de$feature <- rownames(de)
  }
  
  if (is.na(effect_cutoff)) {
    effect_cutoff = effect_line
  }
  if (is.na(stat_cutoff)) {
    stat_cutoff = stat_line
  }
  
  de$statistic <- de[[stat_column]]
  de$effect <- de[[effect_column]]
  print(head(de))
  
  subtitle.str <- ""
  if (!is.na(genotyping_column)) {
    tb <- metadata.df %>% 
      dplyr::rename(genotyping_col = genotyping_column) %>% 
      pull(genotyping_col) %>% 
      table()
    print(tb)
    subtitle.str <- paste0("MUT: ", tb[["MUT"]], ", WT: ", tb[["WT"]], "\n", subtitle_addition)
  } else {
    subtitle.str <- subtitle_addition
  }
  
  de$diffexpressed <- NA_character_
  
  de$feature <- as.character(de$feature)
  
  
  if (length(only_genes)) {
    de$delabel <- NA
    # de$diffexpressed[abs(de$effect) > effect_cutoff & de$statistic < stat_cutoff] <- "Significant"
    for (gene_group in names(only_genes))  {
      current.list <- only_genes[[gene_group]]
      de$diffexpressed[de$feature %in% current.list] <- gene_group
    }
    # paste0(stat_column, " < ", stat_cutoff, " & absolute log2FC > ", effect_line)
    de$delabel[!is.na(de$diffexpressed)] <- de$feature[!is.na(de$diffexpressed)]
    color.vals <- only_genes_colors
    # de$diffexpressed <- factor(de$diffexpressed, levels = c(names(only_genes), "Significant"))
  } else {
    de$diffexpressed[de$effect > effect_cutoff & de$statistic < stat_cutoff] <- "MUT"
    de$diffexpressed[de$effect < -effect_cutoff & de$statistic < stat_cutoff] <- "WT"
    de$delabel <- NA_character_
    de$delabel[!is.na(de$diffexpressed)] <- de$feature[!is.na(de$diffexpressed)]
    color.vals <- c(WT="#333D84", MUT="#BA242A")
    if (!label_genes) {
      color.vals <- c(WT="black", MUT="black")
    }
  }
  
  print(table(de$diffexpressed, useNA = "always"))
  
  ## Split dataframe to order significant points
  de.not_significant <- de %>% dplyr::filter(is.na(diffexpressed) & abs(effect) <= effect_line | statistic >= stat_cutoff)
  de.significant <- de %>% dplyr::filter(is.na(diffexpressed) & abs(effect) > effect_line & statistic < stat_cutoff)
  de.highlighted <- de %>% dplyr::filter(!is.na(diffexpressed))
  
  p1 <- ggplot(de.highlighted, mapping = aes(x=effect, y=-log10(statistic), label = delabel, color = diffexpressed)) + 
    geom_point(data = de.not_significant, color = "grey90") +
    geom_point(data = de.significant, color = "black") +
    geom_point(alpha = 1) +
    geom_vline(xintercept=c(-effect_line, effect_line), col="black", linetype = "longdash") +
    geom_hline(yintercept=-log10(stat_line), col="black", linetype = "longdash") +
    theme_classic() +
    scale_color_manual(values=color.vals, na.value = "grey90") +
    labs(
      title = title, 
      subtitle = subtitle.str,
      x = effect_column,
      y = paste0("-log10(", stat_column, ")")
    ) +
    theme(legend.title = element_blank())
  
  if (label_genes)  {
    p1 <- p1 + geom_text_repel(mapping = aes(segment.linetype = 1, color = diffexpressed), max.overlaps = max.overlaps, fontface = "bold") ## Change to 2 for dashed
  }

  return(p1)
}



# make_volcano_plot_diff_accessibility <- function(df, lfc_thresh = 0.4, pval = 0.05) {
#   df$diffexpressed <- "NO"
#   df$diffexpressed[df$avg_log2FC > lfc_thresh & df$p_val < pval] <- "UP"
#   df$diffexpressed[df$avg_log2FC < -lfc_thresh & df$p_val < pval] <- "DOWN"
#   
#   df$label <- if_else(df$diffexpressed != "NO", df$gene_name, NA_character_)
#   
#   p1 <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), col = diffexpressed, label=label)) + 
#     geom_point() + 
#     theme_minimal() +
#     geom_text_repel(max.overlaps = 20) +
#     scale_color_manual(values=c("blue", "grey80", "red")) +
#     geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), col="red") +
#     geom_hline(yintercept=-log10(pval), col="red") +
#     labs(title  = "Differentially accessible regions")
#   
#   return(p1)
# }
