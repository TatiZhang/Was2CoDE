plot_volcano <- function(df,
                         logFC_variable,
                         pvalue_variable,
                         FCcutoff_quantile = 0.9,
                         padj_cutoff = 0.05,
                         xlim_quantile = c(0.005, 0.995)){
  stopifnot(length(rownames(df)) == nrow(df))
  
  pCutoff <- max(df$p_val[which(df$p_val_adj <= padj_cutoff)])
  FCcutoff <- stats::quantile(abs(df$avg_log2FC), probs = FCcutoff_quantile)
  xlim <- stats::quantile(df$avg_log2FC, probs = xlim_quantile)
  
  plot1 <- EnhancedVolcano::EnhancedVolcano(
    df,
    lab = rownames(df),
    x = logFC_variable,
    y = pvalue_variable,
    xlim = xlim,
    pCutoff = pCutoff,
    FCcutoff = FCcutoff
  )
  
  plot1
}
