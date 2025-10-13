plot_volcano <- function(df,
                         logFC_variable,
                         pvalue_variable,
                         FCcutoff_quantile = 0.9,
                         padj_cutoff = 0.05,
                         xlim_quantile = c(0.005, 0.995)){
  stopifnot(length(rownames(df)) == nrow(df))
  
  pCutoff <- max(de_res$p_val[which(de_res$p_val_adj <= padj_cutoff)])
  FCcutoff <- stats::quantile(abs(de_res$avg_log2FC), probs = FCcutoff_quantile)
  xlim <- stats::quantile(de_res$avg_log2FC, probs = xlim_quantile)
  
  plot1 <- EnhancedVolcano::EnhancedVolcano(
    de_res,
    lab = rownames(de_res),
    x = logFC_variable,
    y = pvalue_variable,
    xlim = xlim,
    pCutoff = pCutoff,
    FCcutoff = FCcutoff
  )
  
  plot1
}
