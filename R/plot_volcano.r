plot_volcano <- function(df,
                         logFC_variable,
                         pvalue_variable,
                         FCcutoff_quantile = 0.9,
                         padj_cutoff = 0.05,
                         xlim_quantile = c(0.005, 0.995)){
  stopifnot(length(rownames(df)) == nrow(df))
  
  padj_vec <- stats::p.adjust(df[,pvalue_variable], method = "BH")
  pCutoff <- max(df[which(padj_vec <= padj_cutoff), pvalue_variable])
  FCcutoff <- stats::quantile(abs(df[,logFC_variable]), probs = FCcutoff_quantile)
  xlim <- stats::quantile(df[,logFC_variable], probs = xlim_quantile)
  
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
