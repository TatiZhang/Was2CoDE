plot_volcano <- function(df,
                         logFC_variable,
                         pvalue_variable,
                         FCcutoff_quantile = 0.9,
                         padj_cutoff = 0.05,
                         xlim_quantile = c(0.005, 0.995),
                         ymax_quantile = 0.995){
  stopifnot(length(rownames(df)) == nrow(df))
  
  # replace NA
  idx <- which(is.na(df[,pvalue_variable]))
  if(length(idx) > 0) df[idx,pvalue_variable] <- 1
  
  padj_vec <- stats::p.adjust(df[,pvalue_variable], 
                              method = "BH")
  gene_idx <- which(padj_vec <= padj_cutoff)
  if(length(gene_idx) > 0){
    pCutoff <- max(df[gene_idx, pvalue_variable])
  } else {
    pCutoff <- 0
  }

  FCcutoff <- stats::quantile(abs(df[,logFC_variable]), 
                              probs = FCcutoff_quantile,
                              na.rm = TRUE)
  xlim <- stats::quantile(df[,logFC_variable], 
                          probs = xlim_quantile,
                          na.rm = TRUE)
  
  ymax <- stats::quantile(-log10(df[,pvalue_variable]), 
                          probs = ymax_quantile)
  df[,pvalue_variable] <- pmax(df[,pvalue_variable], 10^(-ymax))
  ylim <- c(0, ymax)
  
  plot1 <- EnhancedVolcano::EnhancedVolcano(
    df,
    lab = rownames(df),
    x = logFC_variable,
    y = pvalue_variable,
    xlim = xlim,
    ylim = ylim,
    pCutoff = pCutoff,
    FCcutoff = FCcutoff
  )
  
  plot1
}
