#' Plot Signed Log P-values for Two Methods/Datasets
#'
#' Creates a scatter plot comparing signed log p-values between two differential 
#' expression methods/datasets, with optional highlighting of specific genes or pathways.
#' The signed log p-value is calculated as sign(logFC) * -log10(p-value), which
#' allows visualization of both significance and direction of change.
#'
#' @param df A data frame containing gene names and results from two methods.
#'   Must contain columns: Gene, method1_logFC, method1_pval, method2_logFC, method2_pval
#'   where method1 and method2 are replaced with the actual method names.
#' @param method1 Character string specifying the name of the first method.
#'   Used to construct column names (e.g., "DESeq2" looks for "DESeq2_logFC" and "DESeq2_pval").
#' @param method2 Character string specifying the name of the second method.
#'   Used to construct column names (e.g., "edgeR" looks for "edgeR_logFC" and "edgeR_pval").
#' @param bool_capped Logical indicating whether to cap extreme values. Currently unused
#'   but values are automatically capped at 97.5% quantile for visualization.
#' @param highlight_pathways Character vector of gene names to highlight in red.
#'   If NULL (default), genes significant in both methods (FDR < 0.05) will be highlighted.
#'
#' @return A ggplot2 object showing the comparison plot with:
#'   \itemize{
#'     \item Purple shaded regions indicating significance in both methods
#'     \item Dashed lines showing FDR cutoffs
#'     \item Coral line showing the first principal component
#'     \item Red points and labels for highlighted genes
#'     \item Correlation and Fisher test results in the title
#'   }
#'
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @importFrom stats p.adjust cor prcomp fisher.test
#' @importFrom grDevices rgb
#'
#' @export

was2code_plot_signed_logpval <- function(df, 
                                   method1, 
                                   method2,
                                   bool_capped = TRUE,
                                   highlight_pathways = NULL) {
  rownames(df) <- df$Gene
  
  # Remove genes with NA
  tmp <- df[, c("Gene", 
                paste0(method1, "_logFC"),
                paste0(method1, "_pval"),
                paste0(method2, "_logFC"),
                paste0(method2, "_pval"))]
  rm_genes <- which(apply(tmp, 1, function(row) any(is.na(row))))
  if (length(rm_genes) > 0) {
    df <- df[-rm_genes, , drop = FALSE]
  }
  
  # Signed log p-values
  x_val <- sign(df[, paste0(method1, "_logFC")]) * -log10(df[, paste0(method1, "_pval")])
  y_val <- sign(df[, paste0(method2, "_logFC")]) * -log10(df[, paste0(method2, "_pval")])
  names(x_val) <- names(y_val) <- df$Gene
  
  # FDR cutoffs
  x_adj <- p.adjust(df[, paste0(method1, "_pval")], method = "BH")
  y_adj <- p.adjust(df[, paste0(method2, "_pval")], method = "BH")
  names(x_adj) <- names(y_adj) <- df$Gene
  
  # figure out x_fdrcutoff
  if(length(which(x_adj <= 0.05)) >= 1){
    x_fdrcutoff <- max(x_val[which(x_adj > 0.05)])
  } else {
    x_fdrcutoff <- 2*max(abs(x_val))
  }
  x_genes <- names(x_adj)[which(x_adj <= 0.05)]
  
  # figure out y_fdrcutoff
  if(length(which(y_adj <= 0.05)) >= 1){
    y_fdrcutoff <- max(y_val[which(y_adj > 0.05)])
  } else {
    y_fdrcutoff <- 2*max(abs(y_val))
  }
  y_genes <- names(y_adj)[which(y_adj <= 0.05)]
  
  # Cap values
  idx <- intersect(which(!is.na(x_val)), which(!is.na(y_val)))
  x_val <- x_val[idx]; y_val <- y_val[idx]
  x_limit <- c(-1,1) * quantile(abs(x_val), 0.999)
  y_limit <- c(-1,1) * quantile(abs(y_val), 0.999)
  x_val <- pmin(pmax(x_val, 0.975 * x_limit[1]), 0.975 * x_limit[2])
  y_val <- pmin(pmax(y_val, 0.975 * y_limit[1]), 0.975 * y_limit[2])
  
  # PCA
  cor_val <- cor(x_val, y_val)
  pc1_slope <- prcomp(cbind(x_val, y_val), center = FALSE)$rotation[2,1] /
    prcomp(cbind(x_val, y_val), center = FALSE)$rotation[1,1]
  
  # Highlighting DE genes/pathways when highlight_pathways == NULL
  if (is.null(highlight_pathways)) {
    # If no specific pathways provided, highlight genes significant in both methods
    sig_both <- unique(c(x_genes, y_genes))
    genes_to_highlight <- sig_both
  } else {
    # Use provided pathways for highlighting
    genes_to_highlight <- highlight_pathways
  }
  
  ggplot_df <- data.frame(
    gene = names(x_val),
    method1 = x_val,
    method2 = y_val,
    highlight = ifelse(names(x_val) %in% genes_to_highlight, "highlight", "other"),
    label = ifelse(names(x_val) %in% genes_to_highlight, names(x_val), "")
  )
  
  # Fisher enrichment test
  if (!is.null(highlight_pathways)) {
    highlight_set <- intersect(highlight_pathways, ggplot_df$gene)
    sig_genes <- intersect(x_genes, y_genes)
    a <- sum(highlight_set %in% sig_genes)
    b <- sum(highlight_set %in% ggplot_df$gene) - a
    c <- sum(sig_genes %in% ggplot_df$gene) - a
    d <- nrow(ggplot_df) - a - b - c
    fisher_matrix <- matrix(c(a, b, c, d), nrow = 2)
    fisher_pval <- fisher.test(fisher_matrix)$p.value
  } else {
    fisher_pval <- NA
  }
  
  # Plot
  p <- ggplot2::ggplot(ggplot_df, ggplot2::aes(x = method1, y = method2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "gray") +
    ggplot2::geom_rect(data = ggplot_df[1, ],
                       ggplot2::aes(xmin = x_fdrcutoff, xmax = x_limit[2], ymin = y_fdrcutoff, ymax = y_limit[2]),
                       fill = rgb(135, 50, 255, maxColorValue = 255), alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = ggplot_df[1, ],
                       ggplot2::aes(xmin = x_limit[1], xmax = -x_fdrcutoff, ymin = y_limit[1], ymax = -y_fdrcutoff),
                       fill = rgb(135, 50, 255, maxColorValue = 255), alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_point(data = ggplot_df[ggplot_df$highlight == "other", ],
                        ggplot2::aes(color = highlight), alpha = 0.5) +
    ggplot2::geom_point(data = ggplot_df[ggplot_df$highlight == "highlight", ],
                        ggplot2::aes(color = highlight), alpha = 0.8) +
    ggrepel::geom_text_repel(data = ggplot_df[ggplot_df$highlight == "highlight", ],
                             ggplot2::aes(label = label),
                             color = "red", size = 2, show.legend = FALSE) +
    ggplot2::geom_abline(slope = pc1_slope, intercept = 0, color = "coral", linewidth = 1) +
    ggplot2::geom_hline(yintercept = c(-y_fdrcutoff, y_fdrcutoff), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-x_fdrcutoff, x_fdrcutoff), linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("highlight" = "red", "other" = "gray60")) +
    ggplot2::scale_x_continuous(limits = x_limit, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = y_limit, expand = c(0, 0)) +
    ggplot2::labs(
      x = paste0(method1, " (Signed -log10 p-value)"),
      y = paste0(method2, " (Signed -log10 p-value)")
    ) +
    ggplot2::ggtitle(paste0(
      method1, " vs. ", method2,
      " (Cor: ", round(cor_val, 2),
      if (!is.na(fisher_pval)) paste0(", Fisher P: ", signif(fisher_pval, 3)) else "",
      ")\n#genes: ", nrow(ggplot_df),
      ", #highlighted: ", sum(ggplot_df$highlight == "highlight")
    )) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      legend.position = "none"
    )
  
  return(p)
}