#' Build the unified per-gene DE table shared by the four DE wrappers
#'
#' \code{deseq2_helper()}, \code{dreamlet_helper()}, \code{esvd_helper()} and
#' \code{nebula_helper()} each return \code{list(results = <this table>,
#' original = <the method's own object>)}, so a caller can compare methods
#' without knowing each one's column names.
#'
#' \code{padj} is \code{stats::p.adjust(pvalue, "BH")} for every method, over
#' the genes with a non-\code{NA} p-value. This deliberately differs from the
#' \code{padj} DESeq2 reports, which applies BH after its independent
#' filtering; DESeq2's own column is still in \code{original}.
#'
#' @param gene    Character vector of gene names.
#' @param logFC   Case-vs-control log fold change, on the \eqn{\log_2} scale.
#' @param se      Standard error of \code{logFC}, on the same scale.
#' @param pvalue  Two-sided p-value of the method's difference test.
#'
#' @returns A \code{data.frame} with one row per gene, row names \code{gene},
#' and columns \code{gene}, \code{logFC}, \code{se}, \code{pvalue} and
#' \code{padj}.
#' @noRd
.unify_de_results <- function(gene,
                              logFC,
                              se,
                              pvalue){
  stopifnot(length(gene) == length(logFC),
            length(gene) == length(se),
            length(gene) == length(pvalue),
            anyDuplicated(gene) == 0)

  # `p.adjust()` drops NA p-values before counting the tests, so a gene a
  # method declined to test does not inflate everyone else's adjustment
  results_df <- data.frame(gene = as.character(gene),
                           logFC = as.numeric(logFC),
                           se = as.numeric(se),
                           pvalue = as.numeric(pvalue),
                           padj = stats::p.adjust(as.numeric(pvalue),
                                                  method = "BH"),
                           stringsAsFactors = FALSE)
  rownames(results_df) <- results_df$gene

  results_df
}
