#' Run eSVD-DE on a Seurat object
#'
#' A thin wrapper around \code{eSVD2::eSVD_helper()} (eSVD2 >= 1.2.0), kept so
#' the analysis repo's call signature is stable. Everything the wrapper used
#' to do itself (drop individuals with too few cells, reject an underpowered
#' cohort with a warning and \code{NA}) now lives in \code{eSVD2}, which also
#' removes all-zero genes before the fit and reinserts them afterwards with
#' \code{NA} statistics and an FDR of 1, so a caller's
#' \code{sum(fdr_vec < 0.05)} needs no missingness handling.
#'
#' \code{eSVD2::eSVD()} refuses, with an error, what it cannot fit: an
#' all-zero gene, an individual with fewer than 3 cells, or an arm with fewer
#' than 2 individuals. Only \code{eSVD2::eSVD_helper()} turns the first two
#' into a filter and the last into a warning plus \code{NA}, which is why this
#' function no longer calls \code{eSVD()} directly.
#'
#' @param batch_var_prefix  \code{NULL}, or a prefix matching one of
#'                          \code{categorical_vars}.
#' @param case_control_levels  Character vector of length 2: the control level
#'                          first, then the case level.
#' @param case_control_var  Column of \code{seurat_obj@meta.data} holding the
#'                          case-control status; factor or character.
#' @param categorical_vars  Character vector of categorical covariate columns,
#'                          or \code{NULL}.
#' @param id_var            Column of \code{seurat_obj@meta.data} naming each
#'                          cell's individual.
#' @param numerical_vars    Character vector of numerical covariate columns, or
#'                          \code{NULL}.
#' @param seurat_obj        A \code{Seurat} object whose \code{counts} layer
#'                          holds raw counts.
#' @param bool_check_donors When \code{FALSE}, \code{min_cells},
#'                          \code{min_cells_casecontrol}, \code{min_cells_per_id}
#'                          and \code{min_ids} are set to 0, which disables
#'                          them, and the per-individual minimum inside
#'                          \code{eSVD2::eSVD()} is disabled too.
#'                          \code{min_ids_per_arm} is kept as supplied, because
#'                          \code{eSVD2::eSVD()} cannot fit an arm with fewer
#'                          than 2 individuals and would stop with an error;
#'                          keeping the filter turns that into a warning plus
#'                          \code{NA}. All-zero genes are still removed and
#'                          reinserted.
#' @param intermediate_save \code{NULL}, or a file path at which
#'                          \code{eSVD2::eSVD()} saves the object after each
#'                          stage.
#' @param min_cells_casecontrol  Return \code{NA} when either arm has this many
#'                          cells or fewer.
#' @param min_cells_per_id  Drop every individual with fewer than this many
#'                          cells; \code{eSVD2} requires \code{0} or at least 3.
#' @param min_cells         Return \code{NA} when this many cells or fewer
#'                          remain after the drop.
#' @param min_ids           Return \code{NA} when this many individuals or fewer
#'                          remain, pooled across arms.
#' @param min_ids_per_arm   Return \code{NA} when either arm has fewer than this
#'                          many individuals. 2 is the smallest value for which
#'                          the Welch degrees of freedom are defined.
#' @param verbose           Integer; \code{0} is silent.
#' @param ...               Further arguments to \code{eSVD2::eSVD()}, such as
#'                          \code{k}, \code{max_iter}, \code{bool_diet} or
#'                          \code{cap_multiplier}.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{\code{results}}{A \code{data.frame} with one row per gene of
#'   \code{seurat_obj} and columns \code{gene}, \code{logFC}, \code{se}
#'   (both \eqn{\log_2}, from \code{eSVD2::report_results()}), \code{pvalue}
#'   and \code{padj} (BH); the format shared by the four DE wrappers. All-zero
#'   genes have \code{NA} in every statistic, including \code{padj} (the
#'   \code{eSVD} object gives them an FDR of 1). When the cohort was rejected,
#'   every gene's statistics are \code{NA}.}
#'   \item{\code{original}}{Either \code{NA} (the cohort was rejected; the
#'   warning says which filter fired) or the \code{eSVD} object with the added
#'   element \code{gene_status}; see \code{?eSVD2::eSVD_helper}.}
#' }
#' @noRd
esvd_helper <- function(batch_var_prefix, # a variable inside categorical_vars. Can be NULL
                        case_control_levels, # Control and then Case
                        case_control_var,
                        categorical_vars,
                        id_var,
                        numerical_vars,
                        seurat_obj,
                        bool_check_donors = TRUE,
                        intermediate_save = NULL, # NULL or filepath to save intermediary results
                        min_cells_casecontrol = 20,
                        min_cells_per_id = 3,
                        min_cells = 20,
                        min_ids = 4,
                        min_ids_per_arm = 2,
                        verbose = 0,
                        ...){
  if(!requireNamespace("eSVD2", quietly = TRUE)){
    stop("Package 'eSVD2' is required. Install it with ",
         "remotes::install_github('linnykos/eSVD2').")
  }
  # `eSVD_helper()` and the all-zero-gene handling appeared in 1.2.0; an older
  # install has no exported `eSVD_helper` and its `eSVD()` has a different
  # contract, so refuse it up front with the version named.
  esvd2_version <- utils::packageVersion("eSVD2")
  if(esvd2_version < "1.2.0"){
    stop("Package 'eSVD2' >= 1.2.0 is required (installed: ", esvd2_version,
         "). Update it with remotes::install_github('linnykos/eSVD2').")
  }
  stopifnot(length(bool_check_donors) == 1, is.logical(bool_check_donors))

  # `eSVD2::filter_cohort()` disables a filter when its threshold is 0, and
  # `eSVD2::eSVD_helper()` forwards `min_cells_per_id` as `eSVD()`'s
  # `min_cells_per_individual`, so zeroing these four is the off switch.
  # `min_ids_per_arm` is deliberately not zeroed: `eSVD2::eSVD()` refuses an
  # arm with fewer than 2 individuals with an error (Welch df would be 0), so
  # zeroing it would only swap a warning-plus-NA for an error.
  if(!bool_check_donors){
    min_cells <- 0
    min_cells_casecontrol <- 0
    min_cells_per_id <- 0
    min_ids <- 0
  }

  esvd_res <- eSVD2::eSVD_helper(batch_var_prefix = batch_var_prefix,
                                 case_control_levels = case_control_levels,
                                 case_control_var = case_control_var,
                                 categorical_vars = categorical_vars,
                                 id_var = id_var,
                                 numerical_vars = numerical_vars,
                                 seurat_obj = seurat_obj,
                                 intermediate_save = intermediate_save,
                                 min_cells = min_cells,
                                 min_cells_casecontrol = min_cells_casecontrol,
                                 min_cells_per_id = min_cells_per_id,
                                 min_ids = min_ids,
                                 min_ids_per_arm = min_ids_per_arm,
                                 verbose = verbose,
                                 ...)

  # `results` is the format shared by the four DE wrappers; see
  # unify_de_results_claude.R. A rejected cohort still gets a row per gene, so
  # a caller can stack the four methods' tables without a special case.
  if(!inherits(esvd_res, "eSVD")){
    gene_vec <- rownames(seurat_obj)
    na_vec <- rep(NA_real_, length(gene_vec))
    results_df <- .unify_de_results(gene = gene_vec,
                                    logFC = na_vec,
                                    se = na_vec,
                                    pvalue = na_vec)
  } else {
    report_df <- eSVD2::report_results(esvd_res)
    # eSVD2 pads an all-zero gene with p = 1 so its own `fdr_vec` needs no NA
    # handling, but its BH counted only the analyzed genes. Counting the
    # padded genes here would inflate every other gene's `padj`, so they are
    # NA, and `padj` then equals eSVD2's `fdr_vec` on the analyzed genes.
    pvalue_vec <- report_df$pvalue
    pvalue_vec[esvd_res$gene_status[report_df$genes] != "analyzed"] <- NA_real_
    results_df <- .unify_de_results(gene = report_df$genes,
                                    logFC = report_df$logFC,
                                    se = report_df$logFC_se,
                                    pvalue = pvalue_vec)
  }

  list(results = results_df,
       original = esvd_res)
}
