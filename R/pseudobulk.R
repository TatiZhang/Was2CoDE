# Donor-level pseudobulk construction and DESeq2 fitting, shared by the
# DE-method wrappers.
#
# Every wrapper in this package that runs a bulk method on single-cell data
# needs the same four steps: aggregate nuclei to one count vector per donor,
# carry the donor-level numerical covariates across the aggregation, coerce the
# design variables into the factor/scale form the model expects, and drop
# design variables that do not vary.  `build_pseudobulk()` is that pipeline;
# `deseq2_helper()` and `deseq2_tost_helper()` both call it.
#
# The second half of the file is the DESeq2 fit-and-extract pair used by
# `deseq2_tost_helper()`: `.fit_deseq2_mle()` and `.extract_lfc_se()`.  These
# are kept apart from the TOST arithmetic so that file stays about equivalence
# testing, and so any future wrapper needing unshrunken (log2FC, SE) per gene
# can reuse them.  Note `deseq2_helper()` deliberately does NOT use them: it
# calls `DESeq2::DESeq(dds)` with the package defaults, and routing it through
# `.fit_deseq2_mle()` would change its console output and couple it to
# TOST-specific choices.

#' Validate the donor / covariate arguments shared by the DE-method wrappers
#'
#' @param case_control_levels Length-2 character vector, control level first.
#' @param case_control_var Name of the case/control metadata column.
#' @param categorical_vars Character vector of categorical covariates, or `NULL`.
#' @param id_var Name of the donor id metadata column.
#' @param numerical_vars Character vector of numerical covariates, or `NULL`.
#'
#' @return Invisibly `TRUE`; called for its side effect of erroring.
#' @noRd
.check_donor_vars <- function(case_control_levels,
                              case_control_var,
                              categorical_vars,
                              id_var,
                              numerical_vars) {
  stopifnot(all(is.null(categorical_vars)) || (length(unique(categorical_vars)) == length(categorical_vars) && all(is.character(categorical_vars))))
  stopifnot(all(is.null(numerical_vars)) || (length(unique(numerical_vars)) == length(numerical_vars) && all(is.character(numerical_vars))))
  stopifnot(length(case_control_var) == 1,
            is.character(case_control_var),
            length(id_var) == 1,
            is.character(id_var),
            length(case_control_levels) == 2,
            all(is.character(case_control_levels)))

  invisible(TRUE)
}

#' Aggregate a Seurat object to one pseudobulk count vector per donor
#'
#' Summing nuclei per donor before testing is what makes the downstream test
#' donor-level rather than cell-level, and is what avoids pseudoreplication.
#'
#' @param seurat_obj A `Seurat` object with cells from many donors.
#' @param case_control_var,categorical_vars,id_var Metadata column names used
#'   as grouping variables.
#'
#' @return A `Seurat` object with one "cell" per donor x covariate combination,
#'   carrying the original object's variable features.
#' @noRd
.pseudobulk_from_seurat <- function(seurat_obj,
                                    case_control_var,
                                    categorical_vars,
                                    id_var) {
  pseudo_seurat <- Seurat::AggregateExpression(seurat_obj,
                                               assays = "RNA",
                                               return.seurat = TRUE,
                                               group.by = c(id_var, case_control_var, categorical_vars))
  Seurat::VariableFeatures(pseudo_seurat) <- Seurat::VariableFeatures(seurat_obj)

  pseudo_seurat
}

#' Copy donor-constant numerical covariates onto the pseudobulk object
#'
#' Aggregation drops numerical metadata, so it has to be re-attached by hand.
#' Each numerical variable must be constant within a donor -- otherwise it is
#' not a donor-level covariate and cannot survive the aggregation.
#'
#' @param pseudo_seurat Pseudobulk `Seurat` object.
#' @param seurat_obj The original cell-level `Seurat` object.
#' @param numerical_vars Character vector of metadata column names, or `NULL`.
#' @param id_var Donor id column name.
#'
#' @return `pseudo_seurat` with the numerical variables added to `meta.data`.
#' @noRd
.attach_donor_level_numerics <- function(pseudo_seurat,
                                         seurat_obj,
                                         numerical_vars,
                                         id_var) {
  for (variable in numerical_vars) {
    tmp <- rep(NA, length(Seurat::Cells(pseudo_seurat)))
    names(tmp) <- pseudo_seurat@meta.data[, id_var]

    for (person in unique(names(tmp))) {
      idx <- which(seurat_obj@meta.data[, id_var] == person)
      person_idx <- which(names(tmp) == person)
      zz <- seurat_obj@meta.data[idx, variable]
      if (any(abs(diff(range(zz))) >= 1e-4)) {
        stop(paste0("Person (", person, ") violates the variable (", variable, "): ", diff(range(zz))))
      }
      tmp[person_idx] <- mean(zz)
    }

    pseudo_seurat@meta.data[, variable] <- tmp
  }

  pseudo_seurat
}

#' Build the model's colData from the pseudobulk metadata
#'
#' Relevels the case/control factor so the control level is the reference,
#' orders the other factors by decreasing frequency, standardizes the numerical
#' covariates, and drops any factor with no variation across donors (which
#' would otherwise make the design matrix rank deficient).  Note `id_var` is
#' releveled but deliberately not kept: donor is the unit of observation, not a
#' design variable.
#'
#' @param metadata_pseudobulk `data.frame` of pseudobulk metadata.
#' @param case_control_levels Length-2 character vector, control first.
#' @param case_control_var,categorical_vars,id_var,numerical_vars Column names.
#'
#' @return A `data.frame` holding only the design variables that vary.
#' @noRd
.prepare_pseudobulk_metadata <- function(metadata_pseudobulk,
                                         case_control_levels,
                                         case_control_var,
                                         categorical_vars,
                                         id_var,
                                         numerical_vars) {
  metadata_pseudobulk[, case_control_var] <- stats::relevel(factor(metadata_pseudobulk[, case_control_var]),
                                                            ref = case_control_levels[1])

  for (variable in c(categorical_vars, id_var)) {
    tab_vec <- table(metadata_pseudobulk[, variable])
    metadata_pseudobulk[, variable] <- factor(metadata_pseudobulk[, variable],
                                              levels = names(tab_vec)[order(tab_vec, decreasing = TRUE)])
  }

  for (variable in numerical_vars) {
    metadata_pseudobulk[, variable] <- scale(as.numeric(metadata_pseudobulk[, variable]))
  }

  # make sure there's variation among the donors
  metadata_pseudobulk <- metadata_pseudobulk[, c(categorical_vars, numerical_vars, case_control_var), drop = FALSE]
  keep_vars <- c()
  for (j in 1:ncol(metadata_pseudobulk)) {
    if (!is.factor(metadata_pseudobulk[, j])) {
      keep_vars <- c(keep_vars, j)
    } else if (length(unique(metadata_pseudobulk[, j])) > 1) {
      keep_vars <- c(keep_vars, j)
    }
  }

  metadata_pseudobulk[, keep_vars, drop = FALSE]
}

#' Aggregate a Seurat object into a donor-level pseudobulk count matrix and design
#'
#' The shared front half of every pseudobulk DE wrapper in this package: it
#' turns a cell-level `Seurat` object into the `(countData, colData)` pair a
#' bulk model needs, with covariates carried across the aggregation.
#'
#' Only the variable features of `seurat_obj` are returned, matching the
#' behaviour the DE-method wrappers have always had.
#'
#' @param seurat_obj A `Seurat` object.
#' @param case_control_levels Length-2 character vector: control level first,
#'   then case.
#' @param case_control_var Name of the case/control metadata column.
#' @param categorical_vars Character vector of categorical covariates, or `NULL`.
#' @param id_var Name of the donor id metadata column.
#' @param numerical_vars Character vector of donor-constant numerical
#'   covariates, or `NULL`.
#'
#' @return A list with
#'   \describe{
#'     \item{`mat`}{genes x pseudobulk-sample count matrix}
#'     \item{`metadata`}{`data.frame` of design variables, rows aligned to the
#'       columns of `mat`}
#'     \item{`pseudo_seurat`}{the intermediate aggregated `Seurat` object, kept
#'       so callers can reach metadata that was dropped from the design}
#'   }
#' @noRd
build_pseudobulk <- function(seurat_obj,
                             case_control_levels,
                             case_control_var,
                             categorical_vars,
                             id_var,
                             numerical_vars) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required. Install it with install.packages('Seurat').")
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("Package 'SeuratObject' is required. Install it with install.packages('SeuratObject').")
  .check_donor_vars(case_control_levels = case_control_levels,
                    case_control_var = case_control_var,
                    categorical_vars = categorical_vars,
                    id_var = id_var,
                    numerical_vars = numerical_vars)

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = seurat_obj,
                                           case_control_var = case_control_var,
                                           categorical_vars = categorical_vars,
                                           id_var = id_var)
  pseudo_seurat <- .attach_donor_level_numerics(pseudo_seurat = pseudo_seurat,
                                                seurat_obj = seurat_obj,
                                                numerical_vars = numerical_vars,
                                                id_var = id_var)

  mat_pseudobulk <- SeuratObject::LayerData(pseudo_seurat,
                                            layer = "counts",
                                            assay = "RNA",
                                            features = Seurat::VariableFeatures(pseudo_seurat))
  metadata_pseudobulk <- .prepare_pseudobulk_metadata(metadata_pseudobulk = pseudo_seurat@meta.data,
                                                      case_control_levels = case_control_levels,
                                                      case_control_var = case_control_var,
                                                      categorical_vars = categorical_vars,
                                                      id_var = id_var,
                                                      numerical_vars = numerical_vars)

  # Bulk model fitters match countData columns to colData rows by position, so
  # a mismatch here would silently misassign donors rather than error.
  stopifnot(ncol(mat_pseudobulk) == nrow(metadata_pseudobulk),
            all(colnames(mat_pseudobulk) == rownames(metadata_pseudobulk)))

  list(mat = mat_pseudobulk,
       metadata = metadata_pseudobulk,
       pseudo_seurat = pseudo_seurat)
}

# ---------------------------------------------------------------------------
# DESeq2 fitting on the pseudobulk matrix
# ---------------------------------------------------------------------------

#' Fit the DESeq2 negative-binomial GLM with unshrunken coefficients
#'
#' `betaPrior = FALSE` is not optional here: a zero-centered prior on the log2
#' fold change shrinks estimates toward zero, which inflates the count of genes
#' declared "stably expressed".  DESeq2 itself refuses `altHypothesis =
#' "lessAbs"` when `betaPrior = TRUE` for the same reason.
#'
#' The design is `~ <covariates> + <case/control>`, built from the column names
#' of `metadata_pseudobulk` in order.  `.prepare_pseudobulk_metadata()` puts the
#' case/control variable last, which is what makes its coefficient reachable as
#' `"<var>_<case>_vs_<control>"`.
#'
#' @param mat_pseudobulk Integer count matrix, genes x pseudobulk samples.
#' @param metadata_pseudobulk `data.frame` of design variables.
#' @param use_t Passed to [DESeq2::DESeq()]; use a t reference distribution.
#' @param quiet Passed to [DESeq2::DESeq()]; suppress its progress messages.
#'
#' @return A fitted `DESeqDataSet`.
#' @noRd
.fit_deseq2_mle <- function(mat_pseudobulk,
                            metadata_pseudobulk,
                            use_t = FALSE,
                            quiet = FALSE) {
  design_formula <- stats::as.formula(paste0("~ ", paste0(colnames(metadata_pseudobulk), collapse = "+")))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                        colData = metadata_pseudobulk,
                                        design = design_formula)
  DESeq2::DESeq(dds, betaPrior = FALSE, useT = use_t, quiet = quiet)
}

#' Pull the unshrunken log2 fold change and its standard error out of a fit
#'
#' Returns the effect and its uncertainty only -- no p-values.  Whatever test a
#' caller wants (a Wald difference test, a TOST equivalence test, a
#' meta-analytic pool) is built from `(lfc, se)` downstream; for DESeq2's own
#' differential-expression p-values, run [deseq2_helper()].
#'
#' @param dds A fitted `DESeqDataSet`.
#' @param coef_name Name of the coefficient, as in `DESeq2::resultsNames()`.
#' @param cooks_filter Logical; if `TRUE`, genes flagged by DESeq2's Cook's
#'   distance cutoff have their `lfc` and `se` set to `NA`.  An outlier-driven
#'   fit should not be allowed to support an inference either way.
#' @param base_mean_min Genes with `baseMean` below this are set to `NA`.
#'
#' @return A `data.frame` with `base_mean`, `lfc`, `se`, `df` (or `NA` if a
#'   normal reference is in use), and the logical flag `cooks_outlier`.
#' @noRd
.extract_lfc_se <- function(dds,
                            coef_name,
                            cooks_filter = TRUE,
                            base_mean_min = 0) {
  # independentFiltering is switched off because it only ever affects padj,
  # which is not returned here -- and because DESeq2's filter is tuned to
  # maximize rejections of the *difference* null, so it has no business
  # deciding which genes are eligible for some other test.  Filtering is
  # instead an explicit, pre-specifiable baseMean threshold.
  res <- DESeq2::results(dds,
                         name = coef_name,
                         independentFiltering = FALSE,
                         cooksCutoff = cooks_filter)

  lfc <- res$log2FoldChange
  se <- res$lfcSE

  # DESeq2 signals a Cook's-flagged gene by NA-ing its p-value while leaving
  # lfc/se intact, so the flag has to be read off the p-value and propagated
  # by hand.
  cooks_outlier <- is.na(res$pvalue) & !is.na(lfc) & !is.na(se)
  drop_idx <- cooks_outlier | is.na(lfc) | is.na(se) | res$baseMean < base_mean_min
  lfc[drop_idx] <- NA
  se[drop_idx] <- NA

  df_vec <- SummarizedExperiment::mcols(dds)$tDegreesFreedom
  if (is.null(df_vec)) df_vec <- NA_real_

  data.frame(base_mean = res$baseMean,
             lfc = lfc,
             se = se,
             df = df_vec,
             cooks_outlier = cooks_outlier,
             row.names = rownames(res),
             stringsAsFactors = FALSE)
}

