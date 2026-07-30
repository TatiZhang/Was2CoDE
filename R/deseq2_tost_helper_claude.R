# TOST (two one-sided tests) equivalence testing on top of DESeq2.
#
# Design notes, from the project wiki (../Was2CoDE_wiki): see the pages
# [[analysis-tost-for-deg-equivalence]] (Route 1), [[method-tost]],
# [[evolution-equivalence-ci-correspondence]] and [[sesoi]].
#
#   * TOST is a *wrapper*, not a method.  It needs only an estimate, a standard
#     error and a reference distribution.  Here the estimate is DESeq2's
#     pseudobulk log2 fold change and the standard error is its Wald `lfcSE`.
#   * Each of the two one-sided tests is run at the FULL alpha.  TOST is an
#     intersection-union test (Berger & Hsu 1996, Thm 1), so there is no
#     alpha-splitting between them, and the gene-level p-value is
#     max(p_lower, p_upper) (Schuirmann 1987).
#   * The test is built from two explicit one-sided tests, NOT from a
#     "is the 1 - 2*alpha CI inside the bounds?" check.  That recipe is exact
#     only for equal-tailed t-intervals (Berger & Hsu 1996, sec. 5.2).
#   * Fold changes must be UNSHRUNKEN.  `lfcShrink()` / apeglm / ashr pull
#     log2FC toward zero, which is anti-conservative for an equivalence test:
#     it manufactures "stably expressed" genes.  `betaPrior = FALSE` is
#     enforced below and `lfcShrink` is deliberately not called.
#   * The equivalence bound (SESOI) must be pre-specified, so `lfc_threshold`
#     has no default.
#
# The arithmetic here reproduces DESeq2's own `altHypothesis = "lessAbs"`
# branch of `DESeq2::results()`; it is factored out so that (a) it can be unit
# tested without Bioconductor, and (b) the difference test and the equivalence
# test are guaranteed to come from the same (lfc, se) pair.

# ---------------------------------------------------------------------------
# Core TOST arithmetic (no DESeq2 / Seurat dependency)
# ---------------------------------------------------------------------------

#' Upper-tail probability under the normal or t reference distribution
#'
#' @param q Numeric vector of quantiles.
#' @param df `NULL` for the standard normal, or a numeric vector of degrees of
#'   freedom (recycled against `q`) for the t distribution.
#'
#' @return Numeric vector of upper-tail probabilities.
#' @noRd
.tost_upper_tail <- function(q, df = NULL) {
  if (is.null(df)) {
    stats::pnorm(q, lower.tail = FALSE)
  } else {
    stats::pt(q, df = df, lower.tail = FALSE)
  }
}

#' Validate the shared TOST arguments
#'
#' @param lfc Numeric vector of effect estimates.
#' @param se Numeric vector of standard errors.
#' @param lfc_threshold Positive scalar equivalence bound.
#' @param alpha Significance level in (0, 0.5).
#' @param df `NULL` or numeric vector of degrees of freedom.
#'
#' @return Invisibly `TRUE`; called for its side effect of erroring.
#' @noRd
.tost_check_args <- function(lfc, se, lfc_threshold, alpha, df) {
  stopifnot(is.numeric(lfc),
            is.numeric(se),
            length(lfc) == length(se))
  stopifnot(length(lfc_threshold) == 1,
            is.numeric(lfc_threshold),
            is.finite(lfc_threshold),
            lfc_threshold > 0)
  stopifnot(length(alpha) == 1,
            is.numeric(alpha),
            alpha > 0,
            alpha < 0.5)
  stopifnot(is.null(df) || (is.numeric(df) && length(df) %in% c(1, length(lfc))))
  if (any(se[!is.na(se)] <= 0)) stop("All standard errors in 'se' must be positive.")

  invisible(TRUE)
}

#' Two one-sided tests (TOST) from an estimate and its standard error
#'
#' Tests the interval null `H0: |beta| >= lfc_threshold` ("the gene is not
#' stably expressed") against `H1: |beta| < lfc_threshold` ("the gene is stably
#' expressed"), by rejecting both one-sided nulls at the full `alpha`.
#'
#' Also returns the ordinary two-sided Wald p-value for the usual difference
#' test, computed from the same `(lfc, se)` pair so that the two tests are
#' internally consistent, and the per-gene minimum rejectable bound.
#'
#' @param lfc Numeric vector of *unshrunken* log2 fold-change estimates.
#' @param se Numeric vector of standard errors, same length as `lfc`.
#' @param lfc_threshold Positive scalar: the equivalence bound (SESOI) on the
#'   log2 fold-change scale.  Must be pre-specified, never read off the data.
#' @param alpha Level of each one-sided test.  TOST is an intersection-union
#'   test, so this is the full alpha, not alpha / 2.  Used only for
#'   `lfc_threshold_min`; the returned p-values do not depend on it.
#' @param df `NULL` (default) for a standard normal reference, matching
#'   DESeq2's Wald test, or degrees of freedom for a t reference.
#' @param gene_names Optional character vector used as row names.
#'
#' @return A `data.frame` with one row per gene and columns
#'   \describe{
#'     \item{`lfc`, `se`}{the inputs}
#'     \item{`stat_tost`}{DESeq2's `lessAbs` statistic,
#'       `min(max((T - lfc)/se, 0), max((lfc + T)/se, 0))`}
#'     \item{`pvalue_lower`}{one-sided p for `H0: beta <= -lfc_threshold`}
#'     \item{`pvalue_upper`}{one-sided p for `H0: beta >= +lfc_threshold`}
#'     \item{`pvalue_tost`}{`max(pvalue_lower, pvalue_upper)`}
#'     \item{`pvalue_diff`}{ordinary two-sided Wald p for `H0: beta = 0`}
#'     \item{`lfc_threshold_min`}{`|lfc| + q_alpha * se`, the smallest
#'       equivalence bound this gene could have rejected at `alpha`.  This is
#'       the resource-based SESOI: reporting its distribution states what
#'       magnitude of effect the cohort can actually rule out.}
#'   }
#' @noRd
tost_from_estimates <- function(lfc,
                                se,
                                lfc_threshold,
                                alpha = 0.05,
                                df = NULL,
                                gene_names = NULL) {
  .tost_check_args(lfc, se, lfc_threshold, alpha, df)

  # Two one-sided tests, each at the full alpha (intersection-union test).
  #   H0_upper: beta >= +T   rejected when lfc is far enough BELOW +T
  #   H0_lower: beta <= -T   rejected when lfc is far enough ABOVE -T
  z_upper <- (lfc_threshold - lfc) / se
  z_lower <- (lfc + lfc_threshold) / se

  pvalue_upper <- .tost_upper_tail(z_upper, df = df)
  pvalue_lower <- .tost_upper_tail(z_lower, df = df)

  # Schuirmann (1987): the TOST p-value is the larger of the two.
  pvalue_tost <- pmax(pvalue_lower, pvalue_upper)
  stat_tost <- pmin(pmax(z_upper, 0), pmax(z_lower, 0))

  # Ordinary two-sided Wald test, from the same estimate and standard error.
  pvalue_diff <- pmin(1, 2 * .tost_upper_tail(abs(lfc) / se, df = df))

  # Smallest bound this gene could reject: TOST rejects at level alpha iff
  # both z's exceed q_alpha, i.e. iff lfc_threshold > |lfc| + q_alpha * se.
  q_alpha <- if (is.null(df)) {
    stats::qnorm(1 - alpha)
  } else {
    stats::qt(1 - alpha, df = df)
  }
  lfc_threshold_min <- abs(lfc) + q_alpha * se

  res <- data.frame(lfc = lfc,
                    se = se,
                    stat_tost = stat_tost,
                    pvalue_lower = pvalue_lower,
                    pvalue_upper = pvalue_upper,
                    pvalue_tost = pvalue_tost,
                    pvalue_diff = pvalue_diff,
                    lfc_threshold_min = lfc_threshold_min,
                    stringsAsFactors = FALSE)
  if (!is.null(gene_names)) {
    stopifnot(length(gene_names) == length(lfc))
    rownames(res) <- gene_names
  }

  res
}

#' Partition genes into differential / stable / trivial / undetermined
#'
#' Crossing the ordinary difference test with the equivalence test gives four
#' classes rather than two.  The point of the partition is that "few DEGs" is
#' currently reported identically whether the non-significant genes are
#' *provably stable* or merely *undetermined*; only the latter is evidence of
#' being underpowered.
#'
#' Multiplicity: Benjamini-Hochberg is applied separately to the difference
#' p-values and to the TOST p-values.  There is no correction *within* a gene's
#' TOST (intersection-union test).  Note the "stable" count tends to be an
#' under-estimate, since TOST is conservative away from the bound.
#'
#' @param pvalue_diff Numeric vector of two-sided difference p-values.
#' @param pvalue_tost Numeric vector of TOST p-values, same length.
#' @param alpha Level at which both adjusted p-values are thresholded.
#' @param p_adjust_method Passed to [stats::p.adjust()]; default `"BH"`.
#' @param gene_names Optional character vector used as row names.
#'
#' @return A `data.frame` with `padj_diff`, `padj_tost` and a factor `class`
#'   with levels `"differential"`, `"trivial"`, `"stable"`, `"undetermined"`.
#'   Genes with an `NA` in either adjusted p-value get `class = NA`.
#' @noRd
equivalence_partition <- function(pvalue_diff,
                                  pvalue_tost,
                                  alpha = 0.05,
                                  p_adjust_method = "BH",
                                  gene_names = NULL) {
  stopifnot(is.numeric(pvalue_diff),
            is.numeric(pvalue_tost),
            length(pvalue_diff) == length(pvalue_tost))
  stopifnot(length(alpha) == 1, is.numeric(alpha), alpha > 0, alpha < 1)

  padj_diff <- stats::p.adjust(pvalue_diff, method = p_adjust_method)
  padj_tost <- stats::p.adjust(pvalue_tost, method = p_adjust_method)

  is_diff <- padj_diff <= alpha
  is_equiv <- padj_tost <= alpha

  class_vec <- rep(NA_character_, length(pvalue_diff))
  usable <- !is.na(is_diff) & !is.na(is_equiv)
  class_vec[usable & is_diff & !is_equiv] <- "differential"
  class_vec[usable & is_diff & is_equiv] <- "trivial"
  class_vec[usable & !is_diff & is_equiv] <- "stable"
  class_vec[usable & !is_diff & !is_equiv] <- "undetermined"

  res <- data.frame(padj_diff = padj_diff,
                    padj_tost = padj_tost,
                    class = factor(class_vec,
                                   levels = c("differential", "trivial",
                                              "stable", "undetermined")),
                    stringsAsFactors = FALSE)
  if (!is.null(gene_names)) {
    stopifnot(length(gene_names) == length(pvalue_diff))
    rownames(res) <- gene_names
  }

  res
}

# ---------------------------------------------------------------------------
# DESeq2 fitting
# ---------------------------------------------------------------------------

#' Fit the DESeq2 negative-binomial GLM with unshrunken coefficients
#'
#' `betaPrior = FALSE` is not optional here: a zero-centered prior on the log2
#' fold change shrinks estimates toward zero, which inflates the count of genes
#' declared "stably expressed".  DESeq2 itself refuses `altHypothesis =
#' "lessAbs"` when `betaPrior = TRUE` for the same reason.
#'
#' @param mat_pseudobulk Integer count matrix, genes x pseudobulk samples.
#' @param metadata_pseudobulk `data.frame` of design variables.
#' @param use_t Passed to [DESeq2::DESeq()]; use a t reference distribution.
#'
#' @return A fitted `DESeqDataSet`.
#' @noRd
.fit_deseq2_mle <- function(mat_pseudobulk,
                            metadata_pseudobulk,
                            use_t = FALSE) {
  design_formula <- stats::as.formula(paste0("~ ", paste0(colnames(metadata_pseudobulk), collapse = "+")))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                        colData = metadata_pseudobulk,
                                        design = design_formula)
  DESeq2::DESeq(dds, betaPrior = FALSE, useT = use_t, quiet = TRUE)
}

#' Pull the unshrunken log2 fold change and its standard error out of a fit
#'
#' @param dds A fitted `DESeqDataSet`.
#' @param coef_name Name of the coefficient, as in `DESeq2::resultsNames()`.
#' @param cooks_filter Logical; if `TRUE`, genes flagged by DESeq2's Cook's
#'   distance cutoff are set to `NA`.  An outlier-driven fit should not be
#'   allowed to support an equivalence claim.
#' @param base_mean_min Genes with `baseMean` below this are set to `NA`.
#'
#' @return A `data.frame` with `base_mean`, `lfc`, `se`, `pvalue_deseq2`,
#'   `padj_deseq2` and `df` (or `NA` if a normal reference is in use).
#' @noRd
.extract_lfc_se <- function(dds,
                            coef_name,
                            cooks_filter = TRUE,
                            base_mean_min = 0) {
  # independentFiltering is switched off: DESeq2's filter is tuned to maximize
  # rejections of the *difference* null, and reusing it would silently change
  # which genes are eligible to be called stable.  Filtering is instead an
  # explicit, pre-specifiable baseMean threshold.
  res <- DESeq2::results(dds,
                         name = coef_name,
                         independentFiltering = FALSE,
                         cooksCutoff = cooks_filter)

  lfc <- res$log2FoldChange
  se <- res$lfcSE

  # Cook's-flagged genes have their p-value NA'd by DESeq2 but keep an lfc/se;
  # propagate that flag so they cannot be declared stable.
  drop_idx <- is.na(res$pvalue) | is.na(lfc) | is.na(se) | res$baseMean < base_mean_min
  lfc[drop_idx] <- NA
  se[drop_idx] <- NA

  df_vec <- SummarizedExperiment::mcols(dds)$tDegreesFreedom
  if (is.null(df_vec)) df_vec <- NA_real_

  data.frame(base_mean = res$baseMean,
             lfc = lfc,
             se = se,
             pvalue_deseq2 = res$pvalue,
             padj_deseq2 = res$padj,
             df = df_vec,
             row.names = rownames(res),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Top-level wrapper
# ---------------------------------------------------------------------------

#' Donor-level DESeq2 differential expression with a TOST equivalence test
#'
#' Aggregates a `Seurat` object to donor-level pseudobulk, fits the DESeq2
#' negative-binomial GLM with *unshrunken* coefficients, then runs both the
#' ordinary two-sided difference test and a TOST equivalence test against a
#' pre-specified log2 fold-change bound.  Genes are partitioned into
#' `differential` / `trivial` / `stable` / `undetermined`.
#'
#' The signature matches `deseq2_helper()` with three additional arguments;
#' `lfc_threshold` has no default because an equivalence bound read off the
#' data voids the type-I error guarantee.
#'
#' @param case_control_levels Length-2 character vector: control level first,
#'   then case.
#' @param case_control_var Name of the case/control metadata column.
#' @param categorical_vars Character vector of categorical covariates, or
#'   `NULL`.
#' @param id_var Name of the donor id metadata column.
#' @param numerical_vars Character vector of donor-constant numerical
#'   covariates, or `NULL`.
#' @param seurat_obj A `Seurat` object.
#' @param lfc_threshold Positive scalar: the SESOI on the log2 fold-change
#'   scale.  Pre-specify it.
#' @param alpha Level for both tests.  Each of the two one-sided tests runs at
#'   the full `alpha` (intersection-union test).
#' @param use_t Use a t rather than normal reference distribution.  Reasonable
#'   with few donors; `FALSE` matches DESeq2's default Wald test.
#' @param cooks_filter Set genes flagged by Cook's distance to `NA`.
#' @param base_mean_min Drop genes whose `baseMean` is below this.
#' @param p_adjust_method Passed to [stats::p.adjust()].
#'
#' @return A `data.frame` with one row per gene, holding `base_mean`, the
#'   unshrunken `lfc` and `se`, DESeq2's own p-values, the TOST quantities from
#'   [tost_from_estimates()], the adjusted p-values and the `class` factor from
#'   [equivalence_partition()].  The `lfc_threshold` used and the coefficient
#'   name are attached as attributes.
#' @noRd
deseq2_tost_helper <- function(case_control_levels, # Control and then Case
                               case_control_var,
                               categorical_vars,
                               id_var,
                               numerical_vars,
                               seurat_obj,
                               lfc_threshold,
                               alpha = 0.05,
                               use_t = FALSE,
                               cooks_filter = TRUE,
                               base_mean_min = 0,
                               p_adjust_method = "BH") {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Package 'DESeq2' is required. Install it with BiocManager::install('DESeq2').")
  stopifnot(length(lfc_threshold) == 1,
            is.numeric(lfc_threshold),
            is.finite(lfc_threshold),
            lfc_threshold > 0)

  # aggregate to one count vector per donor; see pseudobulk_claude.R
  pseudobulk <- build_pseudobulk(seurat_obj = seurat_obj,
                                 case_control_levels = case_control_levels,
                                 case_control_var = case_control_var,
                                 categorical_vars = categorical_vars,
                                 id_var = id_var,
                                 numerical_vars = numerical_vars)

  dds <- .fit_deseq2_mle(mat_pseudobulk = pseudobulk$mat,
                         metadata_pseudobulk = pseudobulk$metadata,
                         use_t = use_t)

  coef_name <- paste0(case_control_var, "_", case_control_levels[2], "_vs_", case_control_levels[1])
  estimate_df <- .extract_lfc_se(dds = dds,
                                 coef_name = coef_name,
                                 cooks_filter = cooks_filter,
                                 base_mean_min = base_mean_min)

  deseq2_tost(estimate_df = estimate_df,
              lfc_threshold = lfc_threshold,
              alpha = alpha,
              use_t = use_t,
              p_adjust_method = p_adjust_method,
              coef_name = coef_name)
}

#' Run the TOST partition on an already-extracted table of estimates
#'
#' Separated from [deseq2_tost_helper()] so the same partition can be applied
#' to any method that yields a per-gene effect and standard error (dreamlet,
#' NEBULA, eSVD-DE), and so it can be tested without fitting a GLM.
#'
#' @param estimate_df `data.frame` with at least `lfc` and `se` columns, and
#'   optionally `base_mean`, `pvalue_deseq2`, `padj_deseq2`, `df`.
#' @param lfc_threshold,alpha,p_adjust_method See [deseq2_tost_helper()].
#' @param use_t Whether `estimate_df$df` should be used as the reference
#'   distribution's degrees of freedom.
#' @param coef_name Optional string recorded as an attribute.
#'
#' @return A `data.frame`; see [deseq2_tost_helper()].
#' @noRd
deseq2_tost <- function(estimate_df,
                        lfc_threshold,
                        alpha = 0.05,
                        use_t = FALSE,
                        p_adjust_method = "BH",
                        coef_name = NULL) {
  stopifnot(is.data.frame(estimate_df),
            all(c("lfc", "se") %in% colnames(estimate_df)))

  df_vec <- NULL
  if (use_t && "df" %in% colnames(estimate_df) && all(!is.na(estimate_df$df))) {
    df_vec <- estimate_df$df
  }

  # NA lfc/se (dropped genes) cannot go through the TOST arithmetic; run the
  # test on the usable genes and re-expand.
  usable <- !is.na(estimate_df$lfc) & !is.na(estimate_df$se)
  tost_cols <- c("stat_tost", "pvalue_lower", "pvalue_upper", "pvalue_tost",
                 "pvalue_diff", "lfc_threshold_min")
  tost_df <- as.data.frame(matrix(NA_real_,
                                  nrow = nrow(estimate_df),
                                  ncol = length(tost_cols),
                                  dimnames = list(rownames(estimate_df), tost_cols)))
  if (any(usable)) {
    tost_sub <- tost_from_estimates(lfc = estimate_df$lfc[usable],
                                    se = estimate_df$se[usable],
                                    lfc_threshold = lfc_threshold,
                                    alpha = alpha,
                                    df = if (is.null(df_vec)) NULL else df_vec[usable])
    tost_df[usable, ] <- tost_sub[, tost_cols]
  }

  partition_df <- equivalence_partition(pvalue_diff = tost_df$pvalue_diff,
                                        pvalue_tost = tost_df$pvalue_tost,
                                        alpha = alpha,
                                        p_adjust_method = p_adjust_method)

  res <- cbind(estimate_df[, setdiff(colnames(estimate_df), c("lfc", "se")), drop = FALSE],
               estimate_df[, c("lfc", "se"), drop = FALSE],
               tost_df,
               partition_df)
  rownames(res) <- rownames(estimate_df)

  attr(res, "lfc_threshold") <- lfc_threshold
  attr(res, "alpha") <- alpha
  attr(res, "coef_name") <- coef_name

  res
}
