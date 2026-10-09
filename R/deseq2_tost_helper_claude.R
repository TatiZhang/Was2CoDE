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
# This file is about equivalence and nothing else: it returns no
# differential-expression p-values.  Run `deseq2_helper()` for those and join
# the two tables by gene name if the four-way differential / stable /
# "significant but negligible" / undetermined partition is wanted.
#
# The arithmetic here reproduces DESeq2's own `altHypothesis = "lessAbs"`
# branch of `DESeq2::results()`; it is factored out so it can be unit tested
# without Bioconductor and reused for any method that yields a per-gene
# (estimate, standard error) pair.

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
#' Also returns the per-gene minimum rejectable bound.  No difference test is
#' computed here -- see [deseq2_helper()] for that.
#'
#' @param lfc Numeric vector of *unshrunken* log2 fold-change estimates.
#' @param se Numeric vector of standard errors, same length as `lfc`.
#' @param lfc_threshold Positive scalar: the equivalence bound (SESOI) on the
#'   log2 fold-change scale.  Must be pre-specified, never read off the data.
#' @param alpha Level of each one-sided test.  TOST is an intersection-union
#'   test, so this is the full alpha, not alpha / 2.  Used only for
#'   `lfc_threshold_min`; the returned p-values do not depend on it.  This is a
#'   single-test function and knows nothing about the gene set, so the level is
#'   taken at face value -- see the note on `lfc_threshold_min` below.
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
#'     \item{`lfc_threshold_min`}{`|lfc| + q_alpha * se`, the smallest
#'       equivalence bound this gene could have rejected at `alpha`.  This is
#'       the resource-based SESOI: reporting its distribution states what
#'       magnitude of effect the cohort can actually rule out.  **Computed on
#'       the unadjusted, per-test scale**, because this function sees one gene
#'       at a time and multiplicity is not defined without the gene set.
#'       [deseq2_tost()] recomputes it on the multiplicity-corrected scale so
#'       that it agrees with the `equivalent` call.}
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
                    lfc_threshold_min = lfc_threshold_min,
                    stringsAsFactors = FALSE)
  if (!is.null(gene_names)) {
    stopifnot(length(gene_names) == length(lfc))
    rownames(res) <- gene_names
  }

  res
}

#' The raw-p-value threshold that a multiple-testing adjustment actually applied
#'
#' Every method in [stats::p.adjust()] is monotone in the raw p-value, so its
#' rejection set is always `{p <= t}` for some data-dependent `t`.  Recovering
#' that `t` is what lets a per-gene quantity derived from the raw scale (here,
#' `lfc_threshold_min`) be put on the same footing as the adjusted call.  For
#' Benjamini-Hochberg, `t` is the largest rejected p-value, equal to
#' `alpha * k / m` with `k` rejections among `m` tests.
#'
#' The value is exact for the gene set as it stands.  It is *not* a fixed
#' constant of the design: change `lfc_threshold` and every p-value moves, so
#' `k` -- and therefore `t` -- moves with it.
#'
#' @param pvalue Numeric vector of raw p-values (`NA`s ignored).
#' @param alpha The nominal level applied to the adjusted p-values.
#' @param p_adjust_method Passed to [stats::p.adjust()].
#'
#' @return A single numeric threshold on the raw p-value scale.
#' @noRd
.effective_alpha <- function(pvalue, alpha, p_adjust_method = "BH") {
  p <- pvalue[!is.na(pvalue)]
  m <- length(p)
  # no correction, or nothing to correct over: the nominal level is applied
  if (m == 0 || identical(p_adjust_method, "none")) return(alpha)

  padj <- stats::p.adjust(p, method = p_adjust_method)
  rejected <- !is.na(padj) & padj <= alpha
  k <- sum(rejected)

  # Nothing rejected: no observed p-value pins the threshold, so report the
  # level the single most significant gene would have had to beat.
  if (k == 0) return(alpha / m)

  # Benjamini-Hochberg rejects exactly {p <= alpha * k / m}: every rejected
  # p is at most that, and p_(k+1) > (k+1) * alpha / m by construction.
  if (p_adjust_method %in% c("BH", "fdr")) return(alpha * k / m)

  # Any other monotone correction: fall back to the attained level, the largest
  # p that survived.  Still reproduces the rejection set exactly, but it is a
  # realized rather than nominal threshold.
  max(p[rejected])
}

#' Adjust TOST p-values across genes and call the equivalent ones
#'
#' Multiplicity: Benjamini-Hochberg across genes.  There is no correction
#' *within* a gene's TOST -- it is an intersection-union test, so each of the
#' two one-sided tests already runs at the full alpha.  Note the count of
#' equivalent genes tends to be an under-estimate, since TOST is conservative
#' away from the bound.
#'
#' @param pvalue_tost Numeric vector of TOST p-values.
#' @param alpha Level at which the adjusted p-value is thresholded.
#' @param p_adjust_method Passed to [stats::p.adjust()]; default `"BH"`.
#' @param gene_names Optional character vector used as row names.
#'
#' @return A `data.frame` with `padj_tost` and the logical `equivalent`
#'   (`NA` where the p-value is `NA`).
#' @noRd
equivalence_call <- function(pvalue_tost,
                             alpha = 0.05,
                             p_adjust_method = "BH",
                             gene_names = NULL) {
  stopifnot(is.numeric(pvalue_tost))
  stopifnot(length(alpha) == 1, is.numeric(alpha), alpha > 0, alpha < 1)

  padj_tost <- stats::p.adjust(pvalue_tost, method = p_adjust_method)

  res <- data.frame(padj_tost = padj_tost,
                    equivalent = padj_tost <= alpha,
                    stringsAsFactors = FALSE)
  if (!is.null(gene_names)) {
    stopifnot(length(gene_names) == length(pvalue_tost))
    rownames(res) <- gene_names
  }

  res
}

# ---------------------------------------------------------------------------
# Top-level wrapper
# ---------------------------------------------------------------------------

#' Donor-level DESeq2 TOST equivalence test
#'
#' Aggregates a `Seurat` object to donor-level pseudobulk, fits the DESeq2
#' negative-binomial GLM with *unshrunken* coefficients, and tests each gene
#' for equivalence against a pre-specified log2 fold-change bound.
#'
#' This function tests **only** equivalence -- it returns no
#' differential-expression p-values.  For those, run [deseq2_helper()] on the
#' same object and join the two tables by gene name; crossing the two calls
#' recovers the four-way differential / stable / "significant but negligible"
#' / undetermined partition, but that is the caller's business, not this
#' function's.
#'
#' The signature matches `deseq2_helper()` with extra arguments;
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
#' @param alpha Level of the test.  Each of the two one-sided tests runs at the
#'   full `alpha` (intersection-union test).
#' @param use_t Use a t rather than normal reference distribution.  Reasonable
#'   with few donors; `FALSE` matches DESeq2's default Wald test.
#' @param cooks_filter Set genes flagged by Cook's distance to `NA`.
#' @param base_mean_min Drop genes whose `baseMean` is below this.
#' @param p_adjust_method Passed to [stats::p.adjust()].
#'
#' @return A `data.frame` with one row per gene, holding `base_mean`, the
#'   unshrunken `lfc` and `se`, the TOST quantities from
#'   [tost_from_estimates()], and `padj_tost` / `equivalent` from
#'   [equivalence_call()].  `lfc_threshold_min <= lfc_threshold` agrees with
#'   `equivalent` exactly, because the former is put on the
#'   multiplicity-corrected scale.  The `lfc_threshold` and `alpha` used, the
#'   raw-scale threshold the correction actually applied (`alpha_effective`),
#'   and the coefficient name are attached as attributes.
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

  # aggregate to one count vector per donor; see pseudobulk.R
  pseudobulk <- build_pseudobulk(seurat_obj = seurat_obj,
                                 case_control_levels = case_control_levels,
                                 case_control_var = case_control_var,
                                 categorical_vars = categorical_vars,
                                 id_var = id_var,
                                 numerical_vars = numerical_vars)

  dds <- .fit_deseq2_mle(mat_pseudobulk = pseudobulk$mat,
                         metadata_pseudobulk = pseudobulk$metadata,
                         use_t = use_t,
                         quiet = TRUE)

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

#' Run the TOST on an already-extracted table of estimates
#'
#' Separated from [deseq2_tost_helper()] so the same test can be applied to any
#' method that yields a per-gene effect and standard error (dreamlet, NEBULA,
#' eSVD-DE), and so it can be tested without fitting a GLM.
#'
#' @param estimate_df `data.frame` with at least `lfc` and `se` columns, and
#'   optionally `base_mean`, `df`, `cooks_outlier`.
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
                 "lfc_threshold_min")
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

  call_df <- equivalence_call(pvalue_tost = tost_df$pvalue_tost,
                              alpha = alpha,
                              p_adjust_method = p_adjust_method)

  # `lfc_threshold_min` arrives from tost_from_estimates() on the unadjusted
  # per-test scale, which would contradict `equivalent`: a gene whose raw
  # p-value clears alpha but whose adjusted p-value does not would be reported
  # as having a minimum rejectable bound inside `lfc_threshold` while being
  # called not equivalent.  Put it back on the scale the adjustment actually
  # applied, so the two agree exactly.
  alpha_effective <- .effective_alpha(tost_df$pvalue_tost, alpha, p_adjust_method)
  q_effective <- if (is.null(df_vec)) {
    stats::qnorm(1 - alpha_effective)
  } else {
    stats::qt(1 - alpha_effective, df = df_vec)
  }
  tost_df$lfc_threshold_min <- abs(estimate_df$lfc) + q_effective * estimate_df$se

  res <- cbind(estimate_df[, setdiff(colnames(estimate_df), c("lfc", "se")), drop = FALSE],
               estimate_df[, c("lfc", "se"), drop = FALSE],
               tost_df,
               call_df)
  rownames(res) <- rownames(estimate_df)

  attr(res, "lfc_threshold") <- lfc_threshold
  attr(res, "alpha") <- alpha
  attr(res, "alpha_effective") <- alpha_effective
  attr(res, "coef_name") <- coef_name

  res
}
