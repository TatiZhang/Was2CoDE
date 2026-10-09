context("Test deseq2_tost_helper (TOST equivalence testing)")

# ---------------------------------------------------------------------------
# tost_from_estimates: the core arithmetic
# ---------------------------------------------------------------------------

test_that("tost_from_estimates matches the closed-form one-sided p-values", {
  lfc <- c(0, 0.5, -0.5, 3)
  se <- c(1, 0.25, 0.25, 1)
  bound <- 2

  res <- tost_from_estimates(lfc = lfc, se = se, lfc_threshold = bound)

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), length(lfc))

  expect_equal(res$pvalue_upper, stats::pnorm((bound - lfc) / se, lower.tail = FALSE))
  expect_equal(res$pvalue_lower, stats::pnorm((lfc + bound) / se, lower.tail = FALSE))
  expect_equal(res$pvalue_tost, pmax(res$pvalue_lower, res$pvalue_upper))
  expect_true(all(res$pvalue_tost >= 0 & res$pvalue_tost <= 1))
})

test_that("the TOST p-value is symmetric in the sign of the effect", {
  lfc <- c(-1.3, -0.4, 0, 0.4, 1.3)
  se <- rep(0.4, 5)

  res <- tost_from_estimates(lfc = lfc, se = se, lfc_threshold = 1)

  expect_equal(res$pvalue_tost, rev(res$pvalue_tost))
  # and the two one-sided tests swap roles under a sign flip
  expect_equal(res$pvalue_lower, rev(res$pvalue_upper))
})

test_that("TOST is monotone in the bound and in the standard error", {
  # a wider equivalence bound is easier to reject
  p_narrow <- tost_from_estimates(lfc = 0.1, se = 0.3, lfc_threshold = 0.5)$pvalue_tost
  p_wide <- tost_from_estimates(lfc = 0.1, se = 0.3, lfc_threshold = 2.0)$pvalue_tost
  expect_lt(p_wide, p_narrow)

  # a noisier estimate makes equivalence harder to establish
  p_precise <- tost_from_estimates(lfc = 0, se = 0.1, lfc_threshold = 1)$pvalue_tost
  p_noisy <- tost_from_estimates(lfc = 0, se = 1.0, lfc_threshold = 1)$pvalue_tost
  expect_lt(p_precise, p_noisy)

  # an effect far outside the bound is nowhere near equivalent
  expect_gt(tost_from_estimates(lfc = 5, se = 0.1, lfc_threshold = 1)$pvalue_tost, 0.99)
})

test_that("no alpha-splitting: rejection needs both one-sided tests at full alpha", {
  alpha <- 0.05
  # lfc sits just below the upper bound, so only the lower test rejects
  res <- tost_from_estimates(lfc = 0.9, se = 0.3, lfc_threshold = 1)
  expect_lt(res$pvalue_lower, alpha)
  expect_gt(res$pvalue_upper, alpha)
  expect_equal(res$pvalue_tost, res$pvalue_upper)
  expect_gt(res$pvalue_tost, alpha)
})

test_that("lfc_threshold_min is exactly the smallest rejectable bound", {
  alpha <- 0.05
  lfc <- 0.3
  se <- 0.2

  min_bound <- tost_from_estimates(lfc = lfc, se = se,
                                   lfc_threshold = 1, alpha = alpha)$lfc_threshold_min
  expect_equal(min_bound, abs(lfc) + stats::qnorm(1 - alpha) * se)

  # just inside the minimum bound: fails; just outside: rejects
  p_just_below <- tost_from_estimates(lfc = lfc, se = se,
                                      lfc_threshold = min_bound - 1e-6,
                                      alpha = alpha)$pvalue_tost
  p_just_above <- tost_from_estimates(lfc = lfc, se = se,
                                      lfc_threshold = min_bound + 1e-6,
                                      alpha = alpha)$pvalue_tost
  expect_gt(p_just_below, alpha)
  expect_lt(p_just_above, alpha)
})

test_that("supplying df uses a t reference and is more conservative than normal", {
  res_norm <- tost_from_estimates(lfc = 0, se = 0.3, lfc_threshold = 1)
  res_t <- tost_from_estimates(lfc = 0, se = 0.3, lfc_threshold = 1, df = 5)

  expect_equal(res_t$pvalue_upper, stats::pt((1 - 0) / 0.3, df = 5, lower.tail = FALSE))
  expect_gt(res_t$pvalue_tost, res_norm$pvalue_tost)
  expect_gt(res_t$lfc_threshold_min, res_norm$lfc_threshold_min)
})

test_that("tost_from_estimates rejects malformed input", {
  expect_error(tost_from_estimates(lfc = c(0, 1), se = 1, lfc_threshold = 1))
  expect_error(tost_from_estimates(lfc = 0, se = 1, lfc_threshold = 0))
  expect_error(tost_from_estimates(lfc = 0, se = 1, lfc_threshold = -1))
  expect_error(tost_from_estimates(lfc = 0, se = -1, lfc_threshold = 1))
  expect_error(tost_from_estimates(lfc = 0, se = 1, lfc_threshold = 1, alpha = 0.6))
})

# ---------------------------------------------------------------------------
# Calibration: type-I error at the boundary, conservatism away from it, power
# ---------------------------------------------------------------------------

test_that("TOST holds its nominal size at the equivalence boundary", {
  set.seed(1)
  n_sim <- 20000
  alpha <- 0.05
  bound <- 1
  se <- 0.4

  # worst case under H0: the true effect sits exactly on the bound
  lfc_boundary <- stats::rnorm(n_sim, mean = bound, sd = se)
  p_boundary <- tost_from_estimates(lfc = lfc_boundary, se = rep(se, n_sim),
                                    lfc_threshold = bound)$pvalue_tost
  rate_boundary <- mean(p_boundary <= alpha)
  expect_lt(rate_boundary, alpha + 0.01)
  expect_gt(rate_boundary, alpha - 0.01)

  # away from the boundary the test is conservative (Lauzon 2009)
  lfc_far <- stats::rnorm(n_sim, mean = 2 * bound, sd = se)
  rate_far <- mean(tost_from_estimates(lfc = lfc_far, se = rep(se, n_sim),
                                       lfc_threshold = bound)$pvalue_tost <= alpha)
  expect_lt(rate_far, rate_boundary)
  expect_lt(rate_far, 1e-3)
})

test_that("TOST has power to declare equivalence when the effect is truly null", {
  set.seed(2)
  n_sim <- 5000
  alpha <- 0.05
  bound <- 1
  se <- 0.2

  lfc_null <- stats::rnorm(n_sim, mean = 0, sd = se)
  rate <- mean(tost_from_estimates(lfc = lfc_null, se = rep(se, n_sim),
                                   lfc_threshold = bound)$pvalue_tost <= alpha)
  expect_gt(rate, 0.9)
})

test_that("shrinking the fold change toward zero inflates the stable count", {
  # documents the wiki's trap: lfcShrink-style shrinkage is anti-conservative
  # for an equivalence test, so the helper must use unshrunken estimates.
  set.seed(3)
  n_sim <- 20000
  alpha <- 0.05
  bound <- 1
  se <- 0.4

  lfc_boundary <- stats::rnorm(n_sim, mean = bound, sd = se)
  rate_mle <- mean(tost_from_estimates(lfc = lfc_boundary, se = rep(se, n_sim),
                                       lfc_threshold = bound)$pvalue_tost <= alpha)
  rate_shrunk <- mean(tost_from_estimates(lfc = 0.5 * lfc_boundary, se = rep(se, n_sim),
                                          lfc_threshold = bound)$pvalue_tost <= alpha)
  expect_gt(rate_shrunk, rate_mle)
  expect_gt(rate_shrunk, alpha)
})

# ---------------------------------------------------------------------------
# equivalence_call
# ---------------------------------------------------------------------------

test_that("equivalence_call thresholds the adjusted TOST p-value", {
  # with "none" adjustment the thresholding is exact
  pvalue_tost <- c(0.001, 0.049, 0.051, 0.900, NA)

  res <- equivalence_call(pvalue_tost, alpha = 0.05, p_adjust_method = "none")

  expect_equal(res$equivalent, c(TRUE, TRUE, FALSE, FALSE, NA))
  expect_equal(res$padj_tost, pvalue_tost)
})

test_that("equivalence_call applies BH across genes", {
  set.seed(4)
  pvalue_tost <- stats::runif(100)

  res <- equivalence_call(pvalue_tost, alpha = 0.05)

  expect_equal(res$padj_tost, stats::p.adjust(pvalue_tost, method = "BH"))
  # BH is never anti-conservative relative to the raw p-values
  expect_true(all(res$padj_tost >= pvalue_tost))
  expect_equal(res$equivalent, res$padj_tost <= 0.05)
})

test_that("equivalence_call carries gene names through", {
  res <- equivalence_call(c(0.5, 0.5), gene_names = c("GeneA", "GeneB"))
  expect_equal(rownames(res), c("GeneA", "GeneB"))
})

# ---------------------------------------------------------------------------
# Multiplicity: lfc_threshold_min must live on the same scale as `equivalent`
# ---------------------------------------------------------------------------

test_that(".effective_alpha recovers the raw-p threshold the adjustment applied", {
  set.seed(20)
  pvalue <- c(stats::runif(40, 0, 0.01), stats::runif(60, 0, 1))
  alpha <- 0.05

  t_eff <- .effective_alpha(pvalue, alpha, "BH")

  # a p.adjust method's rejection set is always {p <= t}; check that this t
  # reproduces it exactly
  padj <- stats::p.adjust(pvalue, method = "BH")
  expect_equal(pvalue <= t_eff, padj <= alpha)
  # and that it is BH's closed form, alpha * k / m
  expect_equal(t_eff, alpha * sum(padj <= alpha) / length(pvalue))
  # correcting for multiplicity can only tighten the threshold
  expect_lt(t_eff, alpha)
})

test_that(".effective_alpha handles no rejections and no correction", {
  # nothing rejected: fall back to what the most significant gene must beat
  expect_equal(.effective_alpha(rep(0.9, 20), 0.05, "BH"), 0.05 / 20)
  # "none" means the nominal level is already the applied threshold
  expect_equal(.effective_alpha(rep(0.9, 20), 0.05, "none"), 0.05)
  expect_equal(.effective_alpha(numeric(0), 0.05, "BH"), 0.05)
})

test_that("lfc_threshold_min agrees with `equivalent` after multiple testing", {
  # Genes engineered so their TOST p-values straddle the gap between the
  # nominal alpha and the stricter threshold BH actually applies.  Before the
  # correction was carried into lfc_threshold_min, these reported a minimum
  # rejectable bound inside lfc_threshold while being called not equivalent.
  bound <- 0.5
  alpha <- 0.05
  estimate_df <- data.frame(lfc = c(seq(0.28, 0.36, by = 0.005), rep(0.45, 30)),
                            se = 0.1)
  rownames(estimate_df) <- paste0("g", seq_len(nrow(estimate_df)))

  res <- deseq2_tost(estimate_df, lfc_threshold = bound, alpha = alpha)

  # the disputed band is genuinely populated, so this is not a vacuous check
  expect_gt(sum(res$pvalue_tost > attr(res, "alpha_effective") &
                res$pvalue_tost <= alpha), 0)

  expect_equal(res$lfc_threshold_min <= bound, res$equivalent)
  expect_lt(attr(res, "alpha_effective"), alpha)
})

test_that("lfc_threshold_min agrees with `equivalent` under a t reference too", {
  bound <- 0.5
  alpha <- 0.05
  estimate_df <- data.frame(lfc = c(seq(0.20, 0.34, by = 0.005), rep(0.45, 30)),
                            se = 0.1,
                            df = 10)
  rownames(estimate_df) <- paste0("g", seq_len(nrow(estimate_df)))

  res <- deseq2_tost(estimate_df, lfc_threshold = bound, alpha = alpha, use_t = TRUE)

  expect_equal(res$lfc_threshold_min <= bound, res$equivalent)
})

test_that("with no correction, lfc_threshold_min is the unadjusted bound", {
  bound <- 0.5
  alpha <- 0.05
  estimate_df <- data.frame(lfc = c(0.1, 0.3, 0.45), se = 0.1,
                            row.names = c("a", "b", "c"))

  res <- deseq2_tost(estimate_df, lfc_threshold = bound, alpha = alpha,
                     p_adjust_method = "none")

  expect_equal(attr(res, "alpha_effective"), alpha)
  expect_equal(res$lfc_threshold_min,
               abs(estimate_df$lfc) + stats::qnorm(1 - alpha) * estimate_df$se)
  expect_equal(res$lfc_threshold_min <= bound, res$equivalent)
})

# ---------------------------------------------------------------------------
# Agreement with DESeq2's own altHypothesis = "lessAbs"
# ---------------------------------------------------------------------------

test_that("tost_from_estimates reproduces DESeq2's altHypothesis='lessAbs'", {
  skip_if_not_installed("DESeq2")
  set.seed(5)

  dds <- DESeq2::makeExampleDESeqDataSet(n = 300, m = 12, betaSD = 1)
  dds <- DESeq2::DESeq(dds, betaPrior = FALSE, quiet = TRUE)
  bound <- 0.5

  res_mle <- DESeq2::results(dds, name = "condition_B_vs_A",
                             independentFiltering = FALSE, cooksCutoff = FALSE)
  res_deseq2_tost <- DESeq2::results(dds, name = "condition_B_vs_A",
                                     lfcThreshold = bound, altHypothesis = "lessAbs",
                                     independentFiltering = FALSE, cooksCutoff = FALSE)

  keep <- !is.na(res_mle$log2FoldChange) & !is.na(res_mle$lfcSE) & res_mle$lfcSE > 0
  ours <- tost_from_estimates(lfc = res_mle$log2FoldChange[keep],
                              se = res_mle$lfcSE[keep],
                              lfc_threshold = bound)

  expect_gt(sum(keep), 100)
  expect_equal(ours$pvalue_tost, res_deseq2_tost$pvalue[keep], tolerance = 1e-10)
  expect_equal(ours$stat_tost, res_deseq2_tost$stat[keep], tolerance = 1e-10)
  # no difference test is returned -- that is deseq2_helper()'s job
  expect_false("pvalue_diff" %in% colnames(ours))
})

# ---------------------------------------------------------------------------
# End-to-end on a simulated Seurat object
# ---------------------------------------------------------------------------

test_that("deseq2_tost_helper calls the right simulated genes equivalent", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat()
  skip_if(is.null(sim))

  res <- deseq2_tost_helper(case_control_levels = c("control", "case"),
                            case_control_var = "diagnosis",
                            categorical_vars = NULL,
                            id_var = "donor",
                            numerical_vars = NULL,
                            seurat_obj = sim$seurat_obj,
                            lfc_threshold = 0.5,
                            alpha = 0.05)

  expect_true(is.data.frame(res))
  expect_equal(sort(rownames(res)), sort(names(sim$gene_class)))
  expect_true(all(c("lfc", "se", "pvalue_tost", "padj_tost", "equivalent") %in% colnames(res)))
  expect_equal(attr(res, "lfc_threshold"), 0.5)

  # strictly an equivalence test: no differential-expression p-values
  expect_false(any(c("pvalue_diff", "padj_diff", "class",
                     "pvalue_deseq2", "padj_deseq2") %in% colnames(res)))

  gene_class <- sim$gene_class[rownames(res)]
  frac_equivalent <- tapply(res$equivalent, gene_class, mean)

  # well-expressed null genes are positively established as equivalent
  expect_gt(frac_equivalent[["null"]], 0.7)
  # genes with a real 4-fold effect are far outside the bound
  expect_equal(frac_equivalent[["de"]], 0)
  # barely-expressed genes are too imprecise to establish equivalence -- they
  # are NOT equivalent, but for a completely different reason than the DE genes
  expect_equal(frac_equivalent[["noisy"]], 0)

  # ...and that reason is visible in the minimum rejectable bound: noisy genes
  # fail because the bound they could reject is huge, DE genes because |lfc| is
  expect_gt(stats::median(res$lfc_threshold_min[gene_class == "noisy"], na.rm = TRUE),
            stats::median(res$lfc_threshold_min[gene_class == "null"], na.rm = TRUE))
  expect_lt(stats::median(abs(res$lfc[gene_class == "noisy"]), na.rm = TRUE), 1)
  expect_gt(stats::median(abs(res$lfc[gene_class == "de"])), 1.5)
  expect_lt(stats::median(abs(res$lfc[gene_class == "null"])), 0.3)
})

test_that("deseq2_tost_helper accepts a donor-constant numerical covariate", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                           n_gene = c(de = 5, null = 20, noisy = 10), seed = 11)
  skip_if(is.null(sim))

  res <- deseq2_tost_helper(case_control_levels = c("control", "case"),
                            case_control_var = "diagnosis",
                            categorical_vars = NULL,
                            id_var = "donor",
                            numerical_vars = "age",
                            seurat_obj = sim$seurat_obj,
                            lfc_threshold = 1)

  expect_equal(nrow(res), length(sim$gene_class))
  expect_true(all(res$pvalue_tost >= 0 & res$pvalue_tost <= 1, na.rm = TRUE))
})

test_that("deseq2_tost_helper propagates the t reference distribution", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                           n_gene = c(de = 5, null = 20, noisy = 10), seed = 21)
  skip_if(is.null(sim))

  res <- deseq2_tost_helper(case_control_levels = c("control", "case"),
                            case_control_var = "diagnosis",
                            categorical_vars = NULL,
                            id_var = "donor",
                            numerical_vars = NULL,
                            seurat_obj = sim$seurat_obj,
                            lfc_threshold = 0.5,
                            use_t = TRUE)

  # 12 pseudobulk samples minus the intercept and the case/control coefficient
  expect_equal(unique(res$df), 10)
  expect_true(all(!is.na(res$pvalue_tost)))

  # the p-values really do come from the t distribution
  expect_equal(res$pvalue_upper,
               stats::pt((0.5 - res$lfc) / res$se, df = res$df, lower.tail = FALSE))

  gene_class <- sim$gene_class[rownames(res)]
  expect_true(all(res$equivalent[gene_class == "null"]))
  expect_false(any(res$equivalent[gene_class %in% c("de", "noisy")]))
})

test_that("deseq2_tost_helper requires a positive, pre-specified bound", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 20,
                           n_gene = c(de = 2, null = 5, noisy = 3), seed = 12)
  skip_if(is.null(sim))

  expect_error(deseq2_tost_helper(case_control_levels = c("control", "case"),
                                  case_control_var = "diagnosis",
                                  categorical_vars = NULL,
                                  id_var = "donor",
                                  numerical_vars = NULL,
                                  seurat_obj = sim$seurat_obj))
  expect_error(deseq2_tost_helper(case_control_levels = c("control", "case"),
                                  case_control_var = "diagnosis",
                                  categorical_vars = NULL,
                                  id_var = "donor",
                                  numerical_vars = NULL,
                                  seurat_obj = sim$seurat_obj,
                                  lfc_threshold = 0))
})

# ---------------------------------------------------------------------------
# The four-way partition: deseq2_helper() crossed with deseq2_tost()
# ---------------------------------------------------------------------------
#
# Neither function produces the partition on its own -- deseq2_tost_helper()
# tests only equivalence, deseq2_helper() only difference.  Crossing them is
# what separates "few DEGs because the genes are provably stable" from "few
# DEGs because the cohort was underpowered", which is the point of the whole
# exercise (see the wiki's [[analysis-tost-for-deg-equivalence]]).

test_that("deseq2_tost can reuse deseq2_helper's estimates without a second fit", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat()
  skip_if(is.null(sim))

  res_de <- deseq2_helper(case_control_levels = c("control", "case"),
                          case_control_var = "diagnosis",
                          categorical_vars = NULL,
                          id_var = "donor",
                          numerical_vars = NULL,
                          seurat_obj = sim$seurat_obj)$original

  # deseq2_helper() already returns the unshrunken log2FC and its standard
  # error, which is everything the equivalence test needs
  estimate_df <- data.frame(base_mean = res_de$baseMean,
                            lfc = res_de$log2FoldChange,
                            se = res_de$lfcSE,
                            row.names = rownames(res_de))
  res_tost <- deseq2_tost(estimate_df, lfc_threshold = 0.5, alpha = 0.05)

  # This route is exact, not an approximation: independent filtering and the
  # Cook's cutoff only ever touch p-values, never log2FoldChange or lfcSE, so
  # the estimates are bit-identical to the ones deseq2_tost_helper() extracts.
  res_full <- deseq2_tost_helper(case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL,
                                 seurat_obj = sim$seurat_obj,
                                 lfc_threshold = 0.5,
                                 alpha = 0.05)
  res_full <- res_full[rownames(res_tost), ]

  expect_equal(res_tost$lfc, res_full$lfc)
  expect_equal(res_tost$se, res_full$se)
  expect_equal(res_tost$pvalue_tost, res_full$pvalue_tost)
  expect_equal(res_tost$equivalent, res_full$equivalent)
})

test_that("crossing deseq2_helper with deseq2_tost recovers the four-way partition", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  alpha <- 0.05
  sim <- .simulate_donor_seurat()
  skip_if(is.null(sim))

  res_de <- deseq2_helper(case_control_levels = c("control", "case"),
                          case_control_var = "diagnosis",
                          categorical_vars = NULL,
                          id_var = "donor",
                          numerical_vars = NULL,
                          seurat_obj = sim$seurat_obj)$original
  res_tost <- deseq2_tost(data.frame(lfc = res_de$log2FoldChange,
                                     se = res_de$lfcSE,
                                     row.names = rownames(res_de)),
                          lfc_threshold = 0.5,
                          alpha = alpha)

  # the two tables must be joinable gene for gene
  expect_equal(rownames(res_de), rownames(res_tost))

  # DESeq2's independent filtering leaves padj = NA for genes it declined to
  # test.  Treat those as "not differential": a gene was filtered precisely
  # because it had no prospect of significance, so the equivalence test is what
  # decides whether it is stable or merely undetermined.  Dropping them instead
  # would discard exactly the low-information genes the partition exists to
  # classify.
  is_diff <- !is.na(res_de$padj) & res_de$padj <= alpha
  is_equiv <- res_tost$equivalent
  expect_gt(sum(is.na(res_de$padj)), 0)

  partition <- factor(ifelse(is_diff & !is_equiv, "differential",
                      ifelse(is_diff & is_equiv, "significant but negligible",
                      ifelse(!is_diff & is_equiv, "stable", "undetermined"))),
                      levels = c("differential", "significant but negligible",
                                 "stable", "undetermined"))

  gene_class <- sim$gene_class[rownames(res_tost)]
  tab <- table(gene_class, partition)

  # a real 4-fold effect: significant, and not equivalent
  expect_gt(tab["de", "differential"] / sum(gene_class == "de"), 0.9)
  expect_equal(unname(tab["de", "stable"]), 0)

  # well expressed with no true effect: positively established as stable
  expect_gt(tab["null", "stable"] / sum(gene_class == "null"), 0.7)
  expect_equal(unname(tab["null", "undetermined"]), 0)

  # barely expressed with no true effect: neither test can conclude anything.
  # These are the genes that a DEG count alone would misreport as evidence of
  # stability, and they are the only class that is genuinely underpowered.
  expect_equal(unname(tab["noisy", "undetermined"]), sum(gene_class == "noisy"))

  # every gene lands in exactly one class
  expect_equal(sum(tab), length(gene_class))
  expect_false(any(is.na(partition)))
})
