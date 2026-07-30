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
  expect_equal(res$pvalue_diff, pmin(1, 2 * stats::pnorm(abs(lfc) / se, lower.tail = FALSE)))
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
# equivalence_partition
# ---------------------------------------------------------------------------

test_that("equivalence_partition assigns the four classes correctly", {
  # already-adjusted-scale p-values: with Bonferroni-free "none" adjustment the
  # thresholding is exact
  pvalue_diff <- c(0.001, 0.001, 0.900, 0.900, NA)
  pvalue_tost <- c(0.900, 0.001, 0.001, 0.900, 0.001)

  res <- equivalence_partition(pvalue_diff, pvalue_tost,
                               alpha = 0.05, p_adjust_method = "none")

  expect_equal(as.character(res$class),
               c("differential", "trivial", "stable", "undetermined", NA))
  expect_equal(levels(res$class),
               c("differential", "trivial", "stable", "undetermined"))
  expect_equal(res$padj_diff, pvalue_diff)
})

test_that("equivalence_partition applies BH separately to each test", {
  set.seed(4)
  pvalue_diff <- stats::runif(100)
  pvalue_tost <- stats::runif(100)

  res <- equivalence_partition(pvalue_diff, pvalue_tost, alpha = 0.05)

  expect_equal(res$padj_diff, stats::p.adjust(pvalue_diff, method = "BH"))
  expect_equal(res$padj_tost, stats::p.adjust(pvalue_tost, method = "BH"))
  # BH is never anti-conservative relative to the raw p-values
  expect_true(all(res$padj_diff >= pvalue_diff))
})

test_that("equivalence_partition carries gene names through", {
  res <- equivalence_partition(c(0.5, 0.5), c(0.5, 0.5),
                               gene_names = c("GeneA", "GeneB"))
  expect_equal(rownames(res), c("GeneA", "GeneB"))
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
  # and the recomputed difference test matches DESeq2's Wald p-value
  expect_equal(ours$pvalue_diff, res_mle$pvalue[keep], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# End-to-end on a simulated Seurat object
# ---------------------------------------------------------------------------

test_that("deseq2_tost_helper partitions simulated genes into the right classes", {
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
  expect_true(all(c("lfc", "se", "pvalue_tost", "padj_tost", "class") %in% colnames(res)))
  expect_equal(attr(res, "lfc_threshold"), 0.5)

  gene_class <- sim$gene_class[rownames(res)]
  tab <- table(gene_class, res$class)

  # truly differential genes are found, and are never called stable
  expect_gt(tab["de", "differential"] / sum(gene_class == "de"), 0.9)
  expect_equal(unname(tab["de", "stable"]), 0)

  # well-expressed null genes are positively established as stable
  expect_gt(tab["null", "stable"] / sum(gene_class == "null"), 0.7)

  # barely-expressed null genes are undetermined, not stable: this is the
  # distinction the whole partition exists to make
  expect_gt(tab["noisy", "undetermined"] / sum(gene_class == "noisy"), 0.7)

  # the estimated fold changes point the right way
  expect_gt(stats::median(abs(res$lfc[gene_class == "de"])), 1.5)
  expect_lt(stats::median(abs(res$lfc[gene_class == "null"])), 0.3)

  # noisy genes have a much larger minimum rejectable bound than clean ones
  expect_gt(stats::median(res$lfc_threshold_min[gene_class == "noisy"], na.rm = TRUE),
            stats::median(res$lfc_threshold_min[gene_class == "null"], na.rm = TRUE))
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
  expect_true(all(res$class[gene_class == "de"] == "differential"))
  expect_true(all(res$class[gene_class == "noisy"] == "undetermined"))
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
