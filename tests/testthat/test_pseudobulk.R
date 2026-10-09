context("Test pseudobulk construction")

# ---------------------------------------------------------------------------
# .check_donor_vars
# ---------------------------------------------------------------------------

test_that(".check_donor_vars accepts valid argument sets", {
  expect_true(.check_donor_vars(case_control_levels = c("control", "case"),
                                case_control_var = "diagnosis",
                                categorical_vars = NULL,
                                id_var = "donor",
                                numerical_vars = NULL))
  expect_true(.check_donor_vars(case_control_levels = c("control", "case"),
                                case_control_var = "diagnosis",
                                categorical_vars = c("sex", "batch"),
                                id_var = "donor",
                                numerical_vars = "age"))
})

test_that(".check_donor_vars rejects malformed argument sets", {
  ok <- list(case_control_levels = c("control", "case"),
             case_control_var = "diagnosis",
             categorical_vars = NULL,
             id_var = "donor",
             numerical_vars = NULL)

  # three levels instead of two
  expect_error(do.call(.check_donor_vars,
                       utils::modifyList(ok, list(case_control_levels = c("a", "b", "c")))))
  # non-character level names
  expect_error(do.call(.check_donor_vars,
                       utils::modifyList(ok, list(case_control_levels = c(0, 1)))))
  # more than one case/control column
  expect_error(do.call(.check_donor_vars,
                       utils::modifyList(ok, list(case_control_var = c("a", "b")))))
  # duplicated covariate names would duplicate design-matrix columns
  expect_error(do.call(.check_donor_vars,
                       utils::modifyList(ok, list(categorical_vars = c("sex", "sex")))))
  expect_error(do.call(.check_donor_vars,
                       utils::modifyList(ok, list(numerical_vars = c("age", "age")))))
})

# ---------------------------------------------------------------------------
# .pseudobulk_from_seurat
# ---------------------------------------------------------------------------

test_that(".pseudobulk_from_seurat sums each donor's cells exactly", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 40,
                               n_gene = c(de = 3, null = 5, noisy = 2), seed = 30)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = NULL,
                                           id_var = "donor")
  mat <- as.matrix(SeuratObject::LayerData(pseudo_seurat, layer = "counts", assay = "RNA"))

  # one pseudobulk sample per donor
  expect_equal(ncol(mat), length(sim$donor_group))
  expect_equal(nrow(mat), length(sim$gene_class))

  # the aggregated counts are the per-donor row sums of the cell-level counts
  cell_counts <- SeuratObject::LayerData(sim$seurat_obj, layer = "counts")
  donor_of_cell <- sim$seurat_obj@meta.data$donor
  truth <- sapply(names(sim$donor_group), function(donor) {
    Matrix::rowSums(cell_counts[, donor_of_cell == donor, drop = FALSE])
  })
  # pseudobulk column names are "<donor>_<diagnosis>"
  donor_of_column <- sub("_.*$", "", colnames(mat))
  expect_equal(unname(mat), unname(truth[rownames(mat), donor_of_column]))
})

test_that(".pseudobulk_from_seurat carries variable features across aggregation", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 2, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 31)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = NULL,
                                           id_var = "donor")

  expect_equal(Seurat::VariableFeatures(pseudo_seurat),
               Seurat::VariableFeatures(sim$seurat_obj))
})

test_that(".pseudobulk_from_seurat splits on categorical covariates too", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 32)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = "sex",
                                           id_var = "donor")

  # sex is donor-constant, so grouping on it must not create extra samples
  expect_equal(ncol(pseudo_seurat), length(sim$donor_group))
  expect_true("sex" %in% colnames(pseudo_seurat@meta.data))
})

# ---------------------------------------------------------------------------
# .attach_donor_level_numerics
# ---------------------------------------------------------------------------

test_that(".attach_donor_level_numerics copies each donor's value onto its sample", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 4, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 33)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = NULL,
                                           id_var = "donor")
  expect_false("age" %in% colnames(pseudo_seurat@meta.data))

  pseudo_seurat <- .attach_donor_level_numerics(pseudo_seurat = pseudo_seurat,
                                                seurat_obj = sim$seurat_obj,
                                                numerical_vars = "age",
                                                id_var = "donor")

  expect_true("age" %in% colnames(pseudo_seurat@meta.data))
  expect_equal(pseudo_seurat@meta.data$age,
               unname(sim$donor_age[pseudo_seurat@meta.data$donor]))
})

test_that(".attach_donor_level_numerics errors on a covariate that varies within a donor", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 34)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = NULL,
                                           id_var = "donor")

  # n_umi is a cell-level quantity, so it is not a donor-level covariate
  expect_error(.attach_donor_level_numerics(pseudo_seurat = pseudo_seurat,
                                            seurat_obj = sim$seurat_obj,
                                            numerical_vars = "n_umi",
                                            id_var = "donor"),
               "violates the variable")
})

test_that(".attach_donor_level_numerics is a no-op when there are no numerical vars", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 2, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 35)
  skip_if(is.null(sim))

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = sim$seurat_obj,
                                           case_control_var = "diagnosis",
                                           categorical_vars = NULL,
                                           id_var = "donor")
  after <- .attach_donor_level_numerics(pseudo_seurat = pseudo_seurat,
                                        seurat_obj = sim$seurat_obj,
                                        numerical_vars = NULL,
                                        id_var = "donor")

  expect_equal(colnames(after@meta.data), colnames(pseudo_seurat@meta.data))
})

# ---------------------------------------------------------------------------
# .prepare_pseudobulk_metadata
# ---------------------------------------------------------------------------

test_that(".prepare_pseudobulk_metadata sets control as the reference level", {
  metadata <- data.frame(donor = paste0("D", 1:6),
                         diagnosis = c("case", "case", "case", "control", "control", "control"),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = NULL,
                                      id_var = "donor",
                                      numerical_vars = NULL)

  # the reference level is what makes the DESeq2 coefficient "case vs control"
  expect_true(is.factor(res$diagnosis))
  expect_equal(levels(res$diagnosis), c("control", "case"))
  # donor is the unit of observation, not a design variable
  expect_false("donor" %in% colnames(res))
})

test_that(".prepare_pseudobulk_metadata orders categorical levels by decreasing frequency", {
  metadata <- data.frame(donor = paste0("D", 1:6),
                         diagnosis = rep(c("control", "case"), each = 3),
                         batch = c("b1", "b2", "b2", "b2", "b3", "b3"),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = "batch",
                                      id_var = "donor",
                                      numerical_vars = NULL)

  # b2 appears 3x, b3 2x, b1 1x
  expect_equal(levels(res$batch), c("b2", "b3", "b1"))
})

test_that(".prepare_pseudobulk_metadata standardizes numerical covariates", {
  metadata <- data.frame(donor = paste0("D", 1:6),
                         diagnosis = rep(c("control", "case"), each = 3),
                         age = c(60, 65, 70, 75, 80, 85),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = NULL,
                                      id_var = "donor",
                                      numerical_vars = "age")

  expect_equal(mean(as.numeric(res$age)), 0)
  expect_equal(stats::sd(as.numeric(res$age)), 1)
})

test_that(".prepare_pseudobulk_metadata drops factors with no variation", {
  # a covariate constant across donors would make the design rank deficient
  metadata <- data.frame(donor = paste0("D", 1:6),
                         diagnosis = rep(c("control", "case"), each = 3),
                         sex = rep("F", 6),
                         batch = rep(c("b1", "b2"), 3),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = c("sex", "batch"),
                                      id_var = "donor",
                                      numerical_vars = NULL)

  expect_false("sex" %in% colnames(res))
  expect_true(all(c("batch", "diagnosis") %in% colnames(res)))
})

test_that(".prepare_pseudobulk_metadata keeps numerical covariates even if constant", {
  # only factors are checked for variation; a constant numeric passes through
  metadata <- data.frame(donor = paste0("D", 1:4),
                         diagnosis = rep(c("control", "case"), each = 2),
                         age = rep(70, 4),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = NULL,
                                      id_var = "donor",
                                      numerical_vars = "age")

  expect_true("age" %in% colnames(res))
})

test_that(".prepare_pseudobulk_metadata returns columns in design order", {
  metadata <- data.frame(donor = paste0("D", 1:6),
                         diagnosis = rep(c("control", "case"), each = 3),
                         batch = rep(c("b1", "b2", "b3"), 2),
                         age = c(60, 65, 70, 75, 80, 85),
                         stringsAsFactors = FALSE)

  res <- .prepare_pseudobulk_metadata(metadata_pseudobulk = metadata,
                                      case_control_levels = c("control", "case"),
                                      case_control_var = "diagnosis",
                                      categorical_vars = "batch",
                                      id_var = "donor",
                                      numerical_vars = "age")

  # covariates first, case/control last -- the coefficient of interest is the
  # last term of the design formula built from these column names
  expect_equal(colnames(res), c("batch", "age", "diagnosis"))
})

# ---------------------------------------------------------------------------
# build_pseudobulk
# ---------------------------------------------------------------------------

test_that("build_pseudobulk returns an aligned count matrix and design", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 4, n_cell_per_donor = 40,
                               n_gene = c(de = 3, null = 6, noisy = 2), seed = 36)
  skip_if(is.null(sim))

  res <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                          case_control_levels = c("control", "case"),
                          case_control_var = "diagnosis",
                          categorical_vars = "sex",
                          id_var = "donor",
                          numerical_vars = "age")

  expect_named(res, c("mat", "metadata", "pseudo_seurat"))
  expect_equal(ncol(res$mat), nrow(res$metadata))
  expect_equal(colnames(res$mat), rownames(res$metadata))
  expect_equal(nrow(res$mat), length(sim$gene_class))
  expect_equal(sort(rownames(res$mat)), sort(names(sim$gene_class)))

  # counts survive as non-negative whole numbers, as bulk NB models require
  mat <- as.matrix(res$mat)
  expect_true(all(mat >= 0))
  expect_true(all(mat == round(mat)))

  # the design carries the covariates, releveled and scaled
  expect_equal(colnames(res$metadata), c("sex", "age", "diagnosis"))
  expect_equal(levels(res$metadata$diagnosis), c("control", "case"))
  expect_equal(mean(as.numeric(res$metadata$age)), 0)
})

test_that("build_pseudobulk restricts to the variable features", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 2, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 6, noisy = 2), seed = 37)
  skip_if(is.null(sim))

  subset_genes <- names(sim$gene_class)[1:5]
  Seurat::VariableFeatures(sim$seurat_obj) <- subset_genes

  res <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                          case_control_levels = c("control", "case"),
                          case_control_var = "diagnosis",
                          categorical_vars = NULL,
                          id_var = "donor",
                          numerical_vars = NULL)

  expect_equal(sort(rownames(res$mat)), sort(subset_genes))
})

test_that("build_pseudobulk validates its arguments before doing any work", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 2, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 38)
  skip_if(is.null(sim))

  expect_error(build_pseudobulk(seurat_obj = sim$seurat_obj,
                                case_control_levels = "control",
                                case_control_var = "diagnosis",
                                categorical_vars = NULL,
                                id_var = "donor",
                                numerical_vars = NULL))
  expect_error(build_pseudobulk(seurat_obj = sim$seurat_obj,
                                case_control_levels = c("control", "case"),
                                case_control_var = "diagnosis",
                                categorical_vars = c("sex", "sex"),
                                id_var = "donor",
                                numerical_vars = NULL))
})

# ---------------------------------------------------------------------------
# .fit_deseq2_mle
# ---------------------------------------------------------------------------

test_that(".fit_deseq2_mle builds the design from the metadata columns", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 5, n_cell_per_donor = 60,
                                n_gene = c(de = 3, null = 10, noisy = 4), seed = 50)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = "sex",
                                 id_var = "donor",
                                 numerical_vars = "age")
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)

  expect_s4_class(dds, "DESeqDataSet")
  expect_equal(as.character(DESeq2::design(dds)), c("~", "sex + age + diagnosis"))
  # case/control is the last design term, so its coefficient is the one the
  # wrappers name as "<var>_<case>_vs_<control>"
  expect_true("diagnosis_case_vs_control" %in% DESeq2::resultsNames(dds))
})

test_that(".fit_deseq2_mle fits unshrunken coefficients", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 5, n_cell_per_donor = 60,
                                n_gene = c(de = 3, null = 10, noisy = 4), seed = 51)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)

  # betaPrior = FALSE is load-bearing: a zero-centered prior on the log2 fold
  # change would shrink estimates toward zero, which is anti-conservative for
  # the downstream equivalence test.  DESeq2 records the choice as an attribute
  # and refuses altHypothesis = "lessAbs" when it is TRUE.
  expect_false(attr(dds, "betaPrior"))
  expect_silent(DESeq2::results(dds, name = "diagnosis_case_vs_control",
                                lfcThreshold = 1, altHypothesis = "lessAbs"))
})

test_that(".fit_deseq2_mle propagates use_t as residual degrees of freedom", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 5, n_cell_per_donor = 60,
                                n_gene = c(de = 3, null = 10, noisy = 4), seed = 52)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)

  dds_norm <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata, use_t = FALSE)
  dds_t <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata, use_t = TRUE)

  expect_null(SummarizedExperiment::mcols(dds_norm)$tDegreesFreedom)
  # 10 pseudobulk samples minus the intercept and the case/control coefficient
  expect_equal(unique(SummarizedExperiment::mcols(dds_t)$tDegreesFreedom), 8)
})

# ---------------------------------------------------------------------------
# .extract_lfc_se
# ---------------------------------------------------------------------------

test_that(".extract_lfc_se returns DESeq2's unshrunken estimates verbatim", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 80,
                                n_gene = c(de = 4, null = 15, noisy = 6), seed = 53)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)
  est <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                         cooks_filter = FALSE)

  expect_equal(colnames(est),
               c("base_mean", "lfc", "se", "df", "cooks_outlier"))
  expect_equal(rownames(est), rownames(pseudobulk$mat))
  expect_true(all(is.na(est$df)))
  # the effect and its uncertainty only -- no p-values of any kind
  expect_false(any(grepl("pvalue|padj", colnames(est))))

  res <- DESeq2::results(dds, name = "diagnosis_case_vs_control",
                         independentFiltering = FALSE, cooksCutoff = FALSE)
  expect_equal(est$lfc, res$log2FoldChange)
  expect_equal(est$se, res$lfcSE)
  expect_equal(est$base_mean, res$baseMean)
})

test_that(".extract_lfc_se keeps low-expression genes rather than filtering them", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 80,
                                n_gene = c(de = 4, null = 15, noisy = 6), seed = 54)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)
  est <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                         cooks_filter = FALSE)

  # DESeq2's independent filter is tuned to maximize rejections of the
  # *difference* null, so it has no business deciding which genes are eligible
  # for some other test.  With it off, every gene survives with a usable
  # (lfc, se) however lowly expressed -- discarding is the caller's explicit
  # base_mean_min decision, made downstream.
  expect_equal(nrow(est), nrow(pseudobulk$mat))
  low <- est$base_mean < stats::quantile(est$base_mean, 0.25)
  expect_gt(sum(low), 0)
  expect_true(all(!is.na(est$lfc[low])))
  expect_true(all(!is.na(est$se[low])))
})

test_that(".extract_lfc_se NAs out genes below base_mean_min", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 80,
                                n_gene = c(de = 4, null = 15, noisy = 6), seed = 55)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)

  est_all <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                             cooks_filter = FALSE, base_mean_min = 0)
  est_filt <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                              cooks_filter = FALSE, base_mean_min = 100)

  low <- est_all$base_mean < 100
  expect_gt(sum(low), 0)
  expect_true(all(is.na(est_filt$lfc[low])))
  expect_true(all(is.na(est_filt$se[low])))
  # the estimates for retained genes are untouched, and base_mean is always kept
  expect_equal(est_filt$lfc[!low], est_all$lfc[!low])
  expect_equal(est_filt$base_mean, est_all$base_mean)
})

test_that(".extract_lfc_se propagates the Cook's-distance flag onto lfc and se", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 80,
                                n_gene = c(de = 4, null = 15, noisy = 6), seed = 56)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  # plant a single wild count in one donor of one gene to trip Cook's cutoff
  outlier_gene <- rownames(pseudobulk$mat)[1]
  pseudobulk$mat[outlier_gene, 1] <- as.integer(max(pseudobulk$mat) * 500)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata)

  est_on <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                            cooks_filter = TRUE)
  est_off <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                             cooks_filter = FALSE)

  # DESeq2 signals the flag by NA-ing the p-value while leaving lfc/se intact;
  # without the propagation an outlier-driven fit could support an inference
  flagged <- est_on$cooks_outlier
  expect_gt(sum(flagged), 0)
  expect_true(all(is.na(est_on$lfc[flagged])))
  expect_true(all(is.na(est_on$se[flagged])))
  # with the filter off, nothing is flagged and the estimates survive
  expect_false(any(est_off$cooks_outlier))
  expect_true(all(!is.na(est_off$lfc[flagged])))
})

test_that(".extract_lfc_se carries the t degrees of freedom when use_t was set", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 80,
                                n_gene = c(de = 4, null = 15, noisy = 6), seed = 57)
  skip_if(is.null(sim))

  pseudobulk <- build_pseudobulk(seurat_obj = sim$seurat_obj,
                                 case_control_levels = c("control", "case"),
                                 case_control_var = "diagnosis",
                                 categorical_vars = NULL,
                                 id_var = "donor",
                                 numerical_vars = NULL)
  dds <- .fit_deseq2_mle(pseudobulk$mat, pseudobulk$metadata, use_t = TRUE)
  est <- .extract_lfc_se(dds, coef_name = "diagnosis_case_vs_control",
                         cooks_filter = FALSE)

  # 12 pseudobulk samples minus the intercept and the case/control coefficient
  expect_equal(unique(est$df), 10)
})

# ---------------------------------------------------------------------------
# build_pseudobulk on degenerate covariates (2026-10-09)
# ---------------------------------------------------------------------------

## Bug P1: `Seurat::AggregateExpression()` silently drops a grouping variable
## with a single value ("will be ignored"), so a categorical covariate that is
## constant across the cohort (a single-sex cohort, or a batch that is constant
## within one cell type) never reached `.prepare_pseudobulk_metadata()`, which
## is the step meant to drop it from the design. Instead the column subset
## errored with "undefined columns selected". The unit test of
## `.prepare_pseudobulk_metadata()` passed because it built the metadata by
## hand, with the column present.
test_that("build_pseudobulk drops a categorical covariate that is constant across the cohort", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 37)
  skip_if(is.null(sim))
  obj <- sim$seurat_obj
  obj@meta.data$sex <- "F"

  res <- build_pseudobulk(seurat_obj = obj,
                          case_control_levels = c("control", "case"),
                          case_control_var = "diagnosis",
                          categorical_vars = "sex",
                          id_var = "donor",
                          numerical_vars = "age")

  expect_equal(colnames(res$metadata), c("age", "diagnosis"))
  expect_equal(ncol(res$mat), length(sim$donor_group))
})

## A donor in both arms is a metadata error. `group.by` would otherwise split
## that donor into a case sample and a control sample, and the model would
## treat the two halves as independent donors.
test_that("build_pseudobulk rejects a donor that appears in both arms", {
  skip_if_not_installed("Seurat")
  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 20,
                               n_gene = c(de = 2, null = 3, noisy = 1), seed = 38)
  skip_if(is.null(sim))
  obj <- sim$seurat_obj
  flip_idx <- which(obj@meta.data$donor == "D01")[1:5]
  obj@meta.data$diagnosis[flip_idx] <- "case"

  expect_error(build_pseudobulk(seurat_obj = obj,
                                case_control_levels = c("control", "case"),
                                case_control_var = "diagnosis",
                                categorical_vars = NULL,
                                id_var = "donor",
                                numerical_vars = NULL),
               "D01")
})
