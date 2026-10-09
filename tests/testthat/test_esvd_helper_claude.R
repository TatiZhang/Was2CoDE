context("Test esvd_helper")

# Shared fixture: the simulated donors of helper_simulate_pseudobulk.R. Small
# (4 + 4 donors, 30 cells each, 30 genes) and a low-rank, few-iteration fit,
# because every test here is about the wrapper's contract with eSVD2, not
# about the quality of the fit. `diagnosis` is a character column in the
# fixture, which is the form that exposed Bug E1.
.simulate_esvd_seurat <- function(seed = 10, n_donor_per_group = 4,
                                  n_cell_per_donor = 30){
  .simulate_donor_seurat(n_donor_per_group = n_donor_per_group,
                         n_cell_per_donor = n_cell_per_donor,
                         n_gene = c(de = 4, null = 20, noisy = 6),
                         seed = seed)
}

.run_esvd <- function(seurat_obj, max_iter = 3, ...){
  suppressWarnings(suppressMessages(
    esvd_helper(batch_var_prefix = NULL,
                case_control_levels = c("control", "case"),
                case_control_var = "diagnosis",
                categorical_vars = "sex",
                id_var = "donor",
                numerical_vars = "age",
                seurat_obj = seurat_obj,
                k = 2,
                max_iter = max_iter,
                ...)))$original
}

# Rebuild the Seurat object with one gene's counts zeroed, so the object has
# an all-zero gene the way a cell-type or donor subset of a real dataset does.
.zero_one_gene <- function(seurat_obj, gene){
  count_mat <- SeuratObject::LayerData(seurat_obj, layer = "counts")
  count_mat[gene, ] <- 0
  Seurat::CreateSeuratObject(counts = count_mat,
                             meta.data = seurat_obj@meta.data)
}

.skip_unless_esvd2 <- function(){
  skip_if_not_installed("eSVD2", minimum_version = "1.2.0")
  skip_if_not_installed("Seurat")
}

## Bug E1: `esvd_helper()` opened with `droplevels()` on the case/control
## column, which errors on a character column ("no applicable method for
## 'droplevels'"). The fixture's `diagnosis` is character, as it is in the
## analysis repo's downsampling scripts, so this is the default call path.
test_that("esvd_helper accepts a character case/control column and returns an eSVD object", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat()
  skip_if(is.null(sim))
  expect_true(is.character(sim$seurat_obj@meta.data$diagnosis))

  res <- .run_esvd(sim$seurat_obj)

  expect_true(inherits(res, "eSVD"))
  expect_true("gene_status" %in% names(res))
  expect_true(all(res$gene_status == "analyzed"))
  expect_equal(names(res$pvalue_list$fdr_vec), rownames(sim$seurat_obj))
})

## Bug E2: eSVD2 >= 1.2.0 refuses an all-zero gene with an error, and the old
## wrapper called `eSVD()` directly without removing any. A donor subset of a
## real dataset almost always has such genes, so every iteration of the
## downsampling scripts would have landed in their `tryCatch` as "eSVD FAILED".
## The invariant: the all-zero gene is present in the output, flagged, with NA
## statistics and an FDR of exactly 1, so `sum(fdr_vec < 0.05)` is unaffected.
test_that("esvd_helper removes an all-zero gene and reinserts it with NA statistics", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat()
  skip_if(is.null(sim))
  zero_gene <- "noisy-001"
  seurat_obj <- .zero_one_gene(sim$seurat_obj, zero_gene)

  res <- .run_esvd(seurat_obj)

  expect_true(inherits(res, "eSVD"))
  expect_equal(as.character(res$gene_status[zero_gene]), "all_zero")
  expect_equal(sum(res$gene_status == "all_zero"), 1)
  expect_equal(names(res$pvalue_list$fdr_vec), rownames(seurat_obj))
  expect_true(is.na(res$teststat_vec[zero_gene]))
  expect_equal(unname(res$pvalue_list$fdr_vec[zero_gene]), 1)
  expect_false(anyNA(res$pvalue_list$fdr_vec))
  analyzed_vec <- names(res$gene_status)[res$gene_status == "analyzed"]
  expect_false(anyNA(res$teststat_vec[analyzed_vec]))
})

## Bug E3: the old wrapper checked only the pooled donor count (`min_ids`), so
## a cohort of 4 controls and 1 case passed its filter and then errored inside
## eSVD2 ("each arm needs at least 2 individuals"). The contract of every
## helper in this package is: an underpowered cohort is a warning plus NA,
## never an error, because the analysis scripts count NA as "did not run".
test_that("esvd_helper returns NA with a warning when an arm has one donor", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat()
  skip_if(is.null(sim))
  keep_vec <- sim$seurat_obj@meta.data$donor %in% c("D01", "D02", "D03", "D04",
                                                    "D05")
  seurat_obj <- sim$seurat_obj[, keep_vec]
  expect_equal(sort(as.integer(table(seurat_obj@meta.data$diagnosis))),
               c(30L, 120L))

  expect_warning(res <- suppressMessages(
    esvd_helper(batch_var_prefix = NULL,
                case_control_levels = c("control", "case"),
                case_control_var = "diagnosis",
                categorical_vars = "sex",
                id_var = "donor",
                numerical_vars = "age",
                seurat_obj = seurat_obj,
                k = 2,
                max_iter = 3)),
    regexp = "min_ids_per_arm")
  expect_true(length(res$original) == 1 && is.na(res$original))
})

## `bool_check_donors = FALSE` is the only switch the wrapper adds on top of
## eSVD2, and it works by zeroing the thresholds. Pin that it really is an off
## switch: a cohort of 2 + 2 donors fails the default pooled `min_ids = 4`
## (rejection is at `<=`), and runs when the checks are off.
##
## The one threshold it must NOT zero is `min_ids_per_arm`: eSVD2 cannot fit a
## one-donor arm and errors on it, so with the per-arm filter off the wrapper
## would error where every other helper warns and returns NA (code review,
## 2026-10-09). Pin that a 4-vs-1 cohort still comes back NA with the checks
## off.
test_that("bool_check_donors = FALSE disables every cohort filter", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat(n_donor_per_group = 2, n_cell_per_donor = 30)
  skip_if(is.null(sim))
  expect_equal(length(unique(sim$seurat_obj@meta.data$donor)), 4)

  res_checked <- .run_esvd(sim$seurat_obj)
  expect_true(length(res_checked) == 1 && is.na(res_checked))

  res_unchecked <- .run_esvd(sim$seurat_obj, bool_check_donors = FALSE)
  expect_true(inherits(res_unchecked, "eSVD"))
  expect_equal(length(res_unchecked$teststat_vec), nrow(sim$seurat_obj))

  sim_full <- .simulate_esvd_seurat()
  keep_vec <- sim_full$seurat_obj@meta.data$donor %in% c("D01", "D02", "D03",
                                                         "D04", "D05")
  seurat_one_arm <- sim_full$seurat_obj[, keep_vec]
  expect_warning(res_one_arm <- suppressMessages(
    esvd_helper(batch_var_prefix = NULL,
                case_control_levels = c("control", "case"),
                case_control_var = "diagnosis",
                categorical_vars = "sex",
                id_var = "donor",
                numerical_vars = "age",
                seurat_obj = seurat_one_arm,
                bool_check_donors = FALSE,
                k = 2,
                max_iter = 3)),
    regexp = "min_ids_per_arm")
  expect_true(length(res_one_arm$original) == 1 && is.na(res_one_arm$original))
})

## The wrapper must add nothing of its own: the same call through
## `eSVD2::eSVD_helper()` gives the identical object. The fit is deterministic
## in eSVD2 >= 1.2.0 (fixed start vectors), so this is an exact comparison.
test_that("esvd_helper is a pure pass-through to eSVD2::eSVD_helper", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat()
  skip_if(is.null(sim))
  seurat_obj <- .zero_one_gene(sim$seurat_obj, "noisy-002")

  res_wrapper <- .run_esvd(seurat_obj)
  res_direct <- suppressWarnings(suppressMessages(
    eSVD2::eSVD_helper(batch_var_prefix = NULL,
                       case_control_levels = c("control", "case"),
                       case_control_var = "diagnosis",
                       categorical_vars = "sex",
                       id_var = "donor",
                       numerical_vars = "age",
                       seurat_obj = seurat_obj,
                       k = 2,
                       max_iter = 3)))

  expect_equal(res_wrapper$teststat_vec, res_direct$teststat_vec)
  expect_equal(res_wrapper$pvalue_list$fdr_vec, res_direct$pvalue_list$fdr_vec)
  expect_equal(res_wrapper$log2fc_vec, res_direct$log2fc_vec)
  expect_equal(res_wrapper$gene_status, res_direct$gene_status)
})

## Orientation, the same face as the dreamlet sign-flip bug: the contrast is
## fixed by `case_control_levels` (control first, then case), so a simulated
## up-regulated gene must come back with a positive log fold change. The
## simulated DE genes alternate in sign, so a flip cannot hide behind a
## one-directional panel. Two seeds, since one fit could orient by luck.
test_that("esvd_helper orients the log fold change case-vs-control", {
  .skip_unless_esvd2()
  for(seed in c(10, 11)){
    label <- paste0("seed ", seed)
    sim <- .simulate_esvd_seurat(seed = seed)
    skip_if(is.null(sim))

    res <- .run_esvd(sim$seurat_obj, max_iter = 10)
    report_df <- suppressWarnings(eSVD2::report_results(res))
    rownames(report_df) <- report_df$genes
    de_vec <- names(sim$gene_class)[sim$gene_class == "de"]

    expect_equal(sign(report_df[de_vec, "logFC"]),
                 unname(sign(sim$true_lfc[de_vec])),
                 info = label)
    expect_gt(stats::median(abs(report_df[de_vec, "logFC"])), 1)
  }
})

## The unified `results` table re-labels `eSVD2::report_results()`, with one
## deliberate change: eSVD2 pads an all-zero gene with p = 1 and FDR = 1, but
## its BH counted only the analyzed genes. The unified table gives that gene
## NA instead, so plain BH over the non-NA p-values reproduces eSVD2's own
## `fdr_vec` exactly on the analyzed genes.
test_that("esvd_helper's unified results match report_results and eSVD2's FDR", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat()
  skip_if(is.null(sim))
  seurat_obj <- .zero_one_gene(sim$seurat_obj, "noisy-002")

  res <- suppressWarnings(suppressMessages(
    esvd_helper(batch_var_prefix = NULL,
                case_control_levels = c("control", "case"),
                case_control_var = "diagnosis",
                categorical_vars = "sex",
                id_var = "donor",
                numerical_vars = "age",
                seurat_obj = seurat_obj,
                k = 2,
                max_iter = 3)))
  expect_equal(names(res), c("results", "original"))
  expect_true(inherits(res$original, "eSVD"))

  results_df <- res$results
  report_df <- suppressWarnings(eSVD2::report_results(res$original))
  expect_equal(colnames(results_df), c("gene", "logFC", "se", "pvalue", "padj"))
  expect_equal(results_df$gene, rownames(seurat_obj))
  expect_equal(results_df$logFC, report_df$logFC)
  expect_equal(results_df$se, report_df$logFC_se)

  analyzed_vec <- setdiff(results_df$gene, "noisy-002")
  expect_equal(results_df[analyzed_vec, "pvalue"],
               unname(report_df[analyzed_vec, "pvalue"]))
  expect_equal(results_df[analyzed_vec, "padj"],
               unname(res$original$pvalue_list$fdr_vec[analyzed_vec]))
  expect_true(all(is.na(results_df["noisy-002",
                                   c("logFC", "se", "pvalue", "padj")])))
})

## A rejected cohort returns NA from eSVD2. The unified table still has a row
## per gene, all NA, so a caller can stack the four methods without a special
## case for eSVD.
test_that("esvd_helper returns an all-NA unified table for a rejected cohort", {
  .skip_unless_esvd2()
  sim <- .simulate_esvd_seurat(n_donor_per_group = 2, n_cell_per_donor = 30)
  skip_if(is.null(sim))

  res <- suppressWarnings(suppressMessages(
    esvd_helper(batch_var_prefix = NULL,
                case_control_levels = c("control", "case"),
                case_control_var = "diagnosis",
                categorical_vars = "sex",
                id_var = "donor",
                numerical_vars = "age",
                seurat_obj = sim$seurat_obj,
                k = 2,
                max_iter = 3)))

  expect_true(length(res$original) == 1 && is.na(res$original))
  expect_equal(res$results$gene, rownames(sim$seurat_obj))
  expect_true(all(is.na(res$results[, c("logFC", "se", "pvalue", "padj")])))
})
