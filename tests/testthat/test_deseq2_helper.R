context("Test deseq2_helper")

test_that("deseq2_helper returns a DESeq2 results table over the variable features", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 5, null = 20, noisy = 10), seed = 40)
  skip_if(is.null(sim))

  res <- deseq2_helper(case_control_levels = c("control", "case"),
                       case_control_var = "diagnosis",
                       categorical_vars = NULL,
                       id_var = "donor",
                       numerical_vars = NULL,
                       seurat_obj = sim$seurat_obj)$original

  expect_true(all(c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj") %in% colnames(res)))
  expect_equal(sort(rownames(res)), sort(names(sim$gene_class)))

  # the coefficient is oriented case-vs-control, so every DE gene's estimated
  # fold change must carry the sign it was simulated with
  gene_class <- sim$gene_class[rownames(res)]
  is_de <- gene_class == "de"
  de_lfc <- res$log2FoldChange[is_de]
  expect_equal(sign(de_lfc), unname(sign(sim$true_lfc[rownames(res)][is_de])))
  expect_gt(stats::median(abs(de_lfc)), 1.5)
  expect_lt(stats::median(abs(res$log2FoldChange[gene_class == "null"])), 0.3)
})

test_that("deseq2_helper accepts categorical and donor-constant numerical covariates", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 5, null = 20, noisy = 10), seed = 41)
  skip_if(is.null(sim))

  res <- deseq2_helper(case_control_levels = c("control", "case"),
                       case_control_var = "diagnosis",
                       categorical_vars = "sex",
                       id_var = "donor",
                       numerical_vars = "age",
                       seurat_obj = sim$seurat_obj)$original

  expect_equal(nrow(res), length(sim$gene_class))
  expect_true(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
})

test_that("deseq2_helper rejects a covariate that is not donor-constant", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 30,
                                n_gene = c(de = 2, null = 5, noisy = 2), seed = 42)
  skip_if(is.null(sim))

  # n_umi varies cell to cell, so it cannot survive donor aggregation
  expect_error(deseq2_helper(case_control_levels = c("control", "case"),
                             case_control_var = "diagnosis",
                             categorical_vars = NULL,
                             id_var = "donor",
                             numerical_vars = "n_umi",
                             seurat_obj = sim$seurat_obj),
               "violates the variable")
})

## The unified `results` table must be a faithful re-labelling of DESeq2's own
## estimates, with `padj` replaced by plain BH. DESeq2's `padj` is BH after
## independent filtering, so it is NA for genes DESeq2 declined to adjust;
## the unified `padj` is defined for every gene that has a p-value.
test_that("deseq2_helper's unified results match the DESeq2 object", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 5, null = 20, noisy = 10), seed = 40)
  skip_if(is.null(sim))

  res <- deseq2_helper(case_control_levels = c("control", "case"),
                       case_control_var = "diagnosis",
                       categorical_vars = NULL,
                       id_var = "donor",
                       numerical_vars = NULL,
                       seurat_obj = sim$seurat_obj)
  expect_equal(names(res), c("results", "original"))
  expect_true(inherits(res$original, "DESeqResults"))

  results_df <- res$results
  original <- res$original
  expect_equal(colnames(results_df), c("gene", "logFC", "se", "pvalue", "padj"))
  expect_equal(results_df$gene, rownames(original))
  expect_equal(results_df$logFC, original$log2FoldChange)
  expect_equal(results_df$se, original$lfcSE)
  expect_equal(results_df$pvalue, original$pvalue)
  expect_equal(results_df$padj, stats::p.adjust(original$pvalue, method = "BH"))
  expect_equal(is.na(results_df$padj), is.na(results_df$pvalue))
})
