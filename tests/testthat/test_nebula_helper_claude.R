context("Test nebula_helper")

# Shared fixture: the simulated donors of helper_simulate_pseudobulk.R.
# `diagnosis` is a character column in the fixture, which is the form that
# exposed Bug N1. NEBULA fits every gene, so the panel is kept small.
.run_nebula <- function(seurat_obj,
                        case_control_levels = c("control", "case"),
                        categorical_vars = NULL,
                        numerical_vars = NULL){
  res <- NULL
  utils::capture.output(res <- suppressWarnings(suppressMessages(
    nebula_helper(case_control_levels = case_control_levels,
                  case_control_var = "diagnosis",
                  categorical_vars = categorical_vars,
                  id_var = "donor",
                  numerical_vars = numerical_vars,
                  seurat_obj = seurat_obj))))
  summary_df <- res$original$summary
  rownames(summary_df) <- summary_df$gene

  summary_df
}

.simulate_nebula_seurat <- function(seed){
  .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 60,
                         n_gene = c(de = 6, null = 20, noisy = 10),
                         seed = seed)
}

.de_genes_nebula <- function(sim, summary_df){
  intersect(names(sim$gene_class)[sim$gene_class == "de"], rownames(summary_df))
}

## Bug N1: `nebula_helper()` validated `case_control_levels` and then never
## used it. The model matrix took the reference level from the column as
## stored, so a character arm put "case" first alphabetically, the coefficient
## came back as `logFC_<var>control` with every sign flipped, and the
## case-named column the analysis scripts read did not exist. Same root as the
## dreamlet sign flip; here the exposure is the column name.
test_that("nebula_helper names the coefficient by the case level and orients it case-vs-control", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  for(seed in c(40, 41)){
    label <- paste0("seed ", seed)
    sim <- .simulate_nebula_seurat(seed)
    skip_if(is.null(sim))
    expect_true(is.character(sim$seurat_obj@meta.data$diagnosis), info = label)

    summary_df <- .run_nebula(sim$seurat_obj,
                              categorical_vars = "sex",
                              numerical_vars = "age")
    de_vec <- .de_genes_nebula(sim, summary_df)

    expect_true(all(c("logFC_diagnosiscase", "p_diagnosiscase") %in%
                      colnames(summary_df)), info = label)
    expect_false("logFC_diagnosiscontrol" %in% colnames(summary_df), info = label)
    expect_gt(length(de_vec), 3)
    expect_equal(sign(summary_df[de_vec, "logFC_diagnosiscase"]),
                 unname(sign(sim$true_lfc[de_vec])),
                 info = label)
  }
})

## Bug N1, second face: the three storage forms of the same arm must give the
## same column and the same estimates.
test_that("nebula_helper gives the same estimates however the arm is stored", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(42)
  skip_if(is.null(sim))
  obj_character <- sim$seurat_obj
  obj_control_first <- sim$seurat_obj
  obj_control_first@meta.data$diagnosis <- factor(
    obj_control_first@meta.data$diagnosis, levels = c("control", "case"))
  obj_case_first <- sim$seurat_obj
  obj_case_first@meta.data$diagnosis <- factor(
    obj_case_first@meta.data$diagnosis, levels = c("case", "control"))

  res_character <- .run_nebula(obj_character, categorical_vars = "sex")
  res_control_first <- .run_nebula(obj_control_first, categorical_vars = "sex")
  res_case_first <- .run_nebula(obj_case_first, categorical_vars = "sex")

  gene_vec <- sort(rownames(res_character))
  expect_equal(res_control_first[gene_vec, "logFC_diagnosiscase"],
               res_character[gene_vec, "logFC_diagnosiscase"],
               info = "control-first factor vs character")
  expect_equal(res_case_first[gene_vec, "logFC_diagnosiscase"],
               res_character[gene_vec, "logFC_diagnosiscase"],
               info = "case-first factor vs character")
})

## Swapping the levels is a request for control-vs-case: the column is then
## named by the new "case" level and the estimates change sign.
test_that("nebula_helper reverses the sign and renames the column when the levels are swapped", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(43)
  skip_if(is.null(sim))

  res_forward <- .run_nebula(sim$seurat_obj, categorical_vars = "sex")
  res_swapped <- .run_nebula(sim$seurat_obj, categorical_vars = "sex",
                             case_control_levels = c("case", "control"))

  gene_vec <- sort(rownames(res_forward))
  expect_true("logFC_diagnosiscontrol" %in% colnames(res_swapped))
  expect_equal(res_swapped[gene_vec, "logFC_diagnosiscontrol"],
               -res_forward[gene_vec, "logFC_diagnosiscase"])
})

## Bug N2: with no covariates, `pred` has one column, and reordering it with
## `[order_index, ]` dropped it to a vector, so `model.matrix()` failed. The
## helper could not be called without covariates at all.
test_that("nebula_helper runs with no covariates", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(44)
  skip_if(is.null(sim))

  summary_df <- .run_nebula(sim$seurat_obj)
  de_vec <- .de_genes_nebula(sim, summary_df)

  expect_equal(sort(grep("^logFC_", colnames(summary_df), value = TRUE)),
               c("logFC_(Intercept)", "logFC_diagnosiscase"))
  expect_equal(sign(summary_df[de_vec, "logFC_diagnosiscase"]),
               unname(sign(sim$true_lfc[de_vec])))
})

## Bug N3: a categorical covariate constant across the cohort (a single-sex
## cohort; a batch constant within one cell type) has one level, and
## `model.matrix()` refuses it ("contrasts can be applied only to factors with
## 2 or more levels"). `build_pseudobulk()` drops such a covariate from the
## design; this wrapper now does the same.
test_that("nebula_helper drops a categorical covariate that is constant across the cohort", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(45)
  skip_if(is.null(sim))
  obj <- sim$seurat_obj
  obj@meta.data$sex <- "F"

  summary_df <- .run_nebula(obj, categorical_vars = "sex", numerical_vars = "age")

  expect_false(any(grepl("^logFC_sex", colnames(summary_df))))
  expect_true("logFC_age" %in% colnames(summary_df))
  expect_true("logFC_diagnosiscase" %in% colnames(summary_df))
})

## NEBULA is a cell-level model with a donor random effect, so unlike the
## pseudobulk wrappers it may legitimately adjust for a cell-level covariate.
## Only the arm must be donor-constant. Pin both halves of that rule.
test_that("nebula_helper accepts a cell-level numerical covariate but not a donor in both arms", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(46)
  skip_if(is.null(sim))

  summary_df <- .run_nebula(sim$seurat_obj, numerical_vars = "n_umi")
  expect_true("logFC_n_umi" %in% colnames(summary_df))

  obj <- sim$seurat_obj
  flip_idx <- which(obj@meta.data$donor == "D01")[1:20]
  obj@meta.data$diagnosis[flip_idx] <- "case"
  expect_error(.run_nebula(obj), "D01")
})

## A level that is not in the data must stop and name the argument to fix.
test_that("nebula_helper stops when a requested level is absent from the column", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(47)
  skip_if(is.null(sim))

  expect_error(.run_nebula(sim$seurat_obj, case_control_levels = c("control", "Case")),
               "case_control_levels")
  expect_error(.run_nebula(sim$seurat_obj, case_control_levels = c("control", "control")),
               "case_control_levels")
})

## NEBULA reports natural-log coefficients, while DESeq2, dreamlet and eSVD2
## report log2, so the unified `logFC` and `se` are divided by log(2). Pin the
## conversion (a missed conversion would understate every NEBULA fold change
## by 31% beside the other methods) and that the p-value is untouched.
test_that("nebula_helper's unified results are the case coefficient on the log2 scale", {
  skip_if_not_installed("nebula")
  skip_if_not_installed("Seurat")

  sim <- .simulate_nebula_seurat(40)
  skip_if(is.null(sim))

  res <- NULL
  utils::capture.output(res <- suppressWarnings(suppressMessages(
    nebula_helper(case_control_levels = c("control", "case"),
                  case_control_var = "diagnosis",
                  categorical_vars = "sex",
                  id_var = "donor",
                  numerical_vars = NULL,
                  seurat_obj = sim$seurat_obj))))
  expect_equal(names(res), c("results", "original"))

  results_df <- res$results
  summary_df <- res$original$summary
  expect_equal(colnames(results_df), c("gene", "logFC", "se", "pvalue", "padj"))
  expect_equal(results_df$gene, summary_df$gene)
  expect_equal(results_df$logFC * log(2), summary_df$logFC_diagnosiscase)
  expect_equal(results_df$se * log(2), summary_df$se_diagnosiscase)
  expect_equal(results_df$pvalue, summary_df$p_diagnosiscase)
  expect_equal(results_df$padj,
               stats::p.adjust(summary_df$p_diagnosiscase, method = "BH"))

  de_vec <- intersect(names(sim$gene_class)[sim$gene_class == "de"],
                      results_df$gene)
  expect_equal(sign(results_df[de_vec, "logFC"]),
               unname(sign(sim$true_lfc[de_vec])))
})
