context("Test dreamlet_helper")

# Shared fixture: the simulated donors of helper_simulate_pseudobulk.R, with
# the case/control column recoded by the caller. `diagnosis` is a character
# column in the fixture, which is the form that exposed the bug below.
.run_dreamlet <- function(seurat_obj,
                          case_control_levels = c("control", "case"),
                          categorical_vars = NULL,
                          numerical_vars = NULL){
  res <- suppressWarnings(suppressMessages(
    dreamlet_helper(case_control_levels = case_control_levels,
                    case_control_var = "diagnosis",
                    categorical_vars = categorical_vars,
                    id_var = "donor",
                    numerical_vars = numerical_vars,
                    seurat_obj = seurat_obj)))
  res <- as.data.frame(res)
  rownames(res) <- res$ID

  res
}

.de_genes <- function(sim, res){
  intersect(names(sim$gene_class)[sim$gene_class == "de"], rownames(res))
}

## Bug D1: `dreamlet_helper()` did not set the reference level of the
## case/control column and read the coefficient in position 2. A character
## column is ordered alphabetically by the model matrix, so "case" became the
## reference, the coefficient in position 2 was control-vs-case, and every log
## fold change came back with the wrong sign and no error.
##
## The expectation is the sign each gene was simulated with, which does not
## pass through the code under test. The simulated DE genes alternate in sign,
## so a sign flip cannot hide behind a one-directional panel.
test_that("dreamlet_helper orients the log fold change case-vs-control for a character column", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  seed_vec <- c(40, 41, 42)
  for(seed in seed_vec){
    label <- paste0("seed ", seed)
    sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                  n_gene = c(de = 6, null = 20, noisy = 10),
                                  seed = seed)
    skip_if(is.null(sim))
    expect_true(is.character(sim$seurat_obj@meta.data$diagnosis), info = label)

    res <- .run_dreamlet(sim$seurat_obj)
    de_vec <- .de_genes(sim, res)

    expect_gt(length(de_vec), 3)
    expect_equal(sign(res[de_vec, "logFC"]),
                 unname(sign(sim$true_lfc[de_vec])),
                 info = label)
    expect_gt(stats::median(abs(res[de_vec, "logFC"])), 1.5)
  }
})

## Bug D1, second face: the result depended on how the caller happened to
## store the column. The contrast is fixed by `case_control_levels`, so the
## three storage forms of the same data must give the same estimates.
test_that("dreamlet_helper gives the same estimates however the column is stored", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 6, null = 20, noisy = 10),
                                seed = 43)
  skip_if(is.null(sim))

  obj_character <- sim$seurat_obj
  obj_control_first <- sim$seurat_obj
  obj_control_first@meta.data$diagnosis <- factor(
    obj_control_first@meta.data$diagnosis, levels = c("control", "case"))
  obj_case_first <- sim$seurat_obj
  obj_case_first@meta.data$diagnosis <- factor(
    obj_case_first@meta.data$diagnosis, levels = c("case", "control"))

  res_character <- .run_dreamlet(obj_character)
  res_control_first <- .run_dreamlet(obj_control_first)
  res_case_first <- .run_dreamlet(obj_case_first)

  gene_vec <- sort(rownames(res_character))
  expect_equal(sort(rownames(res_control_first)), gene_vec)
  expect_equal(sort(rownames(res_case_first)), gene_vec)
  expect_equal(res_control_first[gene_vec, "logFC"],
               res_character[gene_vec, "logFC"],
               info = "control-first factor vs character")
  expect_equal(res_case_first[gene_vec, "logFC"],
               res_character[gene_vec, "logFC"],
               info = "case-first factor vs character")
})

## Bug D1 with covariates: position 2 stays the case/control term when
## covariates follow it in the formula, so the flip was the same. This pins
## the orientation when the coefficient is selected among several.
test_that("dreamlet_helper orients the log fold change with covariates in the model", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 6, null = 20, noisy = 10),
                                seed = 44)
  skip_if(is.null(sim))

  res <- .run_dreamlet(sim$seurat_obj,
                       categorical_vars = "sex",
                       numerical_vars = "age")
  de_vec <- .de_genes(sim, res)

  expect_equal(sign(res[de_vec, "logFC"]),
               unname(sign(sim$true_lfc[de_vec])))
})

## The contrast follows `case_control_levels`, not the alphabet: with levels
## named so that the case level sorts first ("affected" < "healthy"), and with
## the logical-like coding "FALSE"/"TRUE" that one of the cohorts uses.
test_that("dreamlet_helper follows case_control_levels for other codings of the arm", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  coding_list <- list(c(control = "healthy", case = "affected"),
                      c(control = "FALSE", case = "TRUE"),
                      c(control = "zz_ctrl", case = "aa_case"))
  for(coding in coding_list){
    label <- paste0("coding ", coding[["control"]], " / ", coding[["case"]])
    sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                  n_gene = c(de = 6, null = 20, noisy = 10),
                                  seed = 45)
    skip_if(is.null(sim))
    obj <- sim$seurat_obj
    obj@meta.data$diagnosis <- ifelse(obj@meta.data$diagnosis == "case",
                                      coding[["case"]], coding[["control"]])

    res <- .run_dreamlet(obj, case_control_levels = unname(coding))
    de_vec <- .de_genes(sim, res)

    expect_equal(sign(res[de_vec, "logFC"]),
                 unname(sign(sim$true_lfc[de_vec])),
                 info = label)
  }
})

## Swapping the two levels is the caller asking for control-vs-case, and the
## estimates must change sign with it. This is what makes `case_control_levels`
## the thing that fixes the contrast.
test_that("dreamlet_helper reverses the sign when the two levels are swapped", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 6, n_cell_per_donor = 100,
                                n_gene = c(de = 6, null = 20, noisy = 10),
                                seed = 46)
  skip_if(is.null(sim))

  res_forward <- .run_dreamlet(sim$seurat_obj,
                               case_control_levels = c("control", "case"))
  res_swapped <- .run_dreamlet(sim$seurat_obj,
                               case_control_levels = c("case", "control"))

  gene_vec <- sort(rownames(res_forward))
  expect_equal(res_swapped[gene_vec, "logFC"],
               -res_forward[gene_vec, "logFC"])
})

## A level that is not in the data used to fall through to whatever coefficient
## sat in position 2. It must stop, and name the argument to fix.
test_that("dreamlet_helper stops when a requested level is absent from the column", {
  skip_if_not_installed("dreamlet")
  skip_if_not_installed("Seurat")

  sim <- .simulate_donor_seurat(n_donor_per_group = 3, n_cell_per_donor = 30,
                                n_gene = c(de = 2, null = 5, noisy = 2),
                                seed = 47)
  skip_if(is.null(sim))

  expect_error(.run_dreamlet(sim$seurat_obj,
                             case_control_levels = c("control", "Case")),
               "case_control_levels")
  expect_error(.run_dreamlet(sim$seurat_obj,
                             case_control_levels = c("control", "control")),
               "case_control_levels")
})
