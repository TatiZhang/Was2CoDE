context("Test was2code_dist")

# ─────────────────────────────────────────────────────────────────────────────
# Synthetic inputs for structural and distributional tests.
#
# 4 genes, 4 cases + 4 controls, 100 cells per donor:
#   null_gene     – all donors ~ N(0, 1)
#   location_gene – cases ~ N(5, 1),    controls ~ N(0, 1)         [mean shift]
#   size_gene     – cases ~ N(0, 3),    controls ~ N(0, 0.3)       [variance shift]
#   shape_gene    – cases ~ bimodal ±3 (sd=0.3), controls ~ N(0, 3.015)
#                   [shape shift; sd matched between groups]
# ─────────────────────────────────────────────────────────────────────────────

.make_test_data <- function(seed     = 42,
                            n_case   = 4,
                            n_ctrl   = 4,
                            n_cpd    = 100) {   # cells per donor
  set.seed(seed)
  n_donors  <- n_case + n_ctrl
  n_cells   <- n_donors * n_cpd
  donor_ids <- c(paste0("case", seq_len(n_case)), paste0("ctrl", seq_len(n_ctrl)))
  cell_ids  <- paste0("cell", seq_len(n_cells))

  gene_names <- c("null_gene", "location_gene", "size_gene", "shape_gene")
  mat <- matrix(NA_real_, nrow = 4, ncol = n_cells,
                dimnames = list(gene_names, cell_ids))

  # sd of bimodal ±3, sd=0.3: sqrt(9 + 0.09) = sqrt(9.09) ≈ 3.015
  bimodal_sd <- sqrt(9.09)

  for (i in seq_len(n_donors)) {
    idx     <- seq((i - 1) * n_cpd + 1, i * n_cpd)
    is_case <- i <= n_case
    signs   <- sample(c(-1L, 1L), n_cpd, replace = TRUE)

    mat["null_gene",     idx] <- rnorm(n_cpd, 0, 1)
    mat["location_gene", idx] <- rnorm(n_cpd, if (is_case) 5 else 0, 1)
    mat["size_gene",     idx] <- rnorm(n_cpd, 0, if (is_case) 3 else 0.3)
    mat["shape_gene",    idx] <- if (is_case) rnorm(n_cpd, signs * 3, 0.3)   # bimodal ±3
                                 else         rnorm(n_cpd, 0, bimodal_sd)     # normal, matched sd
  }

  list(
    mat       = mat,
    meta_cell = data.frame(cell_id    = cell_ids,
                           individual = rep(donor_ids, each = n_cpd),
                           stringsAsFactors = FALSE),
    meta_ind  = data.frame(individual = donor_ids,
                           group      = factor(c(rep("case", n_case),
                                                 rep("ctrl", n_ctrl))),
                           stringsAsFactors = FALSE),
    case_ids  = donor_ids[seq_len(n_case)],
    ctrl_ids  = donor_ids[n_case + seq_len(n_ctrl)]
  )
}

td     <- .make_test_data()
result <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                        var2test = "group", ncores = 1)

# Mean absolute value of all entries in a submatrix; diagonal excluded when rows == cols.
.mean_abs <- function(arr, row_ids, col_ids, metric) {
  sub <- arr[row_ids, col_ids, metric]
  if (identical(as.character(row_ids), as.character(col_ids))) diag(sub) <- NA
  mean(abs(sub), na.rm = TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Output structure
# ─────────────────────────────────────────────────────────────────────────────

test_that("output is a named list with one array per gene", {
  expect_type(result, "list")
  expect_length(result, nrow(td$mat))
  expect_named(result, rownames(td$mat))
})

test_that("each gene's array has the correct dimensions and dimnames", {
  n_donors      <- nrow(td$meta_ind)
  expected_3rd  <- c("was2", "location", "size", "shape")

  for (g in result) {
    expect_true(is.array(g))
    expect_equal(length(dim(g)), 3L)
    expect_equal(dim(g), c(n_donors, n_donors, 4L))
    expect_equal(dimnames(g)[[1]], td$meta_ind$individual)
    expect_equal(dimnames(g)[[2]], td$meta_ind$individual)
    expect_equal(dimnames(g)[[3]], expected_3rd)
  }
})

test_that("diagonal entries are zero for every metric", {
  for (g in result) {
    for (metric in c("was2", "location", "size", "shape")) {
      expect_equal(unname(diag(g[,, metric])),
                   rep(0, nrow(td$meta_ind)),
                   label = sprintf("diagonal of %s", metric))
    }
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# 2. Distance properties
# ─────────────────────────────────────────────────────────────────────────────

test_that("was2 and shape are symmetric and non-negative", {
  n <- nrow(td$meta_ind)
  for (g in result) {
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        expect_gte(g[i, j, "was2"],  0)
        expect_gte(g[i, j, "shape"], 0)
        expect_equal(g[i, j, "was2"],  g[j, i, "was2"],  tolerance = 1e-10)
        expect_equal(g[i, j, "shape"], g[j, i, "shape"], tolerance = 1e-10)
      }
    }
  }
})

test_that("location and size are anti-symmetric", {
  n <- nrow(td$meta_ind)
  for (g in result) {
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        expect_equal(g[i, j, "location"], -g[j, i, "location"], tolerance = 1e-10)
        expect_equal(g[i, j, "size"],     -g[j, i, "size"],     tolerance = 1e-10)
      }
    }
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# 3. Input validation errors
# ─────────────────────────────────────────────────────────────────────────────

# Minimal valid inputs reused across validation tests
.min_mat <- function() {
  matrix(rnorm(10), nrow = 2, ncol = 5,
         dimnames = list(c("g1", "g2"), paste0("c", 1:5)))
}
.min_mc <- function() {
  data.frame(cell_id    = paste0("c", 1:5),
             individual = rep(c("d1", "d2"), c(3, 2)),
             stringsAsFactors = FALSE)
}
.min_mi <- function() {
  data.frame(individual = c("d1", "d2"),
             group      = factor(c("case", "ctrl")),
             stringsAsFactors = FALSE)
}
.run <- function(mat = .min_mat(), mc = .min_mc(), mi = .min_mi()) {
  was2code_dist(mat, mc, mi, var2test = "group", ncores = 1)
}

test_that("error when count_input is not a matrix", {
  expect_error(.run(mat = as.data.frame(.min_mat())))
})

test_that("error when count_input has no row names", {
  m <- .min_mat(); rownames(m) <- NULL
  expect_error(.run(mat = m))
})

test_that("error when count_input has no column names", {
  m <- .min_mat(); colnames(m) <- NULL
  expect_error(.run(mat = m))
})

test_that("error when count_input has duplicate row names", {
  m <- .min_mat(); rownames(m) <- c("g1", "g1")
  expect_error(.run(mat = m))
})

test_that("error when count_input has duplicate column names", {
  m <- .min_mat(); colnames(m) <- c("c1", "c1", "c3", "c4", "c5")
  expect_error(.run(mat = m))
})

test_that("error when meta_cell is not a data.frame", {
  expect_error(.run(mc = list(cell_id = paste0("c", 1:5), individual = rep(c("d1","d2"), c(3,2)))))
})

test_that("error when meta_cell is missing required columns", {
  mc <- .min_mc(); mc$cell_id <- NULL
  expect_error(.run(mc = mc))

  mc <- .min_mc(); mc$individual <- NULL
  expect_error(.run(mc = mc))
})

test_that("error when meta_cell cell_id does not match count_input column names", {
  mc <- .min_mc(); mc$cell_id[1] <- "wrong"
  expect_error(.run(mc = mc))
})

test_that("error when meta_cell has duplicate cell IDs", {
  mc <- .min_mc(); mc$cell_id[2] <- mc$cell_id[1]
  expect_error(.run(mc = mc))
})

test_that("error when meta_ind is not a data.frame", {
  expect_error(.run(mi = list(individual = c("d1","d2"), group = c("case","ctrl"))))
})

test_that("error when meta_ind is missing required columns", {
  mi <- .min_mi(); mi$individual <- NULL
  expect_error(.run(mi = mi))

  mi <- .min_mi(); mi$group <- NULL
  expect_error(.run(mi = mi))
})

test_that("error when individual IDs in meta_cell and meta_ind do not match", {
  mi <- .min_mi(); mi$individual <- c("x1", "x2")
  expect_error(.run(mi = mi))
})

test_that("error when meta_ind has duplicate individual IDs", {
  mi <- .min_mi(); mi$individual <- c("d1", "d1")
  expect_error(.run(mi = mi))
})

test_that("error when var2test has fewer than two levels", {
  mi <- .min_mi(); mi$group <- factor(rep("case", 2))
  expect_error(.run(mi = mi))
})

# ─────────────────────────────────────────────────────────────────────────────
# 4. Distributional sensitivity
#
# For each signal gene, the relevant between-group metric must be substantially
# larger than the within-group metric (threshold: 5x for location/size, 3x for
# shape, since bimodal-vs-normal quantile correlation is less extreme).
# ─────────────────────────────────────────────────────────────────────────────

test_that("location differences: between-group location distance >> within-group", {
  g <- result[["location_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "location")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "location") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "location")) / 2

  expect_gt(between, within * 5)
})

test_that("location differences: between-group was2 >> within-group", {
  g <- result[["location_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "was2")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "was2") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "was2")) / 2

  expect_gt(between, within * 5)
})

test_that("variance differences: between-group size distance >> within-group", {
  g <- result[["size_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "size")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "size") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "size")) / 2

  expect_gt(between, within * 5)
})

test_that("variance differences: between-group was2 >> within-group", {
  g <- result[["size_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "was2")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "was2") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "was2")) / 2

  expect_gt(between, within * 5)
})

test_that("shape differences: between-group shape distance >> within-group", {
  g <- result[["shape_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "shape")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "shape") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "shape")) / 2

  expect_gt(between, within * 2)
})

test_that("shape differences: between-group was2 >> within-group", {
  g <- result[["shape_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "was2")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "was2") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "was2")) / 2

  expect_gt(between, within * 2)
})

test_that("metric isolation: location_gene shows large location but small size and shape", {
  g <- result[["location_gene"]]

  loc_between  <- .mean_abs(g, td$case_ids, td$ctrl_ids, "location")
  size_between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "size")

  # Both groups have sd=1, so size should be near zero
  expect_gt(loc_between, size_between * 5)
})

test_that("metric isolation: size_gene shows large size but small location", {
  g <- result[["size_gene"]]

  size_between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "size")
  loc_between  <- .mean_abs(g, td$case_ids, td$ctrl_ids, "location")

  # Both groups are mean-zero, so location should be near zero
  expect_gt(size_between, loc_between * 5)
})

test_that("metric isolation: shape_gene shows large shape but small location and size", {
  g <- result[["shape_gene"]]

  shape_between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "shape")
  loc_between   <- .mean_abs(g, td$case_ids, td$ctrl_ids, "location")
  size_between  <- .mean_abs(g, td$case_ids, td$ctrl_ids, "size")

  # Mean and sd are matched between groups, so location and size should be small
  expect_gt(shape_between, loc_between  * 3)
  expect_gt(shape_between, size_between * 3)
})

test_that("null gene: between-group and within-group was2 are of comparable magnitude", {
  g <- result[["null_gene"]]

  between <- .mean_abs(g, td$case_ids, td$ctrl_ids, "was2")
  within  <- (.mean_abs(g, td$case_ids, td$case_ids, "was2") +
              .mean_abs(g, td$ctrl_ids, td$ctrl_ids, "was2")) / 2

  # Neither should be large in absolute terms (N(0,1) empirical W2 ≈ 0.1–0.3)
  expect_lt(between, 1.0)
  expect_lt(within,  1.0)
  # And between should not dwarf within
  expect_lt(between, within * 5)
})

# ─────────────────────────────────────────────────────────────────────────────
# 5. k parameter
# ─────────────────────────────────────────────────────────────────────────────

test_that("k=NULL (default): no off-diagonal NAs", {
  for (g in result) {
    w2 <- g[,, "was2"]
    off_diag <- w2[row(w2) != col(w2)]
    expect_false(any(is.na(off_diag)),
                 info = "all off-diagonal entries should be non-NA when k=NULL")
  }
})

test_that("k=1: matrix contains NAs (not all pairs compared)", {
  result_k1 <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                             var2test = "group", ncores = 1, k = 1)

  for (g in result_k1) {
    w2 <- g[,, "was2"]
    off_diag <- w2[row(w2) != col(w2)]
    expect_true(any(is.na(off_diag)),
                info = "k=1 should leave some off-diagonal entries as NA")
  }
})

test_that("k=1: fewer filled entries than k=NULL", {
  result_k1 <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                             var2test = "group", ncores = 1, k = 1)

  for (gene in rownames(td$mat)) {
    n_filled_full <- sum(!is.na(result[[gene]][,, "was2"]) &
                           row(result[[gene]][,, "was2"]) != col(result[[gene]][,, "was2"]))
    n_filled_k1   <- sum(!is.na(result_k1[[gene]][,, "was2"]) &
                           row(result_k1[[gene]][,, "was2"]) != col(result_k1[[gene]][,, "was2"]))
    expect_lt(n_filled_k1, n_filled_full,
              label = sprintf("gene %s should have fewer filled pairs with k=1", gene))
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# 6. Reproducibility: ncores=1 and ncores=2 give identical results
# ─────────────────────────────────────────────────────────────────────────────

test_that("results are identical for ncores=1 and ncores=2", {
  res1 <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                        var2test = "group", ncores = 1)
  res2 <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                        var2test = "group", ncores = 2)

  for (gene in names(res1)) {
    expect_equal(res1[[gene]], res2[[gene]], tolerance = 1e-10,
                 label = sprintf("gene %s", gene))
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# 7. Parallelization timing
#
# Uses 20 genes x 8 donors x 200 cells/donor so there is enough per-gene work
# for mclapply to overcome fork overhead.  Skipped on CI and Windows (where
# the parallel path is not taken) and when fewer than 2 cores are available.
# ─────────────────────────────────────────────────────────────────────────────

test_that("ncores=2 is faster than ncores=1 on a moderately large dataset", {
  skip_on_ci()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "fewer than 2 cores available")

  set.seed(99)
  n_donors <- 8L; n_cpd <- 200L; n_genes <- 20L
  donor_ids <- paste0("d", seq_len(n_donors))
  cell_ids  <- paste0("c", seq_len(n_donors * n_cpd))

  mat <- matrix(rnorm(n_genes * n_donors * n_cpd), nrow = n_genes,
                dimnames = list(paste0("g", seq_len(n_genes)), cell_ids))
  mc  <- data.frame(cell_id    = cell_ids,
                    individual = rep(donor_ids, each = n_cpd),
                    stringsAsFactors = FALSE)
  mi  <- data.frame(individual = donor_ids,
                    group      = factor(rep(c("case", "ctrl"), each = n_donors / 2L)),
                    stringsAsFactors = FALSE)

  t_serial   <- system.time(
    was2code_dist(mat, mc, mi, var2test = "group", ncores = 1)
  )[["elapsed"]]
  t_parallel <- system.time(
    was2code_dist(mat, mc, mi, var2test = "group", ncores = 2)
  )[["elapsed"]]

  expect_lt(t_parallel, t_serial * 0.85)
})
