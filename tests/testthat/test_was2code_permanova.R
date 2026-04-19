context("Test was2code_permanova")

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

# Full (no NA) distance matrices with a configurable group effect.
.make_perm_data <- function(n_donors = 20, n_genes = 5, seed = 123,
                            between_dist = 1.5) {
  set.seed(seed)
  donors   <- paste0("ind", seq_len(n_donors))
  group    <- c(rep(0L, n_donors %/% 2), rep(1L, n_donors - n_donors %/% 2))
  meta_ind <- data.frame(individual = donors, group = group,
                         stringsAsFactors = FALSE)

  dist_list <- lapply(seq_len(n_genes), function(i) {
    arr <- array(0, dim = c(n_donors, n_donors, 1),
                 dimnames = list(donors, donors, "was2"))
    for (j in seq_len(n_donors)) for (k in seq(j, n_donors)) {
      d <- abs(rnorm(1, 2)) + ifelse(group[j] == group[k], 0, between_dist)
      arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
    }
    arr
  })
  names(dist_list) <- paste0("gene", seq_len(n_genes))
  list(dist_list = dist_list, meta_ind = meta_ind)
}

# Distance matrices with a controlled NA fraction.
.make_perm_data_na <- function(n_donors = 20, n_genes = 5, seed = 123,
                               na_frac = 0.2, between_dist = 1.5) {
  set.seed(seed)
  donors   <- paste0("ind", seq_len(n_donors))
  group    <- c(rep(0L, n_donors %/% 2), rep(1L, n_donors - n_donors %/% 2))
  meta_ind <- data.frame(individual = donors, group = group,
                         stringsAsFactors = FALSE)

  dist_list <- lapply(seq_len(n_genes), function(i) {
    arr <- array(NA_real_, dim = c(n_donors, n_donors, 1),
                 dimnames = list(donors, donors, "was2"))
    for (j in seq_len(n_donors)) for (k in seq(j, n_donors)) {
      if (runif(1) < na_frac) next
      d <- abs(rnorm(1, 2)) + ifelse(group[j] == group[k], 0, between_dist)
      arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
    }
    arr
  })
  names(dist_list) <- paste0("gene", seq_len(n_genes))
  list(dist_list = dist_list, meta_ind = meta_ind)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Basic functionality
# ─────────────────────────────────────────────────────────────────────────────

test_that("returns a list with pval, F_ob, and F_perm", {
  d      <- .make_perm_data()
  result <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group", n_perm = 99)

  expect_type(result, "list")
  expect_named(result, c("pval", "F_ob", "F_perm"))

  n_genes <- length(d$dist_list)
  expect_type(result$pval, "double")
  expect_length(result$pval, n_genes)
  expect_named(result$pval, names(d$dist_list))
  expect_true(all(result$pval >= 0 & result$pval <= 1))

  expect_type(result$F_ob, "double")
  expect_length(result$F_ob, n_genes)
  expect_named(result$F_ob, names(d$dist_list))

  expect_true(is.matrix(result$F_perm))
  expect_equal(dim(result$F_perm), c(n_genes, 99L))
  expect_equal(rownames(result$F_perm), names(d$dist_list))
})

test_that("works with a single gene", {
  set.seed(1)
  n      <- 10L
  donors <- paste0("ind", seq_len(n))
  group  <- c(rep(0L, 5L), rep(1L, 5L))
  arr    <- array(0, dim = c(n, n, 1), dimnames = list(donors, donors, "was2"))
  for (j in seq_len(n)) for (k in seq(j, n)) {
    d <- if (group[j] == group[k]) 0.1 else 5
    arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
  }
  result <- was2code_permanova(
    list(geneX = arr),
    data.frame(individual = donors, group = group, stringsAsFactors = FALSE),
    var2test = "group", n_perm = 99
  )
  expect_length(result$pval, 1L)
  expect_named(result$pval, "geneX")
  expect_gte(result$pval[["geneX"]], 0)
  expect_lte(result$pval[["geneX"]], 1)
})

test_that("result is reproducible with the same r_seed", {
  d  <- .make_perm_data()
  r1 <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                           n_perm = 99, r_seed = 42)
  r2 <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                           n_perm = 99, r_seed = 42)
  expect_equal(r1$pval,   r2$pval)
  expect_equal(r1$F_ob,   r2$F_ob)
  expect_equal(r1$F_perm, r2$F_perm)
})

# ─────────────────────────────────────────────────────────────────────────────
# 2. NA distance entries
# ─────────────────────────────────────────────────────────────────────────────

test_that("handles partial NA distance matrices without dropping genes", {
  d      <- .make_perm_data_na()
  result <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group", n_perm = 99)

  expect_length(result$pval, length(d$dist_list))
  expect_named(result$pval, names(d$dist_list))
  expect_true(all(result$pval >= 0 & result$pval <= 1 | is.na(result$pval)))
})

test_that("power is retained with a small NA fraction", {
  set.seed(7)
  n      <- 10L
  donors <- paste0("ind", seq_len(n))
  group  <- c(rep(0L, 5L), rep(1L, 5L))
  mi     <- data.frame(individual = donors, group = group, stringsAsFactors = FALSE)

  make_dl <- function(na_frac) {
    lapply(seq_len(3), function(i) {
      arr <- array(NA_real_, dim = c(n, n, 1), dimnames = list(donors, donors, "was2"))
      for (j in seq_len(n)) for (k in seq(j, n)) {
        if (runif(1) < na_frac) next
        d <- if (group[j] == group[k]) abs(rnorm(1, 0.1, 0.01)) else abs(rnorm(1, 5, 0.5))
        arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
      }
      arr
    }) |> stats::setNames(paste0("gene", seq_len(3)))
  }

  res_full <- was2code_permanova(make_dl(0),    mi, var2test = "group", n_perm = 99)
  res_na   <- was2code_permanova(make_dl(0.05), mi, var2test = "group", n_perm = 99)

  expect_true(all(res_full$pval <= 0.05, na.rm = TRUE))
  expect_true(all(res_na$pval   <= 0.05, na.rm = TRUE))
})

# ─────────────────────────────────────────────────────────────────────────────
# 3. Input validation
# ─────────────────────────────────────────────────────────────────────────────

test_that("error when meta_ind is not a data.frame", {
  d <- .make_perm_data()
  expect_error(was2code_permanova(d$dist_list, as.list(d$meta_ind), var2test = "group"))
})

test_that("error when meta_ind is missing required columns", {
  d      <- .make_perm_data()
  mi_bad <- d$meta_ind["individual"]
  expect_error(
    was2code_permanova(d$dist_list, mi_bad, var2test = "group", n_perm = 9),
    "meta_ind is missing columns"
  )
})

test_that("error when individual IDs are not unique", {
  d      <- .make_perm_data()
  mi_dup <- rbind(d$meta_ind, d$meta_ind[1, ])
  expect_error(
    was2code_permanova(d$dist_list, mi_dup, var2test = "group", n_perm = 9),
    "individual IDs in meta_ind are not unique"
  )
})

test_that("error when var2test contains NA values", {
  d <- .make_perm_data()
  d$meta_ind$group[1] <- NA
  expect_error(
    was2code_permanova(d$dist_list, d$meta_ind, var2test = "group", n_perm = 9),
    "var2test contains NA values"
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# 4. Statistical power and calibration
# ─────────────────────────────────────────────────────────────────────────────

test_that("detects strong between-group signal (small p-values)", {
  set.seed(11)
  n      <- 10L
  donors <- paste0("ind", seq_len(n))
  group  <- c(rep(0L, 5L), rep(1L, 5L))
  mi     <- data.frame(individual = donors, group = group, stringsAsFactors = FALSE)

  dist_list <- lapply(seq_len(3), function(i) {
    arr <- array(0, dim = c(n, n, 1), dimnames = list(donors, donors, "was2"))
    for (j in seq_len(n)) for (k in seq(j, n)) {
      d <- if (group[j] == group[k]) abs(rnorm(1, 0.1, 0.01)) else abs(rnorm(1, 5, 0.5))
      arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
    }
    arr
  }) |> stats::setNames(paste0("gene", seq_len(3)))

  result <- was2code_permanova(dist_list, mi, var2test = "group", n_perm = 99)
  expect_true(all(result$pval <= 0.01))
})

test_that("p-values are calibrated under the null and powerful under the alternative", {
  set.seed(99)
  n        <- 10L
  n_trials <- 100L
  donors   <- paste0("ind", seq_len(n))
  group    <- c(rep(0L, 5L), rep(1L, 5L))
  mi       <- data.frame(individual = donors, group = group, stringsAsFactors = FALSE)

  sim_pval <- function(seed, signal) {
    set.seed(seed)
    arr <- array(0, dim = c(n, n, 1), dimnames = list(donors, donors, "was2"))
    for (j in seq_len(n)) for (k in seq(j, n)) {
      d <- if (signal && group[j] != group[k]) abs(rnorm(1, 5, 0.5))
           else                                abs(rnorm(1, 2, 0.5))
      arr[j, k, "was2"] <- arr[k, j, "was2"] <- d
    }
    was2code_permanova(list(gene1 = arr), mi, var2test = "group", n_perm = 99)$pval[[1]]
  }

  seeds      <- sample.int(1000, n_trials)
  alt_pvals  <- sapply(seeds, sim_pval, signal = TRUE)
  null_pvals <- sapply(seeds, sim_pval, signal = FALSE)

  expect_gt(mean(alt_pvals  < 0.05), 0.8)  # high power under the alternative
  expect_lt(mean(null_pvals < 0.05), 0.2)  # calibrated under the null

  ks_alt  <- suppressWarnings(ks.test(alt_pvals,  "punif")$statistic)
  ks_null <- suppressWarnings(ks.test(null_pvals, "punif")$statistic)
  expect_lt(ks_null, ks_alt)               # null p-values are more uniform
})

# ─────────────────────────────────────────────────────────────────────────────
# 5. Parallelization
#
# Correctness: ncores=2 must produce identical output to ncores=1.
# Timing: ncores=2 must be meaningfully faster on a workload with enough
# permutations and donors that fork overhead is amortized.  Skipped on CI and
# Windows (where the parallel path is not taken) and when < 2 cores available.
# ─────────────────────────────────────────────────────────────────────────────

test_that("ncores=2 gives identical results to ncores=1", {
  d  <- .make_perm_data()
  r1 <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                           n_perm = 99, r_seed = 7, ncores = 1)
  r2 <- was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                           n_perm = 99, r_seed = 7, ncores = 2)
  expect_equal(r1$pval,   r2$pval)
  expect_equal(r1$F_ob,   r2$F_ob)
  expect_equal(r1$F_perm, r2$F_perm)
})

test_that("ncores=2 is faster than ncores=1 on a large permutation workload", {
  skip_on_ci()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "fewer than 2 cores available")

  # 20 genes x 30 donors: each permutation is O(genes x donors^2), giving
  # enough per-permutation work for mclapply to overcome fork overhead.
  d <- .make_perm_data(n_donors = 30L, n_genes = 20L, seed = 7)

  t_serial   <- system.time(
    was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                       n_perm = 999L, ncores = 1)
  )[["elapsed"]]
  t_parallel <- system.time(
    was2code_permanova(d$dist_list, d$meta_ind, var2test = "group",
                       n_perm = 999L, ncores = 2)
  )[["elapsed"]]

  expect_lt(t_parallel, t_serial * 0.85)
})

# ─────────────────────────────────────────────────────────────────────────────
# 7. Biology-based tests using was2code_dist canonical dataset
#
# Re-uses .make_test_data() from test_was2code_dist.R (loaded via devtools::load_all).
# The four genes are:
#   null_gene     – no signal between cases and controls
#   location_gene – mean shift (cases vs controls)
#   size_gene     – variance shift
#   shape_gene    – shape shift (bimodal vs normal, sd-matched)
# ─────────────────────────────────────────────────────────────────────────────

.make_bio_data <- function(seed = 42, n_case = 4, n_ctrl = 4, n_cpd = 100) {
  set.seed(seed)
  n_donors  <- n_case + n_ctrl
  n_cells   <- n_donors * n_cpd
  donor_ids <- c(paste0("case", seq_len(n_case)), paste0("ctrl", seq_len(n_ctrl)))
  cell_ids  <- paste0("cell", seq_len(n_cells))
  gene_names <- c("null_gene", "location_gene", "size_gene", "shape_gene")
  mat <- matrix(NA_real_, nrow = 4, ncol = n_cells,
                dimnames = list(gene_names, cell_ids))
  bimodal_sd <- sqrt(9.09)
  for (i in seq_len(n_donors)) {
    idx     <- seq((i - 1) * n_cpd + 1, i * n_cpd)
    is_case <- i <= n_case
    signs   <- sample(c(-1L, 1L), n_cpd, replace = TRUE)
    mat["null_gene",     idx] <- rnorm(n_cpd, 0, 1)
    mat["location_gene", idx] <- rnorm(n_cpd, if (is_case) 5 else 0, 1)
    mat["size_gene",     idx] <- rnorm(n_cpd, 0, if (is_case) 3 else 0.3)
    mat["shape_gene",    idx] <- if (is_case) rnorm(n_cpd, signs * 3, 0.3)
                                 else         rnorm(n_cpd, 0, bimodal_sd)
  }
  td <- list(
    mat       = mat,
    meta_cell = data.frame(cell_id    = cell_ids,
                           individual = rep(donor_ids, each = n_cpd),
                           stringsAsFactors = FALSE),
    meta_ind  = data.frame(individual = donor_ids,
                           group      = factor(c(rep("case", n_case),
                                                 rep("ctrl", n_ctrl))),
                           stringsAsFactors = FALSE)
  )
  dist_res <- was2code_dist(td$mat, td$meta_cell, td$meta_ind,
                            var2test = "group", ncores = 1)
  # was2code_permanova expects arrays with a "was2" slice; extract that slice
  dist_list <- lapply(dist_res, function(arr) arr[,, "was2", drop = FALSE])
  list(dist_list = dist_list, meta_ind = td$meta_ind)
}

bio <- .make_bio_data()

test_that("permanova detects location signal (small p-value)", {
  res <- was2code_permanova(bio$dist_list["location_gene"], bio$meta_ind,
                            var2test = "group", n_perm = 999)
  expect_lt(res$pval[["location_gene"]], 0.05)
})

test_that("permanova detects size signal (small p-value)", {
  res <- was2code_permanova(bio$dist_list["size_gene"], bio$meta_ind,
                            var2test = "group", n_perm = 999)
  expect_lt(res$pval[["size_gene"]], 0.05)
})

test_that("permanova detects shape signal (small p-value)", {
  res <- was2code_permanova(bio$dist_list["shape_gene"], bio$meta_ind,
                            var2test = "group", n_perm = 999)
  expect_lt(res$pval[["shape_gene"]], 0.05)
})

test_that("permanova does not spuriously reject null gene", {
  res <- was2code_permanova(bio$dist_list["null_gene"], bio$meta_ind,
                            var2test = "group", n_perm = 999)
  # With 8 donors and the specific seed, the null gene should not be significant.
  # Use a lenient threshold (0.5) rather than 1.0 to flag implausibly small values.
  expect_gt(res$pval[["null_gene"]], 0.05)
})

test_that("F_ob for signal genes is larger than for the null gene", {
  res <- was2code_permanova(bio$dist_list, bio$meta_ind,
                            var2test = "group", n_perm = 999)
  null_F <- res$F_ob[["null_gene"]]
  expect_gt(res$F_ob[["location_gene"]], null_F)
  expect_gt(res$F_ob[["size_gene"]],     null_F)
  expect_gt(res$F_ob[["shape_gene"]],    null_F)
})

test_that("F_perm dimensions are correct for the full gene set", {
  n_perm <- 199L
  res    <- was2code_permanova(bio$dist_list, bio$meta_ind,
                               var2test = "group", n_perm = n_perm)
  expect_equal(dim(res$F_perm), c(length(bio$dist_list), n_perm))
  expect_equal(rownames(res$F_perm), names(bio$dist_list))
  expect_equal(colnames(res$F_perm), paste0("perm", seq_len(n_perm)))
})
