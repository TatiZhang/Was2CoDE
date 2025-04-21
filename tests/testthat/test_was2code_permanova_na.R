context("Test was2code_permanova_na")
library(testthat)
library(Matrix)
library(foreach)
library(doRNG)

registerDoSEQ()

# Helper to generate test data with optional NA entries
generate_test_data_na <- function(n_samples = 20, n_genes = 5, seed = 123, na_fraction = 0) {
  set.seed(seed)
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, n_samples / 2), rep(1, n_samples / 2))
  batch <- sample(1:3, n_samples, replace = TRUE)
  meta_ind <- data.frame(individual = individuals, group = group, batch = batch)
  
  dist_list <- list()
  for (i in 1:n_genes) {
    dist_mat <- array(NA_real_, dim = c(n_samples, n_samples, 1))
    dimnames(dist_mat) <- list(individuals, individuals, "was2")
    
    for (j in 1:n_samples) {
      for (k in j:n_samples) {
        if (runif(1) < na_fraction) next
        
        base_dist <- abs(rnorm(1, mean = 2))
        group_effect <- ifelse(group[j] == group[k], 0, 1.5)
        val <- base_dist + group_effect
        dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- val
      }
    }
    dist_list[[paste0("gene", i)]] <- dist_mat
  }
  list(dist_list = dist_list, meta_ind = meta_ind)
}

test_that("Basic functionality without covariates (with NAs)", {
  data <- generate_test_data_na(na_fraction = 0.2)
  result <- was2code_permanova_na(data$dist_list, data$meta_ind, var2test = "group", n_perm = 9)
  expect_equal(length(result), length(data$dist_list))
  expect_true(all(result >= 0 & result <= 1 | is.na(result)))
  expect_named(result, names(data$dist_list))
})

test_that("Function works with covariate adjustment (with NAs)", {
  data <- generate_test_data_na(na_fraction = 0.2)
  result <- was2code_permanova_na(data$dist_list, data$meta_ind, var2test = "group", var2adjust = "batch", n_perm = 9)
  expect_equal(length(result), length(data$dist_list))
  expect_true(all(result >= 0 & result <= 1 | is.na(result)))
})

test_that("Function works with residualization (with NAs)", {
  data <- generate_test_data_na(na_fraction = 0.2)
  result <- was2code_permanova_na(data$dist_list, data$meta_ind, var2test = "group", var2adjust = "batch",
                                  residulize_x = TRUE, delta = 0.5, n_perm = 9)
  expect_equal(length(result), length(data$dist_list))
  expect_true(all(result >= 0 & result <= 1 | is.na(result)))
})

test_that("Handles character-type variable correctly (with NAs)", {
  data <- generate_test_data_na(na_fraction = 0.2)
  meta_ind_char <- data$meta_ind
  meta_ind_char$group <- ifelse(meta_ind_char$group == 1, "case", "control")
  result <- was2code_permanova_na(data$dist_list, meta_ind_char, var2test = "group", n_perm = 9)
  expect_equal(length(result), length(data$dist_list))
  expect_true(all(result >= 0 & result <= 1 | is.na(result)))
})

test_that("Errors if metadata is missing required columns", {
  data <- generate_test_data_na()
  meta_ind_missing <- data$meta_ind[, c("individual", "batch")]
  expect_error(was2code_permanova_na(data$dist_list, meta_ind_missing, var2test = "group", n_perm = 9))
})

test_that("Errors if individual IDs are not unique", {
  data <- generate_test_data_na()
  meta_ind_dup <- rbind(data$meta_ind, data$meta_ind[1, ])
  expect_error(was2code_permanova_na(data$dist_list, meta_ind_dup, var2test = "group", n_perm = 9))
})

test_that("Errors if variable to test has NA values", {
  data <- generate_test_data_na()
  meta_ind_na <- data$meta_ind
  meta_ind_na$group[1] <- NA
  expect_error(was2code_permanova_na(data$dist_list, meta_ind_na, var2test = "group", n_perm = 9))
})


test_that("was2code_permanova_na works with only one gene", {
  set.seed(321)
  n_samples <- 10
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, 5), rep(1, 5))
  meta_ind <- data.frame(individual = individuals, group = group)
  
  # One gene only
  dist_mat <- array(NA, dim = c(n_samples, n_samples, 1))
  dimnames(dist_mat) <- list(individuals, individuals, "was2")
  for (j in 1:n_samples) {
    for (k in j:n_samples) {
      if (runif(1) < 0.05) next
      dist <- if (group[j] == group[k]) 0.1 else 5
      dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- dist
    }
  }
  
  dist_list <- list(geneX = dist_mat)
  
  result <- was2code_permanova_na(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    n_perm = 99
  )
  
  expect_equal(length(result), 1)
  expect_named(result, "geneX")
  expect_true(result >= 0 && result <= 1)
})


test_that("Power is retained with small proportion of NAs", {
  set.seed(456)
  n_samples <- 10
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, 5), rep(1, 5))
  meta_ind <- data.frame(individual = individuals, group = group)
  
  dist_list <- list()
  for (i in 1:3) {
    mat <- array(NA, dim = c(n_samples, n_samples, 1))
    dimnames(mat) <- list(individuals, individuals, "was2")
    for (j in 1:n_samples) {
      for (k in j:n_samples) {
        if (runif(1) < 0.05) next
        val <- if (group[j] == group[k]) abs(rnorm(1, 0.1, 0.01)) else abs(rnorm(1, 5, 0.5))
        mat[j, k, "was2"] <- mat[k, j, "was2"] <- val
      }
    }
    dist_list[[paste0("gene", i)]] <- mat
  }
  
  result <- was2code_permanova_na(dist_list, meta_ind, var2test = "group", n_perm = 99)
  expect_equal(length(result), 3)
  expect_true(all(result <= 0.01, na.rm = TRUE))  # still strong signal despite some NAs
})

test_that("p-value distribution is uniform under null and non-uniform under alternative", {
  set.seed(789)
  n_trials <- 100
  n_samples <- 10
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, n_samples / 2), rep(1, n_samples / 2))
  
  simulate_dist_list <- function(seed, with_group_effect = TRUE) {
    set.seed(seed)
    dist_mat <- array(0, dim = c(n_samples, n_samples, 1))
    dimnames(dist_mat) <- list(individuals, individuals, "was2")
    
    for (j in 1:n_samples) {
      for (k in j:n_samples) {
        if (with_group_effect) {
          dist <- if (group[j] == group[k]) {
            abs(rnorm(1, mean = 0.1, sd = 0.01))
          } else {
            abs(rnorm(1, mean = 5, sd = 0.5))
          }
        } else {
          dist <- abs(rnorm(1, mean = 2.0, sd = 0.5))  # same distribution for all
        }
        dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- dist
      }
    }
    
    list(gene1 = dist_mat)
  }
  
  # Generate seeds to ensure reproducibility
  trial_seeds <- sample.int(1000, n_trials)
  
  get_pvals <- function(with_group_effect) {
    sapply(trial_seeds, function(seed) {
      dist_list <- simulate_dist_list(seed = seed, with_group_effect = with_group_effect)
      meta_ind <- data.frame(individual = individuals, group = group)
      pval <- was2code_permanova_na(dist_list, meta_ind, var2test = "group", n_perm = 99)[1]
      
      pval[1]
    })
  }

  alt_pvals <- get_pvals(TRUE)
  null_pvals <- get_pvals(FALSE)
  
  # Perform Kolmogorov-Smirnov test against uniform(0,1)
  ks_alt <- suppressWarnings(ks.test(alt_pvals, "punif")$statistic)
  ks_null <- suppressWarnings(ks.test(null_pvals, "punif")$statistic)
  
  expect_true(mean(alt_pvals < 0.05) > 0.8)  # should have high power
  expect_true(mean(null_pvals < 0.05) < 0.2) # should not reject often under null
  expect_true(ks_null < ks_alt)              # null p-values are more uniform
})

test_that("Identical results when there are no NAs", {
  data <- generate_test_data_na(na_fraction = 0)
  result1 <- was2code_permanova(data$dist_list, data$meta_ind, var2test = "group", n_perm = 99)
  result2 <- was2code_permanova_na(data$dist_list, data$meta_ind, var2test = "group", n_perm = 99)
  expect_equal(result1, result2)
})
