context("Test was2code_permanova")
library(testthat)
library(Matrix)
library(foreach)
library(doRNG)

# Register a sequential backend to avoid warnings about missing parallel backend
registerDoSEQ()

test_that("was2code_permanova function works correctly", {
  # Generate sample data
  set.seed(123)
  n_samples <- 20
  n_genes <- 5
  
  # Create metadata
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, n_samples/2), rep(1, n_samples/2))
  batch <- sample(1:3, n_samples, replace = TRUE)
  meta_ind <- data.frame(
    individual = individuals,
    group = group,
    batch = batch
  )
  
  # Create distance matrices
  dist_list <- list()
  for(i in 1:n_genes) {
    # Create a symmetric distance matrix with was2 dimension
    dist_mat <- array(0, dim = c(n_samples, n_samples, 1))
    dimnames(dist_mat) <- list(individuals, individuals, "was2")
    
    # Fill with random distances
    for(j in 1:n_samples) {
      for(k in j:n_samples) {
        # Generate a distance that's somewhat related to the group
        base_dist <- abs(rnorm(1, mean = 2))
        group_effect <- ifelse(group[j] == group[k], 0, 1.5)
        dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- base_dist + group_effect
      }
    }
    
    gene_name <- paste0("gene", i)
    dist_list[[gene_name]] <- dist_mat
  }
  names(dist_list) <- paste0("gene", 1:n_genes)
  
  # Test case 1: Basic functionality without covariates
  result1 <- was2code_permanova(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    n_perm = 99
  )
  
  # Verify result structure
  expect_equal(length(result1), n_genes)
  expect_true(all(result1 >= 0 & result1 <= 1))
  expect_named(result1, paste0("gene", 1:n_genes))
  
  # Test case 2: With covariate adjustment
  result2 <- was2code_permanova(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    var2adjust = "batch",
    n_perm = 99
  )
  
  expect_equal(length(result2), n_genes)
  expect_true(all(result2 >= 0 & result2 <= 1))
  
  # Test case 3: With residualization
  result3 <- was2code_permanova(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    var2adjust = "batch",
    residulize_x = TRUE,
    n_perm = 99,
    delta = 0.5
  )
  
  expect_equal(length(result3), n_genes)
  expect_true(all(result3 >= 0 & result3 <= 1 | is.na(result3)))
  
  # Test case 4: Input validation - character variable
  meta_ind_char <- meta_ind
  meta_ind_char$group <- ifelse(meta_ind_char$group == 1, "case", "control")
  
  result4 <- was2code_permanova(
    dist_list = dist_list,
    meta_ind = meta_ind_char,
    var2test = "group",
    n_perm = 99
  )
  
  expect_equal(length(result4), n_genes)
  expect_true(all(result4 >= 0 & result4 <= 1))
  
  # # Test case 5: Handle NAs in distance matrices
  # dist_list_with_na <- dist_list
  # for(gene in names(dist_list_with_na)) {
  #   dist_list_with_na[[gene]][1, 2, "was2"] <- dist_list_with_na[[gene]][2, 1, "was2"] <- NA
  # }
  # 
  # result5 <- was2code_permanova(
  #   dist_list = dist_list_with_na,
  #   meta_ind = meta_ind,
  #   var2test = "group",
  #   n_perm = 99
  # )
  # 
  # expect_equal(length(result5), n_genes)
  # expect_true(all(is.na(result5)))
  
  # Test case 6: Error handling - missing columns in metadata
  meta_ind_missing <- meta_ind[, c("individual", "batch")]
  
  expect_error(
    was2code_permanova(
      dist_list = dist_list,
      meta_ind = meta_ind_missing,
      var2test = "group",
      n_perm = 99
    ),
    "names of meta_ind should conttain: individual, group"
  )
  
  # Test case 7: Error handling - non-unique individuals
  meta_ind_dup <- rbind(meta_ind, meta_ind[1, ])
  
  expect_error(
    was2code_permanova(
      dist_list = dist_list,
      meta_ind = meta_ind_dup,
      var2test = "group",
      n_perm = 99
    ),
    "the individual ids in meta_ind are not unique"
  )
  
  # Test case 8: Error handling - NA in var2test
  meta_ind_na <- meta_ind
  meta_ind_na$group[1] <- NA
  
  expect_error(
    was2code_permanova(
      dist_list = dist_list,
      meta_ind = meta_ind_na,
      var2test = "group",
      n_perm = 99
    ),
    "variable to test has NA values"
  )
})

test_that("was2code_permanova detects group differences when signal is strong", {
  set.seed(456)
  n_samples <- 10
  n_genes <- 3
  
  # Create clear group separation
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, n_samples / 2), rep(1, n_samples / 2))
  meta_ind <- data.frame(
    individual = individuals,
    group = group
  )
  
  dist_list <- list()
  for (i in 1:n_genes) {
    dist_mat <- array(0, dim = c(n_samples, n_samples, 1))
    dimnames(dist_mat) <- list(individuals, individuals, "was2")
    
    for (j in 1:n_samples) {
      for (k in j:n_samples) {
        # Within group = low distance; between group = high distance
        if (group[j] == group[k]) {
          dist <- abs(rnorm(1, mean = 0.1, sd = 0.01))  # tight within-group distances
        } else {
          dist <- abs(rnorm(1, mean = 5, sd = 0.5))      # much larger between-group distances
        }
        dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- dist
      }
    }
    
    gene_name <- paste0("gene", i)
    dist_list[[gene_name]] <- dist_mat
  }
  names(dist_list) <- paste0("gene", 1:n_genes)
  
  # Run PERMANOVA
  result_power <- was2code_permanova(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    n_perm = 99
  )
  
  expect_equal(length(result_power), n_genes)
  expect_true(all(result_power >= 0 & result_power <= 1))
  expect_true(all(result_power <= 0.01))  # should have strong signal
})

test_that("p-value distribution is uniform under null and non-uniform under alternative", {
  set.seed(789)
  n_trials <- 100
  n_samples <- 10
  individuals <- paste0("ind", 1:n_samples)
  group <- c(rep(0, n_samples / 2), rep(1, n_samples / 2))
  
  simulate_dist_list <- function(with_group_effect = TRUE) {
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
    
    list(gene1 = dist_mat,
         gene2 = dist_mat)
  }
  
  get_pvals <- function(with_group_effect) {
    replicate(n_trials, {
      dist_list <- simulate_dist_list(with_group_effect)
      meta_ind <- data.frame(individual = individuals, group = group)
      
      pval <- was2code_permanova(
        dist_list = dist_list,
        meta_ind = meta_ind,
        var2test = "group",
        n_perm = 99
      )
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


