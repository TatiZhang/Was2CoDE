context("Test was2code_permanova_na")
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
    dist_mat <- array(NA, dim = c(n_samples, n_samples, 1))
    dimnames(dist_mat) <- list(individuals, individuals, "was2")
    
    # Fill with random distances
    for(j in 1:n_samples) {
      for(k in j:n_samples) {
        bool_idx <- sample(c(TRUE, FALSE), size = 1)
        
        if(bool_idx){
          # Generate a distance that's somewhat related to the group
          base_dist <- abs(rnorm(1, mean = 2))
          group_effect <- ifelse(group[j] == group[k], 0, 1.5)
          dist_mat[j, k, "was2"] <- dist_mat[k, j, "was2"] <- base_dist + group_effect
        }
      }
    }
    
    gene_name <- paste0("gene", i)
    dist_list[[gene_name]] <- dist_mat
  }
  names(dist_list) <- paste0("gene", 1:n_genes)
  
  # Test case 1: Basic functionality without covariates
  result1 <- was2code_permanova_na(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    n_perm = 9
  )
  
  # Verify result structure
  expect_equal(length(result1), n_genes)
  expect_true(all(result1 >= 0 & result1 <= 1))
  expect_named(result1, paste0("gene", 1:n_genes))
  
  # Test case 2: With covariate adjustment
  result2 <- was2code_permanova_na(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    var2adjust = "batch",
    n_perm = 9
  )
  
  expect_equal(length(result2), n_genes)
  expect_true(all(result2 >= 0 & result2 <= 1))
  
  # Test case 3: With residualization
  result3 <- was2code_permanova_na(
    dist_list = dist_list,
    meta_ind = meta_ind,
    var2test = "group",
    var2adjust = "batch",
    residulize_x = TRUE,
    n_perm = 9,
    delta = 0.5
  )
  
  expect_equal(length(result3), n_genes)
  expect_true(all(result3 >= 0 & result3 <= 1 | is.na(result3)))
  
  # Test case 4: Input validation - character variable
  meta_ind_char <- meta_ind
  meta_ind_char$group <- ifelse(meta_ind_char$group == 1, "case", "control")
  
  result4 <- was2code_permanova_na(
    dist_list = dist_list,
    meta_ind = meta_ind_char,
    var2test = "group",
    n_perm = 9
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
    was2code_permanova_na(
      dist_list = dist_list,
      meta_ind = meta_ind_missing,
      var2test = "group",
      n_perm = 99
    )
  )
  
  # Test case 7: Error handling - non-unique individuals
  meta_ind_dup <- rbind(meta_ind, meta_ind[1, ])
  
  expect_error(
    was2code_permanova_na(
      dist_list = dist_list,
      meta_ind = meta_ind_dup,
      var2test = "group",
      n_perm = 99
    )
  )
  
  # Test case 8: Error handling - NA in var2test
  meta_ind_na <- meta_ind
  meta_ind_na$group[1] <- NA
  
  expect_error(
    was2code_permanova_na(
      dist_list = dist_list,
      meta_ind = meta_ind_na,
      var2test = "group",
      n_perm = 99
    )
  )
})
