context("Test was2code_lfc")
library(testthat)
library(Matrix)
library(foreach)
library(doRNG)
library(data.table) 

# Load test data
load("../assets/test_data1.RData")
# test-was2code_lfc.R

test_that("was2code_lfc computes standardized logFC-style values correctly", {
  # Simulate a small distance list with one gene
  set.seed(123)
  gene_name <- "GeneA"
  n_samples <- 4
  dist_components <- c("was2", "location", "size", "shape")
  
  # Create a 4 x 4 x 4 array of random values for each component
  dist_array <- array(runif(n_samples^2 * length(dist_components)), 
                      dim = c(n_samples, n_samples, length(dist_components)),
                      dimnames = list(NULL, NULL, dist_components))
  
  # Assign to a named list
  dist_list <- list()
  dist_list[[gene_name]] <- dist_array
  
  # Define case and control indices
  case_idx <- c(1, 2)
  control_idx <- c(3, 4)
  
  # Run the function
  result <- was2code_lfc(dist_list, case_idx, control_idx)
  
  # Check the output structure
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(dist_components))
  expect_equal(rownames(result), gene_name)
  expect_equal(colnames(result), dist_components)
  
  # Check values are finite or NA (i.e., no unexpected values)
  expect_true(all(is.finite(result) | is.na(result)))
})
