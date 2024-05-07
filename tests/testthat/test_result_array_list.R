context("Test result_array_list")
library(testthat)
test_that("result_array_list outputs correctly", {
  # Create mock data for dist_array_list
  n_gene <- 3
  n_ind <- 2
  gene_ids <- c("Gene1", "Gene2", "Gene3")
  individuals <- c("Ind1", "Ind2")
  mock_dist_array_list <- list()
  metrics <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  # Fill each entry
  for (i in 1:n_gene) {
    array_data <- array(NA, dim = c(n_ind, n_ind, 6),
                        dimnames = list(individuals, individuals, metrics))
    
    for (j in 1:n_ind) {
      for (k in 1:n_ind) {
        for (m in 1:6) {
          array_data[j, k, m] <- runif(1)  # Random float between 0 and 1
        }
      }
    }
    mock_dist_array_list[[gene_ids[i]]] <- array_data}
  meta_ind <- data.frame(individual = individuals)
  rownames(meta_ind) <- individuals
  results <- result_array_list(mock_dist_array_list, meta_ind)
  expect_true(is.list(results))
  expect_length(results, 6)
  expected_names <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  expect_equal(names(results), expected_names)
  # Check each array's dimension and names
  for (result in results) {
    expect_equal(dim(result), c(n_gene, n_ind, n_ind))
    expect_equal(dimnames(result)[[1]], gene_ids)
    expect_equal(dimnames(result)[[2]], individuals)
    expect_equal(dimnames(result)[[3]], individuals)
  }
})
