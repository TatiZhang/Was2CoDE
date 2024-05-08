context("Test ideas_dist_custom")
library(testthat)
library(transport)
library(data.table)
library(doRNG)
# load("tests/assets/test_data1.RData")
load("../assets/test_data1.RData")
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0
# Unit test for ideas_dist_custom
test_that("ideas_dist_custom outputs correctly", {
  # Mock data and parameters for testing
  output <- ideas_dist_custom(
    count_input = count_matrix,
    meta_cell = meta_cell,
    meta_ind = meta_ind,
    var_per_cell = var_per_cell,
    var2test = var2test,
    var2test_type = var2test_type
  )
  
 # 1. Test if output type is a list of six 3-dimensional arrays with correct names
  expect_is(output, "list", info = "Output should be a list.")
  expect_equal(length(output), 6, info = "output list has 6 elements")
  expected_names <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  expect_equal(names(output), expected_names, info = "Output names should match expected names.")
  
  for (i in 1:6) {
    expect_true(is.array(output[[i]]))
    expect_true(is.list(dimnames(output[[i]])))
    expect_true(length(dimnames(output[[i]]))==3)
    # expect_true(length(dimnames(output[[i]][[1]]))== length(rownames(count_matrix)))
    expected_dimnames <- list(
      gene_ids <- rownames(count_matrix),     
      ind_ids <- as.character(meta_ind$individual),
      ind_ids <- as.character(meta_ind$individual)
    )
      expect_true(all(dimnames(output[[i]][[1]])==expected_dimnames[[1]]))
      expect_true(dimnames(output)[[i]][[2]] == expected_dimnames[[2]])
      expect_true(dimnames(output)[[i]][[3]] == expected_dimnames[[3]])
      
    }

 # 2. Test if the each array was computed correctly
    ## Get a sample of known distances/locations/... between the first two individuals is as expected to validate
    ## Extract data for 3 randomly selected genes for the first two individuals
  set.seed(10)
  selected_genes <- sample(nrow(count_matrix), 3)
  # Extract function results for the selected genes
  distances_from_function <- output[[1]][selected_genes, , ]
  location_from_function <- output[[2]][selected_genes, , ]
  location_sign_from_function <- output[[3]][selected_genes, , ]
  size_from_function <- output[[4]][selected_genes, , ]
  size_sign_from_function <- output[[5]][selected_genes, , ]
  shape_from_function <- output[[6]][selected_genes, , ]
  # Manually compute each component
  manual_distances <- numeric(length = 3)
  manual_location <- numeric(length = 3)
  manual_location_sign <- numeric(length = 3)
  manual_size <- numeric(length = 3)
  manual_size_sign <- numeric(length = 3)
  manual_shape <- numeric(length = 3)
  ## Loop through each selected gene
  for (i in seq_along(selected_genes)) {
    gene_index <- selected_genes[i]
    data_ind1 <- sort(count_matrix[gene_index, which(colnames(count_matrix) == meta_ind$individual[1]), drop = FALSE])
    data_ind2 <- sort(count_matrix[gene_index, which(colnames(count_matrix) == meta_ind$individual[2]), drop = FALSE]) # Ensure the data are sorted (necessary for the 1D Wasserstein calculation)
    # Manual calculations
    manual_distances[i] <- transport::wasserstein1d(data_ind1, data_ind2, p = 2) # Compute Wasserstein-2 distance
    manual_location[i] <- (mean(data_ind1) - mean(data_ind2))^2
    manual_location_sign[i] <- mean(data_ind1) - mean(data_ind2)
    manual_size[i] <- (sd(data_ind1) - sd(data_ind2))^2
    manual_size_sign[i] <- sd(data_ind1) - sd(data_ind2)
    manual_shape[i] <- (manual_distance^2 - manual_location - manual_size) / (2 * sd(data_ind1) * sd(data_ind2))
  }
  for (i in seq_along(selected_genes)) {
    expect_equal(distances_from_function[i], manual_distances[i], tolerance = 0.01,
                 info = sprintf("Distances for gene %d do not match.", selected_genes[i]))
    expect_equal(location_from_function[i], manual_location[i], tolerance = 0.01,
                 info = sprintf("Locations for gene %d do not match.", selected_genes[i]))
    expect_equal(location_sign_from_function[i], manual_location_sign[i], tolerance = 0.01,
                 info = sprintf("Location signs for gene %d do not match.", selected_genes[i]))
    expect_equal(size_from_function[i], manual_size[i], tolerance = 0.01,
                 info = sprintf("Sizes for gene %d do not match.", selected_genes[i]))
    expect_equal(size_sign_from_function[i], manual_size_sign[i], tolerance = 0.01,
                 info = sprintf("Size signs for gene %d do not match.", selected_genes[i]))
    expect_equal(shape_from_function[i], manual_shape[i], tolerance = 0.01,
                 info = sprintf("Shapes for gene %d do not match.", selected_genes[i]))
  }
  # Example test to check if diagonal elements of each matrix are zero
  expect_true(all(sapply(output, function(x) all(diag(x[,,1]) == 0))), "Diagonal elements should be zero.")
  
  # Testing error handling
  expect_error(ideas_dist_custom(count_input = "not a matrix", meta_cell, meta_ind, "individual", "var2test", "binary", "NB"),
               "count_matrix is not a matrix")
  })

