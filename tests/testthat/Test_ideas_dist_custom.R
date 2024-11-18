context("Test ideas_dist_custom")
library(testthat)
library(transport)
library(data.table)
library(doRNG)
# load("tests/assets/test_data1.RData")
load("../assets/test_data1.RData")

# Unit test for ideas_dist_custom
test_that("ideas_dist_custom outputs correctly", {
  # Mock data and parameters for testing
  output <- ideas_dist_custom(
    count_input = count_matrix,
    meta_cell = meta_cell,
    meta_ind = meta_ind,
    var_per_cell = var_per_cell,
    var2test = "Study_DesignationCtrl",
    var2test_type = var2test_type
  )
  
 # 1. Test if output type is a list of six 3-dimensional arrays with correct names
  expect_is(output, "list", info = "Output should be a list.")
  expect_equal(length(output), 6, info = "output list has 6 elements")
  expected_names <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  expect_equal(names(output), expected_names, info = "Output names should match expected names.")
  
  for (i in 1:6) {
    gene_ids <- rownames(count_matrix)
    expect_true(is.array(output[[i]]))
    expect_true(is.list(dimnames(output[[i]])))
    expect_true(length(dimnames(output[[i]]))==3)
    expect_true(length(dimnames(output[[i]])[[1]])==length(gene_ids))
    expected_dimnames <- list(
      gene_ids ,     
      meta_ind$individual ,
      meta_ind$individual 
    )
    for (kk in 1:3){
      expect_true(all(dimnames(output)[[i]] == expected_dimnames[[kk]]))
  }
    }

 # 2. Test if the each array was computed correctly
  #   ## Get a sample of known distances/locations/... between the first two individuals is as expected to validate
  #   ## Extract data for 3 randomly selected genes for the first two individuals
  # set.seed(0)
  # selected_genes <- sample(nrow(count_matrix), 3)
  # # Extract function results for the selected genes and the first two individuals
  # distances_from_function <- output[[1]][selected_genes, 1 , 2]
  # location_from_function <- output[[2]][selected_genes, 1, 2]
  # location_sign_from_function <- output[[3]][selected_genes, 1 , 2]
  # size_from_function <- output[[4]][selected_genes, 1 , 2]
  # size_sign_from_function <- output[[5]][selected_genes, 1 , 2]
  # shape_from_function <- output[[6]][selected_genes, 1 , 2]
  # # Manually compute each component
  # manual_distances <- numeric(length = 3)
  # names(manual_distances) <- gene_ids[selected_genes]
  # manual_location <- numeric(length = 3)
  # names(manual_location) <- gene_ids[selected_genes]
  # manual_location_sign <- numeric(length = 3)
  # names(manual_location_sign) <- gene_ids[selected_genes]
  # manual_size <- numeric(length = 3)
  # names(manual_size) <- gene_ids[selected_genes]
  # manual_size_sign <- numeric(length = 3)
  # names(manual_size_sign) <- gene_ids[selected_genes]
  # manual_shape <- numeric(length = 3)
  # names(manual_shape) <- gene_ids[selected_genes]
  # ## Loop through each selected gene
  # for (i in 1:3) {
  #   gene_index <- selected_genes[i]
  #   cols_ind1 <- which(meta_cell$individual == meta_ind$individual[1])
  #   cols_ind2 <- which(meta_cell$individual == meta_ind$individual[2])
  #   data_ind1 <- sort(count_matrix[gene_index, cols_ind1, drop = FALSE])
  #   data_ind2 <- sort(count_matrix[gene_index, cols_ind2, drop = FALSE]) # Ensure the data are sorted (necessary for the 1D Wasserstein calculation)
  #   # Manual calculations
  #   manual_distances[i] <- transport::wasserstein1d(data_ind1, data_ind2, p = 2) # Compute Wasserstein-2 distance
  #   manual_location[i] <- (mean(data_ind1) - mean(data_ind2))^2
  #   manual_location_sign[i] <- mean(data_ind1) - mean(data_ind2)
  #   manual_size[i] <- (sd(data_ind1) - sd(data_ind2))^2
  #   manual_size_sign[i] <- sd(data_ind1) - sd(data_ind2)
  #   quantiles <- seq(0, 1, length.out = 100)
  #   quantiles_1 <- quantile(data_ind1, probs = quantiles)
  #   quantiles_2 <- quantile(data_ind2, probs = quantiles)
  #   quantile_cor_ind1_ind2 <- cor(quantiles_1, quantiles_2)
  #   manual_shape[i] <- abs(2 * sd(data_ind1) * sd(data_ind2) * (1 - quantile_cor_ind1_ind2))
  # }
  # 
  # for (i in seq_along(selected_genes)) {
  #   expect_equal(distances_from_function[i], manual_distances[i], tolerance = 0.01,
  #                info = sprintf("Distances for gene %d do not match.", selected_genes[i]))
  #   expect_equal(location_from_function[i], manual_location[i], tolerance = 0.01,
  #                info = sprintf("Locations for gene %d do not match.", selected_genes[i]))
  #   expect_equal(location_sign_from_function[i], manual_location_sign[i], tolerance = 0.01,
  #                info = sprintf("Location signs for gene %d do not match.", selected_genes[i]))
  #   expect_equal(size_from_function[i], manual_size[i], tolerance = 0.01,
  #                info = sprintf("Sizes for gene %d do not match.", selected_genes[i]))
  #   expect_equal(size_sign_from_function[i], manual_size_sign[i], tolerance = 0.01,
  #                info = sprintf("Size signs for gene %d do not match.", selected_genes[i]))
  #   expect_equal(shape_from_function[i], manual_shape[i], tolerance = 0.01,
  #                info = sprintf("Shapes for gene %d do not match.", selected_genes[i]))
  # }
  # Validate calculations for selected genes and donors
  set.seed(0)
  selected_genes <- sample(nrow(count_matrix), 3)
  for (gene in selected_genes) {
    # Manual calculations
    cols_ind1 <- which(meta_cell$individual == meta_ind$individual[1])
    cols_ind2 <- which(meta_cell$individual == meta_ind$individual[2])
    data_ind1 <- sort(count_matrix[gene, cols_ind1])
    data_ind2 <- sort(count_matrix[gene, cols_ind2])
    
    manual_distance <- transport::wasserstein1d(data_ind1, data_ind2, p = 2)
    expect_equal(output[[1]][gene, 1, 2], manual_distance, tolerance = 0.01)
  }
  })

