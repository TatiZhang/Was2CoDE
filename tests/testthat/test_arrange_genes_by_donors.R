context("Test arrange_genes_by_donors")
library(testthat)
library(transport)
library(data.table)
library(doRNG)

# Load test data
load("../assets/test_data1.RData")

# Ensure count matrix contains only non-negative integers
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0
n_gene <- nrow(count_matrix)

test_that("arrange_genes_by_donors outputs correctly", {
  # Run the function
  dat_res <- arrange_genes_by_donors(count_matrix, meta_ind, meta_cell)
  
  # Basic structure tests
  expect_true(is.list(dat_res), info = "dat_res should be a list")
  expect_equal(length(dat_res), n_gene, tolerance = 1e-8, info = "dat_res length should match number of genes")
  expect_equal(names(dat_res), as.character(rownames(count_matrix)), info = "Names of list entries should match gene ids")
  
  # Test a subset of genes to save time while still validating functionality
  # Use up to 5 genes or all genes if fewer than 5
  test_genes <- sample(1:n_gene, min(5, n_gene))
  
  # Check that each entry in the list is a list itself and contains data for each donor
  for (i in test_genes) {
    gene_data <- dat_res[[i]]
    
    expect_true(is.list(gene_data), 
                info = sprintf("Gene %d data should be a list", i))
    expect_equal(length(gene_data), nrow(meta_ind), tolerance = 1e-8,
                 info = sprintf("For gene %d, list length should match number of donors", i))
    expect_equal(names(gene_data), as.character(meta_ind$individual), 
                 info = sprintf("For gene %d, names should match individuals", i))
    
    # Check data for each donor
    for (j in seq_along(meta_ind$individual)) {
      donor_id <- meta_ind$individual[j]
      donor_cells <- which(meta_cell$individual == donor_id)
      
      expected_data <- count_matrix[i, donor_cells]
      actual_data <- gene_data[[j]]
      
      # Handle potential name differences in expected vs actual data
      if (!is.null(names(expected_data)) && !is.null(names(actual_data))) {
        # Sort both vectors by name to ensure they match
        expected_data <- expected_data[order(names(expected_data))]
        actual_data <- actual_data[order(names(actual_data))]
      }
      
      # Compare lengths
      expect_equal(length(actual_data), length(expected_data), tolerance = 1e-8,
                   info = sprintf("For gene %d and donor %s, data length mismatch", i, donor_id))
      
      # Compare values
      if (length(expected_data) > 0) {
        # Convert to numeric if needed to ensure proper comparison
        expect_equal(as.numeric(actual_data), as.numeric(expected_data), tolerance = 1e-8,
                     info = sprintf("For gene %d and donor %s, data values mismatch", i, donor_id))
      }
    }
  }
})