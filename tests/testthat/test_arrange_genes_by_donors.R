context("Test arrange_genes_by_donors")
library(testthat)
library(transport)
library(data.table)
library(doRNG)
# load("tests/assets/test_data1.RData")
load("../assets/test_data1.RData")
count_matrix <- round(count_matrix)
count_matrix[count_matrix < 0] <- 0
n_gene = nrow(count_matrix)
test_that("arrange_genes_by_donors outputs correctly", {
  dat_res <- arrange_genes_by_donors(count_matrix, meta_ind, meta_cell)
  

  expect_true(is.list(dat_res),"dat_res should be a list")
  expect_true(length(dat_res) == n_gene, "dat_res length should match number of genes")
  expect_equal(names(dat_res), as.character(rownames(count_matrix)), "Names of list entries should match gene ids")
  
  # Check that each entry in the list is a list itself and contains data for each donor
  for (i in seq_along(dat_res)) {
    expect_true(is.list(dat_res[[i]]), "Each gene's data should be a list")
    expect_true(length(dat_res[[i]]) == nrow(meta_ind), "res_ig length should match number of donors")
    expect_equal(names(dat_res[[i]]), as.character(meta_ind$individual), "Names of list entries should match individuals")
    
    for (j in seq_along(dat_res[[i]])) {
      dat_ind <- dat_res[[i]][[j]]
      donor_id = names(dat_res[[i]])[[j]]
      w2use = which(meta_cell$individual == donor_id)  # grab all indexes
      dat_j = count_matrix[i, w2use] #the count matrix rows corresponding to this gene and columns for the donor
      test_result <- all.equal(dat_ind, dat_j, tolerance = 1e-8)
      expect_true(isTRUE(test_result), "Data for individuals should match the expected count matrix slice")
    }
  }
})