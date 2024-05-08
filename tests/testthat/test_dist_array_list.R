context("Test dist_array_list")
library(testthat)

test_that("dist_array_list initializes correctly", {
  # Mock dat_res
  p_genes <- 5
  n_donors <- 3
  gene_ids <- paste("Gene", 1:p_genes, sep="")
  donor_ids <- paste("Donor", 1:n_donors, sep="")
  dat_res <- vector("list", length = p_genes)
  names(dat_res) <- gene_ids
  for (i in 1:p_genes) {
    dat_res[[i]] <- vector("list", length = n_donors)
    names(dat_res[[i]]) <- donor_ids
  }
  array_list <- dist_array_list(dat_res)
  
  expect_true(is.list(array_list))
  expect_equal(length(array_list), p_genes)
  expect_equal(names(array_list), gene_ids)
  
  for (i in seq_along(array_list)) {
    array <- array_list[[i]]
    expect_true(is.array(array))
    expect_equal(dim(array), c(n_donors, n_donors, 6))
    expect_equal(dimnames(array)[[1]], donor_ids)
    expect_equal(dimnames(array)[[2]], donor_ids)
    expect_equal(dimnames(array)[[3]], c("distance",
                                         "location", 
                                         "location_sign", 
                                         "size", 
                                         "size_sign",
                                         "shape"))
    expect_true(all(is.na(array)))
  }
})
