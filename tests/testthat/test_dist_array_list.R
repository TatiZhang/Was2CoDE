context("Test dist_array_list")
library(testthat)

# Create mock data for dat_res and res_ig for 3 genes each with data for 4 individuals
dat_res <- list(
  gene1 = matrix(rnorm(20), ncol = 5),
  gene2 = matrix(rnorm(20), ncol = 5),
  gene3 = matrix(rnorm(20), ncol = 5)
)
res_ig <- list(
  gene1 = list(rnorm(5), rnorm(5), rnorm(5), rnorm(5)),
  gene2 = list(rnorm(5), rnorm(5), rnorm(5), rnorm(5)),
  gene3 = list(rnorm(5), rnorm(5), rnorm(5), rnorm(5))
)

test_that("dist_array_list outputs correct structure and size", {
  results <- dist_array_list(dat_res, res_ig)
  
  expect_true(is.list(results))
  expect_length(results, length(dat_res))
  
  for (result in results) {
    expect_true(is.array(result))
    expect_equal(dim(result), c(4, 4, 6))  # Assuming each gene has data for 4 individuals
  }
})

test_that("dist_array_list has zeros on diagonals", {
  results <- dist_array_list(dat_res, res_ig)
  
  for (result in results) {
    for (i in 1:dim(result)[1]) {
      expect_equal(result[i, i, ], rep(0, 6))
    }
  }
})

test_that("dist_array_list produces symmetric matrices", {
  results <- dist_array_list(dat_res, res_ig)
  
  for (result in results) {
    for (i in 1:dim(result)[1]) {
      for (j in 1:dim(result)[2]) {
        expect_equal(result[i, j, ], result[j, i, ])
      }
    }
  }
})