context("Test divergence")
library(testthat)
library(transport) 
library(Rcpp) 

test_that("divergence function calculates correct values", {
  set.seed(0)
  a <- rnorm(1000,mean=0,sd=1)
  b <- rnorm(1000,mean=2,sd=1)
  
  result <- divergence(a, b)
  
  # Check the length and names of the result
  expect_equal(length(result), 6)
  
  quantiles <- seq(0, 1, length.out = 100)
  quantiles_a <- quantile(a, probs = quantiles)
  quantiles_b <- quantile(b, probs = quantiles)
  quantile_cor_ab <- cor(quantiles_a, quantiles_b)

  # Check specific calculations
  expect_equal(result[1], transport::wasserstein1d(a, b,p=2,wa = NULL, wb = NULL))
  expect_equal(result[2], (mean(a) - mean(b))^2)
  expect_equal(result[3], mean(a) - mean(b))
  expect_equal(result[4], (sd(a) - sd(b))^2)
  expect_equal(result[5], sd(a) - sd(b))
  expect_equal(result[6], abs(2 * sd(a) * sd(b) * (1 - quantile_cor_ab)))
})

test_that("divergence function handles equal inputs", {
  set.seed(0)
  a <- rnorm(1000,mean=0,sd=1)
  b <- a  # Identical to a

  result <- divergence(a, b)
  quantiles <- seq(0, 1, length.out = 100)
  quantiles_a <- quantile(a, probs = quantiles)
  quantiles_b <- quantile(b, probs = quantiles)
  quantile_cor_ab <- cor(quantiles_a, quantiles_b)

  # When inputs are equal, distance and other measures should be zero or neutral
  expect_equal(result[1], 0)
  expect_equal(result[2], 0)
  expect_equal(result[3], 0)
  expect_equal(result[4], 0)
  expect_equal(result[5], 0)
  expect_equal(result[6], 0)
  })
 
test_that("divergence function handles incorrect input types", {
  a <- "not numeric"
  b <- "not numeric"
  expect_error(divergence(a, b))
})


test_that("manual shape calculation is correct", {
  set.seed(0)
  a <- rnorm(1000,mean=0,sd=1)
  b <- rnorm(1000,mean=2,sd=1)
  result <- divergence(a, b)
  
  quantiles <- seq(0, 1, length.out = 1e5)
  quantiles_a <- quantile(a, probs = quantiles)
  quantiles_b <- quantile(b, probs = quantiles)
  quantile_cor_ab <- cor(quantiles_a, quantiles_b)

  shape_manual 		= abs(2 * sd(a) * sd(b) * (1 - quantile_cor_ab))

  result <- divergence(a, b)
  expect_true(abs(result[1]^2 - result[2] - result[4] - shape_manual)<=1e-3)
})
