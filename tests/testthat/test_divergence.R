context("Test divergence")
library(testthat)
library(transport) 
library(Rcpp) 

test_that("divergence function calculates correct values", {
  a <- c(1, 2, 3, 4, 5)
  b <- c(1, 4, 3, 2, 1)
  
  result <- divergence(a, b)
  
  # Check the length and names of the result
  expect_equal(length(result), 6)
  # Check specific calculations
  expect_equal(result[1], transport::wasserstein1d(a, b,p=2,wa = NULL, wb = NULL))
  expect_equal(result[2], (mean(a) - mean(b))^2)
  expect_equal(result[3], mean(a) - mean(b))
  expect_equal(result[4], (sd(a) - sd(b))^2)
  expect_equal(result[5], sd(a) - sd(b))
  expect_equal(result[6], ((transport::wasserstein1d(a, b,p=2,wa = NULL, wb = NULL)^2 - (mean(a) - mean(b))^2 - (sd(a) - sd(b))^2) / (2 * sd(a) * sd(b))))
})

test_that("divergence function handles equal inputs", {
  a <- c(2, 2, 2, 2, 2)
  b <- a  # Identical to a

  result <- divergence(a, b)

  # When inputs are equal, distance and other measures should be zero or neutral
  expect_equal(result[1], 0)
  expect_equal(result[2], 0)
  expect_equal(result[3], 0)
  expect_equal(result[4], 0)
  expect_equal(result[5], 0)
  expect_equal(result[6], ((transport::wasserstein1d(a, b,p=2,wa = NULL, wb = NULL)^2 - (mean(a) - mean(b))^2 - (sd(a) - sd(b))^2) / (2 * sd(a) * sd(b))))
})

 
test_that("divergence function handles incorrect input types", {
  a <- "not numeric"
  b <- "not numeric"
  expect_error(divergence(a, b))
})


test_that("manual shape calculation is correct", {
  set.seed(24)
  a <- rnorm(100,mean=0,sd=1)
  b <- rnorm(100,mean=2,sd=1)
  
  quantiles <- seq(0.1, 0.9, by = 0.1)
  quantiles_a <- quantile(a, probs = quantiles)
  quantiles_b <- quantile(b, probs = quantiles)
  quantile_cor_ab <- cor(quantiles_a, quantiles_b)
  
  mean_a <- mean(a)
  mean_b <- mean(b)
  sd_a <- sd(a)
  sd_b <- sd(b)
  
  cat("Intermediate values:\n")
  cat("mean_a:", mean_a, "\n")
  cat("mean_b:", mean_b, "\n")
  cat("sd_a:", sd_a, "\n")
  cat("sd_b:", sd_b, "\n")
  cat("quantile_cor_ab:", quantile_cor_ab, "\n")
  
  shape_manual 		= abs(2 * sd_a * sd_b * (1 - quantile_cor_ab))

  
  result <- divergence(a, b)
  
  expect_equal(result[6], shape_manual)
})
