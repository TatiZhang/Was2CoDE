context("Test divergence")
library(testthat)
library(transport) 


test_that("divergence function calculates correct values", {
  a <- c(1, 2, 3, 4, 5)
  b <- c(5, 4, 3, 2, 1)
  
  result <- divergence(a, b)
  
  # Check the length and names of the result
  expect_equal(length(result), 6)
  expect_named(result, c("distance", "location", "location_sign", "size", "size_sign", "shape"))
  
  # Check specific calculations
  expect_equal(result["distance"], transport::wasserstein1d(a, b,p=2))
  expect_equal(result["location"], (mean(a) - mean(b))^2)
  expect_equal(result["location_sign"], mean(a) - mean(b))
  expect_equal(result["size"], (sd(a) - sd(b))^2)
  expect_equal(result["size_sign"], sd(a) - sd(b))
  expect_equal(result["shape"], ((result["distance"]^2 - result["location"] - result["size"]) / (2 * sd(a) * sd(b))))
})

test_that("divergence function handles equal inputs", {
  a <- c(2, 2, 2, 2, 2)
  b <- a  # Identical to a
  
  result <- divergence(a, b)
  
  # When inputs are equal, distance and other measures should be zero or neutral
  expect_equal(result["distance"], 0)
  expect_equal(result["location"], 0)
  expect_equal(result["location_sign"], 0)
  expect_equal(result["size"], 0)
  expect_equal(result["size_sign"], 0)
  expect_equal(result["shape"], 0)
})

test_that("divergence function handles empty inputs", {
  a <- numeric(0)
  b <- numeric(0)
  
  # Expect error when input vectors are empty
  expect_error(divergence(a, b))
})

test_that("divergence function handles incorrect input types", {
  a <- "not numeric"
  b <- "not numeric"
  
  # Expect error when input types are incorrect
  expect_error(divergence(a, b))
})
