#' define divergence() using wasserstein1d()
#'
#' @param a Numeric (gene expressions) vector.
#' @param b Numeric (gene expressions) vector.
#' @param p Integer, power parameter for the Wasserstein distance.
#'
#' @return A numeric vector containing six divergence metrics: distance, location difference, 
#' location sign, size difference, size sign, and shape difference.
#' @export
#'
divergence <- function(a, b, p=2) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("Both 'a' and 'b' must be numeric vectors.")
  }
  distance <- transport::wasserstein1d(a, b, p, wa = NULL, wb = NULL)
  location <- (mean(a) - mean(b))^2
  location_sign <- mean(a) - mean(b)
  size <- (sd(a) - sd(b))^2
  size_sign <- sd(a) - sd(b)
  quantiles <- seq(0, 1, length.out = 100)
  quantiles_a <- quantile(a, probs = quantiles)
  quantiles_b <- quantile(b, probs = quantiles)
  quantile_cor_ab <- cor(quantiles_a, quantiles_b)
  shape <- abs(2 * sd(a) * sd(b) * (1 - quantile_cor_ab))
  
  result <- c(
    distance,
    location,
    location_sign,
    size,
    size_sign,
    shape
  )
  

  if (length(result) != 6) {
    stop("Result list from divergence has incorrect length")
  }
  
  return(result)
}
