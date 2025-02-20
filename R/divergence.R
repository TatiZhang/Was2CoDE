#' define divergence() using wasserstein1d()
#'
#' @param a Numeric (gene expressions) vector.
#' @param b Numeric (gene expressions) vector.
#' @param p Integer, power parameter for the Wasserstein distance.
#'
#' @return A numeric vector containing four divergence metrics: distance, location difference, 
#'size difference, and shape difference.
#' @export
#'
divergence <- function(a, b, p=2) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("Both 'a' and 'b' must be numeric vectors.")
  }
  
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  
  if(length(a) == 0 || length(b) == 0) {
    stop("Both 'a' and 'b' have no non-NA values.")
  }
  
  was2 <- transport::wasserstein1d(a, b, p, wa = NULL, wb = NULL)
  location <- mean(a) - mean(b)
  
  sd_a <- stats::sd(a)
  sd_b <- stats::sd(b)
  size <- sd_a - sd_b
  quantiles <- seq(0, 1, length.out = 100)
  quantiles_a <- stats::quantile(a, probs = quantiles)
  quantiles_b <- stats::quantile(b, probs = quantiles)
  quantile_cor_ab <- stats::cor(quantiles_a, quantiles_b)
  shape <- abs(2 * sd_a * sd_b * (1 - quantile_cor_ab))
  
  result <- c(
    was2,
    location,
    size,
    shape
  )
  

  if (length(result) != 4) {
    stop("Result list from divergence has incorrect length")
  }
  
  return(result)
}
