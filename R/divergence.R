#' define divergence() using wasserstein1d()
#'
#' @param a Numeric (gene expressions) vector.
#' @param b Numeric (gene expressions) vector.
#' @param p Integer, power parameter for the Wasserstein distance.
#'
#' @return
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
  shape <- (distance^2 - location - size) / (2 * sd(a) * sd(b))
  
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
