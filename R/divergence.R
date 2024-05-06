#' define divergence() using wasserstein1d()
#'
#' @param a Numeric vector.
#' @param b Numeric vector.
#' @param p Integer, power parameter for the Wasserstein distance.
#' @param wa weights for vector a
#' @param wb weights for vector b
#'
#' @return
#' @export
#'
divergence <- function(a, b, p = 2, wa = NULL, wb = NULL) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("Both 'a' and 'b' must be numeric vectors.")
  }
  if (length(a) != length(b)) {
    stop("'a' and 'b' must be of the same length.")
  }
  distance <- transport::wasserstein1d(a, b, p = p, wa = wa, wb = wb)
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
  names(result) <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  
  if (length(result) != 6) {
    stop("Result list from divergence has incorrect length")
  }
  
  return(result)
}
