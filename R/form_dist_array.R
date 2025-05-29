#' Form a 3D Distance Array from a List of Distance Matrices
#'
#' This internal function constructs a 3D array from a list of distance matrices,
#' where each matrix corresponds to a specific gene and contains pairwise distances
#' between donors.
#'
#' @param dist_array_list A named list of distance matrices. Each element of the list
#'   corresponds to a gene and contains a 3D array where the third dimension includes
#'   the `"was2"` distance.
#'
#' @return A 3D array with dimensions \code{(n_genes, n_donors, n_donors)},
#'   where \code{n_genes} is the number of genes and \code{n_donors} is the number
#'   of donors. The array is indexed by gene and contains the `"was2"` distances
#'   between donors.
#'
#' @examples
#' # Example usage (assuming `dist_list` is a valid input)
#' dist_array <- .form_dist_array(dist_list)
#'
#' @keywords internal
.form_dist_array <- function(dist_array_list) {
  # Extract gene and donor names
  gene_vec <- names(dist_array_list)
  donor_vec <- dimnames(dist_array_list[[1]])[[1]]
  ndonors <- dim(dist_array_list[[1]])[1]
  
  # Initialize an array to store distances
  dist_array <- array(
    NA,
    dim = c(length(gene_vec), length(donor_vec), length(donor_vec)),
    dimnames = list(gene_vec, donor_vec, donor_vec)
  )
  
  # Populate the array with values
  for (gene in gene_vec) {
    dist_array[gene, donor_vec, donor_vec] <- dist_array_list[[gene]][donor_vec, donor_vec, "was2"]
  }
  
  return(dist_array)
}
