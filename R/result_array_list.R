#' Title
#'
#' @param dist_array_list 
#' @param meta_ind 
#'
#' @return
#' @export
#'
#' @examples
result_array_list <- function(dist_array_list, meta_ind) {
  # n_gene, meta_ind, and gene_ids are passed to the function
  n_gene <- length(dist_array_list)
  n_ind <- nrow(meta_ind)  
  result_array_list <- lapply(1:6, function(kk) {
    array(
      dim = c(n_gene, n_ind, n_ind),
      dimnames = list(gene_ids, meta_ind$individual, meta_ind$individual)
    )
  })
  names(result_array_list) <- c("distance", "location", "location_sign", "size", "size_sign", "shape")
  
  for (kk in 1:6) {
    for (i in 1:n_gene) {
      result_array_list[[kk]][i, , ] <- dist_array_list[[i]][, , kk]
    }
  }
  
  #dimensions of the first element and a subset of data for verification
  cat("Dimensions of the first metric array:", dim(result_array_list[[1]]), "\n")
  cat("Data slice from the first metric array for the first gene:\n",result_array_list[[1]][1, 1:2, 1:2],"\n")
  
  return(result_array_list)
}

