#' Generate Result Array List
#'
#' @param dist_array_list 
#' @param meta_ind 
#'
#' @return A list of 3-dimensional arrays, one for each metric (distance, location, location_sign, size, size_sign, shape).
#' Each array has dimensions (number of genes, number of donors, number of donors) and contains the respective metric values.
#'
#' @export
#'
result_array_list <- function(dist_array_list, meta_ind) {
  # n_gene, meta_ind, and gene_ids are passed to the function
  n_gene <- length(dist_array_list)
  gene_ids <- names(dist_array_list)
  
  n_ind <- nrow(meta_ind)  
  result_array_list <- lapply(1:6, function(kk) {
    array(
      dim = c(n_gene, n_ind, n_ind),
      dimnames = list(gene_ids, meta_ind$individual, meta_ind$individual)
    )
  })
  names(result_array_list) <- c("distance", 
                                "location", 
                                "location_sign", 
                                "size", 
                                "size_sign", 
                                "shape")
  
  for (kk in 1:6) {
    for (i in 1:n_gene) {
      result_array_list[[kk]][i, , ] <- dist_array_list[[i]][, , kk]
      dimnames(result_array_list[[kk]]) = list(gene_ids, 
                                              meta_ind$individual, 
                                              meta_ind$individual)
    }
  }

  return(result_array_list)
}

