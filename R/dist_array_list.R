#' To initialize a list of arrays, dist_array_list
#'
#' @param dat_res 
#'
#' @return
#' @export
dist_array_list <- function(dat_res) {
  p <- length(dat_res) # number of genes
  n <- length(dat_res[[1]]) # number of donors
  
  gene_ids <- names(dat_res)
  donor_ids <- names(dat_res[[1]])
  
  array_list <- vector("list", length = p)
  names(array_list) <- gene_ids
  
  array_tmp <- array(NA, dim = c(n, n, 6))
  dimnames(array_tmp)[[1]] <- donor_ids
  dimnames(array_tmp)[[2]] <- donor_ids
  dimnames(array_tmp)[[3]] <- c("distance",
                                "location",
                                "location_sign",
                                "size",
                                "size_sign",
                                "shape")
  
  for(i in 1:p){
    array_list[[i]] <- array_tmp
  }
  
  return(array_list)

}
