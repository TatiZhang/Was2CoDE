#' Compute standardized logFC-style effect sizes from Was2CoDE output
#'
#' @param dist_list A list of distance arrays (from was2code_dist)
#' @param case_idx Integer vector of row indices corresponding to cases
#' @param control_idx Integer vector of row indices corresponding to controls
#'
#' @return A matrix with genes as rows and distance components as columns:
#'         columns = c("was2", "location", "size", "shape")
#' @export
was2code_lfc <- function(dist_list, case_idx, control_idx) {
  dist_components <- c("was2", "location", "size", "shape")
  lfc_mat <- matrix(NA, nrow = length(dist_list), ncol = length(dist_components))
  rownames(lfc_mat) <- names(dist_list)
  colnames(lfc_mat) <- dist_components
  
  # Loop over genes
  for (j in seq_along(dist_list)) {
    dist_gene <- dist_list[[j]]
    
    # loop over components
    for (comp_idx in seq_along(dist_components)) {
      comp <- dist_components[comp_idx]
      tmp <- dist_gene[,,comp]
      diag(tmp) <- NA
      
      numerator <- mean(tmp[case_idx, control_idx], na.rm = TRUE)
      
      sd1 <- sd(tmp[case_idx, case_idx], na.rm = TRUE)
      sd2 <- sd(tmp[control_idx, control_idx], na.rm = TRUE)
      sd1_adj <- sd1 / sqrt(length(case_idx))
      sd2_adj <- sd2 / sqrt(length(control_idx))
      denominator <- sqrt(sd1_adj^2 + sd2_adj^2)
      
      #  Avoid division by 0 or NA
      if (is.finite(numerator) && is.finite(denominator) && denominator > 0) {
        lfc_mat[j, comp_idx] <- numerator / denominator
      } else {
        lfc_mat[j, comp_idx] <- NA
      }
    }
  }
  
  return(lfc_mat)
}

