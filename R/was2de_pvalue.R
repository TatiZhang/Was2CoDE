#' Perform Wilcoxon tests on metrics for Var groups
#'
#' @param dist_list List of distance matrices
#' @param meta_ind Data frame containing metadata information for individuals
#' @param Var dependent variable used for wilcoxon test
#' @return A matrix with rows representing genes and columns representing means of distance metrics for different groups and p-values
#' @export
#'

was2de_pvalue <- function(dist_list, meta_ind, Var,
                          verbose = 0) {
  # Extract gene names from the first distance matrix
  gene_names <- dimnames(dist_list[[1]])[[1]]
  
  # Initialize an empty matrix to store the results
  results_mat <- matrix(nrow = length(gene_names), ncol = 4)
  colnames(results_mat) <- c("mean_dd", "mean_nn", "mean_dn", "p_val")
  rownames(results_mat) <- gene_names
  
  # Extract Var information
  Var <- as.character(meta_ind[[Var]])
  if(is.factor(Var)){
    level_vec <- levels(Var)
  } else {
    level_vec <- unique(Var)
  }
  stopifnot(length(level_vec) == 2)
  
  names(Var) <- meta_ind$individual
  
  # Loop through the distance matrices and categorize the pairs
  for (k in "distance") {
    gene_names <- dimnames(dist_list[[k]])[[1]]
    ind_names <- dimnames(dist_list[[k]])[[2]]
    
    # For each gene
    for (kk in 1:length(gene_names)) {
      if(verbose > 0 & length(gene_names) > 10 & kk %% floor(length(gene_names)/10) == 0) cat('*')
      
      metric_group_dd <- c()
      metric_group_nn <- c()
      metric_group_dn <- c()
      gene <- gene_names[kk]
      
      for (i in 1:(length(ind_names) - 1)) {
        for (j in (i + 1):length(ind_names)) {
          ind_i <- ind_names[i]
          ind_j <- ind_names[j]
          status_i <- Var[ind_i]
          status_j <- Var[ind_j]
          
          if (!is.na(status_i) & !is.na(status_j)) {
            metric <- dist_list[[k]][kk, i, j]
            
            # Group 1: (dementia & dementia)
            if ((status_i == level_vec[2] & status_j == level_vec[2])) {
              metric_group_dd <- c(metric_group_dd, metric)
            }
            # Group 2: (no_dementia & no_dementia)
            if ((status_i == level_vec[1] & status_j == level_vec[1])) {
              metric_group_nn <- c(metric_group_nn, metric)
            }
            # Group 3: (dementia & no_dementia)
            if ((status_i == level_vec[2] & status_j == level_vec[1]) |
                (status_i == level_vec[1] & status_j == level_vec[2])) {
              metric_group_dn <- c(metric_group_dn, metric)
            }
          }
        }
      }

      metric_group_ddnn <- c(metric_group_dd, metric_group_nn)
      
      # remove any NAs
      metric_group_dn <- metric_group_dn[!is.na(metric_group_dn)]
      metric_group_ddnn <- metric_group_ddnn[!is.na(metric_group_ddnn)]
      
      # Perform the Wilcoxon test
      if (length(metric_group_ddnn) > 0 && length(metric_group_dn) > 0) {
        wilcox_test_result <- stats::wilcox.test(metric_group_ddnn, metric_group_dn, alternative = "less")
        p_val <- wilcox_test_result$p.value
      } else {
        p_val <- NA
      }
      
      # Means
      mean_dd <- mean(metric_group_dd, na.rm = TRUE) 
      mean_nn <- mean(metric_group_nn, na.rm = TRUE) 
      mean_dn <- mean(metric_group_dn, na.rm = TRUE) 
      
      # Store the results in the result matrix
      results_mat[gene, ] <- c(mean_dd, mean_nn, mean_dn, p_val)
    }
  }
  
  return(results_mat)
}
