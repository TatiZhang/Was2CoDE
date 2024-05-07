#' Title
#'
#' @param dat_res 
#' @param res_ig 
#'
#' @return
#' @export
dist_array_list <- function(dat_res, res_ig) {
    # Validate input lists
    if (length(dat_res) != length(res_ig)) {
        stop("dat_res and res_ig must have the same number of elements")
    }
  
    # Initialize a list to store the results for each gene
    results <- vector("list", length(dat_res))
    names(results) <- names(dat_res)
  
    # Process each gene's data
    for (i in seq_along(dat_res)) {
        # Retrieve data for this gene from both lists
        data_gene <- dat_res[[i]]
        res_gene <- res_ig[[i]]
    
        # Ensure that data dimensions match
        if (length(data_gene) != length(res_gene)) {
            stop("Data mismatch in length for gene ", names(dat_res)[i])
        }
    
        # Initialize a matrix to store results: rows are individuals, columns are metrics
        results_matrix <- matrix(NA, nrow = length(data_gene), ncol = 2)
        colnames(results_matrix) <- c("Abs Difference", "Squared Difference")
    
        # Calculate the metrics
        for (j in seq_along(data_gene)) {
            abs_diff <- abs(data_gene[j] - res_gene[j])
            sq_diff <- (data_gene[j] - res_gene[j])^2
        
            # Store the results
            results_matrix[j, ] <- c(abs_diff, sq_diff)
        }
    
        # Store the matrix in the results list for this gene
        results[[i]] <- results_matrix
    }
  
    return(results)
}
