#' Arrange Genes by Donors
#'
#' This function arranges gene expression data by donors, 
#' reformatting the count matrix genes and donors. 
#' For each gene in the count matrix, it iterates over 
#' each individual in the metadata and extracts the corresponding expression
#' data, organizing it by donors.
#' 
#' 
#' @param count_matrix A matrix with gene expression data (genes as rows, cells as columns)
#' @param meta_ind A data frame with individual/donor metadata
#' @param meta_cell A data frame with cell metadata
#'
#' @return dat_res, which is a list (of genes) of lists (of donors)
#' @export
arrange_genes_by_donors <- function(count_matrix, meta_ind, meta_cell) {
  n_gene <- nrow(count_matrix)
  gene_ids <- rownames(count_matrix)
  
  # Skip parallel and use sequential processing
  # This version doesn't attempt parallel processing which is causing problems
  foreach::registerDoSEQ()
  
  # Pre-compute the cell indices for each individual to avoid recalculating
  individual_cells <- list()
  for (ind in meta_ind$individual) {
    individual_cells[[as.character(ind)]] <- which(meta_cell$individual == ind)
  }
  
  # Process data sequentially
  dat_res <- list()
  
  # Process each gene sequentially
  for (i_g in 1:n_gene) {
    res_ig <- list()
    
    # For each individual, extract the gene expression data
    for (j in 1:nrow(meta_ind)) {
      ind_j <- meta_ind$individual[j]  # donor's ID
      w2use <- individual_cells[[as.character(ind_j)]]  # Use pre-computed indices
      
      # Extract expression data for this gene and donor
      dat_j <- count_matrix[i_g, w2use]
      res_ig[[j]] <- dat_j
    }
    
    names(res_ig) <- as.character(meta_ind$individual)
    dat_res[[i_g]] <- res_ig
  }
  
  # Set names for the gene-level list
  names(dat_res) <- as.character(gene_ids)
  
  return(dat_res)
}