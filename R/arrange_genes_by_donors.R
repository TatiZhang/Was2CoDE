#' Arrange Genes by Donors
#'
#' This function arranges gene expression data by donors, 
#' reformatting the count matrix genes and donors. 
#' For each gene in the count matrix, it iterates over 
#' each individual in the metadata and extracts the corresponding expression
#' data, organizing it by donors.
#' 
#' 
#' @param count_matrix 
#' @param meta_ind 
#' @param meta_cell 
#'
#' @return dat_res, which is a list (of genes) of lists (of donors)
#' @export
arrange_genes_by_donors <- function(count_matrix, meta_ind, meta_cell) {
  n_gene = nrow(count_matrix)
  gene_ids = rownames(count_matrix)
  dat_res=foreach (i_g = 1:n_gene) %dorng% {
    res_ig = list()
    
    # For loop:
    # For each gene (i_g), the function iterates over each individual 
    # in meta_ind, extracting the corresponding expression data 
    # from the count matrix. 
    # preprocessing to make the data more normally distributed.
    # outputing the residuals that we want
    
    # For each gene (i_g), iterate over each individual in meta_ind
    for (j in 1:nrow(meta_ind)) {
      ind_j = meta_ind$individual[j]  # donor's ID
      w2use = which(meta_cell$individual == ind_j)  # grab all indexes associated with this donor
      # Directly use the count matrix rows corresponding to this gene and columns for the donor
      dat_j = count_matrix[i_g, w2use]
      res_ig[[j]] = dat_j 
    }
    names(res_ig) = as.character(meta_ind$individual)
    res_ig  #output the residual(*)
  }
  names(dat_res) = as.character(gene_ids)
  return(dat_res)
}