#' Compute Wasserstein-2 decomposition for Gene Expression Data
#' Using the divergence function many many many times, 
#' and putting the correct values into the empty 
#' dist_array_list (list of arrays)
#' and then reorganize it to result_array_list
#'
#' @param count_input 
#' @param meta_cell 
#' @param meta_ind 
#' @param var_per_cell 
#' @param var2test 
#' @param var2test_type 
#'
#' @return A list of 3-dimensional arrays, each corresponding to a distance metric (distance, location, location_sign, size, size_sign, shape).
#' Each array has dimensions (number of genes, number of donors, number of donors) and contains the respective metric values.
#'
#' @export
ideas_dist_custom <-
  function(count_input, # the input should be "genes" by "cells"
           meta_cell, meta_ind, var_per_cell, var2test, 
           var2test_type = c("binary", "continuous"), 
           verbose = 0) { 
    # -----------------------------------------------------------------
    # TZ:
    # initializes the var2test_type, validates the structure and content of 
    # meta_cell and meta_ind, which are 2 data frames,
    # checking for necessary columns like cell_id, individual, and var_per_cell
    # ----------------------------------------------------------------- 
    var2test_type   = var2test_type[1]
    
    if(!(is.data.frame(meta_cell))){
      stop("meta_cell should be a data.frame\n")
    }
    
    if(is.data.table(meta_cell)){
      meta_cell = as.data.frame(meta_cell)
    }
    
    if(!(is.data.frame(meta_ind))){
      stop("meta_ind should be a data.frame\n")
    }
    
    if(is.data.table(meta_ind)){
      meta_ind = as.data.frame(meta_ind)
    }
    
    # -----------------------------------------------------------------
    # TZ:
    # validates the structure and content of meta_cell and meta_ind, 
    # checking for necessary columns including cell_id, individual,
    # and the variable per cell (var_per_cell). 
    # ensures that these columns contain unique cell and individual IDs 
    # that match across the metadata and count matrices. 
    # aligning the gene expression data with the corresponding metadata 
    # for each cell and individual.
    # -----------------------------------------------------------------
    
    count_matrix = count_input
    
    if(! is.matrix(count_matrix)){
      stop("count_matrix is not a matrix\n")
    }
    
    
    
    n_cell = ncol(count_matrix)
    n_gene = nrow(count_matrix)
    
    gene_ids = rownames(count_matrix)
    cell_ids = colnames(count_matrix)
    
    if(is.null(gene_ids)){
      stop("count_matrix should have row names for gene ids\n")
    }
    
    if(is.null(cell_ids)){
      stop("count_matrix should have col names for cell ids\n")
    }
    
    if(length(unique(gene_ids)) != n_gene){
      stop("row names of count_matrix (gene ids) are not unique\n")
    }
    
    if(length(unique(cell_ids)) != n_cell){
      stop("col names of count_matrix (cell ids) are not unique\n")
    }
    
    message(sprintf("the count_matrix includes %d genes in %d cells\n", 
                    n_gene, n_cell))
    
    
    if(any(meta_cell$cell_id != colnames(count_matrix))){
      stop("cell_id in meta_cell do not match colnames of count_matrix\n")
    }
    
    # -----------------------------------------------------------------
    # check other aspects of meta_cell
    # -----------------------------------------------------------------    
    
    columns.meta.cell = c("cell_id", "individual", var_per_cell)
    
    if(! all(columns.meta.cell %in% names(meta_cell))){
      str1 = paste(columns.meta.cell, collapse=", ")
      stop(sprintf("names of meta_cell should contain %s\n", str1))
    }
    
    if(length(unique(meta_cell$cell_id)) != nrow(meta_cell)){
      stop("the cell_id's in meta_cell are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # check the input data of meta_ind
    # -----------------------------------------------------------------
    
    columns.meta.ind = c("individual", var2test)
    
    if(! all(columns.meta.ind %in% names(meta_ind))){
      str1 = paste(columns.meta.ind, collapse=", ")
      stop(sprintf("names of meta_ind should contain %s\n", str1))
    }
    
    if(! setequal(meta_cell$individual, meta_ind$individual)){
      stop("the individual ids in meta_cell and meta_ind do not match\n")
    }
    
    if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
      stop("the individual ids in meta_ind are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # TZ:
    # estimates the distribution(distance) for each gene and individual 
    # using Kernel Density Estimation (KDE)
    # -----------------------------------------------------------------
    
    message("estimating distribution for each gene and each individual by kde\n")
    dat_res <- arrange_genes_by_donors(count_matrix, 
                                       meta_ind, 
                                       meta_cell)
    dist_array_list <- dist_array_list(dat_res)
    
    if(verbose > 0) print("Analyzing genes")
    
    for(i in 1:length(dist_array_list)){
      if(verbose > 1 && length(dist_array_list) > 10 && i %% 
         floor(length(dist_array_list) /10) == 0) cat('*')
      if(verbose > 2) print(paste("Gene: ", names(dat_res)[i]))
      
      # Set diagonal elements for each matrix in the 3D array to 0
      for (k in 1:6) {
        diag(dist_array_list[[i]][,,k]) <- 0
        # For loop:
        # For each pair of donors
        # compute the wasserstein distance b/w 2 distributions
        #######
        # this is where I start changing the code
        
        var_levels <- levels(as.factor(meta_ind[[var2test]]))
        if (length(var_levels) != 2) stop("var2test must have exactly two levels.")
        
        level_1 <- var_levels[1]
        level_2 <- var_levels[2]
          
          #v##v#v
          # instead of this:
          # for (j_b in (j_a+1):nrow(meta_ind)) {
          #v##v##v
          # I want to instead:
          # Randomly sample 1 AD donor, and 1 non-AD donor (neither of these two people can be j_a itself)
          # And then, my for (j_b in (these two randomly selected donors))
          #######
          #         
          #         res_b = dat_res[[i]][[j_b]]
          #         
          #         dist_array_list[[i]][j_a, j_b,] = tryCatch(
          #           divergence(res_a, res_b), #calculation
          #           error = function(e) { NA }
          #         )
          #         
          #         dist_array_list[[i]][j_b, j_a,] = dist_array_list[[i]][j_a, j_b,]
          #       }
          #     }
          #     dist_array_list[[i]]
          #   }}
          # 
          # # Calculating the number of NA values in each element of dist_array_list
          # nNA = sapply(dist_array_list, function(x){sum(is.na(c(x)))})
          # table(nNA)
          # Identify AD and non-AD donors, excluding the current donor j_a
        for (j_a in seq_len(nrow(meta_ind))) {
          res_a <- dat_res[[i]][[j_a]]
          
          # Identify AD and non-AD donors, excluding the current donor j_a
          ad_donors <- which(meta_ind[[var2test]] == level_1 & meta_ind$individual != meta_ind$individual[j_a])
          non_ad_donors <- which(meta_ind[[var2test]] == level_2 & meta_ind$individual != meta_ind$individual[j_a])
          
          if (length(ad_donors) > 0 && length(non_ad_donors) > 0) {
            sampled_ad <- sample(ad_donors, 1)
            sampled_non_ad <- sample(non_ad_donors, 1)
            
            # Loop over the randomly sampled donors
            for (j_b in c(sampled_ad, sampled_non_ad)) {
              res_b <- dat_res[[i]][[j_b]]
              
              # Calculate the divergence and update the distance array
              dist_array_list[[i]][j_a, j_b, ] <- tryCatch(
                divergence(res_a, res_b),
                error = function(e) { NA }
              )
              
              # Symmetric assignment
              dist_array_list[[i]][j_b, j_a, ] <- dist_array_list[[i]][j_a, j_b, ]
            }
          }
        }
      }
    }
    
    
    nNA <- sapply(dist_array_list, function(x) { sum(is.na(c(x))) })
    
    if(verbose > 0) print("Finalizing output")
    
    result_array_list <- result_array_list(dist_array_list, 
                                           meta_ind)
    return(result_array_list)
  }
