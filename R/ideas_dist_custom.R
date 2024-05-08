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
#' @param per_cell_adjust 
#'
#' @return
#' @export
ideas_dist_custom <-
  function(count_input, # the input should be "genes" by "cells"
           meta_cell, meta_ind, var_per_cell, var2test, 
           var2test_type = c("binary", "continuous"), 
           per_cell_adjust = c("NB", "both")) { 
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
      
        check_count <- function(v){any(v != round(v) | v < 0)}
        not_count = apply(count_matrix, 1, check_count)
        
        if(any(not_count)){
          str1 = "count_matrix should only include non-negative integers"
          str1 = sprintf("%s, violation in row %d\n", str1, which(not_count)[1])
          stop(str1)
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
      
      for(i in 1:length(dist_array_list)){
     
        # Set diagonal elements for each matrix in the 3D array to 0
        for (k in 1:6) {
          diag(dist_array_list[[i]][,,k]) <- 0
        # For loop:
        # For each pair of donors
        # compute the wasserstein distance b/w 2 distributions
          for (j_a in 1:(nrow(meta_ind)-1)) {
          res_a = dat_res[[i]][[j_a]]
          for (j_b in (j_a+1):nrow(meta_ind)) {
            res_b = dat_res[[i]][[j_b]]

            dist_array_list[[i]][j_a, j_b,] = tryCatch(
              divergence(res_a, res_b), #calculation
              error = function(e) { NA }
            )

            dist_array_list[[i]][j_b, j_a,] = dist_array_list[[i]][j_a, j_b,]
          }
        }
        dist_array_list[[i]]
      }}



    # -----------------------------------------------------------------
    # TZ:
    # conclusion
    # compiles the pairwise distances between individuals for each gene
    # into a 3-dimensional array (dist_array).
    # -----------------------------------------------------------------

    # Calculating the number of NA values in each element of dist_array_list
    nNA = sapply(dist_array_list, function(x){sum(is.na(c(x)))})
    table(nNA)
    
    result_array_list <- result_array_list(dist_array_list, 
                                           meta_ind)
    return(result_array_list)
  }
