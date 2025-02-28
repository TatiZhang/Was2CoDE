# modified from https://github.com/Sun-lab/ideas/blob/main/R/ideas_dist.R
was2code_dist <-
  function(count_input, 
           meta_cell, 
           meta_ind, 
           var_per_cell, 
           var2test) {
    
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
    # check the input data of count_input,when fit_method != dca_direct
    # -----------------------------------------------------------------
    count_matrix = count_input
    
    if(! is.matrix(count_matrix)){
      stop("count_matrix is not a matrix\n")
    }
    
    check_count <- function(v){any(v != round(v) | v < 0)}
    not_count <- apply(count_matrix, 1, check_count)
    
    if(any(not_count)){
      str1 <- "count_matrix should only include non-negative integers"
      str1 <- sprintf("%s, violation in row %d\n", str1, which(not_count)[1])
      stop(str1)
    }
    
    n_cell <- ncol(count_matrix)
    n_gene <- nrow(count_matrix)
    
    gene_ids <- rownames(count_matrix)
    cell_ids <- colnames(count_matrix)
    
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
    
    # -----------------------------------------------------------------
    # check cell_id order of meta_cell, when fit_method != dca_direct
    # -----------------------------------------------------------------
    
    if(any(meta_cell$cell_id != colnames(count_matrix))){
      stop("cell_id in meta_cell do not match colnames of count_matrix\n")
    }
    
    # -----------------------------------------------------------------
    # check other aspects of meta_cell
    # -----------------------------------------------------------------    
    
    columns_meta_cell = c("cell_id", "individual", var_per_cell)
    
    if(! all(columns_meta_cell %in% names(meta_cell))){
      str1 <- paste(columns_meta_cell, collapse=", ")
      stop(sprintf("names of meta_cell should contain %s\n", str1))
    }
    
    if(length(unique(meta_cell$cell_id)) != nrow(meta_cell)){
      stop("the cell_id's in meta_cell are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # check the input data of meta_ind
    # -----------------------------------------------------------------
    
    columns_meta_ind <- c("individual", var2test)
    
    if(! all(columns_meta_ind %in% names(meta_ind))){
      str1 <- paste(columns_meta_ind, collapse=", ")
      stop(sprintf("names of meta_ind should contain %s\n", str1))
    }
    
    if(!setequal(meta_cell$individual, meta_ind$individual)){
      stop("the individual ids in meta_cell and meta_ind do not match\n")
    }
    
    if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
      stop("the individual ids in meta_ind are not unique\n")
    }
    
    if(any(meta_cell[,var_per_cell, drop = FALSE] <= 0.0)){
      str1 <- "the variables listed in 'var_per_cell' will be log transformed,"
      stop(paste(str1, "so they must be positive."))
    }
    
    i_g <- 0 # for debugging reasons
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using kde
    # -----------------------------------------------------------------
    message("estimating distribution for each gene and each individual by kde\n")
    cov_value <- apply(log10(meta_cell[,var_per_cell,drop=FALSE]), 2, stats::median)
    dat_res <- foreach (i_g = 1:n_gene) %dorng% {
      res_ig <- list()
      
      for (j in 1:nrow(meta_ind)) {
        ind_j <- meta_ind$individual[j]
        w2use <- which(meta_cell$individual == ind_j)
        dat_j <- c(count_matrix[i_g, w2use])
        z     <- log10(meta_cell[w2use, var_per_cell, drop = FALSE])
        z     <- as.data.frame(z)
        str_z <- paste(names(z), collapse= "+")
        fmla  <- as.formula(paste("log10(dat_j + 0.5) ~ ", str_z))
        lm_j  <- stats::lm(fmla, data = z)
        base_j <- c(t(lm_j$coefficients) %*% c(1, cov_value))
        res_ig[[j]] <- lm_j$resid + base_j
      }
      names(res_ig) <- as.character(meta_ind$individual)
      res_ig
    }
    
    
    dist_array_list <- foreach (i_g = 1:n_gene) %dorng% {
      res_ig <- dat_res[[i_g]]
      dist_array1 <- array(NA, 
                           dim = c(rep(nrow(meta_ind), 2), 4),
                           dimnames = list(meta_ind$individual, 
                                           meta_ind$individual,
                                           c("was2", "location", "size", "shape")))
      
      for(kk in 1:4) diag(dist_array1[,,kk]) <- 0
      
      for (j_a in 1:(nrow(meta_ind)-1)) {
        res_a <- res_ig[[j_a]]
        
        for (j_b in (j_a+1):nrow(meta_ind)) {
          res_b <- res_ig[[j_b]]
          
          dist_array1[j_a, j_b,] <- tryCatch(
            divergence(res_a, res_b),
            error = function(e) { rep(NA, 4) }
          )
          
          # the row is donor 1, the column is donor 2
          # we need to flip the mean and standard deviation (which are signed, elements 2 and 3)
          dist_array1[j_b, j_a,] <- dist_array1[j_a, j_b,] * c(1,-1,-1,1)
        }
      }
      dist_array1
    }
    names(dist_array_list) <- gene_ids
    
    return(dist_array_list)
  }