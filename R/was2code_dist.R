#' Compute Wasserstein-2 Distance Between Individuals
#'
#' @param count_input A count matrix (genes x cells)
#' @param meta_cell Cell-level metadata including 'individual' and 'cell_id'
#' @param meta_ind Individual-level metadata including 'individual' and test variable
#' @param var_per_cell Variable to normalize per cell (e.g., "nCount_RNA")
#' @param var2test Binary variable to test at the individual level
#' @param ncores Number of cores to use for parallelization
#' @param k Number of opposite-group individuals to compare each person with (default = NULL for all)
#' @return A list of 3D distance arrays (one per gene)
#' @export
#' 

# modified from https://github.com/Sun-lab/ideas/blob/main/R/ideas_dist.R
was2code_dist <-
  function(count_input, 
           meta_cell, 
           meta_ind, 
           var_per_cell, 
           var2test,
           ncores = 2,
           k = NULL,
           verbose = 0) { 
    set.seed(0)
    
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
    # setting up cores
    # -----------------------------------------------------------------
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("doRNG", quietly = TRUE)) {
      stop("Package 'doRNG' is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package 'foreach' is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    
    # # Register parallel backend with the specified number of cores
    # prev_backend <- foreach::getDoParName()
    # on.exit({
    #   if (prev_backend != "sequential") {
    #     foreach::setDefaultCluster(NULL)
    #     if (prev_backend != "doParallel") {
    #       eval(parse(text = paste0("do", prev_backend, "::register", prev_backend, "()")))
    #     }
    #   }
    # })
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    if(!(is.data.frame(meta_cell))){
      stop("meta_cell should be a data.frame\n")
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
    
    if(any(!c("cell_id", "individual") %in% colnames(meta_cell)))
      stop("cell_id and individual must be column names in meta_cell\n")
    
    if(any(!c(var2test, "individual") %in% colnames(meta_ind)))
      stop("var2test and individual must be column names in meta_ind\n")
    
    if(!is.factor(meta_ind[,var2test])){
      message(paste0("Converting '", var2test, "' to a factor"))
      meta_ind[,var2test] <- as.factor(meta_ind[,var2test])
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
    dat_res <- foreach::foreach(i_g = 1:n_gene, .packages = c("stats")) %dorng% {
      res_ig <- list()
      
      for (j in 1:nrow(meta_ind)) {
        ind_j <- meta_ind$individual[j]
        w2use <- which(meta_cell$individual == ind_j)
        dat_j <- c(count_matrix[i_g, w2use])
        z     <- log10(meta_cell[w2use, var_per_cell, drop = FALSE])
        z     <- as.data.frame(z)
        str_z <- paste(names(z), collapse= "+")
        fmla  <- stats::as.formula(paste("log10(dat_j + 0.5) ~ ", str_z))
        lm_j  <- stats::lm(fmla, data = z)
        base_j <- c(t(lm_j$coefficients) %*% c(1, cov_value))
        res_ig[[j]] <- lm_j$resid + base_j
      }
      names(res_ig) <- as.character(meta_ind$individual)
      
      res_ig
    }
    
    if(verbose > 0){
      for(j in 1:length(dat_res[[1]])){
        print(paste0("Residual for donor ", j, " is: ", 
                     paste0(head(dat_res[[1]][[j]], collapse = ","))))
      }
    }
    
    
    # Add group labels (e.g., case vs control)
    group_labels <- meta_ind[,var2test]
    group_levels <- levels(group_labels)
    if (length(group_levels) < 2) {
      stop("The variable for testing must have at least two levels.")
    }
    
    message(sprintf("Found %d factor levels in %s: %s", 
                    length(group_levels), 
                    var2test, 
                    paste(group_levels, collapse=", ")))
    if (any(!group_levels %in% group_labels)) {
      stop("The group labels are not set correctly")
    }
    
    case_inds <- meta_ind$individual[group_labels == group_levels[1]]
    control_inds <- meta_ind$individual[group_labels == group_levels[2]]
    
    # This is the version that can be parallelized
    if(verbose > 1) print("About to start computing divergences")
    if(ncores > 1) {
      if(verbose > 1) print("Starting the parallelized version")
      dist_array_list <- foreach::foreach(i_g = 1:n_gene) %dorng% {
        res_ig <- dat_res[[i_g]]
        dist_array1 <- array(NA, 
                             dim = c(rep(nrow(meta_ind), 2), 4),
                             dimnames = list(meta_ind$individual, 
                                             meta_ind$individual,
                                             c("was2", "location", "size", "shape")))
        
        for(kk in 1:4) diag(dist_array1[,,kk]) <- 0
        
        for (j_a in 1:nrow(meta_ind)) {
          id_a <- meta_ind$individual[j_a]
          label_a <- group_labels[j_a]
          
          # Determine comparison groups
          # 1. Within group comparison
          same_group <- label_a
          same_inds <- meta_ind$individual[group_labels == same_group]
          same_inds <- same_inds[same_inds != id_a]  # Exclude self
          j_a_candidates <- match(same_inds, meta_ind$individual)
          
          # 2. Between group comparison - handle multiple levels
          other_groups <- group_levels[group_levels != label_a]
          opp_inds <- meta_ind$individual[group_labels %in% other_groups]
          j_b_candidates <- match(opp_inds, meta_ind$individual)
          
          # Sample if k is specified
          if (!is.null(k)) {
            # For opposite group
            if (length(j_b_candidates) > k) {
              set.seed(j_a)  # for reproducibility
              j_b_candidates <- sample(j_b_candidates, k)
            }
            
            # For same group
            if (length(j_a_candidates) > k) {
              set.seed(j_a + nrow(meta_ind))  # Different seed for within-group
              j_a_candidates <- sample(j_a_candidates, k)
            }
          }
          
          # Process opposite group comparisons
          for (j_b in j_b_candidates) {
            if (j_b == j_a || is.na(j_b)) next
            
            res_a <- res_ig[[j_a]]
            res_b <- res_ig[[j_b]]
            
            dist_array1[j_a, j_b,] <- tryCatch(
              divergence(res_a, res_b, verbose = verbose),
              error = function(e) { rep(NA, 4) }
            )
            
            dist_array1[j_b, j_a,] <- dist_array1[j_a, j_b,] * c(1, -1, -1, 1)
          }
          
          # Process within-group comparisons
          for (j_b in j_a_candidates) {
            if (j_b == j_a || is.na(j_b)) next
            
            res_a <- res_ig[[j_a]]
            res_b <- res_ig[[j_b]]
            
            dist_array1[j_a, j_b,] <- tryCatch(
              divergence(res_a, res_b, verbose = verbose),
              error = function(e) { rep(NA, 4) }
            )
            
            dist_array1[j_b, j_a,] <- dist_array1[j_a, j_b,] * c(1, -1, -1, 1)
          }
        }
        dist_array1
      }
    } else {
      if(verbose > 1) print("Starting the single-core version")
      
      # this is the single-core version that is purely for debugging purposes
      dist_array_list <- vector("list", length = n_gene)
      for(i_g in 1:n_gene){
        if(verbose > 4) print("asdf")
        if(verbose > 4) print(paste0("Iteration: ", i_g))
        res_ig <- dat_res[[i_g]]
        dist_array1 <- array(NA, 
                             dim = c(rep(nrow(meta_ind), 2), 4),
                             dimnames = list(meta_ind$individual, 
                                             meta_ind$individual,
                                             c("was2", "location", "size", "shape")))
        
        for(kk in 1:4) diag(dist_array1[,,kk]) <- 0
        
        for (j_a in 1:nrow(meta_ind)) {
          if(verbose > 4) print(paste0("Currently working on: j_a=", j_a))
          id_a <- meta_ind$individual[j_a]
          label_a <- group_labels[j_a]
          if(verbose > 4) print(paste0("label_a is: ", label_a))
          
          # Determine comparison groups
          # 1. Within group comparison
          same_group <- label_a
          same_inds <- meta_ind$individual[group_labels == same_group]
          same_inds <- same_inds[same_inds != id_a]  # Exclude self
          j_a_candidates <- match(same_inds, meta_ind$individual)
          
          # 2. Between group comparison - handle multiple levels
          other_groups <- group_levels[group_levels != label_a]
          if(verbose > 4) print(paste0("Other groups are: ", paste(other_groups, collapse = ", ")))
          
          opp_inds <- meta_ind$individual[group_labels %in% other_groups]
          j_b_candidates <- match(opp_inds, meta_ind$individual)
          
          if(verbose > 4) {
            print(paste0("j_b_candidates are: ", paste0(j_b_candidates, collapse = ", ")))
            print(paste0("Length of j_b_candidates: ", length(j_b_candidates)))
            print(paste0("j_a_candidates are: ", paste0(j_a_candidates, collapse = ", ")))
            print(paste0("Length of j_a_candidates: ", length(j_a_candidates)))
          }
          
          # Sample if k is specified
          if (!is.null(k)) {
            # For opposite group
            if (length(j_b_candidates) > k) {
              j_b_candidates <- sample(j_b_candidates, k)
            }
            
            # For same group
            if (length(j_a_candidates) > k) {
              j_a_candidates <- sample(j_a_candidates, k)
            }
          }
          
          # Process opposite group comparisons
          for (j_b in j_b_candidates) {
            if(verbose > 4) print(paste0("Comparing j_a to opposite group: j_b=", j_b))
            if(verbose > 4) print(paste0("j_b is NA: ", is.na(j_b)))
            if (j_b == j_a || is.na(j_b)) next
            
            res_a <- res_ig[[j_a]]
            res_b <- res_ig[[j_b]]
            
            if(verbose > 0){
              print(paste0("Between-group comparison: j_a is: ", j_a, ": j_b is: ", j_b))
              print(paste0("head of res_a: ", paste0(head(res_a), collapse = ", ")))
              print(paste0("head of res_b: ", paste0(head(res_b), collapse = ", ")))
            }
            
            dist_array1[j_a, j_b,] <- tryCatch(
              divergence(res_a, res_b, verbose = verbose),
              error = function(e) { rep(NA, 4) }
            )
            
            dist_array1[j_b, j_a,] <- dist_array1[j_a, j_b,] * c(1, -1, -1, 1)
            if(verbose > 4) print("Completing opposite group comparison")
          }
          
          # Process within-group comparisons
          for (j_b in j_a_candidates) {
            if(verbose > 4) print(paste0("Comparing j_a to within group: j_b=", j_b))
            if (j_b == j_a || is.na(j_b)) next
            
            res_a <- res_ig[[j_a]]
            res_b <- res_ig[[j_b]]
            
            if(verbose > 0){
              print(paste0("Within-group comparison: j_a is: ", j_a, ": j_b is: ", j_b))
              print(paste0("head of res_a: ", paste0(head(res_a), collapse = ", ")))
              print(paste0("head of res_b: ", paste0(head(res_b), collapse = ", ")))
            }
            
            dist_array1[j_a, j_b,] <- tryCatch(
              divergence(res_a, res_b, verbose = verbose),
              error = function(e) { rep(NA, 4) }
            )
            
            dist_array1[j_b, j_a,] <- dist_array1[j_a, j_b,] * c(1, -1, -1, 1)
            if(verbose > 4) print("Completing within-group comparison")
          }
        }
        dist_array_list[[i_g]] <- dist_array1
      }
    }
    
    if(verbose > 4) print("Finishing")

    names(dist_array_list) <- gene_ids
    return(dist_array_list)
  }