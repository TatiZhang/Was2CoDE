# code originally from https://github.com/Sun-lab/ideas/blob/main/R/permanova.R and https://github.com/Sun-lab/ideas/blob/main/R/permanova_utilities.R
was2code_permanova <- function(dist_list, 
                               meta_ind, 
                               var2test, 
                               var2adjust = NULL, 
                               n_perm = 999, 
                               r_seed = 2020, 
                               residulize_x = FALSE, 
                               delta = 0.5,
                               ncores = 1,
                               verbose = 0){
  
  # -----------------------------------------------------------------
  # check dist_array
  # -----------------------------------------------------------------
  n_genes <- length(dist_list)
  
  dist_array <- .form_dist_array(dist_list)
  
  wNA <- which(apply(dist_array, 1, anyNA))
  if(length(wNA) > 0){
    message(sprintf("skip %d gene(s) with NA in the dist_array\n", length(wNA)))
    dist_array <- dist_array[-wNA,,]
  }
  
  # -----------------------------------------------------------------
  # check the input data of meta_ind
  # -----------------------------------------------------------------
  
  if(!(is.data.frame(meta_ind))){
    stop("meta_ind should be a data.frame\n")
  }
  
  columns_meta_ind <- c("individual", var2test, var2adjust)
  
  if(! all(columns_meta_ind %in% names(meta_ind))){
    str1 <- paste(columns_meta_ind, collapse=", ")
    stop("names of meta_ind should conttain: ",str1)
  }
  
  if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
    stop("the individual ids in meta_ind are not unique\n")
  }
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # -----------------------------------------------------------------
  # check the var2test
  # -----------------------------------------------------------------
  
  x <- meta_ind[[var2test]]
  
  if(any(is.na(x))){
    stop("variable to test has NA values.\n")
  }
  
  nUniq_x <- length(unique(x))
  
  message(sprintf("testing for '%s', a %s variable\n", var2test, "binary"))
  # make sure x is a vector of 0 or 1. 
  if(is.character(x)){ x <- as.factor(x) }
  if(is.factor(x)){ x <- as.numeric(x) }
  if(is.numeric(x)){ x <- as.numeric(x == max(x)) }
  
  # -----------------------------------------------------------------
  # start testing
  # -----------------------------------------------------------------
  
  set.seed(r_seed)
  x_perm <- matrix(rep(x, times = n_perm), ncol = n_perm)
  x_perm <- apply(x_perm, 2, sample, size = length(x))
  
  # if there is no cavariate, use standard permanova
  if(is.null(var2adjust)){
    F_ob <- .calc_F_manova(dist_array, label = x)
    F_perm_list <- apply(x_perm, 2, function(perm1){ 
      .calc_F_manova(dist_array, label = perm1)
    }, simplify = FALSE)
    F_perm <- do.call(cbind, F_perm_list)
    
  } else {
    fm1 <- stats::as.formula(paste("~", paste(var2adjust, collapse=" + ")))
    z <- stats::model.matrix(fm1, data = meta_ind)
    
    if(residulize_x){
      resid_perm <- matrix(NA_real_, nrow(meta_ind), n_perm)
      
      # step1: fit
      m1 <- stats::glm(x ~ -1 + z, family = stats::binomial(link="logit"))  # logistic model
      fitted_x <- stats::fitted(m1)
      resid_x <- x - fitted_x
      
      # step2: permutation
      ip <- 0 # index for usable permutations
      id <- 0 # working index
      sd_e <- stats::sd(resid_x) # expected sd
      
      while(ip < n_perm){
        if(verbose == 1 && ip %% floor(n_perm/10) == 0) cat('*')
        if(verbose >= 2) print(paste0("Working on permutation: ", ip))
        
        id <- id + 1
        perm_x <- stats::rbinom(length(fitted_x), 1, prob = fitted_x)
        
        m2 <- tryCatch(stats::glm(perm_x ~ -1 + z, family = stats::binomial(link="logit")), 
                       warning = function(w) { NULL }, 
                       error = function(e) { NULL})
        
        if(! is.null(m2)){
          resid_p_i <- perm_x - fitted(m2) 
          sd_resid_i <- stats::sd(resid_p_i)
          
          if(sd_resid_i > (1-delta)*sd_e & sd_resid_i < (1+delta)*sd_e){
            ip <- ip + 1
            resid_perm[,ip] <- resid_p_i
          }
        }
      }
      
      resid_x <- matrix(resid_x, ncol=1)
      F_ob <- .calc_F_permanovaSZ(dist_array, Rs = resid_x, z = z)
      F_perm  <- .calc_F_permanovaSZ(dist_array, Rs = resid_perm, z = z)
      
    } else {
      F_ob <- .calc_F_permanovaSZ(dist_array, Rs = x, z = z)
      F_perm <- .calc_F_permanovaSZ(dist_array, Rs = x_perm, z = z)
    }
  }
  
  F_perm0 <- apply(F_perm, 2, function(x){return(x-F_ob)} )
  
  if(length(F_ob) == 1){ F_perm0 <- t(F_perm0) }
  
  pval <- rowMeans(cbind(F_perm0, rep(1, nrow(F_perm0))) >= 0)
  names(pval) <- names(dist_list)
    
  if(length(wNA) == 0){
    pval_all <- pval
  }else{
    pval_all <- rep(NA, n_genes)
    pval_all[-wNA] <- pval
  }
  
  return(pval_all)
}

#############

.form_dist_array <- function(dist_list){
  gene_vec <- names(dist_list)
  donor_vec <- dimnames(dist_list[[1]])[[1]]
  
  ndonors <- dim(dist_list[[1]])[1]
  dist_array <- array(NA_real_, 
                      dim = c(length(gene_vec), length(donor_vec), length(donor_vec)),
                      dimnames = list(gene_vec, 
                                      donor_vec,
                                      donor_vec))
  
  for(gene in gene_vec){
    dist_array[gene,donor_vec,donor_vec] <- dist_list[[gene]][donor_vec,donor_vec,"was2"]
  }
  
  return(dist_array)
}

#############

.cal_G <- 
  function(m){# m is a distance matrix
    n <- nrow(m)
    centerM <- diag(n) - 1/n
    G <- -0.5 * centerM %*% (m*m) %*% centerM
    eG <- eigen(G, symmetric = TRUE)
    G <- eG$vector %*% diag(pmax(0, eG$values)) %*% t(eG$vector)
    G
  }

.calc_F_manova <- 
  function(dist_array, label){
    
    uLabel <- unique(label)
    a <- length(uLabel)
    N <- length(label)
    
    n_genes <- dim(dist_array)[1]
    d2 <- dist_array*dist_array
    
    # change the dimension of 3-way array to 2-way array
    dim(d2) <- c(n_genes, N^2)
    sst <- rowSums(d2)/N
    
    # this is an N x N matrix with 0 for between group entries and 
    # 1/group_size for within group entries
    divid_vector <- matrix(0,N,N)
    
    for(ik in 1:a){
      cur_index <- which(label==uLabel[ik])
      divid_vector[cur_index,cur_index] <- 1/length(cur_index)
    }
    divid_vector2 <- c(divid_vector)
    
    ssw <- d2 %*% divid_vector2
    
    Fstat <- ((sst-ssw)*(N-a))/(ssw*(a-1))
    return(Fstat)
  }

.calTrace <- function(G,H){
  sum(diag(H%*%G%*%H))
}

.calc_F_permanovaSZ <- 
  function(dist_array, Rs, z = NULL){
    
    n <- dim(dist_array)[2]
    G <- apply(dist_array, 1, .cal_G)
    dim(G) <- c(n,n,dim(dist_array)[1])
    
    if(is.vector(Rs)){
      Rs <- matrix(Rs, ncol=1)
    } 
    
    i <- 0
    F_stats <- foreach::foreach (i = 1:ncol(Rs), .combine='cbind') %dorng% {
      
      if(is.null(z)){
        xz <- Rs[,i]
      }else{
        xz <- unique(cbind(Rs[,i], z), MARGIN=2)
      }
      
      H <- xz %*% solve(crossprod(xz)) %*% t(xz)
      IH <- diag(nrow(H)) - H
      
      t1 <- apply(G, 3, function(g){.calTrace(g,H)})
      t2 <- apply(G, 3, function(g){.calTrace(g,IH)})
      
      F1 <- t1/t2
    }
    
    return(F_stats)
  }

