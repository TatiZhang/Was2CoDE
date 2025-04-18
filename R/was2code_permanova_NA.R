was2code_permanova_na <- function(dist_list,
                                  meta_ind,
                                  var2test,
                                  var2adjust = NULL,
                                  n_perm = 999,
                                  r_seed = 2020,
                                  residulize_x = FALSE,
                                  delta = 0.5,
                                  verbose = 0) {
  
  ## ------- assemble distance ARRAY (genes × donors × donors) ----------
  dist_array <- .form_dist_array(dist_list)          # your original helper
  
  ## ------- variable to test ------------------------------------------
  x <- meta_ind[,var2test]
  if (anyNA(x)) stop("`var2test` contains NA values.")
  if (is.character(x)) x <- factor(x)
  if (is.factor(x))   x <- as.numeric(x)
  x <- as.numeric(x == max(x))                       # force 0/1
  
  ## ------- permutations ----------------------------------------------
  set.seed(r_seed)
  x_perm <- replicate(n_perm, sample(x), simplify = TRUE)
  
  ## ------- choose the NA‑aware helper -------------------------------
  if (is.null(var2adjust)) {
    
    F_ob   <- .calc_F_manova_na(dist_array, x)
    F_perm <- apply(x_perm, 2, function(col) {
      .calc_F_manova_na(dist_array, col)
    })
    
  } else {   # ---------- covariate‑adjusted path ----------------------
    fm1 <- stats::as.formula(paste("~", paste(var2adjust, collapse=" + ")))
    z <- stats::model.matrix(fm1, data = meta_ind)
    
    if (residulize_x) {
      resid_perm <- matrix(NA_real_, nrow(meta_ind), n_perm)
      
      # step1: fit
      m1 <- stats::glm(x ~ -1 + z, family = stats::binomial(link="logit"))  # logistic model
      fitted_x <- stats::fitted(m1)
      resid_x <- x - fitted_x
      
      # step2: permutation
      ip <- 0 # index for usable permutations
      id <- 0 # working index
      sd_e <- stats::sd(resid_x) # expected sd
      
      while (ip < n_perm) {
        if(verbose == 1 && ip %% floor(n_perm/10) == 0) cat('*')
        if(verbose >= 2) print(paste0("Working on permutation: ", ip))
        
        id <- id + 1
        perm_x <- stats::rbinom(length(fitted_x), 1, prob = fitted_x)
        
        m2 <- tryCatch(stats::glm(perm_x ~ -1 + z, family = stats::binomial(link="logit")), 
                       warning = function(w) { NULL }, 
                       error = function(e) { NULL})
        
        r_i <- p - fitted(glm_p)
        if (abs(sd(r_i) - sd(r_x)) < delta * sd(r_x)) {
          ip <- ip + 1
          r_perm[, ip] <- r_i
        }
      }
      
      F_tmp <- .calc_F_permanovaSZ_na(dist_array, cbind(r_x, r_perm), Z)
      F_ob   <- F_tmp[, 1]
      F_perm <- F_tmp[, -1, drop = FALSE]
      
    } else {
      F_tmp <- .calc_F_permanovaSZ_na(dist_array, cbind(x, x_perm), Z)
      F_ob   <- F_tmp[, 1]
      F_perm <- F_tmp[, -1, drop = FALSE]
    }
  }
  
  ## ------- p‑values ---------------------------------------------------
  F_diff <- sweep(F_perm, 1, F_ob, `-`)
  pval   <- rowMeans(cbind(F_diff, 1) >= 0, na.rm = TRUE)
  names(pval) <- names(dist_list)
  pval
}


############################################
## 1.  Gower‑centering that tolerates NAs ##
############################################

.cal_G_na <- function(m) {
  # m is a (possibly incomplete) distance matrix
  d2 <- m * m
  n  <- nrow(m)
  
  row_means   <- rowMeans(d2, na.rm = TRUE)
  col_means   <- colMeans(d2, na.rm = TRUE)
  grand_mean  <- mean(d2,    na.rm = TRUE)
  
  G <- matrix(0, n, n)
  obs <- !is.na(d2)                 # mask of available pairs
  
  G[obs] <- -0.5 * (d2[obs] - 
                      row_means[row(d2)[obs]] -
                      col_means[col(d2)[obs]] +
                      grand_mean)
  G
}

###############################################
## 2.  PERMANOVA   (no covariates, with NAs) ##
###############################################

.calc_F_manova_na <- function(dist_array, label) {
  
  n_genes <- dim(dist_array)[1]
  N       <- length(label)
  groups  <- split(seq_len(N), label)
  a       <- length(groups)
  
  Fstat <- numeric(n_genes)
  
  for (g in seq_len(n_genes)) {
    d2   <- dist_array[g, , ]^2
    wMat <- !is.na(d2)               # weight / availability matrix
    
    ## ----- total SS -----
    sst  <- sum(d2, na.rm = TRUE) / N        # same scale as Anderson (2001)
    
    ## ----- within‑group SS -----
    SSW  <- 0
    for (ids in groups) {
      if (length(ids) < 2) next
      w_within <- wMat[ids, ids, drop = FALSE]
      d2_within <- d2  [ids, ids, drop = FALSE]
      # mean of available within‑group distances
      SSW <- SSW + sum(d2_within, na.rm = TRUE) /
        max(1, sum(w_within)) * length(ids)
    }
    
    ## ----- pseudo‑F -----
    Fstat[g] <- ((sst - SSW) * (N - a)) / (SSW * (a - 1))
  }
  Fstat
}

#################################################
## 3.  PERMANOVA‑SZ (with covariates, with NAs) ##
#################################################

.calTrace_mask <- function(G, H, mask) {
  # trace(H G H) but only over observed entries
  sum((H %*% G %*% H)[mask])
}

.calc_F_permanovaSZ_na <- function(dist_array, Rs, z = NULL) {
  
  n <- dim(dist_array)[2]
  n_genes <- dim(dist_array)[1]
  
  if (is.vector(Rs)) Rs <- matrix(Rs, ncol = 1)
  
  F_stats <- matrix(NA_real_, n_genes, ncol(Rs))
  
  for (g in seq_len(n_genes)) {
    D  <- dist_array[g, , ]
    G  <- .cal_G_na(D)
    mask <- !is.na(D)                # we’ll reuse this several times
    
    for (perm in seq_len(ncol(Rs))) {
      x <- Rs[, perm]
      X <- if (is.null(z)) x else cbind(x, z)
      H <- X %*% solve(crossprod(X)) %*% t(X)
      IH <- diag(n) - H
      
      t1 <- .calTrace_mask(G, H,  mask)
      t2 <- .calTrace_mask(G, IH, mask)
      
      F_stats[g, perm] <- t1 / t2
    }
  }
  F_stats
}
