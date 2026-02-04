manova <- function(expr_assay, # "SCT"
                   expr_layer, # "data"
                   factor_vars, # has to be strictly factors
                   id_vars, # could be more than one variable
                   seurat_obj,
                   return_sorted = TRUE,
                   verbose = 0){
  
  if(verbose >= 2) print("Extracting epxression matrix")
  expr <- SeuratObject::LayerData(seurat_obj, 
                                  assay = expr_assay, 
                                  layer = expr_layer)
  
  if(length(id_vars) == 1){
    group_key <- factor(seurat_obj@meta.data[,id_vars])
  } else {
    md <- seurat_obj@meta.data[,id_vars]
    group_key <- factor(apply(md, 1, paste, collapse = "||"))
  }
 
  
  mm <- Matrix::sparse.model.matrix(~ 0 + group_key)  # cells x groups
  colnames(mm) <- sub("^group_key", "", colnames(mm))
  n_per_group <- Matrix::colSums(mm)  # group sizes (length = #groups)
  
  # ----------------------------
  # Compute means via sums / n
  # ----------------------------
  if(verbose >= 2) print("Computing donor mean")
  sums <- expr %*% mm                                # genes x groups
  mean_mat <- sweep(as.matrix(sums), 2, n_per_group, "/")  # genes x groups (dense matrix)
  
  # ----------------------------
  # Compute variances (sample variance by default)
  #    Var = (sum(x^2) - n*mean^2)/(n-1)
  # ----------------------------
  if(verbose >= 2) print("Computing donor variance")
  expr_sq <- expr
  expr_sq@x <- expr_sq@x^2                           # sparse-friendly square
  sumsq <- expr_sq %*% mm                            # genes x groups
  sumsq_mat <- as.matrix(sumsq)
  
  var_mat <- sweep(
    sumsq_mat - sweep(mean_mat^2, 2, n_per_group, "*"),
    2,
    (n_per_group - 1),
    "/"
  )
  
  # tiny negative numerical noise -> clamp to 0
  var_mat[var_mat <= 0 & var_mat > -1e-12] <- 0
  var_mat[is.na(var_mat)] <- 0
  
  ########################
  if(verbose >= 2) print("Staring R2 calculation")
  df <- seurat_obj@meta.data[,unique(c(id_vars, factor_vars))]
  df <- unique(df)
  if(nrow(df) != length(levels(group_key))) stop("There is a mismatch between the factors in `id_vars` and the factors that make each cell's covarites unique in `factor_vars`. Likely, you need to include more variables in `id_vars`")
  
  gene_R2 <- matrix(NA, nrow = nrow(mean_mat), ncol = length(factor_vars)+1)
  rownames(gene_R2) <- rownames(mean_mat)
  colnames(gene_R2) <- c(factor_vars, "SST")
  
  for (gene_idx in seq_len(nrow(mean_mat))) {
    if (verbose > 0 && nrow(mean_mat) > 10 && gene_idx %% floor(nrow(mean_mat)/10) == 0) cat("*")
    
    # ----- Celltype R2 (unchanged; computed on ALL lineage-celltype combos) -----
    D_all <- .W2_dist_mat(mean_mat[gene_idx, ], var_mat[gene_idx, ])
    
    res_list <- lapply(1:length(factor_vars), function(j){
      .manova_var_explained(D = D_all, 
                            group = droplevels(df[,factor_vars[j]]))
    })
    
    for(j in 1:length(factor_vars)){
      gene_R2[gene_idx, factor_vars[j]] <- res_list[[j]]$R2
    } 
    gene_R2[gene_idx, "SST"] <-  res_list[[1]]$SST
  }
  
  gene_R2 <- as.data.frame(gene_R2)
  
  if(return_sorted){
    sum_vec <- rowSums(gene_R2[,setdiff(colnames(gene_R2), "SST")])
    gene_R2 <- gene_R2[order(sum_vec, decreasing = TRUE),]
  }
  
  gene_R2
}

###########

.manova_var_explained <- function(D, group) {
  # D: distance matrix (n x n) or dist object, containing UNSQUARED distances
  # group: length-n grouping vector/factor
  
  n <- nrow(D)
  stopifnot(n == ncol(D))
  
  # basic checks (optional but helpful)
  if (anyNA(D)) stop("D has NA values.")
  if (max(abs(D - t(D))) > 1e-8) stop("D must be symmetric.")
  if (max(abs(diag(D))) > 1e-8) stop("diag(D) should be ~0 for a distance matrix.")
  
  # convert one-hot to group labels if needed
  if (length(group) != n) stop("group must have length n.")
  group <- factor(group)
  
  # PERMANOVA uses squared distances
  D2 <- D^2
  
  # total SS
  SST <- sum(D2) / (2 * n)
  
  # within-group SS
  idx_list <- split(seq_len(n), group)
  SSW <- sum(vapply(idx_list, function(ii) {
    nk <- length(ii)
    sum(D2[ii, ii, drop = FALSE]) / (2 * nk)
  }, numeric(1)))
  
  # explained (between / model) SS
  SSA <- SST - SSW
  
  # percent explained
  R2 <- SSA / SST
  
  list(
    SST = SST,
    SSW = SSW,
    SSA = SSA,
    R2 = R2
  )
}

# helper: fast W2 distance matrix for 1D Gaussians
.W2_dist_mat <- function(mu, var) {
  s <- sqrt(var)
  sqrt(outer(mu, mu, "-")^2 + outer(s, s, "-")^2)
}