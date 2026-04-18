dreamlet_helper <- function(case_control_levels, # Control and then Case
                            case_control_var,
                            categorical_vars,
                            id_var,
                            numerical_vars,
                            seurat_obj,
                            min_cells = 5,
                            min_count = 0,
                            min_samples = 4,
                            min_prop = 0){
  
  # check that all the variables in c(case_control_var, categorical_vars, numerical_vars) 
  # are unique within each id_var
  # if not, append that variable to id_var
  vars_check <- c(case_control_var, categorical_vars, numerical_vars)
  meta <- seurat_obj@meta.data
  
  vars_diff <- function(v){
    x <- meta[[v]]
    any(tapply(x, meta[[id_var]], function(vv){
      vv <- vv[!is.na(vv)]
      if(is.numeric(vv)) diff(range(vv)) > 1e-8
      else length(unique(as.character(vv))) > 1
    }), na.rm = TRUE)
  }
  
  vars_append <- vars_check[vapply(vars_check, vars_diff, logical(1))]
  
  if(length(vars_append) > 0){
    augmented_id <- paste0(id_var, "_aug")
    pieces <- c(id_var, vars_append)
    seurat_obj[[augmented_id]] <- apply(meta[, pieces, drop = FALSE], 1, paste, collapse = "_")
    id_var <- augmented_id
  }
  
  seurat_obj$tmp_variable <- rep("tmp", length(Seurat::Cells(seurat_obj)))
  sce <- Seurat::as.SingleCellExperiment(seurat_obj)
  
  pb <- dreamlet::aggregateToPseudoBulk(sce,
                                        assay = "counts",
                                        cluster_id = "tmp_variable",
                                        sample_id = id_var,
                                        verbose = FALSE)
  
  form <- paste0("~ ", case_control_var)
  if(length(categorical_vars) > 0){
    form <- paste0(form, "+", paste(categorical_vars, collapse = "+"))
  }
  if(length(numerical_vars) > 0){
    form <- paste0(form, "+", paste(numerical_vars, collapse = "+"))
  }
  form <- stats::formula(form)
  
  res_proc <- dreamlet::processAssays(pb,
                                      form,
                                      min.cells = min_cells,
                                      min.count = min_count,
                                      min.samples = min_samples,
                                      min.prop = min_prop)
  
  res_dl <- dreamlet::dreamlet(res_proc, form)
  
  res_pvalues <- dreamlet::topTable(res_dl,
                                    coef = paste0(coefNames(res_dl)[2]),
                                    number = Inf) 
  
  return(res_pvalues)
}















