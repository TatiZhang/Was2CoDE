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
  if (!requireNamespace("dreamlet", quietly = TRUE))
    stop("Package 'dreamlet' is required. Install it with BiocManager::install('dreamlet').")
  stopifnot(length(case_control_var) == 1,
            is.character(case_control_var),
            is.character(case_control_levels))

  # The contrast is fixed by `case_control_levels`, never by how the column
  # happens to be stored. A character column is ordered alphabetically by the
  # model matrix ("case" before "control"), which would make the case arm the
  # reference and flip the sign of every log fold change without an error.
  case_control_vec <- as.character(seurat_obj@meta.data[[case_control_var]])
  if(length(case_control_levels) != 2 ||
     case_control_levels[1] == case_control_levels[2] ||
     !all(case_control_levels %in% case_control_vec)){
    stop("`case_control_levels` must be two distinct values that both occur in `",
         case_control_var, "` (control first, then case); received: ",
         paste(case_control_levels, collapse = ", "))
  }
  seurat_obj@meta.data[[case_control_var]] <- stats::relevel(factor(case_control_vec),
                                                             ref = case_control_levels[1])

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
  
  # select the case-vs-control coefficient by name, not by its position among
  # the coefficients
  coef_name <- paste0(case_control_var, case_control_levels[2])
  if(!(coef_name %in% dreamlet::coefNames(res_dl))){
    stop("The coefficient `", coef_name, "` is not in the fitted model, whose ",
         "coefficients are: ", paste(dreamlet::coefNames(res_dl), collapse = ", "))
  }
  res_pvalues <- dreamlet::topTable(res_dl,
                                    coef = coef_name,
                                    number = Inf)
  
  return(res_pvalues)
}















