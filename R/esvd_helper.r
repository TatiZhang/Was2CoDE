esvd_helper <- function(batch_var_prefix, # a variable inside categorical_vars. Can be NULL
                        case_control_levels, # Control and then Case
                        case_control_var,
                        categorical_vars,
                        id_var,
                        numerical_vars,
                        seurat_obj,
                        verbose = 0){
  stopifnot(all(is.null(categorical_vars)) || (length(unique(categorical_vars)) == length(categorical_vars) && all(is.character(categorical_vars))))
  stopifnot(all(is.null(numerical_vars)) || (length(unique(numerical_vars)) == length(numerical_vars) && all(is.character(numerical_vars))))
  stopifnot(length(case_control_var) == 1,
            is.character(case_control_var),
            length(id_var) == 1,
            is.character(id_var),
            length(case_control_levels) == 2,
            all(is.character(case_control_levels)))
  
  # make sure there's an appropriate batch variable
  stopifnot(is.null(batch_var_prefix) || length(grep(batch_var_prefix, categorical_vars) > 1))
  
  # extract count matrix
  mat <- Matrix::t(SeuratObject::LayerData(seurat_obj, 
                                           assay = "RNA", 
                                           layer = "counts"))
  
  if(verbose > 0) print("Processing the covariates")
  if(length(categorical_vars) >= 1){
    categorical_vars_subset <- categorical_vars[sapply(categorical_vars, function(x){
      length(levels(droplevels(seurat_obj@meta.data[,x]))) > 1
    })]
  } else {
    categorical_vars_subset <- NULL
  }
  covariate_dat <- seurat_obj@meta.data[,c(id_var, categorical_vars_subset, numerical_vars), drop = FALSE]
  covariate_df <- data.frame(covariate_dat)
  
  covariate_df[,case_control_var] <- factor(seurat_obj@meta.data[,case_control_var], 
                                            levels = case_control_levels)
  for(variable in setdiff(c(id_var, categorical_vars_subset), case_control_var)){
    covariate_df[,variable] <- factor(covariate_df[,variable], 
                                      levels = names(sort(table(covariate_df[,variable]), 
                                                          decreasing = TRUE)))
  }
  
  covariates <- eSVD2::format_covariates(dat = mat,
                                         covariate_df = covariate_df,
                                         rescale_numeric_variables = numerical_vars)
  
  case_control_variable <- paste0(case_control_var, "_", case_control_levels[2])
  
  ############
  
  if(verbose > 0) print("Initialization")
  eSVD_obj <- eSVD2::initialize_esvd(dat = mat,
                                     covariates = covariates[,-grep(id_var, colnames(covariates))],
                                     case_control_variable = case_control_variable,
                                     bool_intercept = TRUE,
                                     k = 30,
                                     lambda = 0.1,
                                     metadata_case_control = covariates[,case_control_variable],
                                     metadata_individual = covariate_df[,id_var],
                                     verbose = 1)
  
  if(!is.null(batch_var_prefix)){
    omitted_variables <- colnames(eSVD_obj$covariates)[grep(batch_var_prefix, colnames(eSVD_obj$covariates))]
  } else {
    omitted_variables <- NULL
  }
  
  eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Init",
    omitted_variables = c("Log_UMI", omitted_variables)
  )
  
  if(verbose > 0)  print("First fit")
  eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                              l2pen = 0.1,
                              max_iter = 100,
                              offset_variables = setdiff(colnames(eSVD_obj$covariates), case_control_variable),
                              tol = 1e-6,
                              verbose = 1,
                              fit_name = "fit_First",
                              fit_previous = "fit_Init")
  
  eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_First",
    omitted_variables = c("Log_UMI", omitted_variables)
  )
  
  if(verbose > 0) print("Second fit")
  eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                              l2pen = 0.1,
                              max_iter = 100,
                              offset_variables = NULL,
                              tol = 1e-6,
                              verbose = 1,
                              fit_name = "fit_Second",
                              fit_previous = "fit_First")
  
  eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Second",
    omitted_variables = omitted_variables
  )
  
  if(verbose > 0) print("Nuisance estimation")
  eSVD_obj <- eSVD2::estimate_nuisance(input_obj = eSVD_obj,
                                       bool_covariates_as_library = TRUE,
                                       bool_library_includes_interept = TRUE,
                                       bool_use_log = FALSE,
                                       verbose = 1)
  
  if(verbose > 0) print("Computing posterior")
  eSVD_obj <- eSVD2::compute_posterior(input_obj = eSVD_obj,
                                       bool_adjust_covariates = FALSE,
                                       alpha_max = 2*max(mat@x),
                                       bool_covariates_as_library = TRUE,
                                       bool_stabilize_underdispersion = TRUE,
                                       library_min = 0.1,
                                       pseudocount = 0)
  
  if(verbose > 0) print("Computing p-values")
  eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj,
                                            verbose = 1)
  eSVD_obj <- eSVD2::compute_pvalue(input_obj = eSVD_obj)
  
  if(verbose > 0) print("Finalizing")
  eSVD_obj$covariates <- NULL
  eSVD_obj$fit_Init <- NULL
  eSVD_obj$fit_First <- NULL
  eSVD_obj$fit_Second <- NULL
  eSVD_obj$dat <- NULL
  
  eSVD_obj
}