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
  # argument validation shared with the other pseudobulk wrappers; see
  # pseudobulk.R
  .check_donor_vars(case_control_levels = case_control_levels,
                    case_control_var = case_control_var,
                    categorical_vars = categorical_vars,
                    id_var = id_var,
                    numerical_vars = numerical_vars)

  # The contrast is fixed by `case_control_levels`, never by how the column
  # happens to be stored. A character column is ordered alphabetically by the
  # model matrix ("case" before "control"), which would make the case arm the
  # reference and flip the sign of every log fold change without an error.
  case_control_vec <- as.character(seurat_obj@meta.data[[case_control_var]])
  if(case_control_levels[1] == case_control_levels[2] ||
     !all(case_control_levels %in% case_control_vec)){
    stop("`case_control_levels` must be two distinct values that both occur in `",
         case_control_var, "` (control first, then case); received: ",
         paste(case_control_levels, collapse = ", "))
  }
  seurat_obj@meta.data[[case_control_var]] <- stats::relevel(factor(case_control_vec),
                                                             ref = case_control_levels[1])

  # A pseudobulk sample is a donor, so the arm and every numerical covariate
  # must be constant within a donor. This is the same rule `build_pseudobulk()`
  # applies. Splitting a donor on these instead (which this function once did)
  # silently turned a donor in both arms into two donors, and a cell-level
  # numeric into one sample per cell.
  .check_donor_constant_vars(seurat_obj = seurat_obj,
                             id_var = id_var,
                             variables = c(case_control_var, numerical_vars))

  # A categorical covariate that varies within a donor (a donor sequenced in
  # two batches) splits that donor into one sample per level, which is what
  # `build_pseudobulk()` does by grouping on it. dreamlet takes a single
  # `sample_id` column, so the split is done by pasting the levels onto the id.
  meta_df <- seurat_obj@meta.data
  vars_split <- categorical_vars[vapply(categorical_vars, function(variable){
    any(tapply(as.character(meta_df[[variable]]), meta_df[[id_var]], function(x){
      length(unique(x[!is.na(x)])) > 1
    }), na.rm = TRUE)
  }, logical(1))]
  if(length(vars_split) > 0){
    augmented_id <- paste0(id_var, "_aug")
    seurat_obj[[augmented_id]] <- apply(meta_df[, c(id_var, vars_split), drop = FALSE],
                                        1, paste, collapse = "_")
    id_var <- augmented_id
  }

  # dreamlet's `processAssays()` only accepts the object returned by its own
  # `aggregateToPseudoBulk()` (it reads per-sample cell counts and the
  # aggregation parameters from internal slots), so the aggregation cannot go
  # through `build_pseudobulk()`.
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
