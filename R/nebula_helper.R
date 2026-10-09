nebula_helper <- function(case_control_levels, # Control and then Case
                          case_control_var,
                          categorical_vars,
                          id_var,
                          numerical_vars,
                          seurat_obj,
                          verbose = 0){
  if (!requireNamespace("nebula", quietly = TRUE))
    stop("Package 'nebula' is required. Install it with install.packages('nebula').")
  # argument validation shared with the other DE wrappers; see pseudobulk.R
  .check_donor_vars(case_control_levels = case_control_levels,
                    case_control_var = case_control_var,
                    categorical_vars = categorical_vars,
                    id_var = id_var,
                    numerical_vars = numerical_vars)

  # The contrast is fixed by `case_control_levels`, never by how the column
  # happens to be stored. `model.matrix()` takes the reference level from the
  # column as stored, so a character column ("case" before "control") would
  # name the coefficient `<var>control` and flip the sign of every log fold
  # change, and the case-named column a caller reads would not exist.
  case_control_vec <- as.character(seurat_obj@meta.data[[case_control_var]])
  if(case_control_levels[1] == case_control_levels[2] ||
     !all(case_control_levels %in% case_control_vec)){
    stop("`case_control_levels` must be two distinct values that both occur in `",
         case_control_var, "` (control first, then case); received: ",
         paste(case_control_levels, collapse = ", "))
  }
  seurat_obj@meta.data[[case_control_var]] <- stats::relevel(factor(case_control_vec),
                                                             ref = case_control_levels[1])

  # NEBULA is a cell-level model with a donor random effect, so a covariate
  # may vary within a donor (a per-cell QC score is a legitimate adjustment).
  # The arm may not: a donor in both arms is a metadata error, and the random
  # effect would absorb it without complaint.
  .check_donor_constant_vars(seurat_obj = seurat_obj,
                             id_var = id_var,
                             variables = case_control_var)

  # A categorical covariate with one level across the cohort (a single-sex
  # cohort; a batch constant within one cell type) cannot be given contrasts,
  # so drop it from the design, as `build_pseudobulk()` does.
  for(variable in categorical_vars){
    seurat_obj@meta.data[[variable]] <- factor(as.character(seurat_obj@meta.data[[variable]]))
  }
  constant_vars <- categorical_vars[vapply(categorical_vars, function(variable){
    nlevels(seurat_obj@meta.data[[variable]]) < 2
  }, logical(1))]
  if(length(constant_vars) > 0){
    if(verbose > 0){
      print(paste0("Dropping constant categorical covariate(s): ",
                   paste(constant_vars, collapse = ", ")))
    }
    categorical_vars <- setdiff(categorical_vars, constant_vars)
  }
  design_vars <- c(case_control_var, categorical_vars, numerical_vars)

  # `scToNeb()` gained a `verbose` argument only in recent nebula versions,
  # and it guards just the "no assay provided" message, so it is not passed.
  neb_data <- nebula::scToNeb(obj = seurat_obj,
                              assay = "RNA",
                              id = id_var,
                              pred = design_vars,
                              offset = "nCount_RNA")

  # nebula requires the cells of each donor to be contiguous. `drop = FALSE`
  # keeps a one-column `pred` (no covariates) a data.frame.
  order_index <- order(neb_data$id)
  neb_data$count <- neb_data$count[, order_index, drop = FALSE]
  neb_data$id <- neb_data$id[order_index]
  neb_data$pred <- neb_data$pred[order_index, , drop = FALSE]
  neb_data$offset <- neb_data$offset[order_index]

  design_mat <- stats::model.matrix(stats::reformulate(design_vars),
                                    data = neb_data$pred)

  # `cpc = 0` and `mincp = 0` switch off nebula's expression filters so every
  # gene is tested, matching the other wrappers.
  nebula_res <- nebula::nebula(count = neb_data$count,
                               id = neb_data$id,
                               pred = design_mat,
                               offset = neb_data$offset,
                               model = "NBGMM",
                               verbose = verbose > 0,
                               cpc = 0,
                               mincp = 0)

  nebula_res
}
