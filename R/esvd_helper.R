esvd_helper <- function(batch_var_prefix, # a variable inside categorical_vars. Can be NULL
                        case_control_levels, # Control and then Case
                        case_control_var,
                        categorical_vars,
                        id_var,
                        numerical_vars,
                        seurat_obj,
                        bool_check_donors = TRUE,
                        intermediate_save = NULL, # NULL or filepath to save intermediary results
                        min_cells_casecontrol = 20,
                        min_cells_per_id = 3,
                        min_cells = 20,
                        min_ids = 4,
                        verbose = 0,
                        ...){
  if (!requireNamespace("eSVD2", quietly = TRUE))
    stop("Package 'eSVD2' is required. Install it with remotes::install_github('linnykos/eSVD2').")
  
  seurat_obj@meta.data[,case_control_var] <- droplevels(seurat_obj@meta.data[,case_control_var])
  
  if(bool_check_donors){
    # remove any subjects that only have less than min_cells_per_id cells
    indiv_count <- table(seurat_obj@meta.data[,id_var])
    if(any(indiv_count < min_cells_per_id)){
      indiv_rm <- names(indiv_count)[indiv_count < min_cells_per_id]
      
      keep_vec <- !seurat_obj@meta.data[,id_var] %in% indiv_rm
      seurat_obj <- seurat_obj[, keep_vec]
      seurat_obj@meta.data[,id_var] <- droplevels(seurat_obj@meta.data[,id_var])
    }
    # if there are min_cells or less cells, go next
    if(length(Seurat::Cells(seurat_obj)) <= min_cells){
      warning("Not enough cells, returning NA")
      return(NA)
    }
    # if there are min_ids or less donors, go next
    if(length(unique(seurat_obj@meta.data[,id_var])) <= min_ids){
      warning("Not enough donors, returning NA")
      return(NA)
    }
    # if there are min_cells_casecontrol or less cells among either all the cases or controls, go next
    if(any(table(seurat_obj@meta.data[,case_control_var]) <= min_cells_casecontrol)){
      warning("Not enough cell in case or control, returning NA")
      return(NA)
    }
  }
  
  esvd_fn <- utils::getFromNamespace("eSVD", "eSVD2")
  eSVD_obj <- esvd_fn(batch_var_prefix = batch_var_prefix,
                      case_control_levels = case_control_levels,
                      case_control_var = case_control_var,
                      categorical_vars = categorical_vars,
                      id_var = id_var,
                      numerical_vars = numerical_vars,
                      seurat_obj = seurat_obj,
                      intermediate_save = intermediate_save,
                      verbose = verbose,
                      ...)
  
  eSVD_obj
}