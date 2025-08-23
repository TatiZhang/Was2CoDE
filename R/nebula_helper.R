nebula_helper <- function(case_control_levels,
                          case_control_var,
                          categorical_vars,
                          id_var,
                          numerical_vars,
                          seurat_obj){
  neb_data <- nebula::scToNeb(obj = seurat_obj,
                              assay = "RNA",
                              id = id_var,
                              pred = c(case_control_var, categorical_vars, numerical_vars),
                              offset = "nCount_RNA")
  
  order_index <- order(neb_data$id)
  # Reorder the count matrix by the id
  neb_data$count <- neb_data$count[, order_index]
  
  # Reorder the other components of neb_data
  neb_data$id <- neb_data$id[order_index]
  neb_data$pred <- neb_data$pred[order_index, ]
  neb_data$offset <- neb_data$offset[order_index]
  
  df <- eval(parse(text = paste0("stats::model.matrix( ~",
                                 paste0(c(case_control_var, categorical_vars, numerical_vars), collapse = "+"),
                                 ", data = neb_data$pred)")))

  nebula_res <- nebula::nebula(count = neb_data$count,
                               id = neb_data$id,
                               pred = df,
                               offset = neb_data$offset,
                               model = "NBGMM",
                               verbose = TRUE,
                               cpc = 0,
                               mincp = 0)
  
  nebula_res
}
