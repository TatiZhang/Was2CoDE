deseq2_helper <- function(case_control_levels, # Control and then Case
                          case_control_var,
                          categorical_vars,
                          id_var,
                          numerical_vars,
                          seurat_obj){
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Package 'DESeq2' is required. Install it with BiocManager::install('DESeq2').")

  # aggregate to one count vector per donor; see pseudobulk.R
  pseudobulk <- build_pseudobulk(seurat_obj = seurat_obj,
                                 case_control_levels = case_control_levels,
                                 case_control_var = case_control_var,
                                 categorical_vars = categorical_vars,
                                 id_var = id_var,
                                 numerical_vars = numerical_vars)
  # do DESeq2
  dds <- .fit_deseq2_mle(mat_pseudobulk = pseudobulk$mat,
                         metadata_pseudobulk = pseudobulk$metadata,
                         use_t = FALSE,
                         quiet = FALSE)
  # DESeq2::resultsNames(dds)
  deseq2_res <- DESeq2::results(dds,
                                name = paste0(case_control_var, "_", case_control_levels[2], "_vs_", case_control_levels[1]))

  return(deseq2_res)

}
