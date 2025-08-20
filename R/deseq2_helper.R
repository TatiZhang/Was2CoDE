deseq2_helper <- function(case_control_levels,
                          case_control_var,
                          categorical_vars,
                          id_var,
                          numerical_vars,
                          seurat_obj){
  pseudo_seurat <- Seurat::AggregateExpression(seurat_obj, 
                                               assays = "RNA", 
                                               return.seurat = TRUE, 
                                               group.by = c(id_var, case_control_var, categorical_vars))
  Seurat::VariableFeatures(pseudo_seurat) <- Seurat::VariableFeatures(seurat_obj)
  
  # add the numerical variables in manually
  for(variable in numerical_vars){
    tmp <- rep(NA, length(Seurat::Cells(pseudo_seurat)))
    names(tmp) <- pseudo_seurat@meta.data[,id_var]
    
    for(person in unique(names(tmp))){
      idx <- which(seurat_obj@meta.data[,id_var] == person)
      person_idx <- which(names(tmp) == person)
      zz <- seurat_obj@meta.data[idx, variable]
      stopifnot(abs(diff(range(zz))) <= 1e-4)
      tmp[person_idx] <- mean(zz)
    }
    
    pseudo_seurat@meta.data[,variable] <- tmp
  }
  
  mat_pseudobulk <- SeuratObject::LayerData(pseudo_seurat,
                                            layer = "counts",
                                            assay = "RNA",
                                            features = Seurat::VariableFeatures(pseudo_seurat))
  metadata_pseudobulk <- pseudo_seurat@meta.data
  
  metadata_pseudobulk[,case_control_var] <- relevel(factor(metadata_pseudobulk[,case_control_var]), 
                                                    ref = case_control_levels[1])
  
  for(variable in c(categorical_vars, id_var)){
    tab_vec <- table(metadata_pseudobulk[,variable])
    metadata_pseudobulk[,variable] <- factor(metadata_pseudobulk[,variable], 
                                             levels = names(tab_vec)[order(tab_vec, decreasing = TRUE)])
  }
  
  for(variable in numerical_vars){
    metadata_pseudobulk[,variable] <- scale(as.numeric(metadata_pseudobulk[,variable]))
  }
  
  # make sure there's variation among the donors
  metadata_pseudobulk <- metadata_pseudobulk[,c(categorical_vars, numerical_vars, case_control_var)]
  keep_vars <- c()
  for(j in 1:ncol(metadata_pseudobulk)){
    if(!is.factor(metadata_pseudobulk[,j])) {
      keep_vars <- c(keep_vars, j)
    } else if(length(unique(metadata_pseudobulk[,j])) > 1) {
      keep_vars <- c(keep_vars, j)
    }
  }
  metadata_pseudobulk <- metadata_pseudobulk[,keep_vars]
  
  # do DESeq2
  dds <- eval(parse(text = paste0("DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk, colData = metadata_pseudobulk, design = ~ ",
                                  paste0(colnames(metadata_pseudobulk), collapse = "+"),
                                  ")")))
  dds <- DESeq2::DESeq(dds)
  # DESeq2::resultsNames(dds)
  deseq2_res <- DESeq2::results(dds, 
                                name = paste0(case_control_var, "_", case_control_levels[2], "_vs_", case_control_levels[1]))
  
  return(deseq2_res)
  
}