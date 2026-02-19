library(Seurat)
library(SingleCellExperiment)
library(dreamlet)                         

load("~/kzlinlab/projects/microglia_negative_controls/out/kevin/Writeup3/Writeup3_gabitto-microglia_cleaned_only-counts.RData")

categorical_vars <- c("sex", "assay",  "APOE4status")
numerical_vars <- c("PMI","Ageatdeath")
case_control_var <- "ADNC"
case_control_levels <- c("Control", "Case") 
id_var <- "donor_id"
min_cells = 5
min_count = 0
min_samples = 4
min_prop = 0

SEA_AD_dreamlet <- dreamlet_helper(
  case_control_levels = case_control_levels,
  case_control_var = case_control_var,
  categorical_vars = categorical_vars,
  id_var = id_var,
  numerical_vars = numerical_vars,
  seurat_obj = seurat_obj, 
  min_cells = 5,
  min_count = 0,
  min_samples = 4,
  min_prop = 0
)

load("~/kzlinlab/projects/microglia_negative_controls/out/kevin/Writeup3/Writeup3_prater_cleaned_only-counts.RData")

prater_dreamlet <- dreamlet_helper(
  case_control_levels = case_control_levels,
  case_control_var = case_control_var,
  categorical_vars = categorical_vars,
  id_var = id_var,
  numerical_vars = numerical_vars,
  seurat_obj = seurat_obj, 
  min_cells = 5,
  min_count = 0,
  min_samples = 4,
  min_prop = 0
)

# table(seurat_obj$donor_id, seurat_obj$sex)
# table(seurat_obj$donor_id, seurat_obj$assay)
# table(seurat_obj$donor_id, seurat_obj$APOE4status)
# table(seurat_obj$ADNC, seurat_obj$APOE4status)
# We will need some code to automatically detect which variables are needed to "augmnet" the id_var,
#.  so that all the variables in the sce ColData() are constant for the augmented variable

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
  sce <- as.SingleCellExperiment(seurat_obj)
  
  pb <- aggregateToPseudoBulk(sce,
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
  form <- formula(form)
  
  res_proc <- processAssays(pb,
                            form,
                            min.cells = min_cells,
                            min.count = min_count,
                            min.samples = min_samples,
                            min.prop = min_prop)
  
  res_dl <- dreamlet(res_proc, form)
  
  res_pvalues <- topTable(res_dl,
                          coef = paste0(case_control_var, case_control_levels[2]),
                          number = Inf) # see if there's a better way to write this
  
  return(res_pvalues)
}

















