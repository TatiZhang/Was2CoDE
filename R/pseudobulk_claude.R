# Donor-level pseudobulk construction, shared by the DE-method wrappers.
#
# Every wrapper in this package that runs a bulk method on single-cell data
# needs the same four steps: aggregate nuclei to one count vector per donor,
# carry the donor-level numerical covariates across the aggregation, coerce the
# design variables into the factor/scale form the model expects, and drop
# design variables that do not vary.  `build_pseudobulk()` is that pipeline;
# `deseq2_helper()` and `deseq2_tost_helper()` both call it.

#' Validate the donor / covariate arguments shared by the DE-method wrappers
#'
#' @param case_control_levels Length-2 character vector, control level first.
#' @param case_control_var Name of the case/control metadata column.
#' @param categorical_vars Character vector of categorical covariates, or `NULL`.
#' @param id_var Name of the donor id metadata column.
#' @param numerical_vars Character vector of numerical covariates, or `NULL`.
#'
#' @return Invisibly `TRUE`; called for its side effect of erroring.
#' @noRd
.check_donor_vars <- function(case_control_levels,
                              case_control_var,
                              categorical_vars,
                              id_var,
                              numerical_vars) {
  stopifnot(all(is.null(categorical_vars)) || (length(unique(categorical_vars)) == length(categorical_vars) && all(is.character(categorical_vars))))
  stopifnot(all(is.null(numerical_vars)) || (length(unique(numerical_vars)) == length(numerical_vars) && all(is.character(numerical_vars))))
  stopifnot(length(case_control_var) == 1,
            is.character(case_control_var),
            length(id_var) == 1,
            is.character(id_var),
            length(case_control_levels) == 2,
            all(is.character(case_control_levels)))

  invisible(TRUE)
}

#' Aggregate a Seurat object to one pseudobulk count vector per donor
#'
#' Summing nuclei per donor before testing is what makes the downstream test
#' donor-level rather than cell-level, and is what avoids pseudoreplication.
#'
#' @param seurat_obj A `Seurat` object with cells from many donors.
#' @param case_control_var,categorical_vars,id_var Metadata column names used
#'   as grouping variables.
#'
#' @return A `Seurat` object with one "cell" per donor x covariate combination,
#'   carrying the original object's variable features.
#' @noRd
.pseudobulk_from_seurat <- function(seurat_obj,
                                    case_control_var,
                                    categorical_vars,
                                    id_var) {
  pseudo_seurat <- Seurat::AggregateExpression(seurat_obj,
                                               assays = "RNA",
                                               return.seurat = TRUE,
                                               group.by = c(id_var, case_control_var, categorical_vars))
  Seurat::VariableFeatures(pseudo_seurat) <- Seurat::VariableFeatures(seurat_obj)

  pseudo_seurat
}

#' Copy donor-constant numerical covariates onto the pseudobulk object
#'
#' Aggregation drops numerical metadata, so it has to be re-attached by hand.
#' Each numerical variable must be constant within a donor -- otherwise it is
#' not a donor-level covariate and cannot survive the aggregation.
#'
#' @param pseudo_seurat Pseudobulk `Seurat` object.
#' @param seurat_obj The original cell-level `Seurat` object.
#' @param numerical_vars Character vector of metadata column names, or `NULL`.
#' @param id_var Donor id column name.
#'
#' @return `pseudo_seurat` with the numerical variables added to `meta.data`.
#' @noRd
.attach_donor_level_numerics <- function(pseudo_seurat,
                                         seurat_obj,
                                         numerical_vars,
                                         id_var) {
  for (variable in numerical_vars) {
    tmp <- rep(NA, length(Seurat::Cells(pseudo_seurat)))
    names(tmp) <- pseudo_seurat@meta.data[, id_var]

    for (person in unique(names(tmp))) {
      idx <- which(seurat_obj@meta.data[, id_var] == person)
      person_idx <- which(names(tmp) == person)
      zz <- seurat_obj@meta.data[idx, variable]
      if (any(abs(diff(range(zz))) >= 1e-4)) {
        stop(paste0("Person (", person, ") violates the variable (", variable, "): ", diff(range(zz))))
      }
      tmp[person_idx] <- mean(zz)
    }

    pseudo_seurat@meta.data[, variable] <- tmp
  }

  pseudo_seurat
}

#' Build the model's colData from the pseudobulk metadata
#'
#' Relevels the case/control factor so the control level is the reference,
#' orders the other factors by decreasing frequency, standardizes the numerical
#' covariates, and drops any factor with no variation across donors (which
#' would otherwise make the design matrix rank deficient).  Note `id_var` is
#' releveled but deliberately not kept: donor is the unit of observation, not a
#' design variable.
#'
#' @param metadata_pseudobulk `data.frame` of pseudobulk metadata.
#' @param case_control_levels Length-2 character vector, control first.
#' @param case_control_var,categorical_vars,id_var,numerical_vars Column names.
#'
#' @return A `data.frame` holding only the design variables that vary.
#' @noRd
.prepare_pseudobulk_metadata <- function(metadata_pseudobulk,
                                         case_control_levels,
                                         case_control_var,
                                         categorical_vars,
                                         id_var,
                                         numerical_vars) {
  metadata_pseudobulk[, case_control_var] <- stats::relevel(factor(metadata_pseudobulk[, case_control_var]),
                                                            ref = case_control_levels[1])

  for (variable in c(categorical_vars, id_var)) {
    tab_vec <- table(metadata_pseudobulk[, variable])
    metadata_pseudobulk[, variable] <- factor(metadata_pseudobulk[, variable],
                                              levels = names(tab_vec)[order(tab_vec, decreasing = TRUE)])
  }

  for (variable in numerical_vars) {
    metadata_pseudobulk[, variable] <- scale(as.numeric(metadata_pseudobulk[, variable]))
  }

  # make sure there's variation among the donors
  metadata_pseudobulk <- metadata_pseudobulk[, c(categorical_vars, numerical_vars, case_control_var), drop = FALSE]
  keep_vars <- c()
  for (j in 1:ncol(metadata_pseudobulk)) {
    if (!is.factor(metadata_pseudobulk[, j])) {
      keep_vars <- c(keep_vars, j)
    } else if (length(unique(metadata_pseudobulk[, j])) > 1) {
      keep_vars <- c(keep_vars, j)
    }
  }

  metadata_pseudobulk[, keep_vars, drop = FALSE]
}

#' Aggregate a Seurat object into a donor-level pseudobulk count matrix and design
#'
#' The shared front half of every pseudobulk DE wrapper in this package: it
#' turns a cell-level `Seurat` object into the `(countData, colData)` pair a
#' bulk model needs, with covariates carried across the aggregation.
#'
#' Only the variable features of `seurat_obj` are returned, matching the
#' behaviour the DE-method wrappers have always had.
#'
#' @param seurat_obj A `Seurat` object.
#' @param case_control_levels Length-2 character vector: control level first,
#'   then case.
#' @param case_control_var Name of the case/control metadata column.
#' @param categorical_vars Character vector of categorical covariates, or `NULL`.
#' @param id_var Name of the donor id metadata column.
#' @param numerical_vars Character vector of donor-constant numerical
#'   covariates, or `NULL`.
#'
#' @return A list with
#'   \describe{
#'     \item{`mat`}{genes x pseudobulk-sample count matrix}
#'     \item{`metadata`}{`data.frame` of design variables, rows aligned to the
#'       columns of `mat`}
#'     \item{`pseudo_seurat`}{the intermediate aggregated `Seurat` object, kept
#'       so callers can reach metadata that was dropped from the design}
#'   }
#' @noRd
build_pseudobulk <- function(seurat_obj,
                             case_control_levels,
                             case_control_var,
                             categorical_vars,
                             id_var,
                             numerical_vars) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required. Install it with install.packages('Seurat').")
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("Package 'SeuratObject' is required. Install it with install.packages('SeuratObject').")
  .check_donor_vars(case_control_levels = case_control_levels,
                    case_control_var = case_control_var,
                    categorical_vars = categorical_vars,
                    id_var = id_var,
                    numerical_vars = numerical_vars)

  pseudo_seurat <- .pseudobulk_from_seurat(seurat_obj = seurat_obj,
                                           case_control_var = case_control_var,
                                           categorical_vars = categorical_vars,
                                           id_var = id_var)
  pseudo_seurat <- .attach_donor_level_numerics(pseudo_seurat = pseudo_seurat,
                                                seurat_obj = seurat_obj,
                                                numerical_vars = numerical_vars,
                                                id_var = id_var)

  mat_pseudobulk <- SeuratObject::LayerData(pseudo_seurat,
                                            layer = "counts",
                                            assay = "RNA",
                                            features = Seurat::VariableFeatures(pseudo_seurat))
  metadata_pseudobulk <- .prepare_pseudobulk_metadata(metadata_pseudobulk = pseudo_seurat@meta.data,
                                                      case_control_levels = case_control_levels,
                                                      case_control_var = case_control_var,
                                                      categorical_vars = categorical_vars,
                                                      id_var = id_var,
                                                      numerical_vars = numerical_vars)

  # Bulk model fitters match countData columns to colData rows by position, so
  # a mismatch here would silently misassign donors rather than error.
  stopifnot(ncol(mat_pseudobulk) == nrow(metadata_pseudobulk),
            all(colnames(mat_pseudobulk) == rownames(metadata_pseudobulk)))

  list(mat = mat_pseudobulk,
       metadata = metadata_pseudobulk,
       pseudo_seurat = pseudo_seurat)
}
