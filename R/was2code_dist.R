#' Compute Wasserstein-2 Distance Between Individuals
#'
#' @description
#' Computes pairwise Wasserstein-2 distances across donors for each gene.
#' Expression values are assumed to be already denoised and normalized; no
#' transformation or residualization is applied internally.
#'
#' @param count_input A numeric matrix (genes x cells) of preprocessed
#'   expression values. Must have row names (gene IDs) and column names
#'   (cell IDs).
#' @param meta_cell A data.frame with one row per cell. Must contain columns
#'   \code{cell_id} and \code{individual}. Rows must be in the same order as
#'   the columns of \code{count_input}.
#' @param meta_ind A data.frame with one row per donor. Must contain columns
#'   \code{individual} and the column named by \code{var2test}.
#' @param var2test Name of the grouping variable in \code{meta_ind} (must be
#'   a factor or coercible to one, with at least two levels).
#' @param ncores Number of cores for parallelization (default \code{2}).
#' @param k Maximum number of same-group and opposite-group donors to compare
#'   each donor against. \code{NULL} compares against all donors (default
#'   \code{NULL}).
#' @param verbose Integer verbosity level: 0 = silent, 1 = log donor values,
#'   2+ = step-level messages (default \code{0}).
#' @return A named list of length equal to the number of genes. Each element
#'   is a 3D numeric array (donors x donors x 4) with named slices:
#'   \code{was2}, \code{location}, \code{size}, and \code{shape}.
#' @export

# modified from https://github.com/Sun-lab/ideas/blob/main/R/ideas_dist.R
was2code_dist <- function(count_input,
                          meta_cell,
                          meta_ind,
                          var2test,
                          ncores  = 2,
                          k       = NULL,
                          verbose = 0) {
  set.seed(0)

  # ── Coerce data.tables to data.frames ──────────────────────────────────────
  if (!is.data.frame(meta_cell)) stop("meta_cell must be a data.frame")
  if (!is.data.frame(meta_ind))  stop("meta_ind must be a data.frame")
  if (inherits(meta_cell, "data.table")) meta_cell <- as.data.frame(meta_cell)
  if (inherits(meta_ind,  "data.table")) meta_ind  <- as.data.frame(meta_ind)

  # ── Validate count_input ───────────────────────────────────────────────────
  if (!is.matrix(count_input))
    stop("count_input must be a matrix")

  n_gene   <- nrow(count_input)
  n_cell   <- ncol(count_input)
  gene_ids <- rownames(count_input)
  cell_ids <- colnames(count_input)

  if (is.null(gene_ids))                  stop("count_input must have row names (gene IDs)")
  if (is.null(cell_ids))                  stop("count_input must have column names (cell IDs)")
  if (length(unique(gene_ids)) != n_gene) stop("Gene IDs (row names of count_input) are not unique")
  if (length(unique(cell_ids)) != n_cell) stop("Cell IDs (column names of count_input) are not unique")

  message(sprintf("count_input: %d genes x %d cells", n_gene, n_cell))

  # ── Validate meta_cell ─────────────────────────────────────────────────────
  missing_cell_cols <- setdiff(c("cell_id", "individual"), names(meta_cell))
  if (length(missing_cell_cols) > 0)
    stop(sprintf("meta_cell is missing columns: %s", paste(missing_cell_cols, collapse = ", ")))

  if (length(unique(meta_cell$cell_id)) != nrow(meta_cell))
    stop("cell_id values in meta_cell are not unique")

  if (any(meta_cell$cell_id != cell_ids))
    stop("cell_id in meta_cell does not match column names of count_input")

  # ── Validate meta_ind ──────────────────────────────────────────────────────
  required_ind_cols <- c("individual", var2test)
  missing_ind_cols  <- setdiff(required_ind_cols, names(meta_ind))
  if (length(missing_ind_cols) > 0)
    stop(sprintf("meta_ind is missing columns: %s", paste(missing_ind_cols, collapse = ", ")))

  if (length(unique(meta_ind$individual)) != nrow(meta_ind))
    stop("individual IDs in meta_ind are not unique")

  if (!setequal(meta_cell$individual, meta_ind$individual))
    stop("individual IDs in meta_cell and meta_ind do not match")

  if (!is.factor(meta_ind[[var2test]])) {
    message(sprintf("Converting '%s' to a factor", var2test))
    meta_ind[[var2test]] <- as.factor(meta_ind[[var2test]])
  }

  group_labels <- meta_ind[[var2test]]
  group_levels <- levels(group_labels)

  if (length(group_levels) < 2)
    stop("var2test must have at least two factor levels")

  message(sprintf("Group levels in '%s': %s",
                  var2test, paste(group_levels, collapse = ", ")))

  # ── Organise expression values by donor ───────────────────────────────────
  donor_cells <- lapply(meta_ind$individual,
                        function(ind) which(meta_cell$individual == ind))
  names(donor_cells) <- as.character(meta_ind$individual)

  dat_res <- lapply(seq_len(n_gene), function(i_g) {
    lapply(donor_cells, function(w) count_input[i_g, w])
  })

  if (verbose >= 1) {
    for (j in seq_along(dat_res[[1]])) {
      message(sprintf("Donor %d values (head): %s",
                      j, paste(utils::head(dat_res[[1]][[j]]), collapse = ", ")))
    }
  }

  # ── Compute pairwise divergences per gene ──────────────────────────────────
  if (verbose >= 2) message("Computing pairwise divergences...")

  dist_array_list <- if (ncores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_gene), function(i_g) {
      .compute_gene_dist_array(dat_res[[i_g]], meta_ind, group_labels, k, verbose)
    }, mc.cores = ncores)
  } else {
    lapply(seq_len(n_gene), function(i_g) {
      .compute_gene_dist_array(dat_res[[i_g]], meta_ind, group_labels, k, verbose)
    })
  }

  names(dist_array_list) <- gene_ids
  dist_array_list
}


# Build the donor x donor distance array for a single gene.
.compute_gene_dist_array <- function(res_ig, meta_ind, group_labels, k, verbose) {
  n_ind  <- nrow(meta_ind)
  donors <- meta_ind$individual

  dist_array <- array(
    NA_real_,
    dim      = c(n_ind, n_ind, 4),
    dimnames = list(donors, donors, c("was2", "location", "size", "shape"))
  )
  for (d in seq_len(4)) {
    slice <- dist_array[,, d]; diag(slice) <- 0; dist_array[,, d] <- slice
  }

  for (j_a in seq_len(n_ind)) {
    id_a    <- donors[j_a]
    label_a <- group_labels[j_a]

    same_cands  <- match(donors[group_labels == label_a & donors != id_a], donors)
    other_cands <- match(donors[group_labels != label_a], donors)

    if (!is.null(k)) {
      if (length(same_cands)  > k) same_cands  <- sample(same_cands,  k)
      if (length(other_cands) > k) other_cands <- sample(other_cands, k)
    }

    for (j_b in c(same_cands, other_cands)) {
      if (is.na(j_b) || j_b == j_a) next

      d <- tryCatch(
        divergence(res_ig[[j_a]], res_ig[[j_b]], verbose = verbose),
        error = function(e) rep(NA_real_, 4)
      )
      dist_array[j_a, j_b, ] <- d
      dist_array[j_b, j_a, ] <- d * c(1, -1, -1, 1)
    }
  }

  dist_array
}
