#' Run PERMANOVA on Wasserstein-2 Distance Matrices
#'
#' @description
#' Tests whether the grouping variable \code{var2test} explains a significant
#' fraction of variance in the Wasserstein-2 distance matrices produced by
#' \code{\link{was2code_dist}}, using a permutational MANOVA (PERMANOVA).
#' Missing (NA) distance entries — which arise when \code{k} is specified in
#' \code{was2code_dist} — are handled gracefully via an NA-aware Gower-centering
#' approach; no genes are dropped.
#'
#' @param dist_list Named list returned by \code{\link{was2code_dist}}.
#' @param meta_ind A data.frame with one row per donor. Must contain columns
#'   \code{individual} and the column named by \code{var2test}.
#' @param var2test Name of the grouping column in \code{meta_ind} to test.
#' @param n_perm Number of permutations (default \code{999}).
#' @param r_seed Integer random seed for reproducibility (default \code{2020}).
#' @param ncores Number of cores for parallelizing permutations (default \code{1}).
#' @return A list with three elements:
#'   \describe{
#'     \item{pval}{Named numeric vector of permutation p-values, one per gene.}
#'     \item{F_ob}{Named numeric vector of observed pseudo-F statistics, one per gene.}
#'     \item{F_perm}{Numeric matrix (genes × \code{n_perm}) of permuted pseudo-F statistics.
#'       Rows are named by gene; columns by permutation index. Useful for fitting a
#'       parametric null distribution when empirical p-values are too coarse.}
#'   }
#' @export

# code originally from https://github.com/Sun-lab/ideas/blob/main/R/permanova.R
# and https://github.com/Sun-lab/ideas/blob/main/R/permanova_utilities.R
was2code_permanova <- function(dist_list,
                               meta_ind,
                               var2test,
                               n_perm  = 999,
                               r_seed  = 2020,
                               ncores  = 1) {

  # ── Validate meta_ind ──────────────────────────────────────────────────────
  if (!is.data.frame(meta_ind))
    stop("meta_ind must be a data.frame")

  missing_cols <- setdiff(c("individual", var2test), names(meta_ind))
  if (length(missing_cols) > 0)
    stop(sprintf("meta_ind is missing columns: %s", paste(missing_cols, collapse = ", ")))

  if (length(unique(meta_ind$individual)) != nrow(meta_ind))
    stop("individual IDs in meta_ind are not unique")

  # ── Build genes × donors × donors distance array (was2 slice) ─────────────
  dist_array <- .form_dist_array(dist_list)

  n_na <- sum(apply(dist_array, 1, anyNA))
  if (n_na > 0)
    message(sprintf(
      "%d gene(s) have NA distance entries (partial comparisons); using NA-aware PERMANOVA.", n_na))

  # ── Prepare 0/1 group label vector ────────────────────────────────────────
  x <- meta_ind[[var2test]]
  if (anyNA(x)) stop("var2test contains NA values")

  message(sprintf("Testing '%s' with %d permutations", var2test, n_perm))

  if (is.character(x)) x <- as.factor(x)
  if (is.factor(x))   x <- as.numeric(x)
  x <- as.numeric(x == max(x))

  # ── Generate permutation matrix (donors × n_perm) ─────────────────────────
  set.seed(r_seed)
  x_perm <- apply(matrix(rep(x, n_perm), ncol = n_perm), 2, sample, size = length(x))

  # ── Compute observed and permuted F-statistics ─────────────────────────────
  F_ob   <- .calc_F_manova_na(dist_array, label = x)
  F_perm <- do.call(cbind, .perm_lapply(ncores, seq_len(n_perm), function(i) {
    .calc_F_manova_na(dist_array, label = x_perm[, i])
  }))

  # ── p-value: fraction of permutations with F >= observed F ────────────────
  gene_names     <- names(dist_list)
  F_diff         <- sweep(F_perm, 1, F_ob, `-`)
  pval           <- rowMeans(cbind(F_diff, 1) >= 0, na.rm = TRUE)
  names(pval)    <- gene_names
  names(F_ob)    <- gene_names
  rownames(F_perm) <- gene_names
  colnames(F_perm) <- paste0("perm", seq_len(n_perm))

  list(pval = pval, F_ob = F_ob, F_perm = F_perm)
}


# ─── Internal helpers ──────────────────────────────────────────────────────────

# Reshape dist_list into a genes × donors × donors array (was2 slice only).
.form_dist_array <- function(dist_list) {
  gene_names  <- names(dist_list)
  donor_names <- dimnames(dist_list[[1]])[[1]]

  arr <- array(NA_real_,
               dim      = c(length(gene_names), length(donor_names), length(donor_names)),
               dimnames = list(gene_names, donor_names, donor_names))
  for (gene in gene_names)
    arr[gene, , ] <- dist_list[[gene]][donor_names, donor_names, "was2"]
  arr
}

# lapply over indices, optionally parallelized via mclapply (Unix only).
.perm_lapply <- function(ncores, indices, FUN) {
  if (ncores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(indices, FUN, mc.cores = ncores)
  } else {
    lapply(indices, FUN)
  }
}

# Pseudo-F statistic for PERMANOVA (no covariates), NA-aware.
.calc_F_manova_na <- function(dist_array, label) {
  n_genes <- dim(dist_array)[1]
  N       <- length(label)
  groups  <- split(seq_len(N), label)
  a       <- length(groups)

  Fstat <- numeric(n_genes)
  for (g in seq_len(n_genes)) {
    d2    <- dist_array[g, , ]^2
    avail <- !is.na(d2)
    sst   <- sum(d2, na.rm = TRUE) / N

    ssw <- 0
    for (ids in groups) {
      if (length(ids) < 2) next
      d2_w    <- d2[ids, ids, drop = FALSE]
      avail_w <- avail[ids, ids, drop = FALSE]
      ssw <- ssw + sum(d2_w, na.rm = TRUE) / max(1, sum(avail_w)) * length(ids)
    }
    Fstat[g] <- ((sst - ssw) * (N - a)) / (ssw * (a - 1))
  }
  Fstat
}
