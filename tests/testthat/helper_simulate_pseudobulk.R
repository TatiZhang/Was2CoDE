# Shared simulated single-cell data for the pseudobulk / DESeq2 / TOST tests.
# testthat sources helper*.R files before running any test file, so both
# test_pseudobulk.R, test_deseq2_helper.R and test_deseq2_tost_helper_claude.R
# can all use this.

# Donor-level counts with three gene classes:
#   - "de"    : true |log2FC| = 2, well expressed   -> should be differential
#   - "null"  : true log2FC = 0,   well expressed   -> should be provably stable
#   - "noisy" : true log2FC = 0,   barely expressed -> should be undetermined
#
# The DE genes are half up and half down, and are a minority of the panel, so
# that DESeq2's median-of-ratios size factors are not dragged by composition
# bias -- with a large, one-directional DE fraction the case samples get
# systematically down-scaled and every null gene picks up a spurious shift.
#
# Metadata columns: donor (id), diagnosis (case/control), sex (donor-constant
# categorical), age (donor-constant numerical), n_umi (varies within donor, so
# it is NOT a valid donor-level covariate -- used to test that rejection).
.simulate_donor_seurat <- function(n_donor_per_group = 8,
                                   n_cell_per_donor = 200,
                                   n_gene = c(de = 20, null = 60, noisy = 30),
                                   dispersion = 0.02,
                                   seed = 10) {
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  set.seed(seed)

  n_donor <- 2 * n_donor_per_group
  donor_names <- sprintf("D%02d", 1:n_donor)
  donor_group <- rep(c("control", "case"), each = n_donor_per_group)
  names(donor_group) <- donor_names

  gene_class <- rep(names(n_gene), times = n_gene)
  gene_names <- paste0(gene_class, "-",
                       sprintf("%03d", unlist(lapply(n_gene, seq_len), use.names = FALSE)))
  n_gene_total <- length(gene_names)

  # per-cell baseline mean, on the natural scale
  base_mu <- ifelse(gene_class == "noisy", 0.03, 5)
  # true log2 fold change (case vs control), alternating sign among DE genes
  true_lfc <- ifelse(gene_class == "de", 2, 0)
  true_lfc[gene_class == "de"] <- true_lfc[gene_class == "de"] *
    rep(c(1, -1), length.out = sum(gene_class == "de"))

  count_list <- lapply(donor_names, function(donor) {
    is_case <- donor_group[donor] == "case"
    # donor-level random effect: gamma with mean 1 and variance = dispersion
    donor_effect <- stats::rgamma(n_gene_total, shape = 1 / dispersion, rate = 1 / dispersion)
    mu <- base_mu * donor_effect * (2^(true_lfc * is_case))
    mat <- matrix(stats::rpois(n_gene_total * n_cell_per_donor, lambda = rep(mu, times = n_cell_per_donor)),
                  nrow = n_gene_total, ncol = n_cell_per_donor)
    rownames(mat) <- gene_names
    colnames(mat) <- paste0(donor, "c", seq_len(n_cell_per_donor))
    mat
  })

  mat <- Matrix::Matrix(do.call(cbind, count_list), sparse = TRUE)
  metadata <- data.frame(donor = rep(donor_names, each = n_cell_per_donor),
                         diagnosis = rep(donor_group, each = n_cell_per_donor),
                         row.names = colnames(mat),
                         stringsAsFactors = FALSE)

  # donor-constant covariates
  donor_age <- stats::setNames(stats::runif(n_donor, 60, 90), donor_names)
  donor_sex <- stats::setNames(rep(c("F", "M"), length.out = n_donor), donor_names)
  metadata$age <- as.numeric(donor_age[metadata$donor])
  metadata$sex <- as.character(donor_sex[metadata$donor])
  # a cell-level covariate, i.e. one that is NOT donor-constant
  metadata$n_umi <- Matrix::colSums(mat)

  seurat_obj <- Seurat::CreateSeuratObject(counts = mat, meta.data = metadata)
  Seurat::VariableFeatures(seurat_obj) <- gene_names

  list(seurat_obj = seurat_obj,
       gene_class = stats::setNames(gene_class, gene_names),
       true_lfc = stats::setNames(true_lfc, gene_names),
       donor_group = donor_group,
       donor_age = donor_age,
       donor_sex = donor_sex)
}
