# Was2CoDE — CLAUDE.md

## Project overview

Was2CoDE is an R package for individual-level single-cell differential expression analysis. It decomposes the Wasserstein-2 (was2) distance between donor expression distributions into four components: `was2`, `location`, `size`, and `shape`. Built by Tati Zhang; Kevin Lin is the collaborator helping refactor and extend it.

## Key design decisions

### Input data contract
`was2code_dist()` accepts **pre-normalized, denoised** expression matrices. No log10 transformation, no residualization, and no `var_per_cell` parameter. All preprocessing is done upstream by the caller.

### Dependencies
Heavy Bioconductor packages (`DESeq2`, `dreamlet`, `EnhancedVolcano`, `eSVD2`, `nebula`) live in `Suggests`, not `Imports`. Each helper file (`deseq2_helper.R`, `dreamlet_helper.R`, etc.) guards its body with `requireNamespace()`. This avoids GenomeInfoDb startup noise when loading the package.

`doRNG`, `doParallel`, and `foreach` are **not used** — they were removed because `requireNamespace()` does not attach infix operators like `%dorng%`, causing silent NA failures via `tryCatch`. Parallelization uses `parallel::mclapply` / `lapply` directly (same pattern in both `was2code_dist` and `was2code_permanova`).

### Parallelization pattern
Both `was2code_dist` and `was2code_permanova` use the same idiom:
```r
if (ncores > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(..., mc.cores = ncores)
} else {
  lapply(...)
}
```
`ncores = 1` runs sequentially with no cluster setup.

### was2code_permanova
- Merged from two formerly separate functions (`was2code_permanova` + `was2code_permanova_na`).
- No covariate adjustment (`var2adjust`, `residualize_x`, `delta` were removed — data is already preprocessed).
- NA distance entries are handled via NA-aware F-statistic computation (`.calc_F_manova_na`); no genes are dropped.
- Returns a **list** with three elements: `pval` (named numeric vector), `F_ob` (named numeric vector), `F_perm` (genes × n_perm matrix). The `F_ob`/`F_perm` are exposed so callers can fit a parametric null distribution when empirical permutation p-values are too coarse.

## R gotchas found in this codebase

### diag<- on 3D array subsets
`diag(arr[,, d]) <- 0` does **not** write back to `arr` in R. Always extract, modify, and reassign:
```r
slice <- arr[,, d]; diag(slice) <- 0; arr[,, d] <- slice
```

### Array subsetting with drop
`arr[row_ids, col_ids, metric, drop = FALSE]` returns a 3D array; `diag<-` only works on 2D matrices. Drop the `drop = FALSE` argument when you need a matrix slice.

## Testing

Tests use `devtools::load_all()` and `testthat`. Run with:
```r
devtools::load_all(".")
testthat::test_file("tests/testthat/test_was2code_dist.R")
testthat::test_file("tests/testthat/test_was2code_permanova.R")
```

### Test data conventions
`test_was2code_dist.R` generates a canonical 4-gene × 8-donor × 100-cells-per-donor dataset via `.make_test_data()`:
- `null_gene` — all donors ~ N(0, 1); no signal
- `location_gene` — cases ~ N(5, 1), controls ~ N(0, 1)
- `size_gene` — cases ~ N(0, 3), controls ~ N(0, 0.3)
- `shape_gene` — cases bimodal ±3 (sd=0.3), controls ~ N(0, √9.09) [sd-matched]

Shape signal is subtler than location/size at this sample size (8 donors, 100 cells). Thresholds for shape tests use 2× not 5×.

`test_was2code_permanova.R` re-uses this same dataset (via `was2code_dist`) to test that permanova detects signal in the location/size/shape genes and not in the null gene.
