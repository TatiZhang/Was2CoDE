# CLAUDE.md — Was2CoDE R Package

## Workflow Instructions
1. **Always enter plan mode** before starting any non-trivial task.
2. **Use superpower skills** where relevant: `/code-review` and `/verify` for code changes, `/brainstorm` for research directions, `/wiki-query` / `/wiki-ingest` for the project wiki.
3. **After every prompt**, run `/project-state`: refresh the current-state sections of the individual contributor's `CLAUDE_[name].md` in place, and append a dated entry to their `HISTORY_[name].md`. Write only non-obvious things; skip anything already in the code or git history.
4. **This is the method repo, not the analysis or paper repo.** Code here is a general-purpose R package. Dataset-specific analysis lives in `../Was2CoDE_analysis/`; narrative claims live in `../overleaf_tati_subject-de_paper/`. Don't hard-code Prater/Gabitto/Sun/Green specifics into package functions.

## Project Context (High Level)
**Project**: `Was2CoDE` — an R package for individual-level (donor-level) single-cell differential expression. Its namesake method decomposes the Wasserstein-2 (was2) distance between donor expression distributions into four components — `was2`, `location`, `size`, `shape` — and tests them with a PERMANOVA. The package *also* hosts thin wrappers around the four DE methods benchmarked in the companion paper.

**Authors / roles**:
- **Wenjing "Tati" Zhang** (UW Biostat) — package author/maintainer (`wenjiz1@uw.edu`, ORCID 0009-0007-8949-5933)
- **Kevin Z. Lin** (UW Biostatistics) — PI/collaborator, refactoring and extending (`kzlin@uw.edu`)

**Important status note.** The **Was2CoDE method itself was dropped from the companion paper** (decision recorded 2026-07-22 in `../overleaf_tati_subject-de_paper/CLAUDE.md`) to sharpen the narrative. The package remains live and useful for two reasons: (a) it is the reference implementation of the was2 decomposition, potentially a standalone methods contribution; (b) its `*_helper.R` files are the DE-method wrappers (**DESeq2, Dreamlet, NEBULA, eSVD-DE**) that the analysis repo calls. Do not delete method code on the assumption it is unused by the paper.

## Repository Layout
This package lives inside a larger `tati/git/` folder alongside the analysis, paper, and wiki repos:
- `Was2CODE/` — **this repo** (project root; git remote → `github.com/TatiZhang/IdeasCustom`; pkgdown → `tatizhang.github.io/Was2CoDE/`). Note the **directory name, package name, and GitHub repo name all differ** — the repo was renamed from `IdeasCustom` and the remote still carries the old name.
  - `R/` — package source. Core method: `was2code_dist.R`, `was2code_permanova.R`, `was2code_lfc.R`, `divergence.R`, `arrange_genes_by_donors.R`, `result_array_list.R`, `manova.R`. DE-method wrappers: `deseq2_helper.R`, `dreamlet_helper.R`, `nebula_helper.R`, `esvd_helper.R`. Plotting/eval: `plot_signed_logpval.R`, `plot_volcano.R`, `compute_gsea_overlap.R`.
  - `tests/testthat/` — testthat suite (see **Testing** below).
  - `man/`, `NAMESPACE` — roxygen-generated; never hand-edit.
  - `data/` — bundled `.rda`: `housekeeping_hounkpe_df`, `housekeeping_lin_df`, `microglia_prater_df`, `microglia_sun_df`.
  - `vignettes/Was2CoDE.Rmd`, `pkgdown/`, `_pkgdown.yml` — docs site.
  - `raw/` — poster/paper PDFs (git-LFS pointers); excluded from the build via `.Rbuildignore`.
- `../Was2CoDE_analysis/` — analysis-reproducibility repo (`github.com/linnykos/Was2CoDE_analysis`); `code/kevin/Writeup*` and `code/tati/Writeup*`, plus `csv/`, `figures/`, `notes/`. This is the **only** consumer of this package that matters.
- `../overleaf_tati_subject-de_paper/` — the paper (`github.com/linnykos/tati_subject-de_paper`). Read its `CLAUDE.md` / `PAPER_OUTLINE.md` for the current narrative; the paper is mid-revamp and `main.tex` is outdated.
- `../Was2CoDE_wiki/` — the project wiki. See **Wiki** below.
- `../Was2CoDE_paper-summary/` — a *separate* wiki, used only for the paper's Figure-1 embedding of 500+ papers. Don't conflate it with the main wiki.

## Who Is Using This Session?
**Detect the current user** by running: `echo $USER`. This table maps each login to that person's **first-name** context file; it is the source of truth that `/brainstorm` and `/project-state` use to resolve `$USER` to the right filename (so the login `kevinlin` maps to `CLAUDE_kevin.md`, never `CLAUDE_kevinlin.md`).

| Username (login) | Current-state file (first name) | History archive |
|---|---|---|
| `kevinlin` | `CLAUDE_kevin.md` | `HISTORY_kevin.md` |
| _(tati — login TBD)_ | `CLAUDE_tati.md` | `HISTORY_tati.md` |

**File ownership.** Each row above names one person's files, and **only that person's session writes them.** Once `$USER` resolves to a first name, that is the only suffix you may create, edit, append to, rename, or delete — every other collaborator's `CLAUDE_[name].md`, `HISTORY_[name].md`, and `brainstorming_[name].md` is read-only. Read them for context when useful; never modify them, not even to fix a typo or add a cross-reference. This repository is shared, so an edit lands in the owner's working copy immediately and can overwrite state they wrote from their own machine. If a collaborator's file looks wrong or stale, say so instead of editing it. This master `CLAUDE.md` is the exception: it is shared and any collaborator may update it. The single override is the user, in the current turn, *directing you to write* that exact file ("add this to `CLAUDE_tati.md`") — confirm once, then write it. Merely mentioning a collaborator, referring to their file, or being away does not qualify, and attribution does not launder the edit — a change stamped with your name and confined to two lines is still a write to a file you do not own. Record the decision here in the master `CLAUDE.md` and tell the user to contact the owner directly.

After detecting the user, **read that person's `CLAUDE_[name].md` immediately** before doing any other work — it holds the current project state and restores context in ~30 seconds. Do **not** read `HISTORY_[name].md` at startup; it is the append-only session log, consulted only on demand when deep history is needed. If no match is found, ask who the user is and create a new `CLAUDE_[name].md` via `/project-state`. Fill in a collaborator's login the first time they run a session here.

## Wiki
`../Was2CoDE_wiki/` (`/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/tati/git/Was2CoDE_wiki`) is the project's **second brain** (178 pages of ingested sources + synthesis). Pages most relevant to this package: `[[was2code-method]]`, `[[analysis-donors-vs-cells-power]]`, `[[analysis-deg-count-variability]]`, `[[analysis-reproducible-microglia-genes]]`, `[[meta-analysis-de]]`. Run `/wiki-query` to query it and `/wiki-ingest` to add sources; `/brainstorm` and `/project-state` detect it from this `## Wiki` pointer.

## Post-Prompt Update Instructions
After completing each user prompt, run `/project-state`. It will:
- **Refresh in place** the current-state sections of `CLAUDE_[name].md` (Project Status, Key Methodological Details, Open Questions / Next Steps).
- **Append a dated entry at the bottom** of `HISTORY_[name].md` recording new decisions, resolved/open questions, non-obvious code rationale, and empirical findings.

Do NOT record: things already in the code, git history, or reproducible from the code.

## Key Conventions
- New R files drafted by Claude get a `_claude` suffix so the human reviews before integrating.
- Per-person files use **first names** (`CLAUDE_kevin.md`), resolved from logins via the "Who Is Using This Session?" table.
- `man/*.Rd` and `NAMESPACE` are **generated** — edit the roxygen block above the function and re-run `devtools::document()`.
- Every root-level `.md` scaffolding file (`CLAUDE*.md`, `HISTORY*.md`, `brainstorming*.md`) must be listed in `.Rbuildignore`, or `R CMD check` ships it inside the package tarball.

---

# Package Design Decisions

## Input data contract
`was2code_dist()` accepts **pre-normalized, denoised** expression matrices. No log10 transformation, no residualization, and no `var_per_cell` parameter. All preprocessing is done upstream by the caller.

## Dependencies
Heavy Bioconductor packages (`DESeq2`, `dreamlet`, `EnhancedVolcano`, `eSVD2`, `nebula`) live in `Suggests`, not `Imports`. Each helper file (`deseq2_helper.R`, `dreamlet_helper.R`, etc.) guards its body with `requireNamespace()`. This avoids GenomeInfoDb startup noise when loading the package.

`doRNG`, `doParallel`, and `foreach` are **not used** — they were removed because `requireNamespace()` does not attach infix operators like `%dorng%`, causing silent NA failures via `tryCatch`. Parallelization uses `parallel::mclapply` / `lapply` directly (same pattern in both `was2code_dist` and `was2code_permanova`).

## Parallelization pattern
Both `was2code_dist` and `was2code_permanova` use the same idiom:
```r
if (ncores > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(..., mc.cores = ncores)
} else {
  lapply(...)
}
```
`ncores = 1` runs sequentially with no cluster setup.

## was2code_permanova
- Merged from two formerly separate functions (`was2code_permanova` + `was2code_permanova_na`).
- No covariate adjustment (`var2adjust`, `residualize_x`, `delta` were removed — data is already preprocessed).
- NA distance entries are handled via NA-aware F-statistic computation (`.calc_F_manova_na`); no genes are dropped.
- Returns a **list** with three elements: `pval` (named numeric vector), `F_ob` (named numeric vector), `F_perm` (genes × n_perm matrix). The `F_ob`/`F_perm` are exposed so callers can fit a parametric null distribution when empirical permutation p-values are too coarse.

# R gotchas found in this codebase

## diag<- on 3D array subsets
`diag(arr[,, d]) <- 0` does **not** write back to `arr` in R. Always extract, modify, and reassign:
```r
slice <- arr[,, d]; diag(slice) <- 0; arr[,, d] <- slice
```

## Array subsetting with drop
`arr[row_ids, col_ids, metric, drop = FALSE]` returns a 3D array; `diag<-` only works on 2D matrices. Drop the `drop = FALSE` argument when you need a matrix slice.

# Testing

Tests use `devtools::load_all()` and `testthat`. Run with:
```r
devtools::load_all(".")
testthat::test_file("tests/testthat/test_was2code_dist.R")
testthat::test_file("tests/testthat/test_was2code_permanova.R")
```

## Test data conventions
`test_was2code_dist.R` generates a canonical 4-gene × 8-donor × 100-cells-per-donor dataset via `.make_test_data()`:
- `null_gene` — all donors ~ N(0, 1); no signal
- `location_gene` — cases ~ N(5, 1), controls ~ N(0, 1)
- `size_gene` — cases ~ N(0, 3), controls ~ N(0, 0.3)
- `shape_gene` — cases bimodal ±3 (sd=0.3), controls ~ N(0, √9.09) [sd-matched]

Shape signal is subtler than location/size at this sample size (8 donors, 100 cells). Thresholds for shape tests use 2× not 5×.

`test_was2code_permanova.R` re-uses this same dataset (via `was2code_dist`) to test that permanova detects signal in the location/size/shape genes and not in the null gene.
