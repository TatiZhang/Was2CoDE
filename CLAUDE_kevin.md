# CLAUDE_kevin.md — Kevin's Context

> **Current state only.** Every section below is updated *in place* each session — overwrite, don't append. The append-only dated log lives in `HISTORY_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session writes this file and `HISTORY_kevin.md`. Other collaborators may read it for context but must not edit it.

## About Kevin
- Role in project: PI / collaborator. Tati is the package author-maintainer; Kevin refactors, extends, and consumes the package from the analysis repo.
- Background: Biostatistics (UW). Statistical methodology for single-cell / snRNA-seq; developed eSVD-DE (one of the four DE methods wrapped here).
- Email: kzlin@uw.edu

## Environment (local paths)
- Package root: `~/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/tati/git/Was2CODE/` (git remote → `github.com/TatiZhang/IdeasCustom`; pkgdown → `tatizhang.github.io/Was2CoDE/`)
- Analysis repo: `../Was2CoDE_analysis/` (`github.com/linnykos/Was2CoDE_analysis`); `code/kevin/` and `code/tati/`
- Paper: `../overleaf_tati_subject-de_paper/` (`github.com/linnykos/tati_subject-de_paper`)
- Wiki: `../Was2CoDE_wiki/` (178 pages); separate Figure-1-only wiki at `../Was2CoDE_paper-summary/`
- Compute: Bayes server; data outputs under `~/kzlinlab/projects/subject-de/out/`
- Local R: run tests with `devtools::load_all(".")` then `testthat::test_file(...)`

## Project Status (as of 2026-07-30)
- **Scaffolding session.** This repo previously had only a bare technical `CLAUDE.md`. It now carries the full three-file pattern: master `CLAUDE.md` (workflow + layout + who-is-using table + wiki pointer, with the prior technical content preserved below a divider), plus this file and `HISTORY_kevin.md`. `.gitignore` merged with the standard template, `.githooks/` large-file guard added, `.Rbuildignore` extended for the new root `.md` files.
- **Code state** (last substantive commit `6ef6d20`, "MASSIVE update due to claude code"): `was2code_dist` and `was2code_permanova` were heavily refactored — `was2code_permanova_NA.R` and `form_dist_array.R` deleted, the two permanova variants merged into one NA-aware function, `foreach`/`doRNG`/`doParallel` removed in favor of `mclapply`/`lapply`. Working tree is clean.
- **Exported surface (7 functions)**: `arrange_genes_by_donors`, `divergence`, `plot_signed_logpval`, `result_array_list`, `was2code_dist`, `was2code_lfc`, `was2code_permanova`. The `*_helper.R` DE-method wrappers (DESeq2, Dreamlet, NEBULA, eSVD-DE) and `plot_volcano` / `compute_gsea_overlap` are **not exported** — the analysis repo reaches them via `:::` or `load_all()`.
- **The Was2CoDE method is dropped from the companion paper** (2026-07-22 decision, recorded in the paper's `CLAUDE.md`), but the package is still on the paper's critical path because its `*_helper.R` wrappers are what the 4-method benchmark runs through. Package and paper have decoupled: this repo's future is as a standalone methods contribution.

## Key Methodological Details
- **Input contract**: `was2code_dist()` takes pre-normalized, denoised matrices. No log10, no residualization, no `var_per_cell` — all preprocessing is upstream.
- **Decomposition**: was2 distance between donor expression distributions → `was2`, `location`, `size`, `shape`. Location and size are far easier to detect than shape at small n (8 donors × 100 cells → shape tests use a 2× threshold, not 5×).
- **`was2code_permanova` returns a list**, not a vector: `pval`, `F_ob`, `F_perm` (genes × n_perm). `F_ob`/`F_perm` are deliberately exposed so a *parametric* null can be fit when empirical permutation p-values are too coarse — this is the intended path when p-value resolution limits downstream FDR.
- **NA handling**: NA distance entries go through `.calc_F_manova_na` rather than dropping genes.
- **Dependency policy**: heavy Bioconductor packages in `Suggests` with `requireNamespace()` guards, to keep `library(Was2CoDE)` quiet and fast.

## Open Questions / Next Steps
1. Open: **Repo/package/remote name mismatch.** Directory is `Was2CODE`, package is `Was2CoDE`, GitHub remote is still `TatiZhang/IdeasCustom`, and `README.md` still describes "IdeasCustom" with a stale 2024 session-info block. Decide whether to rename the remote (Tati's call — she owns it) and rewrite the README.
2. Open: **Does the package stay Was2CoDE-first or become a DE-wrapper toolkit?** Now that the method is out of the paper, the most-used code here is the four `*_helper.R` wrappers. If those are the real product, they should be exported and documented rather than reached via `:::`.
3. Open: **Standalone methods paper for the was2 decomposition?** Dropping it from the microglia paper leaves the method unpublished. No venue or timeline decided.
4. Open: **Test coverage gaps.** `tests/testthat/` covers `was2code_dist`, `was2code_permanova`, `was2code_lfc`, `divergence`, `arrange_genes_by_donors`, `result_array_list` — nothing covers the four DE-method helpers or the plotting functions.
5. Open: **`git config core.hooksPath .githooks` must be run once per clone** (Kevin's clone is now configured; Tati's is not). Flag it to Tati.
6. Open: `raw/*.pdf` are git-LFS pointers (131–132 bytes on disk). If those PDFs are needed, LFS content must be pulled.
