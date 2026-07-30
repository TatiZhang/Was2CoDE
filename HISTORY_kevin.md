# HISTORY_kevin.md — Kevin's Session Log

> **Append-only, ascending chronological order** (oldest at top, newest at the bottom). Add each session's dated entry to the END of this file. Never read at session startup — consulted only on demand for deep history. Current project state lives in `CLAUDE_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session appends to this file. Other collaborators may read it but must not edit or rewrite entries.

---

### [2026-07-30] (Session 1 — project-setup / project-state scaffolding)
- Ran `/project-setup` and `/project-state` on the `Was2CODE` package repo, which previously had only a bare technical `CLAUDE.md` and none of the three-file collaboration pattern.
- Rewrote the master `CLAUDE.md` onto the standard template: workflow instructions, project context, repository layout (with the three sibling repos), "Who Is Using This Session?" table + file-ownership rules, wiki pointer, post-prompt update instructions, key conventions. All prior technical content (input contract, dependency policy, parallelization idiom, `was2code_permanova` return shape, the two R gotchas, testing conventions) was preserved verbatim below a `---` divider rather than rewritten.
- Created `CLAUDE_kevin.md` and `HISTORY_kevin.md` from the skill templates. Did **not** create `CLAUDE_tati.md` / `HISTORY_tati.md` — per the ownership rule, Tati's files appear when she runs her own session; her row in the who-is-using table is a placeholder with login TBD.
- Merged the standard `.gitignore` template into the existing one (added OS/editor cruft, Python, and LaTeX-artifact blocks; kept the repo's existing `docs`, `pkgdown/.DS_Store`, `*.icloud`, `.cloud` entries). Deliberately did **not** adopt the template's `*.Rproj` ignore — `Was2CoDE.Rproj` is already tracked in this repo, and ignoring it would confuse rather than help.
- Added `.githooks/` (50 MB staged-file guard) and ran `git config core.hooksPath .githooks` for this clone. Non-obvious: git will not run tracked hooks until that config is set, and it is per-clone — Tati must run it after her next pull.
- Added `^CLAUDE_kevin\.md$`, `^HISTORY_kevin\.md$`, `^brainstorming_.*\.md$`, and `^\.githooks$` to `.Rbuildignore`. Non-obvious rationale specific to this repo: unlike the paper/analysis repos, this one is an **R package**, so any root-level file not listed in `.Rbuildignore` gets shipped inside the package tarball and trips `R CMD check`. Future per-person files need the same treatment — noted as a convention in the master `CLAUDE.md`.
- Learned/recorded: the git remote is `github.com/TatiZhang/IdeasCustom`, not `Was2CoDE`. The paper repo's `CLAUDE.md` currently states the remote as `github.com/TatiZhang/Was2CoDE`, which is wrong; the directory name (`Was2CODE`), package name (`Was2CoDE`), pkgdown URL (`tatizhang.github.io/Was2CoDE/`), and remote (`IdeasCustom`) are four different spellings of the same thing. `README.md` still describes the package as "IdeasCustom" with a stale 2024 session-info dump.
- Learned/recorded: only 7 functions are exported; the four DE-method wrappers (`deseq2_helper.R`, `dreamlet_helper.R`, `nebula_helper.R`, `esvd_helper.R`) plus `plot_volcano` and `compute_gsea_overlap` are internal, so the analysis repo must be reaching them via `:::` or `devtools::load_all()`.
- Open: whether the package's center of gravity should shift from the was2 decomposition to the DE-method wrappers, now that the Was2CoDE method was dropped from the paper (2026-07-22) while the wrappers remain on its critical path. Recorded as Open Question 2 in `CLAUDE_kevin.md`.
- Open: the was2 decomposition is now unpublished with no venue — a standalone methods paper is unplanned (Open Question 3).
- Open: no test coverage for any of the four DE-method helpers or the plotting functions (Open Question 4).
- Open: `raw/*.pdf` are git-LFS pointer stubs (131–132 bytes on disk), so those PDFs are not actually present locally.
