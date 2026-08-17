# ProteomicsML 0.2.7

## Documentation

* Added complete roxygen documentation (title, parameter descriptions,
  return value, `\dontrun{}` examples) for all 15 exported functions,
  fixing "Undocumented code objects" and "No examples" NOTEs from
  `R CMD check`.
* Added the `man/` directory (was previously missing from the repo
  entirely, despite `NAMESPACE` declaring exports).
* Removed inaccurate `LazyData: true` from `DESCRIPTION` (no `data/`
  directory exists).

# ProteomicsML 0.2.6

## CI

* `R-CMD-check.yaml` now installs only hard dependencies (Depends/Imports/
  LinkingTo) before running check, instead of forcing `pak` to resolve every
  optional `Suggests` package. Several Suggests are heavy, Bioconductor-only
  packages used solely behind `requireNamespace()` guards at runtime
  (`clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `enrichplot`,
  `EnhancedVolcano`); pre-resolving them in CI was fragile and caused
  R-CMD-check to fail during dependency installation in 0.2.3-0.2.5, before
  the actual check ever ran. `_R_CHECK_FORCE_SUGGESTS_=false` (already set)
  means `R CMD check` itself tolerates their absence.

# ProteomicsML 0.2.5

## Bug fixes

* `EnhancedVolcano` (used in `plot_volcano()`) is also a Bioconductor-only
  package but was missing from `Remotes:`, causing `pak` to fail dependency
  resolution in CI (0.2.4 R-CMD-check still failed). Added
  `bioc::EnhancedVolcano`.

# ProteomicsML 0.2.4

## Bug fixes

* Fixed invalid `Remotes:` syntax in `DESCRIPTION` (`bioc::release/<pkg>` is
  not valid pak/remotes syntax and broke CI dependency installation).
  Corrected to `bioc::<pkg>`, which resolves the Bioconductor release
  matching the installed R version.

# ProteomicsML 0.2.3

## Bug fixes

* Declared all package dependencies actually used by the code in `DESCRIPTION`.
  Previously `dplyr`, `ggplot2`, `tibble`, and `tools` were used unconditionally
  but not listed under `Imports`, and optional packages used for GSEA/enrichment
  and ML workflows (`clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `enrichplot`,
  `Rtsne`, `pls`, `glmnet`, `randomForest`, `pheatmap`, `dynamicTreeCut`,
  `EnhancedVolcano`) were not listed under `Suggests`. This could cause
  `remotes::install_github()` to silently skip needed dependencies and functions
  to fail at runtime with "there is no package called ..." errors.
* Added a `Remotes:` field so Bioconductor-only dependencies
  (`clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `enrichplot`) can be
  resolved automatically.

## Documentation

* README now documents the CSV format `read_wide_proteomics()` expects,
  with a worked example, before the usage examples.
* README now notes that enrichment/GSEA features require Bioconductor
  packages installed via `BiocManager::install()`.

# ProteomicsML 0.2.2

* Initial tracked release.
