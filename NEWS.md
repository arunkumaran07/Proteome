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
