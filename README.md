# ProteomicsML 0.2.2

Modular R package for cancer proteomics:
- Read wide CSVs with auto-detected group row
- Differential expression (ref vs pooled / ref vs target / pairwise)
- Volcano plots
- Reactome GSEA + across-comparison summaries
- PCA, t-SNE, PLS-DA
- LASSO & Random Forest
- Heatmap clustering + module enrichment

## Install
```r
# install.packages("devtools")
devtools::install_github("arunkumaran07/Proteome")
library(ProteomicsML)
run_proteomics_wizard()
```
