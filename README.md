
# ProteomicsML <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->
![R-CMD-check](https://github.com/arunkumaran07/Proteome/actions/workflows/R-CMD-check.yaml/badge.svg)
![pkgdown](https://github.com/arunkumaran07/Proteome/actions/workflows/pkgdown.yaml/badge.svg)
<!-- badges: end -->

Modular R package for cancer proteomics:
- Read wide CSVs with auto-detected group row (DDA/DIA abundance data)
- Differential expression (ref vs pooled / ref vs target / pairwise)
- Volcano plots
- Reactome GSEA + across-comparison summaries
- PCA, t-SNE, PLS-DA
- LASSO & Random Forest
- Heatmap clustering + module enrichment
- - **Interactive analysis wizard** for non-coders

## Installation

```r
# install.packages("remotes")
remotes::install_github("arunkumaran07/Proteome")
library(ProteomicsML)

Example: differential expression
# Demo dataset shipped with the package
csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")

io <- read_wide_proteomics(csv_path)
de <- diff_expr_ttest(io$expr, io$groups, ref = levels(io$groups)[1])

head(de)

Example: launch the wizard
The wizard is an interactive menu that guides you through the most common analyses step by step:

run_proteomics_wizard()

When run, you’ll see a menu in your R console:

=== ProteomicsML: Analysis Wizard ===
1: Differential expression + volcano
2: Reactome summary across comparisons (dot-plots)
3: PCA / t-SNE / PLS-DA
4: Machine learning: LASSO + Random Forest
5: Heatmap + module detection + enrichment
6: Quit

Choose a number to start the workflow you want.
This makes the package accessible to researchers who don’t want to write R code directly.

See full documentation at: https://arunkumaran07.github.io/Proteome/
