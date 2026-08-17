
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
- **Interactive analysis wizard** for non-coders

## Installation

```r
# install.packages("remotes")
remotes::install_github("arunkumaran07/Proteome")
library(ProteomicsML)
```

## Preparing your CSV file

Before running any analysis, your data needs to be a **wide-format CSV**: proteins/genes as rows, samples as columns, with one column holding the protein/gene ID. The package expects a specific layout:

1. **Column 1** — protein/gene identifiers (e.g. gene symbols like `TP53`, `EGFR`).
2. **Remaining columns** — one per sample, with abundance/expression values.
3. **A group row** — the first (or second) data row must contain the group label for each sample (e.g. `Normal` / `Tumour`), in the same column position as that sample's data. `read_wide_proteomics()` auto-detects whether this row is row 1 or row 2.

Example (`example_proteome.csv`, shipped with the package):

```
Gene,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6
Group,Normal,Normal,Normal,Tumour,Tumour,Tumour
TP53,5.1,4.8,5.3,7.2,6.9,7.5
EGFR,2.0,1.8,2.1,3.5,3.7,3.6
BRCA1,8.2,8.1,8.3,6.0,6.2,6.1
MYC,1.1,1.0,1.2,4.5,4.7,4.6
```

Notes:
- Column 1 header can be anything (`Gene`, `Protein`, `ID`, etc.) — it isn't used for matching, only its position (default `id_col = 1`).
- The group row's ID cell (`Group` above) is ignored; only the values in the sample columns are read as group labels.
- Non-numeric or blank sample columns are dropped automatically.
- Save the file as standard CSV (comma-separated, UTF-8). Excel: **File > Save As > CSV UTF-8 (Comma delimited)**.

### Reading the CSV into R

Use `read_wide_proteomics()` to load your file before calling any other function:

```r
# Using your own file:
io <- read_wide_proteomics("path/to/your_proteomics_data.csv")

# Or the demo dataset shipped with the package:
csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
io <- read_wide_proteomics(csv_path)
```

`read_wide_proteomics()` arguments:

| Argument | Default | Description |
|---|---|---|
| `path` | — | Path to the CSV file |
| `id_col` | `1` | Column index containing protein/gene IDs |
| `group_row` | `"auto"` | Which row holds group labels: `"auto"`, `1`, or `2` |

It returns a list:
- `io$expr` — numeric matrix of expression values (proteins x samples)
- `io$groups` — factor of group labels (one per sample)
- `io$samples` — sample names

Once loaded, `io$expr` and `io$groups` are the inputs to every downstream function (`diff_expr_ttest()`, `run_pca()`, `train_lasso()`, `cluster_proteome()`, etc.).

## Example: differential expression

```r
csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")

io <- read_wide_proteomics(csv_path)
de <- diff_expr_ttest(io$expr, io$groups, ref = levels(io$groups)[1])

head(de)
```

## Example: launch the wizard

The wizard is an interactive menu that guides you through the most common analyses step by step — it will prompt you for your CSV path directly, so you don't need to call `read_wide_proteomics()` yourself:

```r
run_proteomics_wizard()
```

When run, you'll see a menu in your R console:

```
=== ProteomicsML: Analysis Wizard ===
1: Differential expression + volcano
2: Reactome summary across comparisons (dot-plots)
3: PCA / t-SNE / PLS-DA
4: Machine learning: LASSO + Random Forest
5: Heatmap + module detection + enrichment
6: Quit
```

Choose a number to start the workflow you want. When prompted for a CSV path, enter the path to your wide-format file as described above. This makes the package accessible to researchers who don't want to write R code directly.

See full documentation at: https://arunkumaran07.github.io/Proteome/
