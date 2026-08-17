# ---- IO ----

#' @importFrom stats p.adjust
NULL

#' Read a wide-format proteomics CSV
#'
#' Reads a "wide" proteomics CSV (proteins/genes in rows, samples in columns)
#' where one column holds protein/gene IDs and one row holds the group label
#' for each sample. The group row location is auto-detected by default. See
#' the package README for the exact expected CSV layout, or use the demo
#' dataset shipped in \code{system.file("extdata", "example_proteome.csv",
#' package = "ProteomicsML")}.
#'
#' @param path Path to the CSV file.
#' @param id_col Integer column index holding protein/gene IDs. Default \code{1}.
#' @param group_row Which row holds the group labels: \code{"auto"} (default),
#'   \code{1}, or \code{2}.
#' @return A list with components:
#'   \item{expr}{Numeric matrix of expression values (proteins x samples).}
#'   \item{groups}{Factor of group labels, one per sample.}
#'   \item{samples}{Character vector of sample names.}
#'   The list also carries a \code{"group_row"} attribute recording which row
#'   was used for group labels.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' io$expr
#' io$groups
#' }
#' @export
read_wide_proteomics <- function(path, id_col = 1, group_row = c("auto",1,2)) {
  group_row <- match.arg(as.character(group_row), c("auto","1","2"))
  df <- utils::read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

  detect_group_row <- function() {
    candidates <- c(1, 2)
    score <- function(r) {
      v <- as.character(df[r, -id_col])
      nonnum <- is.na(suppressWarnings(as.numeric(v)))
      as.numeric(sum(nonnum)) + 2*(tolower(df[r, id_col]) %in% c("group","groups"))
    }
    candidates[which.max(vapply(candidates, score, numeric(1)))]
  }

  gr <- if (group_row == "auto") detect_group_row() else as.integer(group_row)
  groups <- factor(trimws(as.character(df[gr, -id_col])))

  dat <- df[-seq_len(gr), , drop = FALSE]
  ids <- make.unique(as.character(dat[[id_col]]))
  dat <- dat[, -id_col, drop = FALSE]

  m <- apply(dat, 2, function(x) as.numeric(as.character(x)))
  rownames(m) <- ids
  m[!is.finite(m)] <- NA_real_

  keep_cols <- colSums(is.na(m)) < nrow(m)
  m <- m[, keep_cols, drop = FALSE]
  samples <- colnames(m)
  groups  <- droplevels(groups[keep_cols])

  res <- list(expr = m, groups = groups, samples = samples)
  attr(res, "group_row") <- gr
  res
}

# ---- Differential expression ----

#' Differential expression by Welch t-test with BH correction
#'
#' Compares a reference group against all other samples pooled together,
#' using a per-protein Welch t-test on the (assumed log2-scale) expression
#' values, with Benjamini-Hochberg adjustment across proteins.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param ref Reference group label (samples not in \code{ref} are pooled as
#'   the comparison group).
#' @param aliases Named character vector mapping alias labels to their
#'   canonical group name (e.g. \code{c(NB = "Normal")}), applied to \code{ref}
#'   before matching against \code{levels(groups)}.
#' @return A tibble with columns \code{Gene}, \code{log2FC}, \code{pvalue},
#'   \code{padj}, sorted by adjusted p-value.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' de <- diff_expr_ttest(io$expr, io$groups, ref = levels(io$groups)[1])
#' head(de)
#' }
#' @export
diff_expr_ttest <- function(expr, groups, ref, aliases = c(NB = "Normal", Control = "Normal")) {
  stopifnot(is.matrix(expr), length(groups) == ncol(expr))
  groups <- droplevels(groups)
  ref <- trimws(ref)
  if (ref %in% names(aliases)) ref <- aliases[[ref]]
  if (!ref %in% levels(groups)) {
    stop("Reference '", ref, "' not found. Available: ", paste(levels(groups), collapse = ", "))
  }
  grp_ref <- groups == ref
  other   <- !grp_ref
  res <- apply(expr, 1, function(x) {
    if (sum(grp_ref) < 2 || sum(other) < 2) return(c(log2FC = NA, pvalue = NA))
    tt <- stats::t.test(x[other], x[grp_ref], var.equal = FALSE)
    c(
      log2FC = stats::median(x[other], na.rm = TRUE) - stats::median(x[grp_ref], na.rm = TRUE),
      pvalue = tt$p.value
    )
  })
  res <- as.data.frame(t(res))
  res$Gene <- rownames(res)
  res$padj <- p.adjust(res$pvalue, "BH")
  tibble::as_tibble(res[order(res$padj), c("Gene","log2FC","pvalue","padj")])
}

#' Differential expression: reference vs. one specific target group
#'
#' Subsets \code{expr}/\code{groups} to just the \code{ref} and \code{target}
#' groups, then runs \code{\link{diff_expr_ttest}}.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param ref Reference group label.
#' @param target Target group label to compare against \code{ref}.
#' @return A tibble with columns \code{Gene}, \code{log2FC}, \code{pvalue},
#'   \code{padj}, as returned by \code{\link{diff_expr_ttest}}.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' levs <- levels(io$groups)
#' de <- diff_expr_vs(io$expr, io$groups, ref = levs[1], target = levs[2])
#' }
#' @export
diff_expr_vs <- function(expr, groups, ref, target) {
  stopifnot(is.matrix(expr), length(groups) == ncol(expr))
  groups <- droplevels(groups)
  if (!ref %in% levels(groups))  stop("ref not found in groups")
  if (!target %in% levels(groups)) stop("target not found in groups")
  idx <- groups %in% c(ref, target)
  expr2   <- expr[, idx, drop = FALSE]
  groups2 <- droplevels(groups[idx])
  diff_expr_ttest(expr2, groups2, ref = ref)
}

#' Pairwise differential expression across groups
#'
#' If \code{ref} is provided, compares \code{ref} against each other group in
#' turn. If \code{ref} is \code{NULL}, compares every pair of groups.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param ref Optional reference group label. If \code{NULL} (default), all
#'   pairwise combinations of groups are compared.
#' @param global_adjust If \code{TRUE}, also compute a \code{padj_global}
#'   column with BH adjustment applied across all comparisons combined.
#'   Default \code{FALSE}.
#' @return A tibble with columns \code{comparison}, \code{ref}, \code{target},
#'   \code{Gene}, \code{log2FC}, \code{pvalue}, \code{padj} (and
#'   \code{padj_global} if \code{global_adjust = TRUE}), stacked across all
#'   comparisons.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' res <- diff_expr_pairwise(io$expr, io$groups, ref = levels(io$groups)[1])
#' }
#' @export
diff_expr_pairwise <- function(expr, groups, ref = NULL, global_adjust = FALSE) {
  stopifnot(is.matrix(expr), length(groups) == ncol(expr))
  groups <- droplevels(groups)
  levs <- levels(groups)
  build_pairs <- function() {
    if (!is.null(ref)) {
      if (!ref %in% levs) stop("ref not found in groups")
      lapply(setdiff(levs, ref), function(tg) c(ref = ref, target = tg))
    } else {
      combn(levs, 2, simplify = FALSE) |>
        lapply(function(ab) c(ref = ab[1], target = ab[2]))
    }
  }
  pairs <- build_pairs()
  out <- lapply(pairs, function(pr) {
    de <- diff_expr_vs(expr, groups, ref = pr[["ref"]], target = pr[["target"]])
    de$comparison <- paste0(pr[["target"]], "_vs_", pr[["ref"]])
    de$ref <- pr[["ref"]]; de$target <- pr[["target"]]
    de[, c("comparison","ref","target","Gene","log2FC","pvalue","padj")]
  })
  res <- dplyr::bind_rows(out)
  if (isTRUE(global_adjust)) res$padj_global <- p.adjust(res$pvalue, "BH")
  res
}

# ---- Plots ----

#' Volcano plot of differential expression results
#'
#' Plots log2 fold-change against adjusted p-value. Uses
#' \pkg{EnhancedVolcano} if installed, otherwise falls back to a plain
#' \pkg{ggplot2} volcano plot.
#'
#' @param de A differential expression tibble as returned by
#'   \code{\link{diff_expr_ttest}} or \code{\link{diff_expr_vs}}, with columns
#'   \code{Gene}, \code{log2FC}, \code{padj}.
#' @param p_cut Adjusted p-value significance cutoff, used for the plotted
#'   threshold line(s). Default \code{0.05}.
#' @param fc_cut Absolute log2 fold-change cutoff, used for the plotted
#'   threshold line(s). Default \code{1}.
#' @return A ggplot object (or an \pkg{EnhancedVolcano} plot object).
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' de <- diff_expr_ttest(io$expr, io$groups, ref = levels(io$groups)[1])
#' plot_volcano(de)
#' }
#' @export
plot_volcano <- function(de, p_cut = 0.05, fc_cut = 1) {
  if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    EnhancedVolcano::EnhancedVolcano(
      de, lab = de$Gene, x = "log2FC", y = "padj",
      pCutoff = p_cut, FCcutoff = fc_cut, legendPosition = "right"
    )
  } else {
    ggplot2::ggplot(de, ggplot2::aes(x = log2FC, y = -log10(padj))) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = 2) +
      ggplot2::geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Volcano plot", y = "-log10(adj p)", x = "log2 fold-change")
  }
}

# ---- Reactome ----

#' Reactome GSEA from a differential expression table
#'
#' Maps gene symbols to Entrez IDs and runs Reactome gene set enrichment
#' analysis (GSEA) ranked by log2 fold-change. Requires the Bioconductor
#' packages \pkg{ReactomePA}, \pkg{clusterProfiler}, and \pkg{org.Hs.eg.db}
#' (install via \code{BiocManager::install()}); if any are missing, a message
#' is printed and \code{NULL} is returned.
#'
#' @param de A differential expression tibble with columns \code{Gene} and
#'   \code{log2FC}, as returned by \code{\link{diff_expr_ttest}}.
#' @param organism Organism name passed to
#'   \code{ReactomePA::gsePathway()}. Default \code{"human"}.
#' @return A GSEA result object from \code{ReactomePA::gsePathway()}, or
#'   \code{NULL} if the required Bioconductor packages are not installed.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' de <- diff_expr_ttest(io$expr, io$groups, ref = levels(io$groups)[1])
#' gsea <- reactome_gsea_from_de(de)
#' }
#' @export
reactome_gsea_from_de <- function(de, organism = "human") {
  ok <- requireNamespace("ReactomePA", quietly = TRUE) &&
        requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) {
    message("Install ReactomePA, clusterProfiler, org.Hs.eg.db for GSEA.")
    return(NULL)
  }
  r <- de$log2FC; names(r) <- de$Gene
  mapped <- clusterProfiler::bitr(names(r),
                                  fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  m <- merge(de[, c("Gene","log2FC")], mapped, by.x = "Gene", by.y = "SYMBOL")
  ranking <- m$log2FC; names(ranking) <- m$ENTREZID
  ranking <- sort(ranking, decreasing = TRUE)
  ReactomePA::gsePathway(ranking, organism = organism,
                         pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         minGSSize = 10, maxGSSize = 500, verbose = FALSE)
}

#' Summarise a folder of Reactome GSEA result CSVs
#'
#' Reads every CSV in \code{dir} matching \code{pattern} (each expected to be
#' a per-comparison Reactome GSEA export with \code{Description}, \code{NES},
#' and \code{p.adjust} columns), combines them, and builds two dot-plots: the
#' top/bottom pathways by average NES, and the most variable pathways across
#' comparisons.
#'
#' @param dir Directory to search for GSEA CSV files. Default \code{"."}.
#' @param pattern Regular expression matching GSEA CSV filenames. Default
#'   \code{"^ReactomeGSEA_.*\\.csv$"}.
#' @return A list with components:
#'   \item{data}{Combined tibble of all GSEA results read.}
#'   \item{topbottom_plot}{ggplot dot-plot of top/bottom pathways by average NES.}
#'   \item{mostvar_plot}{ggplot dot-plot of the most variable pathways across comparisons.}
#' @examples
#' \dontrun{
#' res <- reactome_across("path/to/gsea_csvs")
#' res$topbottom_plot
#' }
#' @export
reactome_across <- function(dir = ".", pattern = "^ReactomeGSEA_.*\\.csv$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No GSEA CSV files found in: ", dir)
  required_cols <- c("Description", "NES", "p.adjust")
  dfs <- lapply(files, function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(required_cols %in% names(df))) return(NULL)
    tt <- sub("^.*?_(.*)$", "\\1", tools::file_path_sans_ext(basename(f)))
    tibble::tibble(Description = as.character(df$Description),
                   NES = as.numeric(df$NES),
                   p.adjust = as.numeric(df$p.adjust),
                   TumorType = tt)
  })
  x <- dplyr::bind_rows(Filter(Negate(is.null), dfs))
  if (!nrow(x)) stop("No valid GSEA tables found.")

  path_avg <- x |>
    dplyr::group_by(Description) |>
    dplyr::summarise(avg_NES = mean(NES, na.rm = TRUE), .groups = "drop")
  top20 <- head(path_avg[order(-path_avg$avg_NES), ]$Description, 20)
  bot20 <- head(path_avg[order(path_avg$avg_NES), ]$Description, 20)
  sel <- c(top20, bot20)
  plot_df <- x |>
    dplyr::filter(Description %in% sel) |>
    dplyr::mutate(NegLogP = -log10(p.adjust + 1e-10),
                  Description = factor(Description, levels = rev(sel)))
  p1 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(TumorType, Description, size = NegLogP, color = NES)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(name = "-log10(p.adj)") +
    ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() + ggplot2::labs(title = "Top & bottom Reactome pathways")

  path_var <- x |>
    dplyr::group_by(Description) |>
    dplyr::summarise(sd_NES = stats::sd(NES, na.rm = TRUE), .groups = "drop")
  top_var <- head(path_var[order(-path_var$sd_NES), ]$Description, 40)
  plot_df2 <- x |>
    dplyr::filter(Description %in% top_var) |>
    dplyr::mutate(NegLogP = -log10(p.adjust + 1e-10),
                  Description = factor(Description, levels = rev(top_var)))
  p2 <- ggplot2::ggplot(plot_df2,
                        ggplot2::aes(TumorType, Description, size = NegLogP, color = NES)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(name = "-log10(p.adj)") +
    ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() + ggplot2::labs(title = "Most variable Reactome pathways")

  list(data = x, topbottom_plot = p1, mostvar_plot = p2)
}

# ---- Dimensionality reduction & ML ----

#' Principal component analysis of a proteomics matrix
#'
#' Runs \code{stats::prcomp()} on samples (rows = samples) and returns the
#' first two components plotted by group.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Optional factor of group labels, length \code{ncol(expr)},
#'   used to color the plot. Default \code{NULL}.
#' @param scale. Passed to \code{stats::prcomp()}; whether to scale variables
#'   to unit variance before analysis. Default \code{TRUE}.
#' @return A list with components:
#'   \item{scores}{Data frame of PC1/PC2 scores per sample.}
#'   \item{plot}{ggplot scatter plot of PC1 vs PC2.}
#'   \item{pca}{The full \code{prcomp} object.}
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' pca <- run_pca(io$expr, io$groups)
#' pca$plot
#' }
#' @export
run_pca <- function(expr, groups = NULL, scale. = TRUE) {
  pca <- stats::prcomp(t(expr), scale. = scale.)
  df  <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  df$Sample <- rownames(df)
  if (!is.null(groups)) df$Group <- groups
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = Group, label = Sample)) +
    ggplot2::geom_point(size = 3) + ggplot2::theme_minimal()
  list(scores = df, plot = p, pca = pca)
}

#' t-SNE embedding of a proteomics matrix
#'
#' Runs \code{Rtsne::Rtsne()} on samples and returns a 2D embedding plotted by
#' group. Requires the \pkg{Rtsne} package.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Optional factor (or vector) of group labels, length
#'   \code{ncol(expr)}, used to color the plot. Default \code{NULL}.
#' @param perplexity Perplexity parameter passed to \code{Rtsne::Rtsne()}.
#'   Default \code{10}.
#' @param seed Random seed for reproducibility. Default \code{42L}.
#' @return A list with components:
#'   \item{coords}{Data frame of tSNE1/tSNE2 coordinates per sample.}
#'   \item{plot}{ggplot scatter plot of the embedding.}
#'   \item{tsne}{The full \code{Rtsne} result object.}
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' ts <- run_tsne(io$expr, io$groups, perplexity = 1)
#' ts$plot
#' }
#' @export
run_tsne <- function(expr, groups = NULL, perplexity = 10, seed = 42L) {
  stopifnot(requireNamespace("Rtsne", quietly = TRUE))
  set.seed(seed)
  ts <- Rtsne::Rtsne(t(expr), dims = 2, perplexity = perplexity)
  df <- data.frame(tSNE1 = ts$Y[,1], tSNE2 = ts$Y[,2],
                   Sample = colnames(expr),
                   Group  = if (!is.null(groups)) as.character(groups) else NA)
  p <- ggplot2::ggplot(df, ggplot2::aes(tSNE1, tSNE2, color = Group)) +
    ggplot2::geom_point(size = 3) + ggplot2::theme_minimal()
  list(coords = df, plot = p, tsne = ts)
}

#' PLS-DA of a proteomics matrix
#'
#' Runs partial least squares discriminant analysis via
#' \code{pls::plsr()} with a one-hot encoded group response. Requires the
#' \pkg{pls} package.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param ncomp Number of PLS components to fit. Default \code{2}.
#' @return A list with components:
#'   \item{scores}{Data frame of Comp1/Comp2 scores per sample.}
#'   \item{plot}{ggplot scatter plot of the first two components.}
#'   \item{model}{The full \code{plsr} model object.}
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' pls_res <- run_plsda(io$expr, io$groups)
#' pls_res$plot
#' }
#' @export
run_plsda <- function(expr, groups, ncomp = 2) {
  stopifnot(requireNamespace("pls", quietly = TRUE))
  Y <- stats::model.matrix(~ 0 + groups)
  colnames(Y) <- levels(groups)
  fit <- pls::plsr(Y ~ t(expr), ncomp = ncomp, scale = TRUE, method = "simpls")
  scores <- as.data.frame(fit$scores[, 1:2, drop = FALSE]); names(scores) <- c("Comp1","Comp2")
  scores$Group <- groups; scores$Sample <- colnames(expr)
  p <- ggplot2::ggplot(scores, ggplot2::aes(Comp1, Comp2, color = Group)) +
    ggplot2::geom_point(size = 3) + ggplot2::theme_minimal()
  list(scores = scores, plot = p, model = fit)
}

#' Train a multinomial LASSO classifier
#'
#' Fits a cross-validated multinomial LASSO (or elastic net) model via
#' \code{glmnet::cv.glmnet()} at \code{lambda.1se}, and returns the top
#' features per class plus a confusion matrix plot. Requires the
#' \pkg{glmnet} package.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param alpha Elastic net mixing parameter passed to \code{glmnet}
#'   (\code{1} = LASSO, \code{0} = ridge). Default \code{1}.
#' @param seed Random seed for reproducibility. Default \code{42L}.
#' @return A list with components:
#'   \item{model}{The fitted \code{glmnet} model at \code{lambda.1se}.}
#'   \item{best_lambda}{The selected \code{lambda.1se} value.}
#'   \item{top_features}{Data frame of top 10 features per class by
#'     coefficient magnitude.}
#'   \item{confusion_plot}{ggplot heatmap of the training confusion matrix.}
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' fit <- train_lasso(io$expr, io$groups)
#' fit$confusion_plot
#' }
#' @export
train_lasso <- function(expr, groups, alpha = 1, seed = 42L) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))
  x <- t(expr); y <- droplevels(groups)
  set.seed(seed)
  cv <- glmnet::cv.glmnet(x, y, family = "multinomial", alpha = alpha)
  lam <- cv$lambda.1se
  fit <- glmnet::glmnet(x, y, family = "multinomial", alpha = alpha, lambda = lam)
  coefs <- glmnet::coef.glmnet(fit)
  if (is.list(coefs)) {
    top_list <- lapply(names(coefs), function(cl) {
      m <- as.matrix(coefs[[cl]])
      m <- m[rownames(m) != "(Intercept)", , drop = FALSE]
      m <- m[order(abs(m[,1]), decreasing = TRUE), , drop = FALSE]
      head(data.frame(Feature = rownames(m), Coefficient = m[,1], TumorType = cl, row.names = NULL), 10)
    })
    top_df <- dplyr::bind_rows(top_list)
  } else {
    top_df <- tibble::tibble()
  }
  pred <- predict(cv, x, s = lam, type = "class")
  cm <- table(Predicted = as.vector(pred), Actual = y)
  cm_df <- as.data.frame(cm)
  p_cm <- ggplot2::ggplot(cm_df, ggplot2::aes(Actual, Predicted, fill = Freq)) +
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = Freq)) + ggplot2::theme_minimal()
  list(model = fit, best_lambda = lam, top_features = top_df, confusion_plot = p_cm)
}

#' Train a Random Forest classifier
#'
#' Fits a \code{randomForest::randomForest()} classifier and returns variable
#' importance plus a confusion matrix plot. Requires the \pkg{randomForest}
#' package.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param groups Factor of group labels, length \code{ncol(expr)}.
#' @param ntree Number of trees to grow. Default \code{500}.
#' @param mtry Number of variables sampled at each split. Defaults to
#'   \code{floor(sqrt(ncol(x)))} if \code{NULL}.
#' @param seed Random seed for reproducibility. Default \code{42L}.
#' @return A list with components:
#'   \item{model}{The fitted \code{randomForest} object.}
#'   \item{importance}{Data frame of variable importance.}
#'   \item{confusion_plot}{ggplot heatmap of the training confusion matrix.}
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' fit <- train_rf(io$expr, io$groups)
#' fit$confusion_plot
#' }
#' @export
train_rf <- function(expr, groups, ntree = 500, mtry = NULL, seed = 42L) {
  stopifnot(requireNamespace("randomForest", quietly = TRUE))
  set.seed(seed)
  x <- t(expr); y <- droplevels(groups)
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(ncol(x))))
  rf <- randomForest::randomForest(x, y, ntree = ntree, mtry = mtry, importance = TRUE)
  imp <- as.data.frame(randomForest::importance(rf))
  imp$Feature <- rownames(imp)
  pred <- stats::predict(rf, x)
  cm <- table(Predicted = pred, Actual = y)
  cm_df <- as.data.frame(cm)
  p_cm <- ggplot2::ggplot(cm_df, ggplot2::aes(Actual, Predicted, fill = Freq)) +
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = Freq)) + ggplot2::theme_minimal()
  list(model = rf, importance = imp, confusion_plot = p_cm)
}

# ---- Clustering & enrichment ----

#' Z-score rows of a matrix
#'
#' Internal helper: z-scores each row of a matrix (subtract mean, divide by
#' SD), computed by scaling the transpose then transposing back.
#'
#' @param mat A numeric matrix.
#' @return A numeric matrix of the same dimensions as \code{mat}, with rows
#'   z-scored.
#' @keywords internal
zscore_rows <- function(mat) t(scale(t(mat)))

#' Cluster a proteome and detect co-expression modules
#'
#' Selects the most variable proteins, hierarchically clusters samples by
#' z-scored expression, draws a heatmap via \pkg{pheatmap}, and (if
#' \pkg{dynamicTreeCut} is installed) cuts the tree into modules.
#'
#' @param expr Numeric matrix of expression values, proteins x samples.
#' @param topN Number of most-variable proteins to include. Default \code{1000}.
#' @param k_min Minimum module size passed to
#'   \code{dynamicTreeCut::cutreeDynamic()}. Default \code{20}.
#' @return A list with components:
#'   \item{z}{Matrix of z-scored values (samples x features).}
#'   \item{module_df}{Tibble mapping each protein to its detected module.}
#'   A heatmap is also drawn as a side effect via \code{pheatmap::pheatmap()}.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' cl <- cluster_proteome(io$expr, topN = 4, k_min = 2)
#' }
#' @export
cluster_proteome <- function(expr, topN = 1000, k_min = 20) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install 'pheatmap'.")
  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    warning("Install dynamicTreeCut for module detection; returning heatmap only.")
  }
  vars <- apply(expr, 1, stats::var, na.rm = TRUE)
  keep <- order(vars, decreasing = TRUE)[seq_len(min(topN, length(vars)))]
  mat  <- expr[keep, , drop = FALSE]

  z <- zscore_rows(t(mat))  # samples x features
  z[!is.finite(z)] <- 0
  hm_data <- t(z)           # features x samples

  d_rows <- stats::dist(hm_data, method = "euclidean")
  hc_rows <- stats::hclust(d_rows, method = "complete")

  modules <- rep(NA_integer_, nrow(hm_data))
  if (requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    modules <- dynamicTreeCut::cutreeDynamic(dendro = hc_rows,
                                             distM = as.matrix(d_rows),
                                             deepSplit = 2, minClusterSize = k_min)
  }
  module_df <- tibble::tibble(Protein = rownames(hm_data),
                              Module = as.factor(modules))

  pheatmap::pheatmap(hm_data, scale = "row",
                     clustering_distance_rows = d_rows,
                     clustering_method = "complete",
                     show_rownames = FALSE, show_colnames = FALSE)

  list(z = z, module_df = module_df)
}

#' Enrich co-expression modules (GO BP & KEGG)
#'
#' For each module in \code{module_df} with at least 10 mapped genes, runs GO
#' biological process and KEGG enrichment via \pkg{clusterProfiler}. Requires
#' the Bioconductor packages \pkg{clusterProfiler} and \pkg{org.Hs.eg.db}.
#'
#' @param module_df A tibble with columns \code{Protein} and \code{Module},
#'   as returned by \code{\link{cluster_proteome}}.
#' @return A named list of enrichment result objects, with names like
#'   \code{"GO_BP_M1"} and \code{"KEGG_M1"} for each module.
#' @examples
#' \dontrun{
#' csv_path <- system.file("extdata", "example_proteome.csv", package = "ProteomicsML")
#' io <- read_wide_proteomics(csv_path)
#' cl <- cluster_proteome(io$expr, topN = 4, k_min = 2)
#' enr <- enrich_modules(cl$module_df)
#' }
#' @export
enrich_modules <- function(module_df) {
  ok <- requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) stop("Install clusterProfiler and org.Hs.eg.db for module enrichment.")
  gene_map <- clusterProfiler::bitr(module_df$Protein,
                                    fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  mod_entrez <- dplyr::left_join(module_df, gene_map, by = c("Protein" = "SYMBOL")) |>
    dplyr::filter(!is.na(ENTREZID))
  out <- list()
  for (m in stats::na.omit(unique(mod_entrez$Module))) {
    genes <- unique(mod_entrez$ENTREZID[mod_entrez$Module == m])
    if (length(genes) < 10) next
    ego <- clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                     ont = "BP", keyType = "ENTREZID",
                                     pAdjustMethod = "BH", pvalueCutoff = 0.05)
    ekegg <- try(clusterProfiler::enrichKEGG(gene = genes, organism = "hsa",
                                             keyType = "ncbi-geneid", pAdjustMethod = "BH",
                                             pvalueCutoff = 0.05), silent = TRUE)
    out[[paste0("GO_BP_M", m)]]  <- ego
    out[[paste0("KEGG_M", m)]]   <- if (inherits(ekegg, "try-error")) NULL else ekegg
  }
  out
}

# ---- Wizard ----

#' Launch the interactive proteomics analysis wizard
#'
#' Starts an interactive console menu that walks through the package's main
#' workflows (differential expression + volcano, Reactome summaries, PCA /
#' t-SNE / PLS-DA, LASSO + Random Forest, heatmap + module enrichment)
#' without requiring the user to write any R code. Prompts for a CSV path
#' and analysis choices via \code{readline()}/\code{utils::menu()}.
#'
#' @return Invisibly, \code{TRUE}. Called for its interactive side effects.
#' @examples
#' \dontrun{
#' run_proteomics_wizard()
#' }
#' @export
run_proteomics_wizard <- function() {
  cat("\n=== ProteomicsML: Analysis Wizard ===\n")
  repeat {
    choice <- utils::menu(c(
      "Differential expression + volcano",
      "Reactome summary across comparisons (dot-plots)",
      "PCA / t-SNE / PLS-DA",
      "Machine learning: LASSO + Random Forest",
      "Heatmap + module detection + enrichment",
      "Quit"
    ), title = "Select an option")
    if (choice %in% c(0, 6)) break

    if (choice == 1) {
      path <- readline("CSV path (wide format): ")
      io   <- read_wide_proteomics(path, group_row = "auto")
      io$groups <- factor(trimws(io$groups))
      levs <- levels(io$groups)
      cat("Groups detected:", paste(levs, collapse = ", "), "\n")

      ref_i <- utils::menu(levs, title = "Choose the CONTROL/REFERENCE group")
      if (ref_i < 1) { message("Cancelled."); next }
      ref <- levs[ref_i]

      mode <- utils::menu(c(
        "Compare reference vs EACH other group (pairwise)",
        "Compare reference vs ALL OTHERS (pooled)",
        "Compare reference vs ONE chosen group"
      ), title = "Comparison mode")
      if (mode == 1L) {
        res <- diff_expr_pairwise(io$expr, io$groups, ref = ref, global_adjust = FALSE)
        print(dplyr::count(res, comparison))
        first_cmp <- res$comparison[1]
        pdat <- dplyr::filter(res, comparison == first_cmp)
        print(plot_volcano(pdat[, c("Gene","log2FC","pvalue","padj")]))
      } else if (mode == 2L) {
        de <- diff_expr_ttest(io$expr, io$groups, ref = ref)
        print(plot_volcano(de))
      } else if (mode == 3L) {
        targets <- setdiff(levs, ref)
        tgt_i <- utils::menu(targets, title = "Choose TARGET group")
        if (tgt_i < 1) { message("Cancelled."); next }
        target <- targets[tgt_i]
        de <- diff_expr_vs(io$expr, io$groups, ref = ref, target = target)
        print(plot_volcano(de))
      }

    } else if (choice == 2) {
      dir <- readline("Folder with ReactomeGSEA_*.csv: ")
      res <- reactome_across(dir)
      print(res$topbottom_plot); print(res$mostvar_plot)

    } else if (choice == 3) {
      path <- readline("CSV path (wide format): ")
      io   <- read_wide_proteomics(path, group_row = "auto")
      print(run_pca(io$expr, io$groups)$plot)
      if (requireNamespace("Rtsne", quietly = TRUE)) print(run_tsne(io$expr, io$groups)$plot)
      if (requireNamespace("pls", quietly = TRUE))  print(run_plsda(io$expr, io$groups)$plot)

    } else if (choice == 4) {
      path <- readline("CSV path (wide format): ")
      io   <- read_wide_proteomics(path, group_row = "auto")
      print(train_lasso(io$expr, io$groups)$confusion_plot)
      print(train_rf(io$expr, io$groups)$confusion_plot)

    } else if (choice == 5) {
      path <- readline("CSV path (wide format): ")
      io   <- read_wide_proteomics(path, group_row = "auto")
      cl   <- cluster_proteome(io$expr, topN = 1000)
      if (requireNamespace("clusterProfiler", quietly = TRUE)) {
        enr <- enrich_modules(cl$module_df)
        if (length(enr) && requireNamespace("enrichplot", quietly = TRUE)) {
          first <- names(enr)[grepl("^GO_BP", names(enr))][1]
          if (!is.na(first)) print(enrichplot::dotplot(enr[[first]], showCategory = 15))
        }
      }
    }
  }
  invisible(TRUE)
}
