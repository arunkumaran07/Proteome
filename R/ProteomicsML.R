# ---- IO ----

#' Read a "wide" proteomics CSV (genes in rows, samples in columns)
#' Assumes col 1 = gene IDs. First data row (or second) holds group labels.
#' @param path CSV path
#' @param id_col column index of gene/protein IDs (default 1)
#' @param group_row "auto", 1 or 2
#' @return list(expr, groups, samples) with attr(..., "group_row")
#' @export
#' @importFrom stats p.adjust
NULL
devtools::document()
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
#' @param expr matrix [genes x samples]
#' @param groups factor of length ncol(expr)
#' @param ref reference group label
#' @param aliases named vector mapping aliases -> canonical (e.g., NB->Normal)
#' @return tibble(Gene, log2FC, pvalue, padj)
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

#' Differential expression: ref vs specific target
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

#' Pairwise differential expression
#' If ref is provided: ref vs each other group; else all pairwise combos
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

#' Volcano plot (EnhancedVolcano if available, else ggplot fallback)
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

#' Reactome GSEA from a DE table
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

#' Summarise a folder of Reactome GSEA CSVs
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

#' PCA
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

#' t-SNE
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

#' PLS-DA (via pls::plsr with one-hot response)
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

#' Train multinomial LASSO (glmnet) with CV
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

#' Train Random Forest
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

#' Z-score rows (helper)
zscore_rows <- function(mat) t(scale(t(mat)))

#' Cluster proteome & detect modules
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

#' Enrich modules (GO BP & KEGG)
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

#' Interactive proteomics analysis menu
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
