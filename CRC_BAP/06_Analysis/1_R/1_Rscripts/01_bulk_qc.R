# 01_bulk_qc.R
# Bulk expression QC, normalisation, PCA and batch assessment.
# Worker count: configured in run_pipeline.R (one level; BLAS pinned to 1).

run_bulk_qc <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(edgeR)
    library(preprocessCore)
    library(ggplot2)
    library(pheatmap)
  })

  obj <- readRDS(file.path(cfg$paths$processed, "validated_inputs.rds"))
  expr <- obj$expression
  clin <- obj$clinical
  clin <- clin[match(colnames(expr), clin$sample_id), , drop = FALSE]
  group <- factor(clin$group)
  if (length(levels(group)) < 2) group <- factor(rep("all", ncol(expr)))

  is_counts <- all(abs(expr - round(expr)) < 1e-8, na.rm = TRUE) && max(expr, na.rm = TRUE) > 50
  if (is_counts) {
    y <- edgeR::DGEList(counts = expr)
    keep <- edgeR::filterByExpr(y, group = group)
    y <- y[keep, , keep.lib.sizes = FALSE]
    log_expr <- edgeR::cpm(y, log = TRUE, prior.count = 1)
    normalized <- "log2 CPM (edgeR, filterByExpr)"
  } else {
    log_expr <- as.matrix(expr)
    if (max(log_expr, na.rm = TRUE) > 50) log_expr <- log2(log_expr + 1)
    if (min(log_expr, na.rm = TRUE) < 0) log_expr <- log2(log_expr - min(log_expr, na.rm = TRUE) + 1)
    if (anyNA(log_expr)) {
      row_med <- apply(log_expr, 1, median, na.rm = TRUE)
      for (j in seq_len(ncol(log_expr))) log_expr[is.na(log_expr[, j]), j] <- row_med[is.na(log_expr[, j])]
    }
    log_expr <- preprocessCore::normalize.quantiles(log_expr)
    rownames(log_expr) <- rownames(expr)
    colnames(log_expr) <- colnames(expr)
    normalized <- "quantile-normalized log2"
  }

  lib_size <- colSums(expr, na.rm = TRUE)
  detected <- colSums(expr > 0, na.rm = TRUE)
  med <- apply(log_expr, 1, median, na.rm = TRUE)
  cor_med <- apply(log_expr, 2, function(x) stats::cor(x, med, use = "pairwise.complete.obs"))

  pca <- stats::prcomp(t(log_expr), scale. = TRUE)
  pc_n <- min(10, ncol(pca$x))
  pc <- pca$x[, seq_len(pc_n), drop = FALSE]
  md <- tryCatch(
    stats::mahalanobis(pc, colMeans(pc), stats::cov(pc)),
    error = function(e) rowSums(scale(pc)^2)
  )
  outlier <- md > stats::qchisq(0.999, df = pc_n)

  batch <- factor(clin$batch)
  if (nlevels(batch) > 1 && !any(is.na(batch))) {
    log_expr_bc <- limma::removeBatchEffect(log_expr, batch = batch)
    pca_bc <- stats::prcomp(t(log_expr_bc), scale. = TRUE)
    batch_status <- paste("corrected (limma removeBatchEffect;", nlevels(batch), "levels)")
  } else {
    log_expr_bc <- log_expr
    pca_bc <- pca
    batch_status <- "not corrected (one level or missing batch)"
  }

  qc <- data.frame(
    sample_id = colnames(expr),
    group = as.character(clin$group),
    batch = as.character(clin$batch),
    library_size = lib_size,
    detected_features = detected,
    correlation_to_median = cor_med,
    PC1 = pca$x[, 1],
    PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else NA,
    mahalanobis = md,
    outlier_flag = outlier,
    quality_flag = ifelse(outlier, "outlier", "pass"),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(qc, file.path(cfg$paths$qc, "Sample_QC.tsv"), sep = "\t")

  fig_dir <- file.path(cfg$paths$qc, "Figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = clin$group, batch = clin$batch)
  save_plot(function() {
    print(ggplot(pca_df, aes(PC1, PC2, colour = group, shape = batch)) + geom_point(size = 2, alpha = 0.8) + theme_bw() + labs(title = "PCA (uncorrected)"))
  }, file.path(fig_dir, "01_pca_uncorrected"), 7, 5)

  pca_bc_df <- data.frame(PC1 = pca_bc$x[, 1], PC2 = pca_bc$x[, 2], group = clin$group, batch = clin$batch)
  save_plot(function() {
    print(ggplot(pca_bc_df, aes(PC1, PC2, colour = group, shape = batch)) + geom_point(size = 2, alpha = 0.8) + theme_bw() + labs(title = "PCA (batch-corrected for display)"))
  }, file.path(fig_dir, "01_pca_batch_corrected"), 7, 5)

  qc_df <- data.frame(sample_id = qc$sample_id, group = qc$group, library_size = qc$library_size, detected = qc$detected_features)
  save_plot(function() {
    print(ggplot(qc_df, aes(group, library_size, fill = group)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) + theme_bw() + labs(title = "Library size by group", y = "Library size"))
  }, file.path(fig_dir, "01_library_size"), 7, 5)

  hv <- top_variable(log_expr, 500)
  ann_col <- data.frame(group = clin$group, batch = clin$batch, row.names = clin$sample_id)
  save_plot(function() {
    pheatmap::pheatmap(hv, show_rownames = FALSE, show_colnames = FALSE, annotation_col = ann_col, main = "Top 500 variable features")
  }, file.path(fig_dir, "01_top_variable_heatmap"), 8, 7)

  outlier_frac <- mean(outlier)
  decision <- if (outlier_frac > cfg$params$max_outlier_fraction) "PASS_WITH_OUTLIERS" else "PASS"
  log_decision("01_bulk_qc", decision, paste("normalization:", normalized, "; outliers:", sum(outlier), "of", length(outlier), "; batch:", batch_status))

  qc_report <- c(
    "# QC Report — CRC Track A",
    "",
    paste("QC date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Script: `06_Analysis/1_R/1_Rscripts/01_bulk_qc.R`"),
    paste("Samples in:", ncol(expr), "→ after QC:", ncol(log_expr), "(no samples dropped at this step)"),
    paste("Features:", nrow(expr), "→", nrow(log_expr)),
    paste("Normalization:", normalized),
    paste("Batch:", batch_status),
    "",
    "## Per-sample summary",
    paste("Table: `05_QC/Sample_QC.tsv`."),
    paste("Outliers flagged:", sum(outlier), "(", round(100 * outlier_frac, 1), "% )"),
    "",
    "## Decision",
    paste("QC checkpoint:", decision),
    paste("Outliers were flagged, not silently removed. Sensitivity analysis will exclude them."),
    "",
    "## Limitations",
    "Purity and tissue composition were not available at QC time; they are handled in the model/sensitivity analysis."
  )
  writeLines(qc_report, file.path(cfg$paths$qc, "QC_Report.md"))

  saveRDS(
    list(log_expr = log_expr, log_expr_bc = log_expr_bc, pca = pca, pca_bc = pca_bc, qc = qc, clinical = clin, normalized = normalized, is_counts = is_counts),
    file.path(cfg$paths$processed, "processed_bulk.rds")
  )
  invisible(TRUE)
}
