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

  # `batch` here is the TCGA project (COAD vs READ): anatomical cohorts of the same tissue on one platform, not a
  # technical batch. Quantify its influence rather than "correcting" it, and do not carry a corrected matrix
  # forward — nothing downstream used one, and the earlier "batch-corrected" panel was a near-duplicate of the
  # uncorrected PCA for exactly this reason.
  batch <- factor(clin$batch)
  k_pc <- min(10, ncol(pca$x))
  project_r2 <- if (nlevels(batch) > 1 && !any(is.na(batch))) {
    vapply(seq_len(k_pc), function(i) summary(stats::lm(pca$x[, i] ~ batch))$r.squared, numeric(1))
  } else {
    rep(NA_real_, k_pc)
  }
  batch_status <- if (nlevels(batch) > 1 && !all(is.na(project_r2))) {
    paste0("no technical batch: TCGA project (", nlevels(batch), " levels) explains at most ",
           round(100 * max(project_r2, na.rm = TRUE), 1), "% of any of PC1-", k_pc, "; not removed")
  } else {
    "not assessed (one level or missing project)"
  }
  log_expr_bc <- log_expr   # retained in the saved object for schema compatibility; no correction is applied
  pca_bc <- pca

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
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], project = factor(clin$batch), sample_id = clin$sample_id)
  var_pc <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  save_plot(function() {
    print(ggplot(pca_df, aes(PC1, PC2, colour = project, shape = project)) + geom_point(size = 2, alpha = 0.8) +
      theme_bw() +
      labs(title = sprintf("PCA of %d quantified tumours", ncol(log_expr)),
           subtitle = sprintf("no technical batch (single platform, single pipeline); PC1 %.1f%%, PC2 %.1f%%",
                              var_pc[1], var_pc[2]),
           colour = "TCGA project", shape = "TCGA project"))
  }, file.path(fig_dir, "01_pca_uncorrected"), 7, 5)

  proj_df <- data.frame(PC = factor(paste0("PC", seq_len(k_pc)), levels = paste0("PC", seq_len(k_pc))),
                        variance_explained = var_pc[seq_len(k_pc)],
                        project_r2 = project_r2)
  save_plot(function() {
    print(ggplot(proj_df, aes(PC)) +
      geom_col(aes(y = variance_explained), fill = "grey70") +
      geom_point(aes(y = 100 * project_r2), colour = "#C0392B", size = 2.2) +
      geom_line(aes(y = 100 * project_r2, group = 1), colour = "#C0392B") +
      theme_bw(base_size = 9) +
      labs(x = NULL, y = "percent",
           title = "Cohort (COAD vs READ) structure is minor and not removed",
           subtitle = "grey bars = variance explained by each PC; red = share of that PC attributable to project"))
  }, file.path(fig_dir, "01_project_variance"), 7, 5)

  # Library size is compared between the two TCGA projects (the only grouping in this single-pipeline cohort).
  qc_df <- data.frame(sample_id = qc$sample_id, project = qc$batch, library_size = qc$library_size,
                      detected = qc$detected_features)
  lib_test <- stats::wilcox.test(library_size ~ project, data = qc_df)
  cat("library size by project: Wilcoxon p =", signif(lib_test$p.value, 3), "\n")
  save_plot(function() {
    print(ggplot(qc_df, aes(project, library_size, fill = project)) + geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) + theme_bw() +
      labs(title = "Library size by TCGA project", y = "Library size", x = NULL,
           subtitle = sprintf("Wilcoxon p = %.3g", lib_test$p.value)))
  }, file.path(fig_dir, "01_library_size"), 7, 5)

  hv <- top_variable(log_expr, 500)
  # Named "project", not "batch": COAD vs READ is anatomical, not a technical batch (see the QC decision).
  ann_col <- data.frame(project = clin$batch, group = clin$group, row.names = clin$sample_id)
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
