#!/usr/bin/env Rscript
# fix_qc_project_effect.R
# Regenerate the two QC PCA panels that were (correctly) reported as visually identical.
#
# What went wrong: 01_bulk_qc.R treated TCGA-COAD vs TCGA-READ as a "batch" and plotted
# limma::removeBatchEffect() output next to the uncorrected PCA. COAD and READ are anatomical cohorts of the same
# tissue on the same platform, not technical batches, so removing the project term barely moves the cloud — the two
# panels were near-duplicates. The corrected matrix was never used downstream either (02-04 read obj$log_expr).
#
# Fix: panel A becomes the honest PCA (colour = TCGA project), and panel B becomes the quantitative diagnostic that
# the "batch" panel was standing in for — the share of each principal component attributable to the project
# (R^2), next to the variance that component explains. If the project explained little, that is now visible as a
# number rather than as two look-alike scatterplots.
#
# Output: 05_QC/Figures/01_pca_uncorrected.{pdf,jpg}   (panel A)
#         05_QC/Figures/01_project_variance.{pdf,jpg}  (panel B)
# Read-only with respect to 04_Processed_Data.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

obj <- readRDS("04_Processed_Data/processed_bulk.rds")
log_expr <- obj$log_expr
clin <- obj$clinical
clin <- clin[match(colnames(log_expr), clin$sample_id), , drop = FALSE]
project <- factor(clin$batch)
cat("samples:", ncol(log_expr), "| project levels:",
    paste(levels(project), table(project), sep = "=", collapse = ", "), "\n")

save_plot <- function(plot_fun, file_base, width = 7, height = 5) {
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height)
  plot_fun()
  grDevices::dev.off()
  grDevices::jpeg(paste0(file_base, ".jpg"), width = width, height = height, units = "in", res = 300, quality = 95)
  plot_fun()
  grDevices::dev.off()
}

fig_dir <- "05_QC/Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pca <- stats::prcomp(t(log_expr), scale. = TRUE)
stopifnot(length(intersect(rownames(pca$x), rownames(obj$pca$x))) == ncol(log_expr))
r_pc <- stats::cor(pca$x[, 1], obj$pca$x[, 1])
cat("PCA reproduces the stored solution: cor(PC1) =", round(r_pc, 6), "\n")
stopifnot(r_pc > 0.999)

## Panel A — PCA coloured by project, shaped by project
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], project = project,
                     sample_id = rownames(pca$x))
var_pc <- 100 * pca$sdev^2 / sum(pca$sdev^2)
save_plot(function() {
  print(ggplot(pca_df, aes(PC1, PC2, colour = project, shape = project)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_bw() +
    labs(title = "PCA of 624 quantified tumours",
         subtitle = sprintf("no technical batch exists (single platform, single pipeline); PC1 %.1f%%, PC2 %.1f%%",
                            var_pc[1], var_pc[2]),
         colour = "TCGA project", shape = "TCGA project"))
}, file.path(fig_dir, "01_pca_uncorrected"), 7, 5)

## Panel B — project structure per principal component
k <- min(10, ncol(pca$x))
proj <- data.table(
  PC = factor(paste0("PC", seq_len(k)), levels = paste0("PC", seq_len(k))),
  variance_explained = var_pc[seq_len(k)],
  project_r2 = vapply(seq_len(k), function(i) {
    fit <- stats::lm(pca$x[, i] ~ project)
    summary(fit)$r.squared
  }, numeric(1))
)
data.table::fwrite(proj, "05_QC/Project_effect_by_PC.tsv", sep = "\t")
cat("\nproject R^2 per PC:\n"); print(proj, digits = 3)
cat("max project R^2 over PC1-10:", round(max(proj$project_r2), 3), "\n")

save_plot(function() {
  print(ggplot(proj, aes(PC)) +
    geom_col(aes(y = variance_explained), fill = "grey70") +
    geom_point(aes(y = 100 * project_r2), colour = "#C0392B", size = 2.2) +
    geom_line(aes(y = 100 * project_r2, group = 1), colour = "#C0392B") +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "percent",
         title = "Cohort (COAD vs READ) structure is minor and not removed",
         subtitle = "grey bars = variance explained by each PC; red = share of that PC attributable to project"))
}, file.path(fig_dir, "01_project_variance"), 7, 5)

cat("\nwrote 05_QC/Figures/01_pca_uncorrected.* and 01_project_variance.*\n")
