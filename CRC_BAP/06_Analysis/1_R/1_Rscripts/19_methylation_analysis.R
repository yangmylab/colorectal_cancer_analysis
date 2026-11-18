#!/usr/bin/env Rscript
# 19_methylation_analysis.R
# Cross-platform test of the frozen K=2 ecotypes at the DNA-methylation level (TCGA-COAD/READ 450K).
#
# Inputs (built by 16_fetch_tcga_methylation.py + 19_build_gene_methylation.py):
#   03_Raw_Data/TCGA/methylation/gene_beta.tsv.gz   gene-level promoter beta (annotation based)
#   03_Raw_Data/TCGA/methylation/global_beta.tsv    mean beta over all probes per sample
#   06_Analysis/03_Clustering/ecotype_labels.tsv, 03_Raw_Data/clinical.tsv
# Output: 06_Analysis/08_Methylation/*.tsv
#         09_Figures/plots/plot_20 .. plot_24

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(survival)
  library(mclust)
})

METH <- "03_Raw_Data/TCGA/methylation"
OUT <- "06_Analysis/08_Methylation"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source("06_Analysis/1_R/1_Rscripts/utils.R")

save_plot <- function(plot_fun, file_base, width = 7, height = 5) {
  # The composites are laid out at a fixed printed width, so a panel's authored size decides its printed
  # type size: a panel authored at 3.4 in prints its 11 pt text at >=7 pt in a two-column figure.
  # CRC_PANEL_MAX_W caps the authored width (inches; 0 = off) and CRC_PANEL_SCALE multiplies both
  # dimensions, so a re-authored panel set needs no edits to the individual plot calls.
  scale <- as.numeric(Sys.getenv("CRC_PANEL_SCALE", "1"))
  max_w <- as.numeric(Sys.getenv("CRC_PANEL_MAX_W", "0"))
  width <- width * scale
  height <- height * scale
  if (max_w > 0 && width > max_w) {
    shrink <- max_w / width
    width <- width * shrink
    height <- height * shrink
  }
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height)
  plot_fun()
  grDevices::dev.off()
  grDevices::jpeg(paste0(file_base, ".jpg"), width = width, height = height, units = "in", res = 300, quality = 95)
  plot_fun()
  grDevices::dev.off()
}

## 1. Gene-level beta matrix --------------------------------------------------
gene_dt <- fread(file.path(METH, "gene_beta.tsv.gz"), data.table = FALSE, check.names = FALSE)
gene_beta <- as.matrix(gene_dt[, -1, drop = FALSE])
storage.mode(gene_beta) <- "numeric"
rownames(gene_beta) <- gene_dt[[1]]
cat("gene-level beta:", nrow(gene_beta), "genes x", ncol(gene_beta), "samples\n")

glob <- fread(file.path(METH, "global_beta.tsv"), data.table = FALSE)

## 2. Sample annotation ------------------------------------------------------
labels <- fread("06_Analysis/03_Clustering/ecotype_labels.tsv", data.table = FALSE)
clin <- fread("03_Raw_Data/clinical.tsv", data.table = FALSE)
clin$os_time <- clin$os_time_days
meta <- data.frame(sample_id = colnames(gene_beta), stringsAsFactors = FALSE)
meta$ecotype <- labels$ecotype[match(meta$sample_id, labels$sample_id)]
meta <- merge(meta, clin[, c("sample_id", "stage", "age", "sex", "sidedness", "os_time", "os_event")],
              by = "sample_id", all.x = TRUE, sort = FALSE)
meta <- meta[match(colnames(gene_beta), meta$sample_id), , drop = FALSE]
meta$global_beta <- glob$global_beta[match(meta$sample_id, glob$sample_id)]
ok <- !is.na(meta$ecotype)
cat("samples with an RNA ecotype label:", sum(ok), "of", nrow(meta), "\n")

## 3. CIMP score -------------------------------------------------------------
cimp_genes <- c("CACNA1G", "CDKN2A", "CRABP1", "IGF2", "MLH1", "NEUROG1", "RUNX3", "SOCS1")
cimp_present <- intersect(cimp_genes, rownames(gene_beta))
cimp_beta <- gene_beta[cimp_present, , drop = FALSE]
cimp_count <- colSums(cimp_beta > 0.30, na.rm = TRUE)
meta$cimp_genes_methylated <- as.integer(cimp_count)
meta$cimp_status <- ifelse(cimp_count >= 3, "CIMP-high", "CIMP-low")
cat("CIMP panel genes available:", paste(cimp_present, collapse = ", "), "\n")

## 4. Differential methylation between ecotypes (limma, gene level) ----------
mm <- gene_beta
mm[mm < 0.001] <- 0.001
mm[mm > 0.999] <- 0.999
mval <- log2(mm / (1 - mm))
design <- model.matrix(~ 0 + factor(meta$ecotype[ok], levels = c("E1", "E2")))
colnames(design) <- c("E1", "E2")
contr <- limma::makeContrasts(E2 - E1, levels = design)
fit <- eBayes(lmFit(mval[, ok, drop = FALSE], design))
fit <- eBayes(contrasts.fit(lmFit(mval[, ok, drop = FALSE], design), contr))
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
tt$gene <- rownames(tt)
fwrite(tt, file.path(OUT, "differential_methylation_ecotype.tsv"), sep = "\t")
n_sig <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
cat("differentially methylated genes (FDR<0.05, |logFC|>0.5):", n_sig, "of", nrow(tt), "\n")

## 5. Unsupervised methylation clustering ------------------------------------
v <- apply(gene_beta, 1, stats::var, na.rm = TRUE)
top <- gene_beta[order(-v)[seq_len(min(5000, length(v)))], , drop = FALSE]
hc <- hclust(dist(scale(t(top))), method = "ward.D2")
meta$methylation_cluster <- paste0("M", cutree(hc, k = 2))
ari_vs_ecotype <- mclust::adjustedRandIndex(meta$ecotype[ok], meta$methylation_cluster[ok])
cat("ARI (RNA ecotype vs methylation cluster):", round(ari_vs_ecotype, 3), "\n")

## 6. Survival by methylation cluster ----------------------------------------
surv_dat <- meta[!is.na(meta$os_time) & !is.na(meta$os_event) & meta$os_time > 0, , drop = FALSE]
surv_dat$os_time <- as.numeric(surv_dat$os_time)
surv_dat$os_event <- as.integer(surv_dat$os_event)
surv_dat$stage_num <- as.numeric(factor(surv_dat$stage, levels = c("I", "II", "III", "IV")))
km <- NULL
if (nrow(surv_dat) > 50 && sum(surv_dat$os_event) > 20) {
  surv_dat$methylation_cluster <- factor(surv_dat$methylation_cluster)
  f1 <- survival::coxph(survival::Surv(os_time, os_event) ~ methylation_cluster, data = surv_dat)
  f2 <- survival::coxph(survival::Surv(os_time, os_event) ~ methylation_cluster + age + sex + stage_num, data = surv_dat)
  s2 <- summary(f2)
  surv_row <- data.frame(
    analysis = "methylation cluster adjusted (age, sex, stage)",
    n = nrow(surv_dat),
    events = sum(surv_dat$os_event),
    c_index_base = as.numeric(survival::concordance(f1)$concordance),
    c_index_full = as.numeric(survival::concordance(f2)$concordance),
    hr = unname(s2$conf.int[1, 1]),
    lower = unname(s2$conf.int[1, 3]),
    upper = unname(s2$conf.int[1, 4]),
    p = unname(s2$coefficients[1, 5])
  )
  fwrite(surv_row, file.path(OUT, "methylation_survival.tsv"), sep = "\t")
  km <- survival::survfit(survival::Surv(os_time, os_event) ~ methylation_cluster, data = surv_dat)
  print(surv_row)
}

## 7. Summary metrics --------------------------------------------------------
meta_plot <- meta[ok, , drop = FALSE]
wilcox_global <- stats::wilcox.test(global_beta ~ ecotype, data = meta_plot)$p.value
cimp_frac <- round(prop.table(table(meta_plot$ecotype, meta_plot$cimp_status), 1), 3)
summ <- data.frame(
  metric = c("genes_gene_level", "samples", "samples_with_ecotype", "cimp_genes_available",
             "cimp_high", "cimp_low", "dmps_fdr05_lfc05", "ari_rna_vs_methylation",
             "mean_global_beta_E1", "mean_global_beta_E2", "global_beta_wilcox_p"),
  value = c(nrow(gene_beta), ncol(gene_beta), sum(ok), length(cimp_present),
            sum(meta$cimp_status == "CIMP-high"), sum(meta$cimp_status == "CIMP-low"),
            n_sig, round(ari_vs_ecotype, 4),
            round(mean(meta$global_beta[meta$ecotype == "E1"], na.rm = TRUE), 4),
            round(mean(meta$global_beta[meta$ecotype == "E2"], na.rm = TRUE), 4),
            signif(wilcox_global, 3))
)
fwrite(summ, file.path(OUT, "methylation_summary.tsv"), sep = "\t")
fwrite(meta, file.path(OUT, "methylation_sample_scores.tsv"), sep = "\t")
print(summ)
print(cimp_frac)

## 8. Figures -----------------------------------------------------------------
save_plot(function() {
  print(ggplot(meta_plot, aes(ecotype, cimp_genes_methylated, fill = ecotype)) +
    geom_boxplot(outlier.size = 0.8) + geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    scale_fill_manual(values = ECOTYPE_COL) +
    theme_bw() +
    labs(x = "RNA ecotype", y = "CIMP panel genes methylated (of 8)",
         title = "CIMP score by ecotype (TCGA-COAD/READ 450K)"))
}, file.path(FIG, "plot_20"), 6, 5)

save_plot(function() {
  print(ggplot(meta_plot, aes(ecotype, global_beta, fill = ecotype)) +
    geom_boxplot(outlier.size = 0.8) +
    scale_fill_manual(values = ECOTYPE_COL) +
    theme_bw() +
    labs(x = "RNA ecotype", y = "Mean beta over all probes",
         title = sprintf("Global methylation by ecotype (Wilcoxon p=%.2f)", wilcox_global)))
}, file.path(FIG, "plot_21"), 6, 5)

top_genes <- head(intersect(tt$gene[order(tt$P.Value)], rownames(gene_beta)), 40)
hm <- gene_beta[top_genes, ok, drop = FALSE]
ord <- order(meta_plot$ecotype)
save_plot(function() {
  pheatmap::pheatmap(
    hm[, ord, drop = FALSE],
    cluster_cols = FALSE,
    show_colnames = FALSE,
    scale = "row",
    annotation_col = data.frame(ecotype = meta_plot$ecotype[ord], row.names = meta_plot$sample_id[ord]),
    main = "Top differentially methylated genes (E2 vs E1)",
    fontsize_row = 7
  )
}, file.path(FIG, "plot_22"), 9, 8)

if (!is.null(km)) {
  save_plot(function() {
    plot(km, col = unname(ECOTYPE_COL), lwd = 2, xlab = "Days",
         ylab = "Overall survival probability", mark.time = TRUE,
         main = "Methylation cluster and OS")
    legend("topright", legend = levels(surv_dat$methylation_cluster), col = unname(ECOTYPE_COL), lwd = 2, bty = "n")
  }, file.path(FIG, "plot_23"), 7, 5)
}

tab <- as.data.frame(table(ecotype = meta_plot$ecotype, cluster = meta_plot$methylation_cluster))
save_plot(function() {
  print(ggplot(tab, aes(ecotype, cluster, fill = Freq)) +
    geom_tile() + geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "#2E86C1") +
    theme_bw() +
    labs(title = sprintf("RNA ecotype vs methylation cluster (ARI %.2f)", ari_vs_ecotype),
         x = "RNA ecotype", y = "Methylation cluster"))
}, file.path(FIG, "plot_24"), 6, 5)

cat("done\n")
