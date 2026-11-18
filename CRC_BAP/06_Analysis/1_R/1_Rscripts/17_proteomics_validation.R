#!/usr/bin/env Rscript
# 17_proteomics_validation.R
# Cross-platform (protein-level) test of the frozen K=2 ecotypes in CPTAC-COAD.
#
# Input : 03_Raw_Data/CPTAC_COAD/ (cBioPortal datahub: coad_cptac_2019)
# Output: 07_Validation/proteomics_validation.tsv
#         07_Validation/Proteomics_Report.md
#         09_Figures/plots/plot_16 .. plot_19 (+ pdf/jpg)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(pheatmap)
  library(mclust)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

PROC <- "04_Processed_Data"
IN <- "03_Raw_Data/CPTAC_COAD"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

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

## 1. Hallmark gene sets (symbols) -------------------------------------------
# msigdbr cannot fetch offline here, so this falls back to the frozen discovery GMT mapped to symbols.
hl <- hallmark_gene_sets("symbol")
gene_sets <- hl$sets
cat("programme set source:", hl$source, "\n")

## 2. Input matrices ---------------------------------------------------------
read_cbio <- function(path, gene_col, value_transform = FALSE) {
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  mat <- as.matrix(dt[, setdiff(colnames(dt), gene_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  ids <- as.character(dt[[gene_col]])
  # CPTAC protein tables use "SYMBOL|SYMBOL" identifiers
  if (gene_col == "Composite.Element.REF") ids <- sub("\\|.*$", "", ids)
  rownames(mat) <- ids
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  if (value_transform) mat <- log2(mat + 1)
  # collapse duplicated symbols by the row with the largest variance
  if (any(duplicated(rownames(mat)))) {
    v <- apply(mat, 1, function(x) stats::var(x, na.rm = TRUE))
    ord <- order(rownames(mat), -v)
    mat <- mat[ord, , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  mat
}

prot <- read_cbio(file.path(IN, "data_protein_quantification.txt"), "Composite.Element.REF")
mrna <- read_cbio(file.path(IN, "data_mrna_seq_v2_rsem.txt"), "Hugo_Symbol", value_transform = TRUE)
cat("protein:", nrow(prot), "x", ncol(prot), "| mRNA:", nrow(mrna), "x", ncol(mrna), "\n")

score_sets <- function(mat, min_size = 5) {
  z <- t(scale(t(mat)))
  sets <- lapply(gene_sets, function(g) intersect(unique(g), rownames(z)))
  sets <- sets[lengths(sets) >= min_size]
  do.call(cbind, lapply(sets, function(g) {
    if (length(g) == 1) z[g, , drop = TRUE] else colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  }))
}

prot_scores <- score_sets(prot)
mrna_scores <- score_sets(mrna)
common_sets <- intersect(colnames(prot_scores), colnames(mrna_scores))
prot_scores <- prot_scores[, common_sets, drop = FALSE]
mrna_scores <- mrna_scores[, common_sets, drop = FALSE]
cat("programme scores:", nrow(prot_scores), "protein samples x", ncol(prot_scores), "sets\n")

## 3. Frozen classifier from the discovery analysis ---------------------------
ec <- readRDS(file.path(PROC, "ecotypes.rds"))
centroids <- ec$classifier$centroids
features <- ec$classifier$features
shared <- intersect(features, common_sets)
stopifnot(length(shared) >= 20)
cen <- centroids[, match(shared, features), drop = FALSE]
cen_z <- t(scale(t(cen)))
colnames(cen_z) <- shared
rownames(cen_z) <- ec$classifier$levels

assign_ecotype <- function(scores) {
  s <- scores[, shared, drop = FALSE]
  s_z <- scale(s)
  sim <- cor(t(s_z), t(cen_z), use = "pairwise.complete.obs")
  labels <- rownames(cen_z)[apply(sim, 1, which.max)]
  list(labels = labels, similarity = sim)
}

prot_call <- assign_ecotype(prot_scores)
mrna_call <- assign_ecotype(mrna_scores)
common_samples <- intersect(rownames(prot_scores), rownames(mrna_scores))
cat("samples with both layers:", length(common_samples), "\n")

prot_lab <- prot_call$labels[match(common_samples, rownames(prot_scores))]
mrna_lab <- mrna_call$labels[match(common_samples, rownames(mrna_scores))]
agreement <- mean(prot_lab == mrna_lab)
ari <- mclust::adjustedRandIndex(prot_lab, mrna_lab)

## 4. Clinical table + survival ---------------------------------------------
clin <- fread(file.path(IN, "data_clinical_patient.txt"), skip = 4, data.table = FALSE)
sam <- fread(file.path(IN, "data_clinical_sample.txt"), skip = 4, data.table = FALSE)
names(clin)[names(clin) == "PATIENT_ID"] <- "patient_id"
dat <- data.frame(
  sample_id = common_samples,
  ecotype_protein = prot_lab,
  ecotype_mrna = mrna_lab,
  stringsAsFactors = FALSE
)
dat$patient_id <- sam$PATIENT_ID[match(dat$sample_id, sam$SAMPLE_ID)]
dat <- merge(dat, clin, by = "patient_id", all.x = TRUE, sort = FALSE)
dat$os_time <- suppressWarnings(as.numeric(dat[[grep("^OS_MONTHS$", names(dat), value = TRUE)]])) * 30.44
dat$os_event <- ifelse(grepl("DECEASED|DEAD", dat$OS_STATUS, ignore.case = TRUE), 1, 0)
dat$age <- suppressWarnings(as.numeric(dat$AGE))
dat$sex <- dat$SEX

km <- NULL
cox_tab <- NULL
surv_dat <- dat[complete.cases(dat$os_time, dat$os_event) & dat$os_time > 0, , drop = FALSE]
if (nrow(surv_dat) >= 20 && sum(surv_dat$os_event) >= 5) {
  surv_dat$ecotype <- factor(surv_dat$ecotype_protein, levels = c("E1", "E2"))
  fit <- survival::coxph(survival::Surv(os_time, os_event) ~ ecotype, data = surv_dat)
  s <- summary(fit)
  cox_tab <- data.frame(
    cohort = "CPTAC-COAD (protein)",
    n = nrow(surv_dat),
    events = sum(surv_dat$os_event),
    c_index = as.numeric(survival::concordance(fit)$concordance),
    hr = unname(s$conf.int[1, 1]),
    lower = unname(s$conf.int[1, 3]),
    upper = unname(s$conf.int[1, 4]),
    p_hr = unname(s$coefficients[1, 5])
  )
  km <- survival::survfit(survival::Surv(os_time, os_event) ~ ecotype, data = surv_dat)
  lr <- survival::survdiff(survival::Surv(os_time, os_event) ~ ecotype, data = surv_dat)
  cox_tab$logrank_p <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)
}

## 5. Cross-platform concordance per programme ------------------------------
prog_cor <- data.frame(
  programme = shared,
  r_mrna_protein = sapply(shared, function(s) {
    stats::cor(mrna_scores[common_samples, s], prot_scores[common_samples, s], use = "complete.obs")
  })
)
prog_cor <- prog_cor[order(-prog_cor$r_mrna_protein), , drop = FALSE]
fwrite(prog_cor, file.path(OUT, "proteomics_programme_correlations.tsv"), sep = "\t")

## 6. Figures -----------------------------------------------------------------
lab_df <- data.frame(
  sample_id = common_samples,
  ecotype_protein = prot_lab,
  ecotype_mrna = mrna_lab,
  stringsAsFactors = FALSE
)

# plot_16 — protein programme heatmap by ecotype
hm_mat <- t(scale(prot_scores[, shared, drop = FALSE]))
hm_mat <- hm_mat[order(-apply(hm_mat, 1, stats::var)), , drop = FALSE]
hm_mat <- hm_mat[seq_len(min(25, nrow(hm_mat))), , drop = FALSE]
ord <- order(lab_df$ecotype_protein)
save_plot(function() {
  pheatmap::pheatmap(
    hm_mat[, ord, drop = FALSE],
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = data.frame(ecotype = lab_df$ecotype_protein[ord], row.names = lab_df$sample_id[ord]),
    main = "CPTAC-COAD: Hallmark programme scores from protein",
    fontsize_row = 8
  )
}, file.path(FIG, "plot_16"), 9, 7)

# plot_17 — mRNA vs protein programme correlation
save_plot(function() {
  print(ggplot(prog_cor, aes(reorder(programme, r_mrna_protein), r_mrna_protein)) +
    geom_col(fill = "grey35") +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_flip() +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "Pearson r (mRNA vs protein)",
         title = "Programme-level cross-platform concordance (CPTAC-COAD)"))
}, file.path(FIG, "plot_17"), 8, 9)

# plot_18 — KM by protein-assigned ecotype
if (!is.null(km)) {
  save_plot(function() {
    plot(km, col = unname(ECOTYPE_COL), lwd = 2, xlab = "Days",
         ylab = "Overall survival probability", mark.time = TRUE,
         main = "CPTAC-COAD: protein-assigned ecotypes")
    legend("topright", legend = c("E1", "E2"), col = unname(ECOTYPE_COL), lwd = 2, bty = "n")
    legend("bottomleft", bty = "n",
           legend = paste0("log-rank p=", signif(cox_tab$logrank_p, 2), "; HR=",
                           round(cox_tab$hr, 2), " (", round(cox_tab$lower, 2), "-", round(cox_tab$upper, 2), ")"))
  }, file.path(FIG, "plot_18"), 7, 5)
}

# plot_19 — ecotype agreement between RNA and protein
agree_tab <- as.data.frame(table(mRNA = lab_df$ecotype_mrna, Protein = lab_df$ecotype_protein))
save_plot(function() {
  print(ggplot(agree_tab, aes(mRNA, Protein, fill = Freq)) +
    geom_tile() + geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "#2E86C1") +
    theme_bw(base_size = 11) +
    labs(title = sprintf("Ecotype assignment: mRNA vs protein (agreement %.0f%%, ARI %.2f)",
                         100 * agreement, ari), x = "mRNA-assigned", y = "protein-assigned"))
}, file.path(FIG, "plot_19"), 6, 5)

## 7. Tables and report -------------------------------------------------------
res <- data.frame(
  metric = c("n_protein_samples", "n_mrna_samples", "n_both", "n_programmes_scored",
             "ecotype_E1_protein", "ecotype_E2_protein", "agreement_mrna_vs_protein", "ari_mrna_vs_protein"),
  value = c(ncol(prot), ncol(mrna), length(common_samples), length(shared),
            sum(lab_df$ecotype_protein == "E1"), sum(lab_df$ecotype_protein == "E2"),
            round(agreement, 4), round(ari, 4))
)
fwrite(res, file.path(OUT, "proteomics_validation_metrics.tsv"), sep = "\t")
if (!is.null(cox_tab)) {
  fwrite(cox_tab, file.path(OUT, "proteomics_validation.tsv"), sep = "\t")
  print(cox_tab)
}
print(res)
cat("median programme correlation mRNA vs protein:", round(median(prog_cor$r_mrna_protein), 3), "\n")

log_decision("17_proteomics", "PASS",
             paste0("agreement ", round(100 * agreement, 1), "%; ARI ", round(ari, 3),
                    "; median programme r ", round(median(prog_cor$r_mrna_protein), 3),
                    "; gene-set source: ", hl$source,
                    if (!is.null(cox_tab)) paste0("; OS HR ", round(cox_tab$hr, 2), " p ", signif(cox_tab$p_hr, 2)) else ""))
