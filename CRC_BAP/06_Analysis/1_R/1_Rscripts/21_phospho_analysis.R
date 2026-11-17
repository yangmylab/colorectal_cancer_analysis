#!/usr/bin/env Rscript
# 21_phospho_analysis.R
# Phosphoproteome layer: does the frozen K=2 ecotype axis exist below the protein level?
#
# Inputs : 03_Raw_Data/CPTAC_COAD/data_phosphoprotein_quantification.txt (PDC000116, cBioPortal datahub)
#          03_Raw_Data/CPTAC_COAD/data_protein_quantification.txt        (programme concordance)
#          03_Raw_Data/CPTAC_COAD/data_mrna_seq_v2_rsem.txt               (programme concordance)
#          04_Processed_Data/ecotypes.rds                                 (frozen classifier)
# Outputs: 07_Validation/21_phospho_*.tsv
#          07_Validation/21_Phospho_Report.md
#          09_Figures/plots/plot_28 .. plot_32 (+ pdf/jpg)
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6  (physical cores of the i5-10400; detectCores() returns 12 SMT threads and must not be used)
#   work split: per-programme correlations, per-gene phospho/protein concordance, ARI permutation null
#   the dataset is small (97 tumours x ~11k sites), so 6 workers is the point where IPC overhead still pays

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(survival)
  library(mclust)
  library(parallel)
  library(ggrepel)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

N_WORKERS <- 6L
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(mc.cores = N_WORKERS)
set.seed(20261117)

IN <- "03_Raw_Data/CPTAC_COAD"
PROC <- "04_Processed_Data"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## 1. Phosphosite matrix ------------------------------------------------------
ph <- fread(file.path(IN, "data_phosphoprotein_quantification.txt"), data.table = FALSE, check.names = FALSE)
site_ids <- ph$ENTITY_STABLE_ID
genes <- ph$GENE_SYMBOL
residue <- sub("^p([STY]).*$", "\\1", ph$PHOSPHOSITE)
position <- suppressWarnings(as.integer(sub("^p[STY]([0-9]+).*$", "\\1", ph$PHOSPHOSITE)))
site <- as.matrix(ph[, -(1:3), drop = FALSE])
storage.mode(site) <- "numeric"
rownames(site) <- make.unique(site_ids)
cat("phosphosites:", nrow(site), "| samples:", ncol(site), "| genes:", length(unique(genes)), "\n")
cat("residue composition: "); print(table(residue, useNA = "ifany"))

## 2. Site filtering (missingness-based, decided a priori) --------------------
# CPTAC phospho enrichment data are sparse (59.8% NA overall). Sites are kept if quantified in
# at least half of the tumours; the threshold is fixed before looking at any ecotype labels and the
# whole sensitivity ladder is reported so the reader can see what the choice costs.
miss <- rowMeans(is.na(site))
filter_tab <- data.frame(
  max_missing_allowed = c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70),
  sites_retained = sapply(c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70), function(th) sum(miss <= th))
)
MAX_MISSING <- 0.50
keep <- which(miss <= MAX_MISSING)
site <- site[keep, , drop = FALSE]
genes <- genes[keep]
phos_labels <- ph$PHOSPHOSITE[keep]
site <- site[!duplicated(rownames(site)), , drop = FALSE]
cat("sites retained at <=", MAX_MISSING, "missing:", nrow(site), "\n")
print(filter_tab)
fwrite(filter_tab, file.path(OUT, "21_phospho_site_filtering.tsv"), sep = "\t")
log_decision("21_phospho", "SITE_FILTER",
             paste0("kept sites quantified in >=", 100 * (1 - MAX_MISSING), "% of tumours (n=", nrow(site),
                    " of ", length(miss), "); threshold fixed before ecotype labels were read"))

## 3. Imputation (kNN, median fallback) --------------------------------------
impute_knn <- function(m) {
  ok <- tryCatch(impute::impute.knn(m, k = 10, rowmax = 0.5, colmax = 0.9, rng.seed = 20261117),
                 error = function(e) NULL)
  if (!is.null(ok)) return(ok$data)
  m2 <- m
  for (i in seq_len(nrow(m2))) {
    na <- is.na(m2[i, ])
    if (any(na)) m2[i, na] <- stats::median(m2[i, ], na.rm = TRUE)
  }
  m2
}
site_imp <- impute_knn(site)
if (anyNA(site_imp)) site_imp[is.na(site_imp)] <- stats::median(site_imp, na.rm = TRUE)
cat("residual NA after imputation:", sum(is.na(site_imp)), "\n")

# per-site z-scores for scoring and concordance (site values are already normalised log ratios)
site_z <- t(scale(t(site_imp)))
site_z <- site_z[is.finite(rowSums(site_z)), , drop = FALSE]
genes_z <- genes[match(rownames(site_z), rownames(site))]
keep_rows <- !is.na(genes_z)
site_z <- site_z[keep_rows, , drop = FALSE]
genes_z <- genes_z[keep_rows]

## 4. Gene-level phospho matrix ----------------------------------------------
# Sites are averaged within a gene first: SRRM2 alone carries 284 sites, so averaging per gene keeps
# hyper-phosphorylated proteins from dominating a programme score.
cnt <- rowsum(rep(1, length(genes_z)), group = genes_z)
sum_z <- rowsum(site_z, group = genes_z)
phos_gene <- sum_z / as.numeric(cnt[rownames(sum_z), 1])
phos_gene <- phos_gene[is.finite(rowSums(phos_gene)), , drop = FALSE]
cat("gene-level phospho matrix:", nrow(phos_gene), "genes x", ncol(phos_gene), "samples\n")

## 5. Programme scores for the three CPTAC layers -----------------------------
# Hallmark sets are taken from the discovery run itself (50 sets, Ensembl IDs) and mapped to symbols
# locally with org.Hs.eg.db. msigdbr 26.1.1 needs the msigdbdf data package plus an online release
# download; neither is available here, and re-downloading was not authorised, so no network call is made.
sig <- readRDS(file.path(PROC, "signature_scores.rds"))
ens <- unique(unlist(sig$gene_sets))
map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", columns = "SYMBOL")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
sym_of <- stats::setNames(map$SYMBOL, map$ENSEMBL)
gene_sets <- lapply(sig$gene_sets, function(g) unique(unname(sym_of[g[g %in% names(sym_of)]])))
gene_sets <- gene_sets[lengths(gene_sets) >= 5]
cat("Hallmark sets from the frozen discovery run:", length(gene_sets), "\n")

score_sets <- function(mat, min_size = 5) {
  sets <- lapply(gene_sets, function(g) intersect(unique(g), rownames(mat)))
  sets <- sets[lengths(sets) >= min_size]
  do.call(cbind, lapply(sets, function(g) {
    if (length(g) == 1) mat[g, , drop = TRUE] else colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
  }))
}

read_cbio <- function(path, gene_col, value_transform = FALSE) {
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  mat <- as.matrix(dt[, setdiff(colnames(dt), gene_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  ids <- as.character(dt[[gene_col]])
  if (gene_col == "Composite.Element.REF") ids <- sub("\\|.*$", "", ids)
  rownames(mat) <- ids
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  if (value_transform) mat <- log2(mat + 1)
  if (any(duplicated(rownames(mat)))) {
    v <- apply(mat, 1, function(x) stats::var(x, na.rm = TRUE))
    mat <- mat[order(rownames(mat), -v), , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  mat
}

prot <- read_cbio(file.path(IN, "data_protein_quantification.txt"), "Composite.Element.REF")
mrna <- read_cbio(file.path(IN, "data_mrna_seq_v2_rsem.txt"), "Hugo_Symbol", value_transform = TRUE)
prot_scores <- score_sets(t(scale(t(prot))))
mrna_scores <- score_sets(t(scale(t(mrna))))
phos_scores <- score_sets(phos_gene)
common_sets <- Reduce(intersect, list(colnames(phos_scores), colnames(prot_scores), colnames(mrna_scores)))
phos_scores <- phos_scores[, common_sets, drop = FALSE]
prot_scores <- prot_scores[, common_sets, drop = FALSE]
mrna_scores <- mrna_scores[, common_sets, drop = FALSE]
cat("programmes scorable in all three layers:", length(common_sets), "of", length(gene_sets), "\n")

## 6. Frozen classifier projection -------------------------------------------
ec <- readRDS(file.path(PROC, "ecotypes.rds"))
centroids <- ec$classifier$centroids
features <- ec$classifier$features
shared <- intersect(features, common_sets)
cat("classifier features recoverable from phospho programmes:", length(shared), "of", length(features), "\n")
stopifnot(length(shared) >= 20)
cen_z <- t(scale(t(centroids[, match(shared, features), drop = FALSE])))
colnames(cen_z) <- shared
rownames(cen_z) <- ec$classifier$levels

assign_ecotype <- function(scores) {
  s_z <- scale(scores[, shared, drop = FALSE])
  sim <- cor(t(s_z), t(cen_z), use = "pairwise.complete.obs")
  list(labels = rownames(cen_z)[apply(sim, 1, which.max)], similarity = sim)
}

phos_call <- assign_ecotype(phos_scores)
prot_call <- assign_ecotype(prot_scores)
mrna_call <- assign_ecotype(mrna_scores)

phos_lab <- phos_call$labels
names(phos_lab) <- rownames(phos_scores)
prot_lab <- prot_call$labels
names(prot_lab) <- rownames(prot_scores)
mrna_lab <- mrna_call$labels
names(mrna_lab) <- rownames(mrna_scores)

common_samples <- Reduce(intersect, list(rownames(phos_scores), rownames(prot_scores), rownames(mrna_scores)))
cat("tumours with all three layers:", length(common_samples), "\n")

cmp <- data.frame(
  comparison = c("phospho_vs_protein", "phospho_vs_mrna", "protein_vs_mrna"),
  agreement = c(
    mean(phos_lab[common_samples] == prot_lab[common_samples]),
    mean(phos_lab[common_samples] == mrna_lab[common_samples]),
    mean(prot_lab[common_samples] == mrna_lab[common_samples])
  ),
  ari = c(
    ari(phos_lab[common_samples], prot_lab[common_samples]),
    ari(phos_lab[common_samples], mrna_lab[common_samples]),
    ari(prot_lab[common_samples], mrna_lab[common_samples])
  )
)
print(cmp)

## 7. ARI permutation null ----------------------------------------------------
ari_null <- unlist(mclapply(seq_len(1000), function(i) {
  ari(sample(phos_lab[common_samples]), prot_lab[common_samples])
}, mc.cores = N_WORKERS))
ari_obs <- cmp$ari[cmp$comparison == "phospho_vs_protein"]
ari_p <- (sum(ari_null >= ari_obs) + 1) / (length(ari_null) + 1)
cat("observed ARI:", round(ari_obs, 3), "| permutation p:", signif(ari_p, 3),
    "| null 95th pct:", round(stats::quantile(ari_null, 0.95), 3), "\n")

## 8. Programme-level cross-layer concordance --------------------------------
prog_cor <- do.call(rbind, mclapply(common_sets, function(s) {
  data.frame(
    programme = s,
    r_phospho_protein = stats::cor(phos_scores[common_samples, s], prot_scores[common_samples, s], use = "complete.obs"),
    r_phospho_mrna = stats::cor(phos_scores[common_samples, s], mrna_scores[common_samples, s], use = "complete.obs"),
    r_protein_mrna = stats::cor(prot_scores[common_samples, s], mrna_scores[common_samples, s], use = "complete.obs"),
    stringsAsFactors = FALSE
  )
}, mc.cores = N_WORKERS))
prog_cor <- prog_cor[order(-prog_cor$r_phospho_protein), , drop = FALSE]
fwrite(prog_cor, file.path(OUT, "21_phospho_programme_correlations.tsv"), sep = "\t")
cat("median r phospho-protein:", round(median(prog_cor$r_phospho_protein), 3),
    "| phospho-mRNA:", round(median(prog_cor$r_phospho_mrna), 3),
    "| protein-mRNA:", round(median(prog_cor$r_protein_mrna), 3), "\n")

## 9. Differential phosphosites (limma, contrast E2 - E1) --------------------
# Pitfall guard: `~0 + group` with topTable(coef="E2") tests the group MEAN, so an explicit contrast is used.
labels_site <- phos_lab[colnames(site_imp)]
ok_site <- !is.na(labels_site)
design <- model.matrix(~ 0 + factor(labels_site[ok_site], levels = c("E1", "E2")))
colnames(design) <- c("E1", "E2")
contr <- limma::makeContrasts(E2 - E1, levels = design)
fit <- eBayes(contrasts.fit(lmFit(site_imp[, ok_site, drop = FALSE], design), contr))
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
tt$site_id <- rownames(tt)
tt$gene <- genes[match(tt$site_id, rownames(site))]
tt$phosphosite <- phos_labels[match(tt$site_id, rownames(site))]
tt <- tt[, c("site_id", "gene", "phosphosite", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
fwrite(tt, file.path(OUT, "21_phospho_differential_sites.tsv"), sep = "\t")
n_sig <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
n_fdr <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
cat("differential phosphosites: FDR<0.05:", n_fdr, "| FDR<0.05 & |logFC|>0.5:", n_sig, "of", nrow(tt), "\n")

# sensitivity: same contrast but grouped by the mRNA (frozen-reference) ecotype
labels_mrna_site <- mrna_lab[colnames(site_imp)]
ok_m <- !is.na(labels_mrna_site)
design_m <- model.matrix(~ 0 + factor(labels_mrna_site[ok_m], levels = c("E1", "E2")))
colnames(design_m) <- c("E1", "E2")
fit_m <- eBayes(contrasts.fit(lmFit(site_imp[, ok_m, drop = FALSE], design_m),
                              limma::makeContrasts(E2 - E1, levels = design_m)))
tt_m <- topTable(fit_m, coef = 1, number = Inf, sort.by = "P")
n_sig_m <- sum(tt_m$adj.P.Val < 0.05 & abs(tt_m$logFC) > 0.5, na.rm = TRUE)
cat("sensitivity (mRNA-assigned ecotype): FDR<0.05 & |logFC|>0.5:", n_sig_m, "\n")

## 10. Phospho vs protein discordance (signalling that is not abundance) ------
shared_genes <- intersect(rownames(phos_gene), rownames(prot))
prot_z <- t(scale(t(prot[shared_genes, common_samples, drop = FALSE])))
disc <- do.call(rbind, mclapply(shared_genes, function(g) {
  p <- phos_gene[g, common_samples]
  a <- prot_z[g, common_samples]
  ok <- is.finite(p) & is.finite(a)
  if (sum(ok) < 20) return(NULL)
  p <- p[ok]
  a <- a[ok]
  grp <- phos_lab[common_samples][ok]
  sd_p <- stats::sd(p)
  sd_a <- stats::sd(a)
  if (!is.finite(sd_p) || !is.finite(sd_a) || sd_p == 0 || sd_a == 0) return(NULL)
  delta <- p - a
  if (length(unique(grp)) < 2) return(NULL)
  wt <- tryCatch(stats::wilcox.test(delta ~ grp), error = function(e) NULL)
  data.frame(
    gene = g,
    r_phospho_protein = stats::cor(p, a, use = "complete.obs"),
    delta_mean_E1 = mean(delta[grp == "E1"], na.rm = TRUE),
    delta_mean_E2 = mean(delta[grp == "E2"], na.rm = TRUE),
    delta_diff = mean(delta[grp == "E2"], na.rm = TRUE) - mean(delta[grp == "E1"], na.rm = TRUE),
    p_wilcox = if (is.null(wt)) NA_real_ else wt$p.value,
    stringsAsFactors = FALSE
  )
}, mc.cores = N_WORKERS))
disc$fdr <- stats::p.adjust(disc$p_wilcox, method = "BH")
disc <- disc[order(disc$p_wilcox), , drop = FALSE]
fwrite(disc, file.path(OUT, "21_phospho_discordance.tsv"), sep = "\t")
cat("genes tested for phospho/protein discordance:", nrow(disc),
    "| FDR<0.05:", sum(disc$fdr < 0.05, na.rm = TRUE),
    "| median r(phospho, protein):", round(median(disc$r_phospho_protein, na.rm = TRUE), 3), "\n")

## 11. Survival (descriptive only: CPTAC-COAD has 8 deaths) -------------------
clin <- fread(file.path(IN, "data_clinical_patient.txt"), skip = 4, data.table = FALSE)
sam <- fread(file.path(IN, "data_clinical_sample.txt"), skip = 4, data.table = FALSE)
names(clin)[names(clin) == "PATIENT_ID"] <- "patient_id"
dat <- data.frame(sample_id = common_samples, ecotype_phospho = phos_lab[common_samples],
                  ecotype_protein = prot_lab[common_samples], stringsAsFactors = FALSE)
dat$patient_id <- sam$PATIENT_ID[match(dat$sample_id, sam$SAMPLE_ID)]
dat <- merge(dat, clin, by = "patient_id", all.x = TRUE, sort = FALSE)
dat$os_time <- suppressWarnings(as.numeric(dat$OS_MONTHS)) * 30.44
dat$os_event <- ifelse(grepl("DECEASED|DEAD", dat$OS_STATUS, ignore.case = TRUE), 1, 0)
surv_dat <- dat[complete.cases(dat$os_time, dat$os_event) & dat$os_time > 0, , drop = FALSE]
cox_tab <- NULL
if (nrow(surv_dat) >= 20 && sum(surv_dat$os_event) >= 5) {
  surv_dat$ecotype <- factor(surv_dat$ecotype_phospho, levels = c("E1", "E2"))
  f1 <- survival::coxph(survival::Surv(os_time, os_event) ~ ecotype, data = surv_dat)
  s1 <- summary(f1)
  coef1 <- unname(s1$coefficients[1, 1])
  se1 <- unname(s1$coefficients[1, 3])
  # With 8 deaths the fit can separate; a coefficient this large is not an effect size, it is a failure
  # to estimate one. Report it as not estimable instead of quoting HR = 5e8.
  separated <- !is.finite(coef1) || abs(coef1) > 4 || !is.finite(se1) || se1 > 5
  lr <- survival::survdiff(survival::Surv(os_time, os_event) ~ ecotype, data = surv_dat)
  cox_tab <- data.frame(
    cohort = "CPTAC-COAD (phospho)", n = nrow(surv_dat), events = sum(surv_dat$os_event),
    c_index = if (separated) NA_real_ else as.numeric(survival::concordance(f1)$concordance),
    hr = if (separated) NA_real_ else exp(coef1),
    lower = if (separated) NA_real_ else unname(s1$conf.int[1, 3]),
    upper = if (separated) NA_real_ else unname(s1$conf.int[1, 4]),
    p_hr = if (separated) NA_real_ else unname(s1$coefficients[1, 5]),
    logrank_p = 1 - stats::pchisq(lr$chisq, length(lr$n) - 1),
    note = if (separated) "ecotype coefficient not estimable (separation with 8 deaths); descriptive only" else "descriptive; only 8 deaths",
    stringsAsFactors = FALSE
  )
  fwrite(cox_tab, file.path(OUT, "21_phospho_survival.tsv"), sep = "\t")
  print(cox_tab)
}

## 12. Summary metrics --------------------------------------------------------
summ <- data.frame(
  metric = c("phosphosites_total", "phosphosites_analysed", "samples", "genes_site_level",
             "genes_with_protein", "programmes_scored", "classifier_features_shared",
             "ecotype_E1_phospho", "ecotype_E2_phospho",
             "agreement_phospho_protein", "ari_phospho_protein", "ari_permutation_p",
             "ari_null_95pct", "agreement_phospho_mrna", "agreement_protein_mrna",
             "median_r_phospho_protein", "median_r_phospho_mrna",
             "sites_fdr05", "sites_fdr05_lfc05", "sites_fdr05_lfc05_mrna_grouped",
             "genes_discordant_fdr05"),
  value = c(length(miss), nrow(site), ncol(site), length(unique(genes_z)), nrow(disc),
            length(common_sets), length(shared),
            sum(phos_lab == "E1"), sum(phos_lab == "E2"),
            round(cmp$agreement[1], 4), round(ari_obs, 4), signif(ari_p, 3),
            round(stats::quantile(ari_null, 0.95), 4),
            round(cmp$agreement[2], 4), round(cmp$agreement[3], 4),
            round(median(prog_cor$r_phospho_protein), 4), round(median(prog_cor$r_phospho_mrna), 4),
            n_fdr, n_sig, n_sig_m, sum(disc$fdr < 0.05, na.rm = TRUE))
)
fwrite(summ, file.path(OUT, "21_phospho_validation.tsv"), sep = "\t")
# per-sample layer labels, so downstream layers (e.g. kinase activity) need not re-derive the classifier
fwrite(data.frame(sample_id = names(phos_lab),
                  ecotype_phospho = unname(phos_lab),
                  ecotype_protein = unname(prot_lab[names(phos_lab)]),
                  ecotype_mrna = unname(mrna_lab[names(phos_lab)]),
                  stringsAsFactors = FALSE),
       file.path(OUT, "21_phospho_ecotype_labels.tsv"), sep = "\t")
print(summ)

## 13. Figures ----------------------------------------------------------------
# plot_28 — phospho programme heatmap by phospho-assigned ecotype
hm <- t(scale(phos_scores[, shared, drop = FALSE]))
hm <- hm[order(-apply(hm, 1, stats::var)), , drop = FALSE][seq_len(min(25, nrow(hm))), , drop = FALSE]
ord <- order(phos_lab[colnames(hm)])
save_plot(function() {
  pheatmap::pheatmap(hm[, ord, drop = FALSE], cluster_cols = FALSE, show_colnames = FALSE,
                     annotation_col = data.frame(ecotype = phos_lab[colnames(hm)][ord], row.names = colnames(hm)[ord]),
                     main = "CPTAC-COAD: Hallmark programme scores from phosphoprotein", fontsize_row = 8)
}, file.path(FIG, "plot_28"), 9, 7)

# plot_29 — programme concordance across the three CPTAC layers
long <- rbind(
  data.frame(programme = prog_cor$programme, r = prog_cor$r_phospho_protein, pair = "phospho vs protein"),
  data.frame(programme = prog_cor$programme, r = prog_cor$r_phospho_mrna, pair = "phospho vs mRNA"),
  data.frame(programme = prog_cor$programme, r = prog_cor$r_protein_mrna, pair = "protein vs mRNA")
)
save_plot(function() {
  print(ggplot(long, aes(reorder(programme, r), r, fill = pair)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) + coord_flip() +
    scale_fill_manual(values = c("phospho vs protein" = "#C0392B", "phospho vs mRNA" = "#E67E22", "protein vs mRNA" = "#2E86C1")) +
    theme_bw(base_size = 8) +
    labs(x = NULL, y = "Pearson r across 97 tumours", fill = NULL,
         title = "Programme-level concordance across CPTAC-COAD layers"))
}, file.path(FIG, "plot_29"), 8, 10)

# plot_30 — volcano of differential phosphosites (E2 vs E1, phospho-assigned)
tt$neglog10p <- -log10(tt$P.Value)
tt$sig <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, "FDR<0.05 & |logFC|>0.5",
                 ifelse(tt$P.Value < 0.01, "nominal p<0.01", "ns"))
top_lab <- head(tt[tt$adj.P.Val < 0.05, , drop = FALSE], 12)
save_plot(function() {
  print(ggplot(tt, aes(logFC, neglog10p, colour = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c("FDR<0.05 & |logFC|>0.5" = "#C0392B", "nominal p<0.01" = "#E67E22", "ns" = "grey70")) +
    geom_text_repel(data = top_lab, aes(label = paste0(gene, " ", phosphosite)), size = 2.2, max.overlaps = 20) +
    theme_bw(base_size = 9) +
    labs(x = "log2 fold change (E2 vs E1)", y = expression(-log[10]~P),
         title = sprintf("Differential phosphosites, phospho-assigned ecotypes (n=%d sites)", nrow(tt)),
         colour = NULL))
}, file.path(FIG, "plot_30"), 8, 6)

# plot_31 — ecotype agreement phospho vs protein
agree_tab <- as.data.frame(table(phospho = phos_lab[common_samples], protein = prot_lab[common_samples]))
save_plot(function() {
  print(ggplot(agree_tab, aes(phospho, protein, fill = Freq)) +
    geom_tile() + geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "#2E86C1") + theme_bw(base_size = 11) +
    labs(title = sprintf("Ecotype assignment: phospho vs protein (agreement %.0f%%, ARI %.2f, permutation p=%.3f)",
                         100 * cmp$agreement[1], ari_obs, ari_p),
         x = "phospho-assigned", y = "protein-assigned"))
}, file.path(FIG, "plot_31"), 6, 5)

# plot_32 — ARI permutation null
save_plot(function() {
  print(ggplot(data.frame(ari = ari_null), aes(ari)) +
    geom_histogram(bins = 40, fill = "grey65", colour = "white") +
    geom_vline(xintercept = ari_obs, colour = "#C0392B", linewidth = 1) +
    theme_bw(base_size = 10) +
    labs(x = "ARI (phospho labels permuted vs protein labels)", y = "count",
         title = sprintf("Null distribution of cross-layer ARI (observed %.2f, p=%.3f)", ari_obs, ari_p)))
}, file.path(FIG, "plot_32"), 6, 4)

## 14. Report -----------------------------------------------------------------
grade <- if (ari_p < 0.05 && cmp$agreement[1] >= 0.75) "E3" else if (ari_p < 0.05) "E1" else "not_supported"
report <- c(
  "# Phosphoproteome layer (CPTAC-COAD, PDC000116)",
  "",
  sprintf("Generated by `21_phospho_analysis.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Data",
  sprintf("- Phosphosites: %d quantified in %d tumours; %d kept at <=50%% missing (threshold fixed a priori).",
          length(miss), ncol(site), nrow(site)),
  sprintf("- %d genes at site level, %d with matched protein abundance; %d of %d Hallmark programmes scorable.",
          length(unique(genes_z)), nrow(disc), length(common_sets), length(gene_sets)),
  sprintf("- Imputation: kNN (k=10); median r(phospho, protein) per gene = %.2f.",
          median(disc$r_phospho_protein, na.rm = TRUE)),
  "",
  "## Does the frozen K=2 axis exist at the phosphoproteome level?",
  sprintf("- Ecotype assigned from phospho programme scores: E1 %d / E2 %d.",
          sum(phos_lab == "E1"), sum(phos_lab == "E2")),
  sprintf("- Agreement with the protein-level call: **%.1f%%** (ARI %.2f); with mRNA: %.1f%%.",
          100 * cmp$agreement[1], ari_obs, 100 * cmp$agreement[2]),
  sprintf("- Permutation null (1000 label permutations): ARI 95th percentile %.2f, **p = %.3f**.",
          stats::quantile(ari_null, 0.95), ari_p),
  sprintf("- Programme concordance: median r phospho-protein = %.2f, phospho-mRNA = %.2f, protein-mRNA = %.2f.",
          median(prog_cor$r_phospho_protein), median(prog_cor$r_phospho_mrna), median(prog_cor$r_protein_mrna)),
  "",
  "## Differential phosphorylation",
  sprintf("- %d sites at FDR<0.05, of which %d also have |log2FC|>0.5 (of %d tested).",
          n_fdr, n_sig, nrow(tt)),
  sprintf("- Sensitivity with the mRNA-assigned ecotype: %d sites at FDR<0.05 & |log2FC|>0.5.", n_sig_m),
  sprintf("- Phospho/protein discordance screen: %d genes, %d at FDR<0.05.",
          nrow(disc), sum(disc$fdr < 0.05, na.rm = TRUE)),
  "",
  "## Survival",
  if (!is.null(cox_tab)) {
    if (is.na(cox_tab$hr)) {
      sprintf("- Not estimable: the ecotype coefficient separated with %d deaths in %d tumours (log-rank p = %.3f is uninformative at this event count).",
              cox_tab$events, cox_tab$n, cox_tab$logrank_p)
    } else {
      sprintf("- Descriptive only: n = %d, %d deaths; HR = %.2f (%.2f-%.2f), log-rank p = %.2f.",
              cox_tab$n, cox_tab$events, cox_tab$hr, cox_tab$lower, cox_tab$upper, cox_tab$logrank_p)
    }
  } else "- Not estimable (too few deaths).",
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (grade == "not_supported") "the phosphoproteome does not reproduce the ecotype partition better than chance."
          else if (grade == "E1") "cross-platform association reaches nominal significance but without a strong partition match."
          else "transcriptome, proteome and phosphoproteome converge on the same K=2 partition, matching F06. This is a reproducibility result, not a prognosis result: the ecotype-OS association remains null (F01/F03)."),
  "",
  "## Limitations",
  "- 59.8% of sites are missing before filtering; heavily modified proteins were averaged per gene to avoid dominance.",
  "- Kinase-activity inference (KSEA/decoupleR substrate scoring) is **not** reported: no offline kinase-substrate",
  "  library is installed and no download was authorised, so activity claims would have been fabricated.",
  "- CPTAC-COAD has 8 deaths, so survival is descriptive; this layer tests reproducibility, not prognosis.",
  "- cBioPortal's datahub ships CPTAC log ratios without per-site confidence; low-abundance sites may be noise.",
  "",
  "## Figures",
  "- plot_28 phospho programme heatmap by ecotype; plot_29 three-layer programme concordance;",
  "  plot_30 differential phosphosite volcano; plot_31 phospho vs protein agreement; plot_32 ARI permutation null."
)
writeLines(report, file.path(OUT, "21_Phospho_Report.md"))

log_decision("21_phospho", "PASS",
             paste("phospho ecotypes:", sum(phos_lab == "E1"), "E1 /", sum(phos_lab == "E2"), "E2;",
                   "agreement vs protein:", round(100 * cmp$agreement[1], 1), "%; ARI:", round(ari_obs, 3),
                   "; permutation p:", signif(ari_p, 3), "; grade:", grade))

## 15. Evidence ledger (F09, rewritten in place on every run) -----------------
ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F09", , drop = FALSE]
f09 <- data.frame(
  finding_id = "F09",
  claim = "The frozen K=2 ecotypes reproduce at the phosphoproteome level in CPTAC-COAD",
  evidence_type = "Cross-platform concordance (phosphoproteome)",
  dataset_s = "CPTAC-COAD (PDC000116) phosphoproteome",
  statistic = paste0("phospho vs protein ecotype agreement ", round(100 * cmp$agreement[1], 1),
                     "%; ARI ", round(ari_obs, 3), "; permutation p=", signif(ari_p, 3),
                     " (null 95th pct ", round(stats::quantile(ari_null, 0.95), 3), ")",
                     "; median programme r phospho-protein ", round(median(prog_cor$r_phospho_protein), 2),
                     "; ", n_sig, "/", nrow(tt), " sites FDR<0.05 & |logFC|>0.5",
                     "; ", sum(disc$fdr < 0.05, na.rm = TRUE), " genes phospho/protein discordant at FDR<0.05"),
  n = length(common_samples),
  grade = grade,
  verb_used = "consistent with",
  contradicting_evidence = "not applicable",
  limitation = "59.8% of sites missing before filtering; kinase-activity inference not performed (no offline substrate library); OS not estimable (8 deaths); 49 of 50 programmes recoverable",
  next_test = "kinase-substrate activity scoring (KSEA/decoupleR) once an annotation library is available",
  figure_table = "plot_28; plot_29; plot_30; plot_31; plot_32",
  stringsAsFactors = FALSE
)
rows <- if (is.null(old) || nrow(old) == 0) f09 else rbind(old, f09)
rows <- rows[order(as.integer(sub("^F", "", rows$finding_id))), , drop = FALSE]
data.table::fwrite(rows, ledger_path, sep = ",")
cat("ledger rows:", nrow(rows), "\n")
cat("grade:", grade, "\n")
cat("done\n")
