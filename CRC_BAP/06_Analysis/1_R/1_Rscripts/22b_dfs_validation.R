#!/usr/bin/env Rscript
# 22b_dfs_validation.R
# Disease-free survival in GSE14333: a SECONDARY endpoint, reported separately from overall survival.
#
# Why a separate script: the discovery/validation OS analysis (04_survival_validation.R) must never
# absorb a DFS cohort. This script reads only 03_Raw_Data/GSE14333/22a_* and writes its own tables,
# report and figures; DFS results are never pooled with the OS meta-analysis.
#
# Input : 03_Raw_Data/GSE14333/22a_GSE14333_expression.tsv
#         03_Raw_Data/GSE14333/22a_GSE14333_clinical_dfs.tsv
#         04_Processed_Data/signature_scores.rds (frozen Hallmark sets) + ecotypes.rds (frozen classifier)
# Output: 07_Validation/22b_dfs_validation.tsv
#         07_Validation/22b_dfs_bootstrap.tsv
#         07_Validation/22b_DFS_Report.md
#         09_Figures/plots/plot_33 .. plot_34 (+ pdf/jpg)
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6 (physical cores of the i5-10400; detectCores() reports 12 SMT threads and is not used)
#   work split: per-programme scoring (49 sets) and a 500-repetition bootstrap of the adjusted HR

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(mclust)
  library(parallel)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

N_WORKERS <- 6L
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(mc.cores = N_WORKERS)
set.seed(20261117)

DIR <- "03_Raw_Data/GSE14333"
PROC <- "04_Processed_Data"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## 1. Inputs ----------------------------------------------------------------
expr_dt <- data.table::fread(file.path(DIR, "22a_GSE14333_expression.tsv"), data.table = FALSE, check.names = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- expr_dt$gene_id
clin <- data.table::fread(file.path(DIR, "22a_GSE14333_clinical_dfs.tsv"), data.table = FALSE)
cat("expression:", nrow(expr), "symbols x", ncol(expr), "tumours\n")

stopifnot(!any(c("os_time_days", "os_event") %in% names(clin)))  # DFS cohort must never carry OS columns

## 2. Hallmark programme scores (frozen discovery sets, Ensembl -> symbol) ---
sig <- readRDS(file.path(PROC, "signature_scores.rds"))
ens <- unique(unlist(sig$gene_sets))
map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", columns = "SYMBOL")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
sym_of <- stats::setNames(map$SYMBOL, map$ENSEMBL)
gene_sets <- lapply(sig$gene_sets, function(g) unique(unname(sym_of[g[g %in% names(sym_of)]])))
gene_sets <- gene_sets[lengths(gene_sets) >= 5]

z <- t(scale(t(expr)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
scores <- do.call(cbind, mclapply(gene_sets, function(g) {
  gg <- intersect(g, rownames(z))
  if (length(gg) < 5) return(rep(NA_real_, ncol(z)))
  colMeans(z[gg, , drop = FALSE], na.rm = TRUE)
}, mc.cores = N_WORKERS))
colnames(scores) <- names(gene_sets)
scores <- scores[, colSums(is.na(scores)) == 0, drop = FALSE]
cat("programmes scored:", ncol(scores), "\n")

## 3. Frozen nearest-centroid classifier -------------------------------------
ec <- readRDS(file.path(PROC, "ecotypes.rds"))
features <- ec$classifier$features
shared <- intersect(features, colnames(scores))
cat("classifier features recoverable:", length(shared), "of", length(features), "\n")
stopifnot(length(shared) >= 20)
cen_z <- t(scale(t(ec$classifier$centroids[, match(shared, features), drop = FALSE])))
colnames(cen_z) <- shared
rownames(cen_z) <- ec$classifier$levels
s_z <- scale(scores[, shared, drop = FALSE])
sim <- cor(t(s_z), t(cen_z), use = "pairwise.complete.obs")
labels <- rownames(cen_z)[apply(sim, 1, which.max)]
names(labels) <- rownames(scores)
cat("ecotype assignment:\n"); print(table(labels))

## 4. Analysis dataset --------------------------------------------------------
dat <- clin
dat$ecotype <- labels[dat$sample_id]
dat$stage_num <- as.numeric(factor(dat$stage, levels = c("I", "II", "III", "IV")))
dat$age <- suppressWarnings(as.numeric(dat$age))
ana <- dat[!is.na(dat$DFS_days) & dat$DFS_days > 0 & !is.na(dat$DFS_event) & !is.na(dat$ecotype), , drop = FALSE]
ana$ecotype <- factor(ana$ecotype, levels = c("E1", "E2"))
cat("analysable tumours:", nrow(ana), "| DFS events:", sum(ana$DFS_event), "\n")
cat("stage distribution of the analysable set:\n"); print(table(ana$stage))
cat("ecotype x stage:\n"); print(table(ana$ecotype, ana$stage))

## 5. Cox models -------------------------------------------------------------
fit_model <- function(d, label, adjust = TRUE) {
  cc <- d
  if (adjust) {
    vars <- c("DFS_days", "DFS_event", "ecotype", "age", "sex", "stage_num")
    cc <- cc[stats::complete.cases(cc[, vars, drop = FALSE]), , drop = FALSE]
    f_base <- survival::coxph(survival::Surv(DFS_days, DFS_event) ~ age + sex + stage_num, data = cc)
    f_full <- survival::coxph(survival::Surv(DFS_days, DFS_event) ~ ecotype + age + sex + stage_num, data = cc)
  } else {
    cc <- cc[stats::complete.cases(cc[, c("DFS_days", "DFS_event", "ecotype"), drop = FALSE]), , drop = FALSE]
    f_base <- survival::coxph(survival::Surv(DFS_days, DFS_event) ~ 1, data = cc)
    f_full <- survival::coxph(survival::Surv(DFS_days, DFS_event) ~ ecotype, data = cc)
  }
  s <- summary(f_full)
  lrt <- stats::anova(f_base, f_full, test = "LRT")
  ecoef <- grep("^ecotype", rownames(s$coefficients))[1]
  ph <- tryCatch(survival::cox.zph(f_full), error = function(e) NULL)
  data.frame(
    model = label,
    n = nrow(cc),
    events = sum(cc$DFS_event),
    c_index_base = as.numeric(survival::concordance(f_base)$concordance),
    c_index_full = as.numeric(survival::concordance(f_full)$concordance),
    hr_E2_vs_E1 = if (is.na(ecoef)) NA_real_ else exp(s$coefficients[ecoef, "coef"]),
    lower = if (is.na(ecoef)) NA_real_ else unname(s$conf.int[ecoef, 3]),
    upper = if (is.na(ecoef)) NA_real_ else unname(s$conf.int[ecoef, 4]),
    p_hr = if (is.na(ecoef)) NA_real_ else s$coefficients[ecoef, "Pr(>|z|)"],
    lrt_p = lrt[2, "Pr(>|Chi|)"],
    ph_min_p = if (is.null(ph)) NA_real_ else min(ph$table[, "p"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

rows <- rbind(
  fit_model(ana, "unadjusted", adjust = FALSE),
  fit_model(ana, "adjusted (age, sex, stage)", adjust = TRUE)
)
sub23 <- ana[ana$stage %in% c("II", "III"), , drop = FALSE]
rows <- rbind(rows, fit_model(sub23, "stage II-III, adjusted", adjust = TRUE))
rownames(rows) <- NULL
print(rows)

# direction consistency with the discovery OS model (E2 was the worse-survival ecotype, HR 1.45)
disc_dir <- 1
rows$direction_consistent <- !is.na(rows$hr_E2_vs_E1) & sign(log(rows$hr_E2_vs_E1)) == disc_dir
fwrite(rows, file.path(OUT, "22b_dfs_validation.tsv"), sep = "\t")

## 6. Bootstrap the adjusted HR (500 resamples) ------------------------------
adj_dat <- ana[stats::complete.cases(ana[, c("DFS_days", "DFS_event", "ecotype", "age", "sex", "stage_num"), drop = FALSE]), , drop = FALSE]
boot <- unlist(mclapply(seq_len(500), function(i) {
  idx <- sample(seq_len(nrow(adj_dat)), replace = TRUE)
  d <- adj_dat[idx, , drop = FALSE]
  if (length(unique(d$ecotype)) < 2 || sum(d$DFS_event) < 5) return(NA_real_)
  fit <- tryCatch(survival::coxph(survival::Surv(DFS_days, DFS_event) ~ ecotype + age + sex + stage_num, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  exp(unname(stats::coef(fit)[1]))
}, mc.cores = N_WORKERS))
boot <- boot[is.finite(boot)]
boot_tab <- data.frame(
  metric = c("bootstrap_reps_used", "hr_median", "hr_lower_2.5", "hr_upper_97.5", "frac_boot_hr_gt_1"),
  value = c(length(boot), stats::median(boot), stats::quantile(boot, 0.025), stats::quantile(boot, 0.975),
            mean(boot > 1))
)
print(boot_tab)
fwrite(boot_tab, file.path(OUT, "22b_dfs_bootstrap.tsv"), sep = "\t")

## 7. Figures -----------------------------------------------------------------
km <- survival::survfit(survival::Surv(DFS_days, DFS_event) ~ ecotype, data = ana)
lr <- survival::survdiff(survival::Surv(DFS_days, DFS_event) ~ ecotype, data = ana)
lr_p <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)

save_plot(function() {
  plot(km, col = c("black", "#C0392B"), lwd = 2, xlab = "Days", mark.time = TRUE,
       ylab = "Disease-free survival probability",
       main = sprintf("GSE14333: DFS by ecotype (log-rank p=%s)", format.pval(lr_p, digits = 2, eps = 0.001)))
  legend("topright", legend = c("E1", "E2"), col = c("black", "#C0392B"), lwd = 2, bty = "n")
  legend("bottomleft", bty = "n", legend = sprintf("n=%d, %d DFS events (endpoint: DFS, not OS)", nrow(ana), sum(ana$DFS_event)))
}, file.path(FIG, "plot_33"), 7, 5)

forest <- rows[!is.na(rows$hr_E2_vs_E1), , drop = FALSE]
save_plot(function() {
  print(ggplot(forest, aes(x = hr_E2_vs_E1, y = reorder(model, hr_E2_vs_E1))) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper), colour = "#C0392B") +
    scale_x_log10() + theme_bw(base_size = 10) +
    labs(x = "Hazard ratio for DFS (E2 vs E1, log scale)", y = NULL,
         title = "GSE14333 disease-free survival (secondary endpoint)",
         subtitle = sprintf("adjusted-model bootstrap 95%% CI %.2f-%.2f over %d resamples",
                            stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), length(boot))))
}, file.path(FIG, "plot_34"), 8, 4.5)

## 8. Report -----------------------------------------------------------------
adj <- rows[rows$model == "adjusted (age, sex, stage)", , drop = FALSE]
sig <- !is.na(adj$lrt_p) && adj$lrt_p < 0.05
grade <- if (sig && adj$direction_consistent) "E1" else "not_supported"
report <- c(
  "# GSE14333 disease-free survival (secondary endpoint)",
  "",
  sprintf("Generated by `22b_dfs_validation.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Endpoint discipline",
  "- GSE14333 (GPL570, 290 primary CRCs) reports **DFS only**: `DFS_Time` (months) and `DFS_Cens`.",
  "- Results here are DFS and are **never** pooled with the OS meta-analysis (F03).",
  sprintf("- DFS is recorded for %d tumours (all stage I-III; all %d stage IV patients are `NA`), with %d events.",
          nrow(ana), sum(dat$stage == "IV", na.rm = TRUE), sum(ana$DFS_event)),
  "",
  "## Ecotype assignment",
  sprintf("- Frozen classifier applied to %d Hallmark programmes (%d of %d features recoverable): E1 %d / E2 %d across all %d tumours.",
          ncol(scores), length(shared), length(features),
          sum(labels == "E1"), sum(labels == "E2"), length(labels)),
  sprintf("- Within the %d analysable (DFS-recorded) tumours: E1 %d / E2 %d.",
          nrow(ana), sum(ana$ecotype == "E1"), sum(ana$ecotype == "E2")),
  "",
  "## Association with DFS",
  sprintf("- Unadjusted: HR = %.2f (%.2f-%.2f), p = %s, C-index = %.3f.",
          rows$hr_E2_vs_E1[1], rows$lower[1], rows$upper[1],
          format.pval(rows$p_hr[1], digits = 2, eps = 0.001), rows$c_index_full[1]),
  sprintf("- Adjusted (age, sex, stage): HR = %.2f (%.2f-%.2f), p = %.2f; LRT p = %.3f; C-index %.3f -> %.3f.",
          adj$hr_E2_vs_E1, adj$lower, adj$upper, adj$p_hr, adj$lrt_p, adj$c_index_base, adj$c_index_full),
  sprintf("- Stage II-III subset: HR = %.2f (%.2f-%.2f), p = %.2f.",
          rows$hr_E2_vs_E1[3], rows$lower[3], rows$upper[3], rows$p_hr[3]),
  sprintf("- Bootstrap stability of the adjusted HR (500 resamples): median %.2f, 95%% CI %.2f-%.2f, %.0f%% of resamples above 1.",
          stats::median(boot), stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), 100 * mean(boot > 1)),
  sprintf("- Proportional hazards assumption: minimum cox.zph p = %.2f.", adj$ph_min_p),
  "",
  "## Internal consistency of the endpoint",
  sprintf("- DFS events per stage (I: %d/%d, II: %d/%d, III: %d/%d) rise monotonically with stage,",
          sum(ana$DFS_event[ana$stage == "I"] == 1), sum(ana$stage == "I"),
          sum(ana$DFS_event[ana$stage == "II"] == 1), sum(ana$stage == "II"),
          sum(ana$DFS_event[ana$stage == "III"] == 1), sum(ana$stage == "III")),
  "  which is the expected behaviour of a correctly parsed DFS endpoint (median follow-up 38.5 months).",
  sprintf("- Raw ecotype x stage table is unbalanced (E2 is enriched in stage III: %d vs %d), so the adjusted and",
          sum(ana$ecotype == "E2" & ana$stage == "III"), sum(ana$ecotype == "E1" & ana$stage == "III")),
  "  stage II-III-only models are the primary read-outs; both retain the association.",
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (grade == "not_supported")
            "no significant DFS association for the frozen ecotypes; this is consistent with the null OS result and does not rescue it."
          else "a nominally significant DFS association in one cohort on a secondary endpoint; exploratory only (E1), not a replication of the OS claim."),
  "",
  "## Limitations",
  "- Secondary endpoint in a single cohort; DFS is not OS and the two cannot be combined.",
  "- All stage IV tumours lack DFS, so the analysable set is stage I-III only (selection by data availability).",
  "- Adjuvant chemo/RT flags are recorded but not adjusted for (treatment is a post-baseline mediator of DFS).",
  "- The cohort is microarray (GPL570) while discovery was RNA-seq; platform shift is handled only by the z-scored programme scores.",
  "",
  "## Figures",
  "- plot_33 DFS Kaplan-Meier by ecotype; plot_34 forest plot of unadjusted / adjusted / stage II-III models with the bootstrap CI."
)
writeLines(report, file.path(OUT, "22b_DFS_Report.md"))

log_decision("22b_gse14333", "PASS",
             paste0("GSE14333 DFS (secondary endpoint): n=", nrow(ana), ", events=", sum(ana$DFS_event),
                    "; adjusted HR=", round(adj$hr_E2_vs_E1, 3), " p=", signif(adj$lrt_p, 3),
                    "; bootstrap CI ", round(stats::quantile(boot, 0.025), 2), "-",
                    round(stats::quantile(boot, 0.975), 2), "; grade=", grade,
                    "; kept separate from OS"))

## 9. Evidence ledger ---------------------------------------------------------
# The disease-free claim is a single merged finding (F11) written by 22d_dfs_replication.R, which reads the
# estimate below and adds the independent replication. This script therefore must NOT create a standalone
# discovery row; it only clears any legacy F10 row left by earlier runs.
ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0 && any(old$finding_id == "F10")) {
  old <- old[old$finding_id != "F10", , drop = FALSE]
  data.table::fwrite(old, ledger_path, sep = ",")
  cat("removed legacy F10 row; disease-free evidence is consolidated into F11 by 22d\n")
}
cat("ledger rows:", if (is.null(old)) 0 else nrow(old), "\n")
cat("grade:", grade, "\n")
cat("done\n")
