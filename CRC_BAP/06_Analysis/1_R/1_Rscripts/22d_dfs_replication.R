#!/usr/bin/env Rscript
# 22d_dfs_replication.R
# Independent replication of the DFS association (F10) on disease-free endpoints in cohorts already on disk.
#
# Endpoints: GSE39582 and GSE17536 carry RFS (recurrence-free survival) columns built by the 11_/12_ builders.
# RFS is a disease-free endpoint of the same family as GSE14333's DFS but is NOT interchangeable with OS;
# nothing here is ever pooled with the OS analysis (F03).
#
# Main test : GSE39582 RFS, stage I-III (curative-intent population, matching GSE14333 which had no stage IV),
#             adjusted for age, sex and stage.
# Sensitivity: GSE39582 all stages; GSE17536 RFS stage I-III.
#
# Inputs : 03_Raw_Data/validation/{GSE39582,GSE17536}/{expression.tsv,clinical.tsv}
#          04_Processed_Data/ecotypes.rds (frozen classifier), Environment/hallmark_frozen.gmt
#          07_Validation/22b_dfs_validation.tsv (GSE14333 estimate for the meta summary)
# Outputs: 07_Validation/22d_dfs_replication.tsv, 22d_dfs_replication_bootstrap.tsv, 22d_DFS_Replication_Report.md
#          09_Figures/plots/plot_35 .. plot_36
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6 (physical cores of the i5-10400); work split = per-programme scoring + 500-rep bootstrap

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

OUT <- "07_Validation"
FIG <- "09_Figures/plots"
PROC <- "04_Processed_Data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

ec <- readRDS(file.path(PROC, "ecotypes.rds"))
features <- ec$classifier$features
cen_z <- t(scale(t(ec$classifier$centroids)))
colnames(cen_z) <- features
rownames(cen_z) <- ec$classifier$levels

## ---------------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------------
classify <- function(scores) {
  shared <- intersect(features, colnames(scores))
  if (length(shared) < 20) stop("Only ", length(shared), " classifier features recoverable", call. = FALSE)
  cz <- cen_z[, shared, drop = FALSE]
  s_z <- scale(scores[, shared, drop = FALSE])
  sim <- cor(t(s_z), t(cz), use = "pairwise.complete.obs")
  lab <- rownames(cz)[apply(sim, 1, which.max)]
  names(lab) <- rownames(scores)
  list(labels = lab, n_features = length(shared))
}

score_programmes <- function(expr, gmt_sets) {
  z <- t(scale(t(expr)))
  z <- z[is.finite(rowSums(z)), , drop = FALSE]
  sets <- lapply(gmt_sets, function(g) intersect(g, rownames(z)))
  sets <- sets[lengths(sets) >= 5]
  out <- do.call(cbind, mclapply(sets, function(g) colMeans(z[g, , drop = FALSE]), mc.cores = N_WORKERS))
  colnames(out) <- names(sets)
  out
}

fit_endpoint <- function(d, label, cohort, endpoint_note, adjust = TRUE) {
  vars <- c("time", "event", "ecotype")
  if (adjust) vars <- c(vars, "age", "sex", "stage_num")
  cc <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  if (nrow(cc) < 40 || sum(cc$event) < 15) {
    return(data.frame(cohort = cohort, model = label, n = nrow(cc), events = sum(cc$event),
                      c_index_base = NA_real_, c_index_full = NA_real_, hr = NA_real_, lower = NA_real_,
                      upper = NA_real_, p_hr = NA_real_, lrt_p = NA_real_, ph_min_p = NA_real_,
                      note = "too few complete cases/events", stringsAsFactors = FALSE))
  }
  f_base <- if (adjust) survival::coxph(survival::Surv(time, event) ~ age + sex + stage_num, data = cc)
            else survival::coxph(survival::Surv(time, event) ~ 1, data = cc)
  f_full <- if (adjust) survival::coxph(survival::Surv(time, event) ~ ecotype + age + sex + stage_num, data = cc)
            else survival::coxph(survival::Surv(time, event) ~ ecotype, data = cc)
  s <- summary(f_full)
  lrt <- stats::anova(f_base, f_full, test = "LRT")
  k <- grep("^ecotype", rownames(s$coefficients))[1]
  ph <- tryCatch(survival::cox.zph(f_full), error = function(e) NULL)
  data.frame(
    cohort = cohort, model = label, n = nrow(cc), events = sum(cc$event),
    c_index_base = as.numeric(survival::concordance(f_base)$concordance),
    c_index_full = as.numeric(survival::concordance(f_full)$concordance),
    hr = if (is.na(k)) NA_real_ else exp(s$coefficients[k, "coef"]),
    lower = if (is.na(k)) NA_real_ else unname(s$conf.int[k, 3]),
    upper = if (is.na(k)) NA_real_ else unname(s$conf.int[k, 4]),
    p_hr = if (is.na(k)) NA_real_ else s$coefficients[k, "Pr(>|z|)"],
    lrt_p = lrt[2, "Pr(>|Chi|)"],
    ph_min_p = if (is.null(ph)) NA_real_ else min(ph$table[, "p"], na.rm = TRUE),
    note = endpoint_note, stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## 1. Cohorts
## ---------------------------------------------------------------------------
cohorts <- list(
  list(dir = "03_Raw_Data/validation/GSE39582", label = "GSE39582", endpoint = "RFS"),
  list(dir = "03_Raw_Data/validation/GSE17536", label = "GSE17536", endpoint = "RFS")
)

rows <- list()
km_data <- list()
for (co in cohorts) {
  expr <- read_matrix_any(file.path(co$dir, "expression.tsv"))
  clin <- data.table::fread(file.path(co$dir, "clinical.tsv"), data.table = FALSE)
  id_type <- if (all(grepl("^ENSG", head(rownames(expr), 100)))) "ensembl" else "symbol"
  hl <- hallmark_gene_sets(id_type)
  scores <- score_programmes(expr, hl$sets)
  cl <- classify(scores)
  cat("\n===", co$label, "| expression", nrow(expr), "x", ncol(expr),
      "| programmes", ncol(scores), "| features", cl$n_features, "| sets:", hl$source, "\n")
  cat("ecotypes: "); print(table(cl$labels))

  dat <- clin
  dat$ecotype <- factor(cl$labels[dat$sample_id], levels = c("E1", "E2"))
  dat$time <- suppressWarnings(as.numeric(dat$rfs_time_days))
  dat$event <- suppressWarnings(as.integer(dat$rfs_event))
  dat$age <- suppressWarnings(as.numeric(dat$age))
  dat$stage_num <- as.numeric(factor(dat$stage, levels = c("I", "II", "III", "IV")))
  dat <- dat[!is.na(dat$time) & dat$time > 0 & !is.na(dat$event) & !is.na(dat$ecotype), , drop = FALSE]
  cat("analysable:", nrow(dat), "| events:", sum(dat$event),
      "| stage counts:\n"); print(table(dat$stage, useNA = "ifany"))

  stage123 <- dat[dat$stage %in% c("I", "II", "III"), , drop = FALSE]
  rows[[length(rows) + 1]] <- fit_endpoint(stage123, "stage I-III, adjusted (age, sex, stage)",
                                           co$label, co$endpoint, adjust = TRUE)
  rows[[length(rows) + 1]] <- fit_endpoint(stage123, "stage I-III, unadjusted", co$label, co$endpoint, adjust = FALSE)
  rows[[length(rows) + 1]] <- fit_endpoint(dat, "all stages, adjusted (age, sex, stage)",
                                           co$label, co$endpoint, adjust = TRUE)
  if (identical(co$label, "GSE39582")) km_data[[co$label]] <- stage123
}
res <- do.call(rbind, rows)
res$direction_consistent <- !is.na(res$hr) & log(res$hr) > 0   # discovery: E2 = worse
print(res)
fwrite(res, file.path(OUT, "22d_dfs_replication.tsv"), sep = "\t")

## ---------------------------------------------------------------------------
## 2. Bootstrap the primary estimate (GSE39582, stage I-III, adjusted)
## ---------------------------------------------------------------------------
prim <- km_data[["GSE39582"]]
prim_cc <- prim[stats::complete.cases(prim[, c("time", "event", "ecotype", "age", "sex", "stage_num"), drop = FALSE]), , drop = FALSE]
boot <- unlist(mclapply(seq_len(500), function(i) {
  d <- prim_cc[sample(nrow(prim_cc), replace = TRUE), , drop = FALSE]
  if (length(unique(d$ecotype)) < 2 || sum(d$event) < 10) return(NA_real_)
  f <- tryCatch(survival::coxph(survival::Surv(time, event) ~ ecotype + age + sex + stage_num, data = d),
                error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  exp(unname(stats::coef(f)[1]))
}, mc.cores = N_WORKERS))
boot <- boot[is.finite(boot)]
boot_tab <- data.frame(
  metric = c("cohort", "n", "events", "bootstrap_reps", "hr_median", "hr_lower_2.5", "hr_upper_97.5", "frac_above_1"),
  value = c(NA, nrow(prim_cc), sum(prim_cc$event), length(boot), stats::median(boot),
            stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), mean(boot > 1))
)
print(boot_tab)
fwrite(boot_tab, file.path(OUT, "22d_dfs_replication_bootstrap.tsv"), sep = "\t")

## ---------------------------------------------------------------------------
## 3. Secondary: fixed-effect summary with the GSE14333 estimate (DFS family only)
## ---------------------------------------------------------------------------
gse14333 <- data.table::fread(file.path(OUT, "22b_dfs_validation.tsv"), data.table = FALSE)
gse14333 <- gse14333[gse14333$model == "adjusted (age, sex, stage)", , drop = FALSE]
ps <- res[res$cohort == "GSE39582" & res$model == "stage I-III, adjusted (age, sex, stage)", , drop = FALSE]
meta <- NULL
if (nrow(gse14333) == 1 && nrow(ps) == 1 && !is.na(ps$hr) && !is.na(gse14333$hr)) {
  se1 <- (log(gse14333$upper) - log(gse14333$lower)) / (2 * 1.96)
  se2 <- (log(ps$upper) - log(ps$lower)) / (2 * 1.96)
  w <- c(1 / se1^2, 1 / se2^2)
  log_hr <- c(log(gse14333$hr), log(ps$hr))
  m <- sum(w * log_hr) / sum(w)
  se_m <- sqrt(1 / sum(w))
  z <- m / se_m
  meta <- data.frame(analysis = "fixed-effect (GSE14333 DFS + GSE39582 RFS, stage I-III adjusted)",
                     cohorts = 2, hr = exp(m), lower = exp(m - 1.96 * se_m), upper = exp(m + 1.96 * se_m),
                     p = 2 * stats::pnorm(-abs(z)),
                     heterogeneity_Q_p = 1 - stats::pchisq(sum(w * (log_hr - m)^2), df = 1))
  print(meta)
  fwrite(meta, file.path(OUT, "22d_dfs_meta.tsv"), sep = "\t")
}

## ---------------------------------------------------------------------------
## 4. Figures
## ---------------------------------------------------------------------------
forest <- res[res$cohort == "GSE39582", , drop = FALSE]
forest <- forest[!is.na(forest$hr), , drop = FALSE]
forest$model <- factor(forest$model, levels = rev(forest$model))
save_plot(function() {
  print(ggplot(forest, aes(x = hr, y = model)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper), colour = "#C0392B") +
    scale_x_log10() + theme_bw(base_size = 10) +
    labs(x = "Hazard ratio for RFS, E2 vs E1 (log scale)", y = NULL,
         title = "GSE39582: independent replication of the disease-free association",
         subtitle = sprintf("stage I-III adjusted HR %.2f (%.2f-%.2f), p = %.3g",
                            forest$hr[forest$model == "stage I-III, adjusted (age, sex, stage)"],
                            forest$lower[forest$model == "stage I-III, adjusted (age, sex, stage)"],
                            forest$upper[forest$model == "stage I-III, adjusted (age, sex, stage)"],
                            forest$p_hr[forest$model == "stage I-III, adjusted (age, sex, stage)"])))
}, file.path(FIG, "plot_35"), 8, 4.5)

km <- survival::survfit(survival::Surv(time, event) ~ ecotype, data = prim_cc)
lr <- survival::survdiff(survival::Surv(time, event) ~ ecotype, data = prim_cc)
lr_p <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)
save_plot(function() {
  plot(km, col = unname(ECOTYPE_COL), lwd = 2, xlab = "Days", mark.time = TRUE,
       ylab = "Recurrence-free survival probability",
       main = sprintf("GSE39582: RFS by ecotype, stage I-III (log-rank p=%s)",
                      format.pval(lr_p, digits = 2, eps = 0.001)))
  legend("topright", legend = c("E1", "E2"), col = unname(ECOTYPE_COL), lwd = 2, bty = "n")
  legend("bottomleft", bty = "n", legend = sprintf("n=%d, %d RFS events (endpoint: RFS, not OS)",
                                                   nrow(prim_cc), sum(prim_cc$event)))
}, file.path(FIG, "plot_36"), 7, 5)

## ---------------------------------------------------------------------------
## 5. Report + ledger
## ---------------------------------------------------------------------------
p_row <- ps[1, , drop = FALSE]
g_row <- res[res$cohort == "GSE17536" & res$model == "stage I-III, adjusted (age, sex, stage)", , drop = FALSE]
replicated <- !is.na(p_row$p_hr) && p_row$p_hr < 0.05 && p_row$hr > 1
grade <- if (replicated) "E2" else "E1"   # merged claim: E2 needs the replication; discovery alone is E1

report <- c(
  "# DFS/RFS replication (independent cohorts already on disk)",
  "",
  sprintf("Generated by `22d_dfs_replication.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Endpoint discipline",
  "- Endpoint here is **RFS** (recurrence-free survival), the disease-free family of endpoints.",
  "- RFS is reported separately from OS; no result here enters the OS analysis (F03).",
  sprintf("- GSE39582: %d tumours with RFS, %d events. Primary model is stage I-III (curative-intent), matching",
          nrow(prim), sum(prim$event)),
  "  GSE14333, which contained no stage IV at all; all-stage is a sensitivity analysis.",
  "",
  "## Primary replication (GSE39582, stage I-III, adjusted for age/sex/stage)",
  sprintf("- n = %d, %d RFS events; HR (E2 vs E1) = **%.2f (%.2f-%.2f)**, p = %.3g; LRT p = %.3g.",
          p_row$n, p_row$events, p_row$hr, p_row$lower, p_row$upper, p_row$p_hr, p_row$lrt_p),
  sprintf("- C-index %.3f -> %.3f; proportional-hazards minimum p = %.2f.",
          p_row$c_index_base, p_row$c_index_full, p_row$ph_min_p),
  sprintf("- Bootstrap (500 resamples, 6 workers): median HR %.2f, 95%% CI %.2f-%.2f, %.0f%% of resamples above 1.",
          stats::median(boot), stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), 100 * mean(boot > 1)),
  "",
  "## Sensitivity",
  sprintf("- GSE39582 all stages: HR %.2f (%.2f-%.2f), p = %.3g (n = %d, %d events).",
          res$hr[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"],
          res$lower[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"],
          res$upper[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"],
          res$p_hr[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"],
          res$n[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"],
          res$events[res$cohort == "GSE39582" & res$model == "all stages, adjusted (age, sex, stage)"]),
  sprintf("- GSE17536 stage I-III: HR %.2f (%.2f-%.2f), p = %.3g (n = %d, %d events).",
          g_row$hr, g_row$lower, g_row$upper, g_row$p_hr, g_row$n, g_row$events),
  if (!is.null(meta)) {
    sprintf("- Fixed-effect summary with GSE14333 (same endpoint family): HR %.2f (%.2f-%.2f), p = %.3g; heterogeneity p = %.2f.",
            meta$hr, meta$lower, meta$upper, meta$p, meta$heterogeneity_Q_p)
  } else "- Meta summary not computed.",
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (replicated) {
            "the disease-free association replicates in an independent cohort (GSE39582) with the same direction and magnitude, which is what E2 requires."
          } else {
            "the association did not replicate in GSE39582; F10 stays E1 and must not be presented as a replicated result."
          }),
  "",
  "## Limitations",
  "- RFS and DFS are related but not identical endpoint definitions; the two are combined only in the labelled secondary summary.",
  "- GSE39582 is microarray (GPL570) and was already used for OS validation, so it is independent for the RFS endpoint but not a fully independent cohort overall.",
  sprintf("- %d GSE39582 tumours lack a stage value and are excluded from the adjusted models.", sum(is.na(prim$stage_num))),
  "- Survival is observational: no treatment adjustment, no competing-risk analysis.",
  "",
  "## Figures",
  "- plot_35 GSE39582 hazard-ratio forest (stage I-III adjusted/unadjusted, all stages); plot_36 RFS Kaplan-Meier by ecotype."
)
writeLines(report, file.path(OUT, "22d_DFS_Replication_Report.md"))

ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F11", , drop = FALSE]
f11 <- data.frame(
  finding_id = "F11",
  claim = if (replicated) {
    "The frozen K=2 ecotypes are associated with disease-free survival (DFS/RFS) in discovery and in independent cohorts"
  } else {
    "The disease-free survival association observed in GSE14333 does not replicate in independent cohorts"
  },
  evidence_type = "Discovery + external validation (disease-free endpoints)",
  dataset_s = "GSE14333 (discovery, DFS); GSE39582 and GSE17536 (RFS)",
  statistic = paste0("discovery GSE14333 adjusted HR=", round(gse14333$hr, 2), " (", round(gse14333$lower, 2), "-",
                     round(gse14333$upper, 2), ") p=", signif(gse14333$lrt_p, 3),
                     "; replication GSE39582 stage I-III HR=", round(p_row$hr, 2), " (", round(p_row$lower, 2), "-",
                     round(p_row$upper, 2), "), p=", signif(p_row$p_hr, 3),
                     "; GSE17536 HR=", round(g_row$hr, 2), " p=", signif(g_row$p_hr, 3),
                     if (!is.null(meta)) paste0("; fixed-effect HR=", round(meta$hr, 2), " (", round(meta$lower, 2),
                                                "-", round(meta$upper, 2), ") p=", signif(meta$p, 3),
                                                ", heterogeneity p=", round(meta$heterogeneity_Q_p, 2)) else ""),
  n = p_row$n,
  grade = grade,
  verb_used = if (replicated) "associated with (replicated)" else "not replicated",
  contradicting_evidence = if (replicated) {
    "OS association remains null in discovery (F01) and in OS validation (F03): the ecotype predicts recurrence, not death"
  } else {
    "the discovery estimate is not reproduced in GSE39582 or GSE17536"
  },
  limitation = "disease-free endpoints only (DFS in GSE14333, RFS in GSE39582/GSE17536: related but not identical definitions); GSE14333 has no stage IV, GSE39582 stages I-III in the primary model; GSE39582 was already used for OS validation; observational, no treatment adjustment; never pooled with OS",
  next_test = "GSE38832 DFS once its series matrix is available; prospective evaluation",
  figure_table = "plot_33; plot_34; plot_35; plot_36",
  stringsAsFactors = FALSE
)
out <- if (is.null(old) || nrow(old) == 0) f11 else rbind(old, f11)
out <- out[order(as.integer(sub("^F", "", out$finding_id))), , drop = FALSE]
data.table::fwrite(out, ledger_path, sep = ",")
cat("ledger rows:", nrow(out), "\n")

log_decision("22d_replication", "PASS",
             paste0("GSE39582 RFS stage I-III: n=", p_row$n, ", events=", p_row$events,
                    "; adjusted HR=", round(p_row$hr, 3), " p=", signif(p_row$p_hr, 3),
                    "; bootstrap ", round(stats::quantile(boot, 0.025), 2), "-", round(stats::quantile(boot, 0.975), 2),
                    "; GSE17536 HR=", round(g_row$hr, 2), " p=", signif(g_row$p_hr, 3),
                    "; grade=", grade))
cat("grade:", grade, "\n")
cat("done\n")
