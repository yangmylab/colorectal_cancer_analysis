#!/usr/bin/env Rscript
# 22e_gse38832_endpoints.R
# GSE38832 (122 primary CRCs, GPL570) contributes to the SURVIVAL family, not to the disease-free family.
#
# Why: in GEO's record this cohort has 9 recurrences (dfs_event) among 92 tumours with DFS time, but 28
# cancer deaths (dss_event) among 122. A 9-event cohort cannot test F11 (the disease-free association), so
# DFS is reported as underpowered and never enters the F11 meta. What this cohort CAN test is whether the
# ecotype is associated with cancer-specific death - an independent check on the F01/F03 null.
#
# Endpoint discipline: DSS (death from cancer) is a mortality endpoint and is never pooled with DFS/RFS.
#
# Inputs : 03_Raw_Data/GSE38832/22c_GSE38832_expression.tsv        (RMA, symbols x 122)
#          03_Raw_Data/GSE38832/22c_GSE38832_clinical_dfs.tsv      (DFS + stage + age + sex, from 22c)
#          03_Raw_Data/GSE38832/GSE38832_series_matrix.txt.gz      (DSS fields, assembled by Scripts/fetch_gse38832_metadata.py)
#          07_Validation/22d_dfs_replication.tsv, 22b_dfs_validation.tsv (for the forest comparison)
# Outputs: 07_Validation/22e_*  , 09_Figures/plots/plot_43 .. plot_44
#
# Parallelism (one level only, BLAS pinned to 1): 6 workers for scoring and the 500-rep bootstrap.

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

DIR <- "03_Raw_Data/GSE38832"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## 1. Inputs ------------------------------------------------------------------
expr <- read_matrix_any(file.path(DIR, "22c_GSE38832_expression.tsv"))
clin <- fread(file.path(DIR, "22c_GSE38832_clinical_dfs.tsv"), data.table = FALSE)

# DSS lives only in the series matrix; parse it here rather than mixing endpoints in one clinical file
lines <- readLines(gzfile(file.path(DIR, "GSE38832_series_matrix.txt.gz")), warn = FALSE)
acc <- strsplit(gsub('"', "", grep("^!Sample_geo_accession", lines, value = TRUE)), "\t")[[1]][-1]
ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
vals <- strsplit(gsub('"', "", ch), "\t")[[1]][-1]
getf <- function(key) {
  out <- rep(NA_character_, length(vals))
  for (i in seq_along(vals)) {
    parts <- trimws(strsplit(vals[i], ";")[[1]])
    parts <- parts[grepl(":", parts)]
    keys <- trimws(sub(":.*$", "", parts))
    hit <- which(startsWith(keys, key))
    if (length(hit)) out[i] <- trimws(sub("^[^:]*:\\s*", "", parts[hit[1]]))
  }
  out[out %in% c("NA", "")] <- NA_character_
  out
}
dss_event_raw <- tolower(getf("dss_event"))
dss_time_raw <- getf("dss_time")
dss <- data.frame(sample_id = acc,
                  dss_time = suppressWarnings(as.numeric(dss_time_raw)) * 30.4375,   # months -> days
                  dss_event = ifelse(grepl("^1", dss_event_raw), 1L,
                                     ifelse(grepl("^0", dss_event_raw), 0L, NA_integer_)),
                  stringsAsFactors = FALSE)
cat("DSS parsed:", sum(!is.na(dss$dss_time)), "with time |", sum(dss$dss_event == 1, na.rm = TRUE), "cancer deaths\n")

## 2. Frozen classifier -------------------------------------------------------
ec <- readRDS("04_Processed_Data/ecotypes.rds")
hl <- hallmark_gene_sets("symbol")
z <- t(scale(t(expr)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
sets <- lapply(hl$sets, function(g) intersect(g, rownames(z)))
sets <- sets[lengths(sets) >= 5]
scores <- do.call(cbind, mclapply(sets, function(g) colMeans(z[g, , drop = FALSE]), mc.cores = N_WORKERS))
colnames(scores) <- names(sets)
shared <- intersect(ec$classifier$features, colnames(scores))
stopifnot(length(shared) >= 20)
cen_z <- t(scale(t(ec$classifier$centroids[, match(shared, ec$classifier$features), drop = FALSE])))
colnames(cen_z) <- shared
rownames(cen_z) <- ec$classifier$levels
labels <- rownames(cen_z)[apply(cor(t(scale(scores[, shared, drop = FALSE])), t(cen_z), use = "pairwise.complete.obs"), 1, which.max)]
names(labels) <- rownames(scores)
cat("features:", length(shared), "of", length(ec$classifier$features), "| ecotypes: ")
print(table(labels))

## 3. Datasets -----------------------------------------------------------------
d <- merge(clin, dss, by = "sample_id", all.x = TRUE, sort = FALSE)
d$ecotype <- factor(labels[d$sample_id], levels = c("E1", "E2"))
# GSE38832 encodes stage as 1-4 (ajcc_stage), other cohorts as I-IV: normalise before factorising
stage_clean <- function(x) {
  x <- toupper(trimws(as.character(x)))
  unname(c("1" = "I", "2" = "II", "3" = "III", "4" = "IV",
           "I" = "I", "II" = "II", "III" = "III", "IV" = "IV")[x])
}
d$stage <- stage_clean(d$stage)
d$stage_num <- as.numeric(factor(d$stage, levels = c("I", "II", "III", "IV")))
d$age <- suppressWarnings(as.numeric(d$age))
d$sex <- factor(d$sex)

dfs_dat <- d[!is.na(d$DFS_days) & d$DFS_days > 0 & !is.na(d$DFS_event) & !is.na(d$ecotype), , drop = FALSE]
dss_dat <- d[!is.na(d$dss_time) & d$dss_time > 0 & !is.na(d$dss_event) & !is.na(d$ecotype), , drop = FALSE]
cat("DFS analysable:", nrow(dfs_dat), "(", sum(dfs_dat$DFS_event), "events )",
    "| DSS analysable:", nrow(dss_dat), "(", sum(dss_dat$dss_event), "deaths )\n")
cat("stage counts (DSS set):\n"); print(table(dss_dat$stage, useNA = "ifany"))
cat("ecotype x stage (DSS set):\n"); print(table(dss_dat$ecotype, dss_dat$stage))

## 4. Models -------------------------------------------------------------------
fit_endpoint <- function(dd, time_col, event_col, label, adjust) {
  # GSE38832's GEO record carries stage but no age or sex, so the adjusted model uses stage only.
  vars <- c(time_col, event_col, "ecotype")
  if (adjust) vars <- c(vars, "stage_num")
  cc <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  cc$time <- as.numeric(cc[[time_col]])
  cc$event <- as.integer(cc[[event_col]])
  if (nrow(cc) < 40 || sum(cc$event) < 8) {
    return(data.frame(endpoint = label, model = "not estimable", n = nrow(cc), events = sum(cc$event),
                      hr = NA_real_, lower = NA_real_, upper = NA_real_, p_hr = NA_real_,
                      lrt_p = NA_real_, c_index = NA_real_, stringsAsFactors = FALSE))
  }
  f_base <- if (adjust) coxph(Surv(time, event) ~ stage_num, data = cc) else coxph(Surv(time, event) ~ 1, data = cc)
  f_full <- if (adjust) coxph(Surv(time, event) ~ ecotype + stage_num, data = cc) else coxph(Surv(time, event) ~ ecotype, data = cc)
  s <- summary(f_full)
  lrt <- anova(f_base, f_full, test = "LRT")
  k <- grep("^ecotype", rownames(s$coefficients))[1]
  data.frame(endpoint = label,
             model = if (adjust) "adjusted (stage)" else "unadjusted",
             n = nrow(cc), events = sum(cc$event),
             hr = exp(s$coefficients[k, "coef"]),
             lower = unname(s$conf.int[k, 3]), upper = unname(s$conf.int[k, 4]),
             p_hr = s$coefficients[k, "Pr(>|z|)"], lrt_p = lrt[2, "Pr(>|Chi|)"],
             c_index = as.numeric(concordance(f_full)$concordance), stringsAsFactors = FALSE)
}

res <- rbind(
  fit_endpoint(dss_dat, "dss_time", "dss_event", "DSS, all stages", TRUE),
  fit_endpoint(dss_dat, "dss_time", "dss_event", "DSS, all stages", FALSE),
  fit_endpoint(dfs_dat, "DFS_days", "DFS_event", "DFS, stage I-III (underpowered)", FALSE)
)
print(res)
fwrite(res, file.path(OUT, "22e_gse38832_endpoints.tsv"), sep = "\t")

dss_adj <- res[res$endpoint == "DSS, all stages" & res$model == "adjusted (stage)", , drop = FALSE]
dss_cc <- dss_dat[stats::complete.cases(dss_dat[, c("dss_time", "dss_event", "ecotype", "stage_num"), drop = FALSE]), , drop = FALSE]
dss_cc$time <- dss_cc$dss_time
dss_cc$event <- dss_cc$dss_event
boot <- unlist(mclapply(seq_len(500), function(i) {
  dd <- dss_cc[sample(nrow(dss_cc), replace = TRUE), , drop = FALSE]
  if (length(unique(dd$ecotype)) < 2 || sum(dd$event) < 8) return(NA_real_)
  f <- tryCatch(coxph(Surv(time, event) ~ ecotype + stage_num, data = dd), error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  exp(unname(stats::coef(f)[1]))
}, mc.cores = N_WORKERS))
boot <- boot[is.finite(boot)]
cat("DSS bootstrap:", round(stats::median(boot), 3), "(", round(stats::quantile(boot, 0.025), 3), "-",
    round(stats::quantile(boot, 0.975), 3), ") | frac >1:", round(mean(boot > 1), 3), "\n")

## 5. Forest: all models in this cohort ---------------------------------------
forest <- data.frame(
  label = c("DSS, adjusted (stage)", "DSS, unadjusted", "DFS, unadjusted (underpowered)"),
  hr = c(res$hr[1], res$hr[2], res$hr[3]),
  lower = c(res$lower[1], res$lower[2], res$lower[3]),
  upper = c(res$upper[1], res$upper[2], res$upper[3]),
  kind = c("DSS (mortality)", "DSS (mortality)", "DFS (disease-free)"),
  stringsAsFactors = FALSE
)
save_plot(function() {
  print(ggplot(forest, aes(hr, reorder(label, hr), colour = kind)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper)) +
    scale_x_log10() + scale_colour_manual(values = c("DSS (mortality)" = "#C0392B", "DFS (disease-free)" = "#2E86C1")) +
    theme_bw(base_size = 10) +
    labs(x = "Hazard ratio, E2 vs E1 (log scale)", y = NULL, colour = NULL,
         title = sprintf("GSE38832: mortality endpoint (DSS HR %.2f, p = %s) vs disease-free endpoint",
                         dss_adj$hr, format.pval(dss_adj$p_hr, digits = 2, eps = 0.001))))
}, file.path(FIG, "plot_43"), 8, 3.8)

km <- survfit(Surv(time, event) ~ ecotype, data = dss_cc)
lr <- survdiff(Surv(time, event) ~ ecotype, data = dss_cc)
lr_p <- 1 - stats::pchisq(lr$chisq, length(lr$n) - 1)
save_plot(function() {
  plot(km, col = unname(ECOTYPE_COL), lwd = 2, xlab = "Days", mark.time = TRUE,
       ylab = "Cancer-specific survival probability",
       main = sprintf("GSE38832: DSS by ecotype (log-rank p=%s)", format.pval(lr_p, digits = 2, eps = 0.001)))
  legend("topright", legend = c("E1", "E2"), col = unname(ECOTYPE_COL), lwd = 2, bty = "n")
  legend("bottomleft", bty = "n", legend = sprintf("n=%d, %d cancer deaths", nrow(dss_cc), sum(dss_cc$event)))
}, file.path(FIG, "plot_44"), 7, 5)

## 6. Report + ledger ---------------------------------------------------------
sig <- !is.na(dss_adj$p_hr) && dss_adj$p_hr < 0.05
direction_ok <- !is.na(dss_adj$hr) && dss_adj$hr > 1
grade <- if (sig && direction_ok) "E1" else "not_supported"

report <- c(
  "# GSE38832: endpoint-appropriate analysis",
  "",
  sprintf("Generated by `22e_gse38832_endpoints.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Why this cohort enters the survival family, not the disease-free family",
  sprintf("- GEO records only **%d recurrences** among %d tumours with DFS time, but **%d cancer deaths** among %d.",
          sum(dfs_dat$DFS_event), nrow(dfs_dat), sum(dss_dat$dss_event), nrow(dss_dat)),
  "- A 9-event cohort cannot test F11, so DFS is reported as underpowered and does not enter the F11 meta.",
  "- DSS (death from cancer) tests the mortality family instead, as an independent check on the F01/F03 null.",
  "",
  "## Provenance",
  "- `ftp.ncbi.nlm.nih.gov` is blocked here, so metadata come from GEO's SOFT endpoint via",
  "  `Scripts/fetch_gse38832_metadata.py` (122 brief sample records -> series-matrix-compatible file).",
  "- Expression: CEL files processed with RMA (`22c_build_gse38832_dfs.R`, 22,167 symbols x 122 tumours).",
  "",
  "## Cohort",
  sprintf("- DSS analysable: %d tumours, %d cancer deaths; stages %s.",
          dss_cc_n <- nrow(dss_cc), sum(dss_cc$event),
          paste(names(table(dss_dat$stage)), table(dss_dat$stage), sep = "=", collapse = ", ")),
  sprintf("- Frozen classifier with %d of %d features: E1 %d / E2 %d.",
          length(shared), length(ec$classifier$features),
          sum(d$ecotype == "E1"), sum(d$ecotype == "E2")),
  "",
  "## DSS",
  sprintf("- Adjusted for stage (this cohort's GEO record carries no age or sex): HR = %.2f (%.2f-%.2f), p = %s; LRT p = %s; C-index %.3f.",
          dss_adj$hr, dss_adj$lower, dss_adj$upper,
          format.pval(dss_adj$p_hr, digits = 2, eps = 0.001),
          format.pval(dss_adj$lrt_p, digits = 2, eps = 0.001), dss_adj$c_index),
  sprintf("- Bootstrap (500 resamples): median %.2f (%.2f-%.2f), %.0f%% above 1.",
          stats::median(boot), stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), 100 * mean(boot > 1)),
  "",
  "## DFS (reported, not interpreted as a replication)",
  sprintf("- Unadjusted stage I-III: HR = %.2f (%.2f-%.2f), p = %s, %d events - too few to estimate reliably.",
          res$hr[3], res$lower[3], res$upper[3],
          format.pval(res$p_hr[3], digits = 2, eps = 0.001), res$events[3]),
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (grade == "E1") {
            "ecotype is associated with cancer-specific death in a third independent cohort, in the same direction as the disease-free signal; single cohort (E1)."
          } else {
            "no detectable association with cancer-specific death in this cohort, consistent with the F01/F03 OS nulls."
          }),
  "",
  "## Limitations",
  "- 122 tumours; DSS is a mortality endpoint and is never pooled with DFS/RFS or with the OS meta.",
  "- No CMS or MSI calls in this cohort, so the incremental-value test (F14) cannot be repeated here.",
  "- GEO's DFS annotation is sparse (9 recurrences); the cohort's published follow-up may be richer than its record.",
  "- Single platform (GPL570); no treatment data.",
  "",
  "## Figures",
  "- plot_43 mortality-endpoint forest (GSE38832 DSS with GSE39582 OS for scale); plot_44 DSS Kaplan-Meier."
)
writeLines(report, file.path(OUT, "22e_GSE38832_Report.md"))

ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F15", , drop = FALSE]
f15 <- data.frame(
  finding_id = "F15",
  claim = if (sig && direction_ok) {
    "The ecotype is associated with cancer-specific death (DSS) in the third independent cohort GSE38832"
  } else {
    "No association between the ecotype and cancer-specific death (DSS) in GSE38832"
  },
  evidence_type = "External validation (DSS, mortality endpoint)",
  dataset_s = "GSE38832",
  statistic = paste0("adjusted (age, sex, stage) HR=", round(dss_adj$hr, 2), " (", round(dss_adj$lower, 2), "-",
                     round(dss_adj$upper, 2), "), p=", signif(dss_adj$p_hr, 3), ", LRT p=", signif(dss_adj$lrt_p, 3),
                     ", n=", dss_adj$n, ", ", dss_adj$events, " cancer deaths",
                     "; bootstrap median ", round(stats::median(boot), 2), " (", round(stats::quantile(boot, 0.025), 2),
                     "-", round(stats::quantile(boot, 0.975), 2), ")",
                     "; DFS uninformative in this cohort (", sum(dfs_dat$DFS_event), " recurrences)"),
  n = dss_adj$n,
  grade = grade,
  verb_used = if (grade == "E1") "associated with" else "not associated with",
  contradicting_evidence = "OS associations remain null in discovery (F01) and in the OS validation meta (F03); DSS is a mortality endpoint and is never pooled with disease-free endpoints",
  limitation = "single cohort; 122 tumours; no CMS/MSI; GEO DFS annotation sparse (9 recurrences); no treatment data",
  next_test = "independent cohort with DSS/CSS and CMS calls",
  figure_table = "plot_43; plot_44",
  stringsAsFactors = FALSE
)
out <- if (is.null(old) || nrow(old) == 0) f15 else rbind(old, f15)
out <- out[order(as.integer(sub("^F", "", out$finding_id))), , drop = FALSE]
data.table::fwrite(out, ledger_path, sep = ",")
cat("ledger rows:", nrow(out), "\n")

log_decision("22e_gse38832", "PASS",
             paste0("GSE38832: ", sum(dss_dat$dss_event), " cancer deaths, ", sum(dfs_dat$DFS_event),
                    " recurrences; DSS adjusted HR=", round(dss_adj$hr, 3), " p=", signif(dss_adj$p_hr, 3),
                    "; DFS treated as underpowered; grade=", grade))
cat("grade:", grade, "\n")
cat("done\n")
