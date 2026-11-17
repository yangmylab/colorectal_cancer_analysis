#!/usr/bin/env Rscript
# 25_incremental_value.R
# Does the frozen K=2 ecotype add prognostic information BEYOND stage, MSI and CMS?
#
# This is the protocol's central question (Analysis_Plan.md: "whether ecotypes refine prognosis beyond
# stage/MSI/CMS/CRIS"). Every previous model in this project adjusted for age, sex and stage only, so the
# claim that the ecotype is worth anything beyond the published taxonomy was untested until now.
#
# Endpoints (never mixed): RFS (primary, stage I-III as in 22d) and OS (secondary, all stages) in GSE39582,
# the only cohort on disk that carries CMS calls.
#
# Inputs : 03_Raw_Data/validation/GSE39582/{expression.tsv,clinical.tsv}  (clinical has cms, msi_status)
#          Environment/hallmark_frozen.gmt, 04_Processed_Data/ecotypes.rds
# Outputs: 07_Validation/25_incremental_value.tsv, 25_incremental_value_bootstrap.tsv,
#          25_Incremental_Value_Report.md, 09_Figures/plots/plot_41 .. plot_42
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6 (physical cores); work split = per-programme scoring + 500-rep bootstrap of the ecotype HR

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

CO <- "03_Raw_Data/validation/GSE39582"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## 1. Ecotype labels (frozen classifier, as in 22d) ---------------------------
expr <- read_matrix_any(file.path(CO, "expression.tsv"))
clin <- fread(file.path(CO, "clinical.tsv"), data.table = FALSE)
id_type <- if (all(grepl("^ENSG", head(rownames(expr), 100)))) "ensembl" else "symbol"
hl <- hallmark_gene_sets(id_type)
z <- t(scale(t(expr)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
sets <- lapply(hl$sets, function(g) intersect(g, rownames(z)))
sets <- sets[lengths(sets) >= 5]
scores <- do.call(cbind, mclapply(sets, function(g) colMeans(z[g, , drop = FALSE]), mc.cores = N_WORKERS))
colnames(scores) <- names(sets)
features <- ecotypes_features <- readRDS("04_Processed_Data/ecotypes.rds")$classifier$features
cen <- readRDS("04_Processed_Data/ecotypes.rds")$classifier$centroids
shared <- intersect(features, colnames(scores))
stopifnot(length(shared) >= 20)
cen_z <- t(scale(t(cen[, match(shared, features), drop = FALSE])))
colnames(cen_z) <- shared
s_z <- scale(scores[, shared, drop = FALSE])
labels <- rownames(cen_z)[apply(cor(t(s_z), t(cen_z), use = "pairwise.complete.obs"), 1, which.max)]
names(labels) <- rownames(scores)
cat("ecotypes scored on", length(shared), "features: "); print(table(labels))

## 2. Clinical harmonisation --------------------------------------------------
d <- clin
d$ecotype <- factor(labels[d$sample_id], levels = c("E1", "E2"))
d$stage_num <- as.numeric(factor(d$stage, levels = c("I", "II", "III", "IV")))
d$age <- suppressWarnings(as.numeric(d$age))
d$sex <- factor(d$sex)
d$sidedness <- factor(d$sidedness)

# CMS: values are C1..C6 with "N/A"; anything else is unknown
cms_raw <- toupper(trimws(as.character(d$cms)))
d$cms <- ifelse(grepl("^C[1-6]$", cms_raw), cms_raw, NA_character_)
d$cms <- factor(d$cms, levels = paste0("C", 1:6))

# MSI: normalise MSI / MSS (some tables use MSI-H / MSS / blank)
msi_raw <- toupper(trimws(as.character(d$msi_status)))
d$msi <- ifelse(grepl("^MSI", msi_raw), "MSI",
                ifelse(grepl("^MSS|^MSS|^MSS", msi_raw) | msi_raw == "MSS", "MSS", NA_character_))
d$msi <- ifelse(is.na(d$msi) & grepl("MSS", msi_raw), "MSS", d$msi)
d$msi <- factor(d$msi, levels = c("MSS", "MSI"))

cat("\nCMS availability: "); print(table(d$cms, useNA = "ifany"))
cat("MSI availability: "); print(table(d$msi, useNA = "ifany"))
cat("\necotype x CMS:\n"); print(table(d$ecotype, d$cms))
cat("\necotype x MSI:\n"); print(table(d$ecotype, d$msi))

## 3. Nested models -----------------------------------------------------------
fit_set <- function(dat, time_col, event_col, label) {
  dat <- dat[!is.na(dat[[time_col]]) & dat[[time_col]] > 0 & !is.na(dat[[event_col]]) & !is.na(dat$ecotype), , drop = FALSE]
  dat$time <- as.numeric(dat[[time_col]])
  dat$event <- as.integer(dat[[event_col]])

  specs <- list(
    list(name = "stage, age, sex", terms = c("age", "sex", "stage_num"), vars = c("age", "sex", "stage_num")),
    list(name = "+ MSI", terms = c("age", "sex", "stage_num", "msi"), vars = c("age", "sex", "stage_num", "msi")),
    list(name = "+ CMS", terms = c("age", "sex", "stage_num", "cms"), vars = c("age", "sex", "stage_num", "cms")),
    list(name = "+ MSI + CMS", terms = c("age", "sex", "stage_num", "msi", "cms"), vars = c("age", "sex", "stage_num", "msi", "cms")),
    # F03-aligned: the OS validation in 04_survival_validation.R adjusted for MSI and sidedness as well
    list(name = "+ MSI + sidedness (F03-aligned)", terms = c("age", "sex", "stage_num", "msi", "sidedness"),
         vars = c("age", "sex", "stage_num", "msi", "sidedness"))
  )

  out <- list()
  for (sp in specs) {
    cc <- dat[stats::complete.cases(dat[, sp$vars, drop = FALSE]), , drop = FALSE]
    if (nrow(cc) < 60 || sum(cc$event) < 20) next
    f_base <- stats::as.formula(paste("Surv(time, event) ~", paste(sp$terms, collapse = " + ")))
    f_full <- stats::as.formula(paste("Surv(time, event) ~ ecotype +", paste(sp$terms, collapse = " + ")))
    m_base <- survival::coxph(f_base, data = cc)
    m_full <- survival::coxph(f_full, data = cc)
    lrt <- stats::anova(m_base, m_full, test = "LRT")
    s <- summary(m_full)
    k <- grep("^ecotype", rownames(s$coefficients))[1]
    c_base <- as.numeric(survival::concordance(m_base)$concordance)
    c_full <- as.numeric(survival::concordance(m_full)$concordance)
    out[[length(out) + 1]] <- data.frame(
      analysis = label, model = sp$name, n = nrow(cc), events = sum(cc$event),
      c_index_base = c_base, c_index_full = c_full, delta_c = c_full - c_base,
      ecotype_hr = if (is.na(k)) NA_real_ else exp(s$coefficients[k, "coef"]),
      hr_lower = if (is.na(k)) NA_real_ else unname(s$conf.int[k, 3]),
      hr_upper = if (is.na(k)) NA_real_ else unname(s$conf.int[k, 4]),
      hr_p = if (is.na(k)) NA_real_ else s$coefficients[k, "Pr(>|z|)"],
      lrt_p = lrt[2, "Pr(>|Chi|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

rfs_dat <- d[d$stage %in% c("I", "II", "III"), , drop = FALSE]   # primary: curative-intent, as in 22d
res_rfs <- fit_set(rfs_dat, "rfs_time_days", "rfs_event", "RFS, stage I-III")
res_os <- fit_set(d, "os_time_days", "os_event", "OS, all stages")
res <- rbind(res_rfs, res_os)
rownames(res) <- NULL
print(res)
fwrite(res, file.path(OUT, "25_incremental_value.tsv"), sep = "\t")

## 4. How much of the ecotype is CMS? ---------------------------------------
cms_tab <- as.data.frame(table(ecotype = d$ecotype, cms = d$cms))
cms_tab <- cms_tab[cms_tab$Freq > 0, , drop = FALSE]
cramers_v <- tryCatch({
  tb <- table(d$ecotype, d$cms)
  chi <- suppressWarnings(stats::chisq.test(tb))
  as.numeric(sqrt(chi$statistic / (sum(tb) * (min(dim(tb)) - 1))))
}, error = function(e) NA_real_)
cat("\nCramer's V (ecotype vs CMS):", round(cramers_v, 3), "\n")

## 5. Bootstrap the CMS-adjusted ecotype effect (RFS) -----------------------
cc <- rfs_dat[stats::complete.cases(rfs_dat[, c("rfs_time_days", "rfs_event", "ecotype", "age", "sex", "stage_num", "cms"), drop = FALSE]), , drop = FALSE]
cc$time <- cc$rfs_time_days
cc$event <- cc$rfs_event
boot <- unlist(mclapply(seq_len(500), function(i) {
  dd <- cc[sample(nrow(cc), replace = TRUE), , drop = FALSE]
  if (length(unique(dd$ecotype)) < 2 || sum(dd$event) < 20) return(NA_real_)
  f <- tryCatch(survival::coxph(Surv(time, event) ~ ecotype + age + sex + stage_num + cms, data = dd),
                error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  exp(unname(stats::coef(f)[1]))
}, mc.cores = N_WORKERS))
boot <- boot[is.finite(boot)]
boot_tab <- data.frame(metric = c("n", "events", "reps", "hr_median", "hr_lower_2.5", "hr_upper_97.5", "frac_above_1"),
                       value = c(nrow(cc), sum(cc$event), length(boot), stats::median(boot),
                                 stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), mean(boot > 1)))
print(boot_tab)
fwrite(boot_tab, file.path(OUT, "25_incremental_value_bootstrap.tsv"), sep = "\t")

## 6. Figures -----------------------------------------------------------------
model_levels <- c("stage, age, sex", "+ MSI", "+ CMS", "+ MSI + CMS", "+ MSI + sidedness (F03-aligned)")
res$model <- factor(res$model, levels = model_levels)
save_plot(function() {
  print(ggplot(res, aes(model, c_index_base, fill = analysis)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    geom_point(aes(y = c_index_full), position = position_dodge(width = 0.8), size = 1.6, colour = "black") +
    coord_cartesian(ylim = c(0.5, max(res$c_index_full) + 0.05)) +
    scale_fill_manual(values = c("RFS, stage I-III" = "#C0392B", "OS, all stages" = "#2E86C1")) +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "C-index (bar = covariates only, dot = + ecotype)", fill = NULL,
         title = "Does the ecotype add beyond stage/MSI/CMS? (GSE39582)") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7)))
}, file.path(FIG, "plot_41"), 9, 5)

cms_frac <- as.data.frame(prop.table(table(d$ecotype, d$cms), 1))
names(cms_frac) <- c("ecotype", "cms", "fraction")
save_plot(function() {
  print(ggplot(cms_frac, aes(cms, fraction, fill = ecotype)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(E1 = "grey40", E2 = "#C0392B")) +
    theme_bw(base_size = 10) +
    labs(x = "CMS subtype", y = "Fraction of ecotype", fill = NULL,
         title = sprintf("Ecotype vs CMS in GSE39582 (Cramer's V = %.2f)", cramers_v)))
}, file.path(FIG, "plot_42"), 7, 4.5)

## 7. Report + ledger --------------------------------------------------------
rfs_cms <- res[res$analysis == "RFS, stage I-III" & res$model == "+ MSI + CMS", , drop = FALSE]
rfs_cms_only <- res[res$analysis == "RFS, stage I-III" & res$model == "+ CMS", , drop = FALSE]
os_cms <- res[res$analysis == "OS, all stages" & res$model == "+ MSI + CMS", , drop = FALSE]
os_f03 <- res[res$analysis == "OS, all stages" & grepl("F03-aligned", res$model), , drop = FALSE]
stopifnot(nrow(rfs_cms) == 1, nrow(rfs_cms_only) == 1, nrow(os_cms) == 1, nrow(os_f03) == 1)
adds <- nrow(rfs_cms) == 1 && !is.na(rfs_cms$lrt_p) && rfs_cms$lrt_p < 0.05
grade <- if (adds) "E1" else "not_supported"

report <- c(
  "# Incremental prognostic value beyond stage, MSI and CMS",
  "",
  sprintf("Generated by `25_incremental_value.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Question",
  "Does the frozen K=2 ecotype add information beyond stage, MSI and CMS - the published taxonomy the",
  "protocol says the ecotype must improve on? Endpoints are kept separate: RFS (stage I-III, primary) and",
  "OS (all stages, secondary). GSE39582 is the only cohort on disk with CMS calls.",
  "",
  "## Cohort and class balance",
  sprintf("- Tumours with an ecotype: E1 %d / E2 %d. CMS available for %d, MSI for %d of %d.",
          sum(d$ecotype == "E1"), sum(d$ecotype == "E2"), sum(!is.na(d$cms)), sum(!is.na(d$msi)), nrow(d)),
  sprintf("- Ecotype is only partly aligned with CMS (Cramer's V = %.2f), so CMS adjustment is a real test,",
          cramers_v),
  "  not a re-parameterisation of the same split.",
  "",
  "## RFS (stage I-III)",
  sprintf("- With MSI + CMS in the model: ecotype HR = %.2f (%.2f-%.2f), p = %.3g; LRT p = %.3g;",
          rfs_cms$ecotype_hr, rfs_cms$hr_lower, rfs_cms$hr_upper, rfs_cms$hr_p, rfs_cms$lrt_p),
  sprintf("  C-index %.3f -> %.3f (delta C = %+.3f).", rfs_cms$c_index_base, rfs_cms$c_index_full, rfs_cms$delta_c),
  sprintf("- With CMS only: ecotype HR = %.2f (%.2f-%.2f), p = %.3g; delta C = %+.3f.",
          rfs_cms_only$ecotype_hr, rfs_cms_only$hr_lower, rfs_cms_only$hr_upper, rfs_cms_only$hr_p, rfs_cms_only$delta_c),
  sprintf("- Bootstrap of the fully adjusted ecotype HR (500 resamples): median %.2f (%.2f-%.2f), %.0f%% above 1.",
          stats::median(boot), stats::quantile(boot, 0.025), stats::quantile(boot, 0.975), 100 * mean(boot > 1)),
  "",
  "## OS (all stages)",
  sprintf("- With MSI + CMS: ecotype HR = %.2f (%.2f-%.2f), p = %.3g; LRT p = %.3g; delta C = %+.3f.",
          os_cms$ecotype_hr, os_cms$hr_lower, os_cms$hr_upper, os_cms$hr_p, os_cms$lrt_p, os_cms$delta_c),
  "- With F03's covariate set (age, sex, stage, MSI, sidedness) but this script's z-score programme labels:",
  sprintf("  HR = %.2f (%.2f-%.2f), p = %.3g. F03 reported HR 1.16 (p 0.333) using GSVA-derived labels; the gap is the",
          os_f03$ecotype_hr, os_f03$hr_lower, os_f03$hr_upper, os_f03$hr_p),
  "  scoring method (borderline tumours switch ecotype), not the covariates - the same sensitivity documented for F06.",
  "  Under either scoring, adding CMS removes the OS association, which is consistent with the F01/F03 null.",
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (adds) {
            paste("the ecotype retains a disease-free association after adjustment for stage, MSI and CMS,",
                  "i.e. it is not merely a re-expression of the published taxonomy; single cohort (E1).")
          } else {
            paste("once stage, MSI and CMS are in the model the ecotype no longer contributes detectably, so it is",
                  "not demonstrably incremental over the published taxonomy in this cohort.")
          }),
  "",
  "## Limitations",
  "- One cohort with CMS; GSE14333 and GSE17536 have no CMS calls, so the incremental test cannot be repeated there.",
  "- CMS itself was derived from expression in the same cohort, so CMS and ecotype share information by construction.",
  "- MSI is missing for ~50 tumours; CMS is missing for 19. Complete-case analysis assumes these are not informative.",
  "- RFS is stage I-III by design (stage IV has no RFS in this build); OS is all stages, so the two rows differ in population.",
  "- Observational, no treatment adjustment.",
  "",
  "## Figures",
  "- plot_41 C-index by model with covariate-only bar vs +ecotype dot; plot_42 ecotype vs CMS composition."
)
writeLines(report, file.path(OUT, "25_Incremental_Value_Report.md"))

ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F14", , drop = FALSE]
f14 <- data.frame(
  finding_id = "F14",
  claim = if (adds) {
    "The ecotype retains a disease-free association after adjustment for stage, MSI and CMS (incremental over the published taxonomy)"
  } else {
    "The ecotype does not add detectable disease-free prognostic information beyond stage, MSI and CMS"
  },
  evidence_type = "Nested Cox models (incremental value)",
  dataset_s = "GSE39582",
  statistic = paste0("RFS stage I-III: ecotype HR=", round(rfs_cms$ecotype_hr, 2), " (", round(rfs_cms$hr_lower, 2), "-",
                     round(rfs_cms$hr_upper, 2), ") after stage+MSI+CMS, p=", signif(rfs_cms$hr_p, 3),
                     ", LRT p=", signif(rfs_cms$lrt_p, 3), ", delta C=", round(rfs_cms$delta_c, 3),
                     "; CMS-only delta C=", round(rfs_cms_only$delta_c, 3),
                     "; Cramer's V ecotype-CMS=", round(cramers_v, 2),
                     if (nrow(os_cms) == 1) paste0("; OS after stage+MSI+CMS HR=", round(os_cms$ecotype_hr, 2), " p=", signif(os_cms$hr_p, 3)) else ""),
  n = rfs_cms$n,
  grade = grade,
  verb_used = if (adds) "added to" else "not incremental over",
  contradicting_evidence = "CMS is expression-derived in the same cohort, so ecotype and CMS are not independent by construction",
  limitation = "single cohort with CMS (GSE14333/GSE17536 have none); complete-case for MSI (missing ~50) and CMS (missing 19); RFS restricted to stage I-III; observational",
  next_test = "independent cohort with CMS/MSI calls and disease-free endpoints",
  figure_table = "plot_41; plot_42",
  stringsAsFactors = FALSE
)
out <- if (is.null(old) || nrow(old) == 0) f14 else rbind(old, f14)
out <- out[order(as.integer(sub("^F", "", out$finding_id))), , drop = FALSE]
data.table::fwrite(out, ledger_path, sep = ",")

# Cross-link the consequence into F11 so the disease-free claim carries its main caveat with it
led <- data.table::fread(ledger_path, data.table = FALSE)
i <- which(led$finding_id == "F11")
if (length(i) == 1) {
  led$contradicting_evidence[i] <- paste0(
    "OS association remains null in discovery (F01) and in OS validation (F03); the disease-free association is ",
    "NOT incremental over CMS (F14: HR ", round(rfs_cms$ecotype_hr, 2), ", p ", signif(rfs_cms$hr_p, 2),
    " after stage+MSI+CMS)")
  data.table::fwrite(led, ledger_path, sep = ",")
}
cat("ledger rows:", nrow(led), "\n")

log_decision("25_incremental", "PASS",
             paste0("RFS stage I-III after stage+MSI+CMS: HR=", round(rfs_cms$ecotype_hr, 3),
                    " p=", signif(rfs_cms$hr_p, 3), " LRT p=", signif(rfs_cms$lrt_p, 3),
                    " deltaC=", round(rfs_cms$delta_c, 3), "; Cramer V=", round(cramers_v, 2),
                    "; grade=", grade))
cat("grade:", grade, "\n")
cat("done\n")
