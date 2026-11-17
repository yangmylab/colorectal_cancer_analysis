# 07_robustness.R
# K sensitivity and outlier-exclusion sensitivity for the primary ecotype-survival result.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(SNFtool)
  library(yaml)
})

cfg <- yaml::read_yaml("config/pipeline_config.yaml")
ec <- readRDS(file.path(cfg$paths$processed, "ecotypes.rds"))
surv <- readRDS(file.path(cfg$paths$processed, "survival_results.rds"))
qc <- data.table::fread(file.path(cfg$paths$qc, "Sample_QC.tsv"), data.table = FALSE)
dat0 <- surv$dat

covs <- c("age", "sex", "stage_num", "msi_status", "sidedness")
covs <- covs[covs %in% names(dat0)]
covs <- covs[vapply(dat0[covs], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
for (cv in covs) if (is.character(dat0[[cv]])) dat0[[cv]] <- factor(dat0[[cv]])

fit_k <- function(K) {
  lab <- SNFtool::spectralClustering(ec$W, K, type = 3)
  names(lab) <- rownames(ec$combined)
  dat <- dat0
  dat$ecotype <- factor(paste0("E", lab[as.character(dat$sample_id)]))
  dat <- dat[!is.na(dat$ecotype), , drop = FALSE]
  f_base <- as.formula(paste("Surv(time, event) ~", if (length(covs) > 0) paste(covs, collapse = " + ") else "1"))
  f_full <- as.formula(paste("Surv(time, event) ~ ecotype", if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""))
  base <- survival::coxph(f_base, data = dat, x = TRUE)
  full <- survival::coxph(f_full, data = dat, x = TRUE)
  a <- tryCatch(stats::anova(base, full, test = "LRT"), error = function(e) NULL)
  p <- if (!is.null(a) && "Pr(>|Chi|)" %in% names(as.data.frame(a))) as.data.frame(a)[2, "Pr(>|Chi|)"] else NA_real_
  data.frame(
    K = K,
    n = nrow(dat),
    events = sum(dat$event),
    c_index_base = as.numeric(survival::concordance(base)$concordance),
    c_index_full = as.numeric(survival::concordance(full)$concordance),
    lrt_p = p,
    min_group_fraction = min(table(dat$ecotype)) / nrow(dat),
    stringsAsFactors = FALSE
  )
}

K_tab <- do.call(rbind, lapply(2:5, fit_k))
dir.create(file.path("08_Robustness"), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(K_tab, file.path("08_Robustness", "K_sensitivity.tsv"), sep = "\t")

# Outlier-exclusion sensitivity for the primary ecotype labels
outliers <- qc$sample_id[qc$outlier_flag]
dat_out <- dat0[!dat0$sample_id %in% outliers, , drop = FALSE]
dat_out$ecotype <- factor(as.character(dat_out$ecotype))
f_base <- as.formula(paste("Surv(time, event) ~", if (length(covs) > 0) paste(covs, collapse = " + ") else "1"))
f_full <- as.formula(paste("Surv(time, event) ~ ecotype", if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""))
b <- survival::coxph(f_base, data = dat_out, x = TRUE)
f <- survival::coxph(f_full, data = dat_out, x = TRUE)
a <- tryCatch(stats::anova(b, f, test = "LRT"), error = function(e) NULL)
p_out <- if (!is.null(a) && "Pr(>|Chi|)" %in% names(as.data.frame(a))) as.data.frame(a)[2, "Pr(>|Chi|)"] else NA_real_
out_tab <- data.frame(
  analysis = paste0("primary K=", ec$K, " excluding QC outliers"),
  n = nrow(dat_out),
  events = sum(dat_out$event),
  c_index_base = as.numeric(survival::concordance(b)$concordance),
  c_index_full = as.numeric(survival::concordance(f)$concordance),
  lrt_p = p_out,
  stringsAsFactors = FALSE
)
data.table::fwrite(out_tab, file.path("08_Robustness", "outlier_sensitivity.tsv"), sep = "\t")

primary <- K_tab[K_tab$K == ec$K, , drop = FALSE]
best <- K_tab[which.max(K_tab$c_index_full), , drop = FALSE]
delta <- best$c_index_full - primary$c_index_full
decision <- if (delta > 0.01 && !is.na(best$lrt_p) && best$lrt_p < 0.05) {
  paste("A finer K (", best$K, ") improves C-index by ", round(delta, 3), " and reaches LRT p < 0.05; report as exploratory and consider a pre-specified follow-up.", sep = "")
} else {
  paste0("No finer K improves C-index by >0.01 with LRT p < 0.05; retain the pre-specified stable K=", ec$K, " solution as primary.")
}

report <- c(
  "# Robustness and sensitivity report — CRC Track A",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Primary result",
  paste("- K =", primary$K, ", n =", primary$n, ", events =", primary$events),
  paste("- C-index base/full:", round(primary$c_index_base, 3), "/", round(primary$c_index_full, 3)),
  paste("- LRT p:", signif(primary$lrt_p, 3)),
  "",
  "## K sensitivity (same outcome, same covariates)",
  "",
  "| K | n | events | C-index base | C-index full | LRT p | min group fraction |",
  "|---|---|---|---|---|---|---|",
  apply(K_tab, 1, function(r) paste0("| ", r[["K"]], " | ", r[["n"]], " | ", r[["events"]], " | ", round(as.numeric(r[["c_index_base"]]), 3), " | ", round(as.numeric(r[["c_index_full"]]), 3), " | ", signif(as.numeric(r[["lrt_p"]]), 3), " | ", round(as.numeric(r[["min_group_fraction"]]), 3), " |")),
  "",
  "## Outlier-exclusion sensitivity",
  paste("- Excluding", length(outliers), "QC outliers: n =", out_tab$n, ", events =", out_tab$events),
  paste("- C-index base/full:", round(out_tab$c_index_base, 3), "/", round(out_tab$c_index_full, 3), "; LRT p:", signif(out_tab$lrt_p, 3)),
  "",
  "## Decision",
  decision,
  "",
  "## Interpretation",
  paste0("The primary K=", ec$K, " solution is retained because it is the selected stability solution in the discovery data. Finer partitions are reported as exploratory sensitivity analyses."),
  "",
  "## Limitations",
  "- TCGA-COAD/READ only; OS not CSS; MSI and purity unavailable.",
  "- No independent external validation yet.",
  "- Hallmark programmes only; no protein, methylation, single-cell or spatial layers in this run."
)
writeLines(report, file.path("08_Robustness", "Robustness_Report.md"))

cat("K sensitivity:\n")
print(K_tab)
cat("\nOutlier sensitivity:\n")
print(out_tab)
cat("\nDecision:", decision, "\n")
