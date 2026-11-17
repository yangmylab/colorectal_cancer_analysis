#!/usr/bin/env Rscript
# 99_synthetic_dryrun.R
# Generate a synthetic multi-view cohort, run the full pipeline, and verify all modules execute.
# SYNTHETIC DATA ONLY — this is a code/decision-logic test, not biological evidence.

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20261117)

base_dir <- "99_Synthetic_DryRun"
data_dir <- file.path(base_dir, "data")
val_dir <- file.path(data_dir, "validation", "SYNTH_VAL")
for (d in c(base_dir, data_dir, val_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(base_dir, "DRYRUN_FAILED.md"))
unlink(file.path(base_dir, "logs", "PIPELINE_FAILED.md"))

n <- 120
p <- 2000
K <- 4
ecotype_true <- sample(seq_len(K), n, replace = TRUE, prob = c(0.30, 0.30, 0.25, 0.15))
names(ecotype_true) <- paste0("S", seq_len(n))

base <- stats::rnorm(p, mean = 6, sd = 1.5)
effect <- matrix(0, p, K)
for (k in seq_len(K)) {
  idx <- sample(p, 150)
  effect[idx, k] <- stats::rnorm(150, mean = 1.2, sd = 0.4)
}

expr <- base + effect[, ecotype_true] + matrix(stats::rnorm(n * p, 0, 0.6), p, n)
rownames(expr) <- paste0("G", seq_len(p))
colnames(expr) <- names(ecotype_true)
expr <- pmax(expr, 0)

# Signatures: four ecotype marker sets + eight random sets
sig_sets <- list()
for (k in seq_len(K)) {
  idx <- order(effect[, k], decreasing = TRUE)[seq_len(100)]
  sig_sets[[paste0("SYNTH_ECOTYPE_", k)]] <- rownames(expr)[idx]
}
for (j in seq_len(8)) {
  sig_sets[[paste0("SYNTH_RANDOM_", j)]] <- sample(rownames(expr), 50)
}
gmt_lines <- vapply(names(sig_sets), function(nm) paste(c(nm, "synthetic", sig_sets[[nm]]), collapse = "\t"), character(1))
writeLines(gmt_lines, file.path(data_dir, "signatures.gmt"))

# Methylation view
pm <- 300
meth_effect <- matrix(0, pm, K)
for (k in seq_len(K)) meth_effect[sample(pm, 60), k] <- stats::rnorm(60, 0.9, 0.3)
meth <- plogis(stats::rnorm(pm, 0, 1) + meth_effect[, ecotype_true] + matrix(stats::rnorm(pm * n, 0, 0.4), pm, n))
rownames(meth) <- paste0("M", seq_len(pm))
colnames(meth) <- colnames(expr)

# Protein view
pp <- 300
prot_effect <- matrix(0, pp, K)
for (k in seq_len(K)) prot_effect[sample(pp, 60), k] <- stats::rnorm(60, 1.0, 0.3)
prot <- 8 + prot_effect[, ecotype_true] + matrix(stats::rnorm(pp * n, 0, 0.5), pp, n)
rownames(prot) <- paste0("P", seq_len(pp))
colnames(prot) <- colnames(expr)

stage <- sample(c("I", "II", "III", "IV"), n, replace = TRUE, prob = c(0.2, 0.35, 0.3, 0.15))
stage_num <- as.numeric(factor(stage, levels = c("I", "II", "III", "IV")))
age <- round(stats::runif(n, 40, 85))
sex <- sample(c("M", "F"), n, replace = TRUE)
msi <- sample(c("MSS", "MSI"), n, replace = TRUE, prob = c(0.85, 0.15))
side <- sample(c("left", "right"), n, replace = TRUE)
batch <- sample(c("B1", "B2", "B3"), n, replace = TRUE)
purity <- round(stats::runif(n, 0.35, 0.9), 2)

lp <- 1.0 * (ecotype_true == 4) + 0.6 * (ecotype_true == 2) + 0.35 * stage_num + 0.02 * (age - 60)
make_surv <- function(scale_h) {
  h <- scale_h * exp(lp)
  te <- stats::rexp(n, h)
  tc <- stats::runif(n, 365, 3000)
  list(time = pmin(te, tc), event = as.integer(te <= tc))
}
os <- make_surv(0.0009)
css <- make_surv(0.0011)
rfs <- make_surv(0.0015)

clinical <- data.frame(
  sample_id = colnames(expr),
  patient_id = paste0("P", seq_len(n)),
  group = "tumour",
  disease_status = "CRC",
  treatment_status = "na",
  timepoint = "baseline",
  tissue = "colon",
  cell_type = "epithelial",
  technology = "RNA-seq",
  platform = "synthetic",
  library_prep = "synthetic",
  repository_id = "SYNTH",
  stage = stage,
  msi_status = msi,
  sidedness = side,
  age = age,
  sex = sex,
  batch = batch,
  purity = purity,
  os_time_days = round(os$time),
  os_event = os$event,
  css_time_days = round(css$time),
  css_event = css$event,
  rfs_time_days = round(rfs$time),
  rfs_event = rfs$event,
  true_ecotype = paste0("E", ecotype_true),
  stringsAsFactors = FALSE
)

write_tsv <- function(mat, path, id_col) {
  df <- data.frame(id = rownames(mat), mat, check.names = FALSE)
  colnames(df)[1] <- id_col
  data.table::fwrite(df, path, sep = "\t")
}

write_tsv(expr, file.path(data_dir, "expression.tsv"), "gene_id")
write_tsv(meth, file.path(data_dir, "methylation.tsv"), "feature_id")
write_tsv(prot, file.path(data_dir, "protein.tsv"), "feature_id")
data.table::fwrite(clinical, file.path(data_dir, "clinical.tsv"), sep = "\t")

# Validation cohort
nv <- 100
ev <- sample(seq_len(K), nv, replace = TRUE, prob = c(0.30, 0.30, 0.25, 0.15))
expr_v <- base + effect[, ev] + matrix(stats::rnorm(nv * p, 0, 0.6), p, nv)
rownames(expr_v) <- rownames(expr)
colnames(expr_v) <- paste0("V", seq_len(nv))
expr_v <- pmax(expr_v, 0)
stage_v <- sample(c("I", "II", "III", "IV"), nv, replace = TRUE, prob = c(0.2, 0.35, 0.3, 0.15))
stage_num_v <- as.numeric(factor(stage_v, levels = c("I", "II", "III", "IV")))
age_v <- round(stats::runif(nv, 40, 85))
lp_v <- 1.0 * (ev == 4) + 0.6 * (ev == 2) + 0.35 * stage_num_v + 0.02 * (age_v - 60)
h_v <- 0.0009 * exp(lp_v)
te_v <- stats::rexp(nv, h_v)
tc_v <- stats::runif(nv, 365, 3000)
clinical_v <- data.frame(
  sample_id = colnames(expr_v),
  patient_id = paste0("VP", seq_len(nv)),
  group = "tumour",
  disease_status = "CRC",
  tissue = "colon",
  stage = stage_v,
  msi_status = sample(c("MSS", "MSI"), nv, replace = TRUE, prob = c(0.85, 0.15)),
  sidedness = sample(c("left", "right"), nv, replace = TRUE),
  age = age_v,
  sex = sample(c("M", "F"), nv, replace = TRUE),
  batch = "VAL",
  purity = round(stats::runif(nv, 0.35, 0.9), 2),
  os_time_days = round(pmin(te_v, tc_v)),
  os_event = as.integer(te_v <= tc_v),
  css_time_days = round(pmin(te_v, tc_v)),
  css_event = as.integer(te_v <= tc_v),
  rfs_time_days = round(pmin(te_v, tc_v)),
  rfs_event = as.integer(te_v <= tc_v),
  true_ecotype = paste0("E", ev),
  stringsAsFactors = FALSE
)
write_tsv(expr_v, file.path(val_dir, "expression.tsv"), "gene_id")
data.table::fwrite(clinical_v, file.path(val_dir, "clinical.tsv"), sep = "\t")

writeLines("SYNTHETIC DATA READY — not biological evidence", file.path(data_dir, "READY"))

writeLines(
  c(
    "# Synthetic dry run",
    "",
    "**SYNTHETIC PIPELINE VALIDATION — NOT BIOLOGICAL EVIDENCE.**",
    "",
    paste("Samples:", n, "Features:", p, "True ecotypes:", K),
    "Generated: expression, methylation, protein, signatures, clinical, validation cohort.",
    "The pipeline is executed by `99_synthetic_dryrun.R`."
  ),
  file.path(base_dir, "README.md")
)

status <- system2(
  "Rscript",
  c("06_Analysis/1_R/1_Rscripts/run_pipeline.R", "--config", "config/pipeline_config_synthetic.yaml", "--synthetic")
)

if (status != 0) {
  writeLines(paste("Synthetic pipeline failed with status", status), file.path(base_dir, "DRYRUN_FAILED.md"))
  quit(status = status)
}

report_src <- file.path(base_dir, "logs", "DRYRUN_REPORT.md")
if (file.exists(report_src)) file.copy(report_src, file.path(base_dir, "DRYRUN_REPORT.md"), overwrite = TRUE)
writeLines("Synthetic dry run completed; see DRYRUN_REPORT.md", file.path(base_dir, "DRYRUN_COMPLETE.md"))
