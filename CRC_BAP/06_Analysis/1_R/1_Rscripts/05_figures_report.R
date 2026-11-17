# 05_figures_report.R
# Assemble summary figures, decision log and final/dry-run report.

run_report <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(survival)
  })

  bulk <- readRDS(file.path(cfg$paths$processed, "processed_bulk.rds"))
  sig <- readRDS(file.path(cfg$paths$processed, "signature_scores.rds"))
  ec <- readRDS(file.path(cfg$paths$processed, "ecotypes.rds"))
  surv <- readRDS(file.path(cfg$paths$processed, "survival_results.rds"))

  fig_dir <- cfg$paths$figures
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  tab_dir <- cfg$paths$tables
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  # Summary figure: stability + KM
  stab <- ec$stability
  save_plot(function() {
    print(ggplot(stab, aes(K, mean_ari)) + geom_line() + geom_point() + geom_hline(yintercept = cfg$params$stability_ari_min, linetype = 2) + theme_bw() + labs(x = "K", y = "Bootstrap ARI", title = "Ecotype stability"))
  }, file.path(fig_dir, "Figure_1_stability"), 6, 4)

  km <- surv$km
  save_plot(function() {
    plot(km, col = seq_along(levels(surv$dat$ecotype)), lwd = 2, xlab = "Days", ylab = "Survival probability", main = "Kaplan-Meier by ecotype")
    legend("topright", legend = levels(surv$dat$ecotype), col = seq_along(levels(surv$dat$ecotype)), lwd = 2, bty = "n")
  }, file.path(fig_dir, "Figure_2_km"), 7, 5)

  # Tables
  data.table::fwrite(ec$stability, file.path(tab_dir, "Table_stability.tsv"), sep = "\t")
  data.table::fwrite(surv$perf, file.path(tab_dir, "Table_survival_performance.tsv"), sep = "\t")
  data.table::fwrite(surv$validation, file.path(tab_dir, "Table_validation.tsv"), sep = "\t")

  # Decision log markdown
  dec_path <- get0("DECISION_LOG_PATH", envir = .GlobalEnv, ifnotfound = file.path(cfg$paths$logs, "DECISION_LOG.tsv"))
  dec <- if (file.exists(dec_path)) data.table::fread(dec_path, data.table = FALSE) else data.frame()
  dec_md <- c("# Decision Log", "", paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "")
  if (nrow(dec) > 0) {
    dec_md <- c(dec_md, "| time | step | decision | rationale |", "|---|---|---|---|", apply(dec, 1, function(r) paste0("| ", r[["timestamp"]], " | ", r[["step"]], " | ", r[["decision"]], " | ", r[["rationale"]], " |")))
  } else {
    dec_md <- c(dec_md, "No decisions recorded.")
  }
  writeLines(dec_md, file.path(cfg$paths$logs, "DECISION_LOG.md"))

  if (synthetic) {
    report <- c(
      "# Synthetic dry-run report",
      "",
      "**SYNTHETIC PIPELINE VALIDATION — NOT BIOLOGICAL EVIDENCE.**",
      "",
      "This run used simulated data solely to verify that every module executes and that the decision logic works.",
      "",
      "## Results",
      paste("- Samples:", ncol(bulk$log_expr)),
      paste("- Features:", nrow(bulk$log_expr)),
      paste("- Signatures scored:", nrow(sig$scores)),
      paste("- Selected K:", ec$K),
      paste("- Best bootstrap ARI:", round(max(ec$stability$mean_ari, na.rm = TRUE), 3)),
      paste("- Classifier CV accuracy:", round(ec$classifier$cv_accuracy, 3)),
      paste("- Survival endpoint:", surv$perf$endpoint),
      paste("- Events:", surv$perf$events),
      paste("- C-index (full):", round(surv$perf$c_index_full, 3)),
      paste("- C-index (base):", round(surv$perf$c_index_base, 3)),
      paste("- LRT p:", signif(surv$perf$lrt_p, 3)),
      paste("- Validation cohorts:", nrow(surv$validation)),
      "",
      "## Decision",
      "The pipeline ran end-to-end. Synthetic outputs are isolated under `99_Synthetic_DryRun/` and must never be reported as findings.",
      "",
      "## What remains blocked",
      "The real analysis is blocked at G1 because no CRC datasets are present in `03_Raw_Data/`."
    )
    writeLines(report, file.path(cfg$paths$logs, "DRYRUN_REPORT.md"))
  } else {
    ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
    ledger <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else data.frame()
    report <- c(
      "# Final Report — CRC Track A",
      "",
      paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      "## 1. Research question and hypothesis",
      "See `00_Project_Design/Analysis_Plan.md` §1–2.",
      "",
      "## 2. Datasets and provenance",
      "See `01_Dataset_Discovery/Dataset_Inventory.csv` and `Dataset_Selection_Report.md`.",
      "",
      "## 3. Sample annotation",
      paste("Samples:", ncol(bulk$log_expr), "; donors:", length(unique(bulk$clinical$patient_id))),
      "G1 annotation gate: see `02_Metadata/Metadata_Audit_Report.md`.",
      "",
      "## 4. Design and statistical plan",
      "See `00_Project_Design/Analysis_Plan.md` §4.",
      "",
      "## 5. Quality control",
      paste("Outliers flagged:", sum(bulk$qc$outlier_flag), "of", nrow(bulk$qc)),
      "See `05_QC/QC_Report.md`.",
      "",
      "## 6. Methods",
      "Scripts: `06_Analysis/1_R/1_Rscripts/`; logs: `Logs/`; environment: `Environment/`.",
      "",
      "## 7. Main findings",
      if (nrow(ledger) > 0) "See `13_Reports/Findings_Ledger.csv`." else "No findings yet (analysis incomplete).",
      "",
      "## 8. Biological interpretation",
      "Not written until findings are graded and independently validated.",
      "",
      "## 9. Independent validation",
      "See `07_Validation/Validation_Report.md`.",
      "",
      "## 10. Robustness",
      "See `08_Robustness/Robustness_Report.md`.",
      "",
      "## 11. Limitations",
      "Retrospective public data; computational-only; no wet-lab validation.",
      "",
      "## 12. Reproducibility package",
      "Scripts, configs and decision log are in this project. A Zenodo DOI is pending.",
      "",
      "## 13. Proposed follow-up experiments",
      "IHC/CODEX validation of ecotype markers; organoid perturbation of prioritized targets.",
      "",
      "## 14. Reviewer scorecard",
      "| Class | Issue | Evidence | Corrective action | Status |",
      "|---|---|---|---|---|",
      "| Critical | None yet | — | — | — |",
      "| Major | No wet-lab validation | design | state limitation | open |",
      "| Minor | Public-data heterogeneity | design | sensitivity analyses | open |",
      "",
      "**Verdict:** analysis incomplete — G1 gate must pass before claims.",
      "",
      "## 15. Evidence-grade distribution",
      "No E1–E4 claims yet."
    )
    writeLines(report, file.path("13_Reports", "Final_Report.md"))
  }

  handoff <- c(
    "# Handoff",
    "",
    paste("Updated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## State",
    if (synthetic) "- Synthetic dry run completed; real analysis blocked at G1 (no local CRC data)." else "- Real analysis pipeline completed; inspect Findings_Ledger and validation report.",
    "",
    "## Next action",
    if (synthetic) "- Download Tier 1 datasets into `03_Raw_Data/`, create `READY`, then run `bash run_all.sh`." else "- Review findings, robustness and validation before manuscript claims.",
    "",
    "## Key files",
    "- `00_Project_Design/Analysis_Plan.md`",
    "- `02_Metadata/Sample_Annotation.csv`",
    "- `13_Reports/Findings_Ledger.csv`",
    "- `Logs/DECISION_LOG.tsv`"
  )
  handoff_path <- if (synthetic) file.path(dirname(cfg$paths$raw), "handoff.md") else "handoff.md"
  writeLines(handoff, handoff_path)

  log_decision("05_report", "PASS", if (synthetic) "synthetic dry-run report written" else "final report written")
  invisible(TRUE)
}
