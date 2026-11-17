#!/usr/bin/env Rscript
# run_pipeline.R
# Orchestrates the CRC Track A BAP pipeline with explicit decision logging.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

config_path <- get_arg("--config", "config/pipeline_config.yaml")
synthetic <- "--synthetic" %in% args

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

cfg <- yaml::read_yaml(config_path)
root <- normalizePath(file.path(dirname(config_path), ".."))
setwd(root)

for (d in unique(unlist(cfg$paths))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

DECISION_LOG_PATH <- file.path(cfg$paths$logs, "DECISION_LOG.tsv")
assign("DECISION_LOG_PATH", DECISION_LOG_PATH, envir = .GlobalEnv)
if (file.exists(DECISION_LOG_PATH)) file.remove(DECISION_LOG_PATH)

set.seed(cfg$params$seed)
options(mc.cores = cfg$params$n_workers)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

source("06_Analysis/1_R/1_Rscripts/utils.R")
source("06_Analysis/1_R/1_Rscripts/00_validate_inputs.R")
source("06_Analysis/1_R/1_Rscripts/01_bulk_qc.R")
source("06_Analysis/1_R/1_Rscripts/02_signatures_deconvolution.R")
source("06_Analysis/1_R/1_Rscripts/03_snf_ecotypes.R")
source("06_Analysis/1_R/1_Rscripts/04_survival_validation.R")
source("06_Analysis/1_R/1_Rscripts/05_figures_report.R")

run_step <- function(name, fun) {
  message("\n=== ", name, " ===")
  tryCatch(
    fun(),
    error = function(e) {
      msg <- paste("PIPELINE FAILED at", name, ":", conditionMessage(e))
      writeLines(msg, file.path(cfg$paths$logs, "PIPELINE_FAILED.md"))
      log_decision(name, "FAIL", conditionMessage(e))
      stop(msg, call. = FALSE)
    }
  )
}

run_step("00_validate_inputs", function() run_validate_inputs(cfg, synthetic))
run_step("01_bulk_qc", function() run_bulk_qc(cfg, synthetic))
run_step("02_signatures", function() run_signatures(cfg, synthetic))
run_step("03_ecotypes", function() run_ecotypes(cfg, synthetic))
run_step("04_survival", function() run_survival_validation(cfg, synthetic))
run_step("05_report", function() run_report(cfg, synthetic))

writeLines(
  c(
    paste("PIPELINE COMPLETE", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("config:", config_path),
    paste("synthetic:", synthetic)
  ),
  file.path(cfg$paths$logs, "PIPELINE_COMPLETE.md")
)

message("\nPipeline complete. Logs: ", cfg$paths$logs)
