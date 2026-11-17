# 00_validate_inputs.R
# Validate inputs and enforce the BAP G1 sample-annotation gate.
# Real data: stops if G1 is blocked. Synthetic dry run: exercises the same checks.

run_validate_inputs <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
  })

  root <- getwd()
  for (d in unique(unlist(cfg$paths))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(cfg$paths$logs, "00_validate_inputs.log")
  con <- file(log_file, open = "wt")
  sink(con, split = TRUE)
  on.exit({ sink(); close(con) }, add = TRUE)

  cat("== 00_validate_inputs ==\n")
  cat("synthetic:", synthetic, "\n")
  cat("project root:", root, "\n")

  expr_path <- cfg$inputs$expression
  clin_path <- cfg$inputs$clinical
  missing_files <- c(expr_path, clin_path)[!file.exists(c(expr_path, clin_path))]
  if (length(missing_files) > 0) {
    msg <- paste(
      "GATE G1 — BLOCKED",
      paste("Blocking issue: missing input files:", paste(missing_files, collapse = ", ")),
      "Evidence: no local CRC data",
      "Samples affected: all planned samples",
      "Why it matters: sample and group assignment cannot be verified",
      "Proposed reading (NOT applied): download Tier 1 datasets into 03_Raw_Data/",
      "Decision needed: confirm/correct/supply the missing data",
      sep = "\n"
    )
    writeLines(msg, file.path(cfg$paths$logs, "G1_BLOCKED.md"))
    cat(msg, "\n")
    log_decision("00_validate", "G1_BLOCKED", paste("missing:", paste(missing_files, collapse = ", ")))
    stop(msg, call. = FALSE)
  }

  expr <- read_matrix_any(expr_path)
  clin <- data.table::fread(clin_path, data.table = FALSE)
  required <- c("sample_id", "patient_id", "tissue", "disease_status")
  missing_cols <- setdiff(required, colnames(clin))
  if (length(missing_cols) > 0) {
    msg <- paste("GATE G1 — BLOCKED", paste("Blocking issue: clinical.tsv missing columns:", paste(missing_cols, collapse = ", ")), sep = "\n")
    writeLines(msg, file.path(cfg$paths$logs, "G1_BLOCKED.md"))
    log_decision("00_validate", "G1_BLOCKED", paste("missing columns:", paste(missing_cols, collapse = ", ")))
    stop(msg, call. = FALSE)
  }

  if (!"group" %in% colnames(clin)) clin$group <- clin$disease_status
  if (!"treatment_status" %in% colnames(clin)) clin$treatment_status <- "na"
  if (!"timepoint" %in% colnames(clin)) clin$timepoint <- "baseline"
  if (!"cell_type" %in% colnames(clin)) clin$cell_type <- "epithelial"
  if (!"technology" %in% colnames(clin)) clin$technology <- "RNA-seq"
  if (!"platform" %in% colnames(clin)) clin$platform <- "unknown"
  if (!"library_prep" %in% colnames(clin)) clin$library_prep <- "unknown"
  if (!"repository_id" %in% colnames(clin)) clin$repository_id <- NA_character_
  if (!"batch" %in% colnames(clin)) clin$batch <- "batch_1"

  samples_expr <- colnames(expr)
  samples_clin <- as.character(clin$sample_id)
  common <- intersect(samples_expr, samples_clin)
  if (length(common) < cfg$params$min_samples) {
    msg <- paste("GATE G1 — BLOCKED", paste("Blocking issue: only", length(common), "samples match between expression and clinical; minimum is", cfg$params$min_samples), sep = "\n")
    writeLines(msg, file.path(cfg$paths$logs, "G1_BLOCKED.md"))
    log_decision("00_validate", "G1_BLOCKED", paste("matched samples:", length(common)))
    stop(msg, call. = FALSE)
  }

  clin <- clin[match(common, samples_clin), , drop = FALSE]
  expr <- expr[, common, drop = FALSE]

  ann <- data.frame(
    sample_id = clin$sample_id,
    repository_id = clin$repository_id,
    patient_id = clin$patient_id,
    group = clin$group,
    group_source = "clinical.tsv",
    disease_status = clin$disease_status,
    treatment_status = clin$treatment_status,
    timepoint = clin$timepoint,
    tissue = clin$tissue,
    cell_type = clin$cell_type,
    technology = clin$technology,
    age = if ("age" %in% colnames(clin)) clin$age else NA,
    sex = if ("sex" %in% colnames(clin)) clin$sex else NA,
    batch = clin$batch,
    platform = clin$platform,
    library_prep = clin$library_prep,
    biological_replicate_id = clin$patient_id,
    technical_replicate_id = NA_character_,
    ambiguity_flag = FALSE,
    ambiguity_note = "",
    quality_flag = "pass",
    stringsAsFactors = FALSE
  )

  critical <- c("sample_id", "patient_id", "group", "disease_status", "tissue")
  for (cc in critical) {
    bad <- is.na(ann[[cc]]) | ann[[cc]] == ""
    ann$ambiguity_flag[bad] <- TRUE
    ann$ambiguity_note[bad] <- paste0("missing critical field: ", cc)
  }

  data.table::fwrite(ann, file.path(cfg$paths$metadata, "Sample_Annotation.csv"))
  amb <- ann[ann$ambiguity_flag, c("sample_id", "repository_id", "ambiguity_note"), drop = FALSE]
  if (nrow(amb) > 0) {
    amb$conflict_field <- "critical_field"
    amb$values_found <- amb$ambiguity_note
    amb$sources <- "clinical.tsv"
    amb$why_it_matters <- "Cannot verify sample/group assignment"
    amb$status <- "open"
    amb$resolution <- "confirm/correct metadata"
    amb$resolved_by <- "user"
    amb$resolved_date <- ""
    amb <- amb[, c("sample_id", "repository_id", "conflict_field", "values_found", "sources", "why_it_matters", "status", "resolution", "resolved_by", "resolved_date")]
  } else {
    amb <- data.frame(sample_id = character(), repository_id = character(), conflict_field = character(), values_found = character(), sources = character(), why_it_matters = character(), status = character(), resolution = character(), resolved_by = character(), resolved_date = character())
  }
  data.table::fwrite(amb, file.path(cfg$paths$metadata, "Ambiguous_Samples.csv"))

  bap_py <- "/home/meng/.reasonix/skills/bioanalyst/scripts/bap.py"
  if (file.exists(bap_py)) {
    out <- system2("python3", c(bap_py, "check-annotation", file.path(cfg$paths$metadata, "Sample_Annotation.csv")), stdout = TRUE, stderr = TRUE)
    cat(paste(out, collapse = "\n"), "\n")
    st <- attr(out, "status")
    if (!is.null(st) && st != 0) {
      if (synthetic) {
        log_decision("00_validate", "SYNTHETIC_ANNOTATION_FAIL", "check-annotation failed; synthetic run continues for code validation")
      } else {
        log_decision("00_validate", "G1_BLOCKED", "bap.py check-annotation failed")
        stop("GATE G1 — BLOCKED: bap.py check-annotation failed", call. = FALSE)
      }
    }
    if (synthetic) {
      log_decision("00_validate", "SYNTHETIC_G1", "synthetic annotation checked; project G1 remains BLOCKED (no real data)")
    } else {
      gout <- system2("python3", c(bap_py, "gate", root, "--gate", "g1"), stdout = TRUE, stderr = TRUE)
      cat(paste(gout, collapse = "\n"), "\n")
      gst <- attr(gout, "status")
      if (!is.null(gst) && gst == 2) {
        log_decision("00_validate", "G1_BLOCKED", "bap.py gate g1 returned 2")
        stop("GATE G1 — BLOCKED: bap.py gate g1 returned 2", call. = FALSE)
      }
    }
  } else {
    log_decision("00_validate", "WARN", "bap.py not found; annotation schema checked in R only")
  }

  saveRDS(list(expression = expr, clinical = clin, annotation = ann), file.path(cfg$paths$processed, "validated_inputs.rds"))
  log_decision("00_validate", "PASS", paste("samples:", ncol(expr), "genes:", nrow(expr)))
  cat("G1 validated:", ncol(expr), "samples,", nrow(expr), "features\n")
  invisible(TRUE)
}
