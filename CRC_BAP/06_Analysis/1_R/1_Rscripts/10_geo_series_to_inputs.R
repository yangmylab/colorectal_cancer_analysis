#!/usr/bin/env Rscript
# 10_geo_series_to_inputs.R
# Convert a locally downloaded GEO series matrix into the CRC_BAP input format.
# This script does NOT download data and does NOT invent survival values.
#
# Usage:
#   Rscript 10_geo_series_to_inputs.R --series-matrix <file> --outdir <dir>

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

series_path <- get_arg("--series-matrix")
outdir <- get_arg("--outdir", "03_Raw_Data")

if (is.null(series_path) || !file.exists(series_path)) {
  stop("Provide --series-matrix <existing file>", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

g <- GEOquery::getGEO(filename = series_path)
expr <- Biobase::exprs(g)
pd <- Biobase::pData(g)

expr_df <- data.frame(gene_id = rownames(expr), expr, check.names = FALSE)
data.table::fwrite(expr_df, file.path(outdir, "expression.tsv"), sep = "\t")

pick_col <- function(patterns) {
  hits <- unique(unlist(lapply(patterns, function(p) grep(p, colnames(pd), ignore.case = TRUE, value = TRUE))))
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

clean_value <- function(x) {
  if (length(x) == 0 || is.na(x)) return(rep(NA_character_, nrow(pd)))
  x <- as.character(x)
  sub("^[^:]+:\\s*", "", x)
}

col_sample <- if ("geo_accession" %in% colnames(pd)) "geo_accession" else rownames(pd)
col_patient <- pick_col(c("patient", "subject", "individual", "donor"))
col_tissue <- pick_col(c("tissue", "source_name_ch1", "organ"))
col_disease <- pick_col(c("disease", "diagnosis", "condition", "group"))
col_stage <- pick_col(c("stage", "tnm"))
col_msi <- pick_col(c("msi", "microsatellite", "mmr"))
col_side <- pick_col(c("sidedness", "side", "location"))
col_age <- pick_col(c("age"))
col_sex <- pick_col(c("sex", "gender"))
col_batch <- pick_col(c("batch", "plate", "scan"))

sample_id <- if (col_sample == "geo_accession") as.character(pd$geo_accession) else as.character(col_sample)
patient_id <- if (!is.na(col_patient)) clean_value(pd[[col_patient]]) else sample_id
tissue <- if (!is.na(col_tissue)) clean_value(pd[[col_tissue]]) else NA_character_
disease <- if (!is.na(col_disease)) clean_value(pd[[col_disease]]) else NA_character_
stage <- if (!is.na(col_stage)) clean_value(pd[[col_stage]]) else NA_character_
msi <- if (!is.na(col_msi)) clean_value(pd[[col_msi]]) else NA_character_
side <- if (!is.na(col_side)) clean_value(pd[[col_side]]) else NA_character_
age <- if (!is.na(col_age)) suppressWarnings(as.numeric(clean_value(pd[[col_age]]))) else NA_real_
sex <- if (!is.na(col_sex)) clean_value(pd[[col_sex]]) else NA_character_
batch <- if (!is.na(col_batch)) clean_value(pd[[col_batch]]) else "batch_1"

group <- ifelse(grepl("normal|healthy|control|adjacent", disease, ignore.case = TRUE), "control",
  ifelse(grepl("tumour|tumor|cancer|carcinoma|adenoma", disease, ignore.case = TRUE), "tumour", disease)
)

clinical <- data.frame(
  sample_id = sample_id,
  patient_id = patient_id,
  group = group,
  group_source = "GEO series matrix",
  disease_status = disease,
  treatment_status = "na",
  timepoint = "baseline",
  tissue = tissue,
  cell_type = "epithelial",
  technology = "GEO",
  platform = if ("platform_id" %in% colnames(pd)) as.character(pd$platform_id) else "GEO",
  library_prep = "unknown",
  repository_id = if ("geo_accession" %in% colnames(pd)) as.character(pd$geo_accession) else NA_character_,
  stage = stage,
  msi_status = msi,
  sidedness = side,
  age = age,
  sex = sex,
  batch = batch,
  purity = NA_real_,
  os_time_days = NA_real_,
  os_event = NA_integer_,
  css_time_days = NA_real_,
  css_event = NA_integer_,
  rfs_time_days = NA_real_,
  rfs_event = NA_integer_,
  stringsAsFactors = FALSE
)

data.table::fwrite(clinical, file.path(outdir, "clinical_available.tsv"), sep = "\t")
data.table::fwrite(clinical, file.path(outdir, "clinical_template.tsv"), sep = "\t")

missing <- c(
  if (anyNA(clinical$tissue) || all(clinical$tissue == "")) "tissue",
  if (anyNA(clinical$disease_status) || all(clinical$disease_status == "")) "disease_status",
  if (anyNA(clinical$group) || all(clinical$group == "")) "group",
  if (all(is.na(clinical$os_event))) "os_time_days/os_event",
  if (all(is.na(clinical$css_event))) "css_time_days/css_event",
  if (all(is.na(clinical$rfs_event))) "rfs_time_days/rfs_event"
)

message("Wrote: ", file.path(outdir, "expression.tsv"))
message("Wrote: ", file.path(outdir, "clinical_available.tsv"))
message("Wrote: ", file.path(outdir, "clinical_template.tsv"))
if (length(missing) > 0) {
  message("Fill before running: ", paste(missing, collapse = ", "))
} else {
  message("All required fields appear present; review and save as clinical.tsv.")
}
message("Do not fabricate survival values. Use the publication supplementary table.")
