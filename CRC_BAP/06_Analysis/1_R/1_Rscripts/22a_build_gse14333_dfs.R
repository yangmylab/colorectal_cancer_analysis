#!/usr/bin/env Rscript
# 22a_build_gse14333_dfs.R
# Build GSE14333 (Jorissen et al., GPL570) as a DFS-endpoint cohort.
#
# Endpoint discipline: GSE14333 carries DFS only, so this builder writes DFS columns and NOTHING else
# (no os_*, css_* or rfs_* columns). The OS validation path in 04_survival_validation.R reads
# 03_Raw_Data/validation/<cohort>/ and therefore can never pick this cohort up by accident.
#
# Input : 03_Raw_Data/GSE14333/GSE14333_series_matrix.txt.gz
# Output: 03_Raw_Data/GSE14333/22a_GSE14333_expression.tsv   symbols x tumours, log2 MAS5/quantile
#         03_Raw_Data/GSE14333/22a_GSE14333_clinical_dfs.tsv  explicit DFS_* columns only
#
# Parallelism: none needed here (probe mapping is a single vectorised AnnotationDbi join over 54,674
# probes, ~10 s); heavy steps run 6 workers elsewhere, BLAS pinned to 1 thread per worker.

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

DIR <- "03_Raw_Data/GSE14333"
series_path <- file.path(DIR, "GSE14333_series_matrix.txt.gz")
if (!file.exists(series_path)) stop("Missing ", series_path, call. = FALSE)

## 1. Series matrix: metadata + table block ----------------------------------
lines <- readLines(gzfile(series_path), warn = FALSE)
tab_start <- grep("^!series_matrix_table_begin", lines)
if (length(tab_start) != 1) stop("Cannot locate the expression table block", call. = FALSE)
n_probes <- length(lines) - tab_start - 3
expr_dt <- data.table::fread(series_path, skip = tab_start, nrows = n_probes, data.table = FALSE, check.names = FALSE)
cat("expression block:", nrow(expr_dt), "probes x", ncol(expr_dt) - 1, "samples\n")

ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
if (length(ch) != 1) stop("Expected exactly one characteristics row", call. = FALSE)
vals <- strsplit(gsub('"', "", ch), "\t")[[1]][-1]
acc <- strsplit(gsub('"', "", grep("^!Sample_geo_accession", lines, value = TRUE)), "\t")[[1]][-1]
title <- strsplit(gsub('"', "", grep("^!Sample_title", lines, value = TRUE)), "\t")[[1]][-1]
cat("samples:", length(vals), "| geo accessions:", length(acc), "\n")
stopifnot(length(vals) == ncol(expr_dt) - 1)

getf <- function(key) {
  hit <- grepl(paste0(key, ":"), vals)
  out <- rep(NA_character_, length(vals))
  out[hit] <- sub(paste0(".*", key, ": *(.*?);.*$"), "\\1", paste0(vals[hit], ";"))
  out <- trimws(out)
  out[out %in% c("NA", "")] <- NA_character_
  out
}

## 2. Probe -> symbol (keep the highest-expressed probe per gene) -------------
mat <- as.matrix(expr_dt[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
probe_mean <- rowMeans(mat, na.rm = TRUE)
map <- AnnotationDbi::select(hgu133plus2.db, keys = expr_dt[[1]],
                             columns = c("SYMBOL", "ENSEMBL"), keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
map$probe_mean <- probe_mean[match(map$PROBEID, expr_dt[[1]])]
map <- map[order(-map$probe_mean), , drop = FALSE]
map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
expr_gene <- mat[match(map$PROBEID, expr_dt[[1]]), , drop = FALSE]
rownames(expr_gene) <- map$SYMBOL
colnames(expr_gene) <- acc
cat("gene-level expression:", nrow(expr_gene), "symbols x", ncol(expr_gene), "tumours\n")
cat("value range:", paste(round(range(expr_gene, na.rm = TRUE), 2), collapse = " .. "),
    "(already log2; no transform applied)\n")

## 3. Clinical: DFS only -----------------------------------------------------
duke <- getf("DukesStage")
stage <- c(A = "I", B = "II", C = "III", D = "IV")[duke]
loc <- getf("Location")
dfs_months <- suppressWarnings(as.numeric(getf("DFS_Time")))
dfs_cens <- suppressWarnings(as.integer(getf("DFS_Cens")))
age <- suppressWarnings(as.numeric(getf("Age_Diag")))
sex_raw <- getf("Gender")

clinical <- data.frame(
  sample_id = acc,
  patient_id = title,
  platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
  stage = unname(stage),
  duke = duke,
  sidedness = ifelse(loc == "Right", "right",
                     ifelse(loc == "Left", "left",
                            ifelse(loc == "Rectum", "rectum",
                                   ifelse(loc == "Colon", "colon", NA_character_)))),
  age = age,
  sex = ifelse(grepl("^M", sex_raw, ignore.case = TRUE), "M",
               ifelse(grepl("^F", sex_raw, ignore.case = TRUE), "F", NA_character_)),
  adjuvant_rt = getf("AdjXRT"),
  adjuvant_chemo = getf("AdjCTX"),
  DFS_time_months = dfs_months,
  DFS_days = dfs_months * 30.4375,
  DFS_event = ifelse(is.na(dfs_cens), NA_integer_, ifelse(dfs_cens == 1L, 0L, 1L)),
  stringsAsFactors = FALSE
)
clinical <- clinical[match(colnames(expr_gene), clinical$sample_id), , drop = FALSE]
stopifnot(!any(c("os_time_days", "os_event", "css_time_days", "rfs_time_days") %in% names(clinical)))

## 4. Write ------------------------------------------------------------------
expr_out <- data.frame(gene_id = rownames(expr_gene), expr_gene, check.names = FALSE)
data.table::fwrite(expr_out, file.path(DIR, "22a_GSE14333_expression.tsv"), sep = "\t")
data.table::fwrite(clinical, file.path(DIR, "22a_GSE14333_clinical_dfs.tsv"), sep = "\t")

cat("\nDFS availability:", sum(!is.na(clinical$DFS_days)), "of", nrow(clinical), "tumours\n")
cat("DFS events:", sum(clinical$DFS_event == 1, na.rm = TRUE),
    "| censored:", sum(clinical$DFS_event == 0, na.rm = TRUE), "\n")
cat("DFS follow-up (months) median:", round(stats::median(clinical$DFS_time_months, na.rm = TRUE), 2),
    "| range:", paste(round(range(clinical$DFS_time_months, na.rm = TRUE), 2), collapse = " .. "), "\n")
cat("stage:\n"); print(table(clinical$stage, useNA = "ifany"))
cat("sidedness:\n"); print(table(clinical$sidedness, useNA = "ifany"))
cat("DFS events by stage:\n"); print(table(clinical$stage, clinical$DFS_event, useNA = "ifany"))

log_decision("22a_gse14333", "BUILD",
             paste0("GSE14333 DFS cohort built: ", nrow(expr_gene), " symbols x ", ncol(expr_gene),
                    " tumours; DFS available for ", sum(!is.na(clinical$DFS_days)),
                    " (", sum(clinical$DFS_event == 1, na.rm = TRUE), " events); OS columns deliberately absent"))
cat("done\n")
