#!/usr/bin/env Rscript
# 22c_build_gse38832_dfs.R
# Build GSE38832 (Tripathi et al., GPL570, 122 primary CRCs) as an INDEPENDENT DFS cohort.
#
# Endpoint discipline: identical to 22a — only DFS/RFS columns are written, and an assertion forbids
# os_*/css_* columns so the OS validation path (03_Raw_Data/validation/) can never absorb this cohort.
#
# Inputs : 03_Raw_Data/GSE38832/GSE38832_RAW.tar              (122 CEL.gz, verified 2026-11-17)
#          03_Raw_Data/GSE38832/GSE38832_series_matrix.txt.gz (metadata; REQUIRED for outcomes)
# Outputs: 03_Raw_Data/GSE38832/22c_GSE38832_expression.tsv   symbols x tumours (RMA, log2)
#          03_Raw_Data/GSE38832/22c_GSE38832_clinical_dfs.tsv DFS/RFS columns only
#          03_Raw_Data/GSE38832/22c_GSE38832_rma.rds          cached ExpressionSet (gitignored)
#
# Note on parallelism: affy::rma has no parallel backend, so CEL -> expression is single-threaded
# (BLAS pinned to 1); the parallel work happens in 22d (6 workers, one level).

suppressPackageStartupMessages({
  library(data.table)
  library(affy)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(hgu133plus2cdf)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

DIR <- "03_Raw_Data/GSE38832"
tar_path <- file.path(DIR, "GSE38832_RAW.tar")
series_path <- file.path(DIR, "GSE38832_series_matrix.txt.gz")
cel_dir <- file.path(DIR, "cel")
if (!file.exists(tar_path)) stop("Missing ", tar_path, call. = FALSE)

## 1. Extract CEL files (plain tar, no compression) ---------------------------
dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
cels <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", ignore.case = TRUE, full.names = TRUE)
if (length(cels) < 100) {
  cat("extracting CEL files ...\n")
  status <- utils::untar(tar_path, exdir = cel_dir, tar = "internal")
  if (status != 0) status <- system2("tar", c("-xf", shQuote(tar_path), "-C", shQuote(cel_dir)))
  cels <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", ignore.case = TRUE, full.names = TRUE)
}
cat("CEL files:", length(cels), "\n")
if (length(cels) == 0) stop("No CEL files extracted", call. = FALSE)

## 2. RMA (background-correct, quantile-normalise, log2) ---------------------
rma_path <- file.path(DIR, "22c_GSE38832_rma.rds")
if (file.exists(rma_path)) {
  eset <- readRDS(rma_path)
  cat("loaded cached RMA matrix\n")
} else {
  ab <- affy::ReadAffy(filenames = basename(cels), celfile.path = cel_dir)
  cat("arrays read:", length(ab), "| platform:", annotation(ab), "\n")
  eset <- affy::rma(ab)
  saveRDS(eset, rma_path)
}
expr <- Biobase::exprs(eset)
colnames(expr) <- sub("_.*$", "", sub("\\.CEL\\.gz$", "", colnames(expr), ignore.case = TRUE))
cat("RMA matrix:", nrow(expr), "probes x", ncol(expr), "arrays | range",
    paste(round(range(expr, na.rm = TRUE), 2), collapse = " .. "), "\n")

## 3. Probe -> symbol (highest mean probe per gene) ---------------------------
probe_mean <- rowMeans(expr, na.rm = TRUE)
map <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr),
                             columns = c("SYMBOL", "ENSEMBL"), keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", , drop = FALSE]
map$probe_mean <- probe_mean[match(map$PROBEID, rownames(expr))]
map <- map[order(-map$probe_mean), , drop = FALSE]
map <- map[!duplicated(map$SYMBOL), , drop = FALSE]
expr_gene <- expr[match(map$PROBEID, rownames(expr)), , drop = FALSE]
rownames(expr_gene) <- map$SYMBOL
cat("gene-level matrix:", nrow(expr_gene), "symbols x", ncol(expr_gene), "tumours\n")
expr_out <- data.frame(gene_id = rownames(expr_gene), expr_gene, check.names = FALSE)
data.table::fwrite(expr_out, file.path(DIR, "22c_GSE38832_expression.tsv"), sep = "\t")

## 4. Series matrix -> clinical ---------------------------------------------
if (!file.exists(series_path)) {
  cat("\nexpression built, but metadata is missing:", series_path, "\n")
  cat("Download (NCBI is unreachable from this host):\n")
  cat("  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38832/matrix/GSE38832_series_matrix.txt.gz\n")
  cat("into", DIR, "and re-run this script. Without outcomes the cohort cannot be analysed.\n")
  log_decision("22c_gse38832", "BLOCKED",
               paste0("expression built (", nrow(expr_gene), " symbols x ", ncol(expr_gene),
                      " arrays); series matrix metadata missing"))
  quit(status = 0)
}

lines <- readLines(gzfile(series_path), warn = FALSE)
acc <- strsplit(gsub('"', "", grep("^!Sample_geo_accession", lines, value = TRUE)), "\t")[[1]][-1]
ch_rows <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
cat("samples:", length(acc), "| characteristic rows:", length(ch_rows), "\n")

# characteristics may arrive as several key/value rows or one packed row; handle both
per_sample <- lapply(seq_along(acc), function(i) list())
for (row in ch_rows) {
  vals <- strsplit(gsub('"', "", row), "\t")[[1]][-1]
  if (length(vals) != length(acc)) next
  for (i in seq_along(vals)) {
    kv <- trimws(strsplit(vals[i], ";")[[1]])
    kv <- kv[grepl(":", kv)]
    for (p in kv) {
      k <- trimws(sub(":.*$", "", p))
      v <- trimws(sub("^[^:]*:\\s*", "", p))
      per_sample[[i]][[k]] <- v
    }
  }
}
keys <- sort(unique(unlist(lapply(per_sample, names))))
cat("discovered characteristic keys:\n"); print(keys)
fwrite(data.frame(key = keys), file.path(DIR, "22c_GSE38832_characteristic_keys.tsv"), sep = "\t")

pick <- function(pattern) {
  k <- grep(pattern, keys, ignore.case = TRUE, value = TRUE)
  if (length(k) == 0) return(rep(NA_character_, length(acc)))
  vapply(per_sample, function(x) if (!is.null(x[[k[1]]])) x[[k[1]]] else NA_character_, character(1))
}

time_raw <- pick("(dfs|rfs|recurrence|relapse).*(time|month|follow)|(time|month).*(dfs|rfs|recurrence|relapse)")
event_raw <- pick("(dfs|rfs|recurrence|relapse).*(event|status|cens)|(event|status|cens).*(dfs|rfs|recurrence|relapse)")
time_val <- suppressWarnings(as.numeric(time_raw))
event_txt <- tolower(event_raw)
event_val <- ifelse(is.na(event_txt), NA_integer_,
                    ifelse(grepl("^(1|yes|true|event|recurrence|relapse|dead|deceased)", event_txt), 1L,
                           ifelse(grepl("^(0|no|false|none|censor)", event_txt), 0L, NA_integer_)))

stage_raw <- pick("stage|dukes|tnm")
sex_raw <- pick("gender|sex")
age_raw <- pick("age")
cat("\nparsed DFS time range (raw units):",
    paste(round(range(time_val, na.rm = TRUE), 2), collapse = " .. "),
    "| median:", round(stats::median(time_val, na.rm = TRUE), 2), "\n")
cat("event table:\n"); print(table(event_val, useNA = "ifany"))

# unit decision: values around 10-150 are months (same convention as GSE14333/GSE17536);
# values around 1-15 are years. Decided from the observed median and logged.
med <- stats::median(time_val, na.rm = TRUE)
unit <- if (is.na(med)) NA_character_ else if (med > 20) "months" else "years"
factor_days <- if (identical(unit, "months")) 30.4375 else 365.25
cat("unit decision:", unit, "\n")

clinical <- data.frame(
  sample_id = acc,
  platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
  stage = stage_raw,
  age = suppressWarnings(as.numeric(age_raw)),
  sex = ifelse(grepl("^m", sex_raw, ignore.case = TRUE), "M",
               ifelse(grepl("^f", sex_raw, ignore.case = TRUE), "F", NA_character_)),
  DFS_time_raw = time_val,
  DFS_time_unit = unit,
  DFS_days = time_val * factor_days,
  DFS_event = event_val,
  stringsAsFactors = FALSE
)
clinical <- clinical[match(colnames(expr_gene), clinical$sample_id), , drop = FALSE]
stopifnot(!any(c("os_time_days", "os_event", "css_time_days", "rfs_time_days") %in% names(clinical)))

data.table::fwrite(clinical, file.path(DIR, "22c_GSE38832_clinical_dfs.tsv"), sep = "\t")
cat("\nDFS available:", sum(!is.na(clinical$DFS_days)), "of", nrow(clinical),
    "| events:", sum(clinical$DFS_event == 1, na.rm = TRUE), "\n")

log_decision("22c_gse38832", "BUILD",
             paste0("GSE38832 DFS cohort built: ", nrow(expr_gene), " symbols x ", ncol(expr_gene),
                    " arrays (RMA); DFS available for ", sum(!is.na(clinical$DFS_days)),
                    " (", sum(clinical$DFS_event == 1, na.rm = TRUE), " events); unit=", unit,
                    "; OS columns absent by assertion"))
cat("done\n")
