#!/usr/bin/env Rscript
# 12_build_gse17536_validation.R
# Build GSE17536 (Moffitt, GPL570) as an independent validation cohort.
# Times in the GEO characteristics are months; converted to days.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(data.table)
})

series_path <- "03_Raw_Data/GSE17536/GSE17536_series_matrix.txt.gz"
outdir <- "03_Raw_Data/validation/GSE17536"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(series_path)) stop("Missing ", series_path, call. = FALSE)

g <- GEOquery::getGEO(filename = series_path)
expr <- Biobase::exprs(g)
pd <- Biobase::pData(g)

map <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr), columns = c("ENSEMBL", "SYMBOL"), keytype = "PROBEID")
map <- map[!is.na(map$ENSEMBL) & map$ENSEMBL != "", , drop = FALSE]
probe_mean <- rowMeans(expr, na.rm = TRUE)
map$probe_mean <- probe_mean[match(map$PROBEID, names(probe_mean))]
map <- map[order(-map$probe_mean), , drop = FALSE]
map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
expr_gene <- expr[map$PROBEID, , drop = FALSE]
rownames(expr_gene) <- map$ENSEMBL
colnames(expr_gene) <- rownames(pd)

char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
get_field <- function(pattern) {
  vals <- rep(NA_character_, nrow(pd))
  for (cc in char_cols) {
    x <- as.character(pd[[cc]])
    idx <- grepl(pattern, x, ignore.case = TRUE)
    vals[idx] <- sub("^[^:]+:\\s*", "", x[idx])
  }
  vals
}

event01 <- function(x) {
  ifelse(grepl("no death|no recurrence|no event|alive", x, ignore.case = TRUE), 0L,
    ifelse(grepl("death|recurrence|event", x, ignore.case = TRUE), 1L, NA_integer_))
}

sample_id <- rownames(pd)
age <- suppressWarnings(as.numeric(get_field("^age")))
sex_raw <- get_field("gender|^sex")
stage_raw <- suppressWarnings(as.integer(get_field("ajcc_stage|tnm\\.stage|^stage")))
os_event <- event01(get_field("overall_event"))
os_time <- suppressWarnings(as.numeric(get_field("overall survival follow-up time")))
css_event <- event01(get_field("dss_event"))
css_time <- suppressWarnings(as.numeric(get_field("dss_time")))
rfs_event <- event01(get_field("dfs_event"))
rfs_time <- suppressWarnings(as.numeric(get_field("dfs_time")))

clinical <- data.frame(
  sample_id = sample_id,
  patient_id = as.character(pd$title),
  group = "tumour",
  group_source = "GEO series matrix",
  disease_status = "CRC",
  treatment_status = NA_character_,
  timepoint = "baseline",
  tissue = "colon",
  cell_type = "epithelial",
  technology = "microarray",
  platform = "Affymetrix HG-U133 Plus 2.0",
  library_prep = "unknown",
  repository_id = as.character(pd$geo_accession),
  stage = ifelse(is.na(stage_raw) | stage_raw < 1 | stage_raw > 4, NA_character_, c("I", "II", "III", "IV")[stage_raw]),
  msi_status = NA_character_,
  sidedness = NA_character_,
  age = age,
  sex = ifelse(grepl("^M", sex_raw, ignore.case = TRUE), "M", ifelse(grepl("^F", sex_raw, ignore.case = TRUE), "F", NA_character_)),
  batch = "GSE17536",
  purity = NA_real_,
  os_time_days = os_time * 30.4375,
  os_event = os_event,
  css_time_days = css_time * 30.4375,
  css_event = css_event,
  rfs_time_days = rfs_time * 30.4375,
  rfs_event = rfs_event,
  stringsAsFactors = FALSE
)

clinical <- clinical[match(colnames(expr_gene), clinical$sample_id), , drop = FALSE]
expr_out <- data.frame(gene_id = rownames(expr_gene), expr_gene, check.names = FALSE)
data.table::fwrite(expr_out, file.path(outdir, "expression.tsv"), sep = "\t")
data.table::fwrite(clinical, file.path(outdir, "clinical.tsv"), sep = "\t")

cat("GSE17536 validation built\n")
cat("expression:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")
cat("OS events:", sum(clinical$os_event, na.rm = TRUE), "OS missing:", sum(is.na(clinical$os_time_days)), "\n")
cat("CSS events:", sum(clinical$css_event, na.rm = TRUE), "RFS events:", sum(clinical$rfs_event, na.rm = TRUE), "\n")
cat("stage missing:", sum(is.na(clinical$stage)), "\n")
