#!/usr/bin/env Rscript
# 11_build_gse39582_validation.R
# Build the GSE39582 external-validation inputs from the local series matrix.
# Probes are mapped to Ensembl via hgu133plus2.db; clinical fields are parsed from
# GEO characteristics. Survival is derived from os.event/os.delay (months).

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(data.table)
})

series_path <- "03_Raw_Data/GSE39582/GSE39582_series_matrix.txt.gz"
outdir <- "03_Raw_Data/validation/GSE39582"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(series_path)) stop("Missing ", series_path, call. = FALSE)

g <- GEOquery::getGEO(filename = series_path)
expr <- Biobase::exprs(g)
pd <- Biobase::pData(g)

# Probe -> Ensembl mapping; keep the highest-mean probe per Ensembl
map <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr), columns = c("ENSEMBL", "SYMBOL"), keytype = "PROBEID")
map <- map[!is.na(map$ENSEMBL) & map$ENSEMBL != "", , drop = FALSE]
probe_mean <- rowMeans(expr, na.rm = TRUE)
map$probe_mean <- probe_mean[match(map$PROBEID, names(probe_mean))]
map <- map[order(-map$probe_mean), , drop = FALSE]
map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
expr_gene <- expr[map$PROBEID, , drop = FALSE]
rownames(expr_gene) <- map$ENSEMBL
colnames(expr_gene) <- pd$title

# Parse characteristics fields
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

sample_id <- as.character(pd$title)
dataset <- get_field("^dataset")
sex_raw <- get_field("^Sex")
age_raw <- get_field("age\\.at\\.diagnosis")
stage_raw <- get_field("tnm\\.stage")
location <- get_field("tumor\\.location")
mmr <- get_field("mmr\\.status")
os_event <- suppressWarnings(as.integer(get_field("os\\.event")))
os_delay <- suppressWarnings(as.numeric(get_field("os\\.delay")))
rfs_event <- suppressWarnings(as.integer(get_field("rfs\\.event")))
rfs_delay <- suppressWarnings(as.numeric(get_field("rfs\\.delay")))
chemo <- get_field("chemotherapy\\.adjuvant")
cms <- get_field("cit\\.molecularsubtype")
source <- as.character(pd$source_name_ch1)

stage_num <- suppressWarnings(as.integer(stage_raw))
stage <- ifelse(is.na(stage_num) | stage_num < 1 | stage_num > 4, NA_character_, c("I", "II", "III", "IV")[stage_num])
sex <- ifelse(grepl("^M", sex_raw, ignore.case = TRUE), "M", ifelse(grepl("^F", sex_raw, ignore.case = TRUE), "F", NA_character_))
msi <- ifelse(grepl("pMMR", mmr, ignore.case = TRUE), "MSS", ifelse(grepl("dMMR", mmr, ignore.case = TRUE), "MSI", NA_character_))
side <- ifelse(grepl("proximal", location, ignore.case = TRUE), "right", ifelse(grepl("distal", location, ignore.case = TRUE), "left", NA_character_))
group <- ifelse(grepl("normal|healthy", source, ignore.case = TRUE), "control", "tumour")

clinical <- data.frame(
  sample_id = sample_id,
  patient_id = sample_id,
  group = group,
  group_source = "GEO series matrix",
  disease_status = "CRC",
  treatment_status = chemo,
  timepoint = "baseline",
  tissue = "colon",
  cell_type = "epithelial",
  technology = "microarray",
  platform = "Affymetrix HG-U133 Plus 2.0",
  library_prep = "unknown",
  repository_id = as.character(pd$geo_accession),
  stage = stage,
  msi_status = msi,
  sidedness = side,
  age = suppressWarnings(as.numeric(age_raw)),
  sex = sex,
  batch = ifelse(is.na(dataset), "GSE39582", dataset),
  purity = NA_real_,
  os_time_days = os_delay * 30.4375,
  os_event = os_event,
  css_time_days = NA_real_,
  css_event = NA_integer_,
  rfs_time_days = rfs_delay * 30.4375,
  rfs_event = rfs_event,
  cms = cms,
  stringsAsFactors = FALSE
)

# Keep samples with expression columns
clinical <- clinical[match(colnames(expr_gene), clinical$sample_id), , drop = FALSE]
clinical <- clinical[!is.na(clinical$sample_id), , drop = FALSE]
expr_gene <- expr_gene[, clinical$sample_id, drop = FALSE]

expr_out <- data.frame(gene_id = rownames(expr_gene), expr_gene, check.names = FALSE)
data.table::fwrite(expr_out, file.path(outdir, "expression.tsv"), sep = "\t")
data.table::fwrite(clinical, file.path(outdir, "clinical.tsv"), sep = "\t")

cat("GSE39582 validation built\n")
cat("expression:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")
cat("tumours:", sum(clinical$group == "tumour"), "controls:", sum(clinical$group == "control"), "\n")
cat("OS events:", sum(clinical$os_event, na.rm = TRUE), "OS time missing:", sum(is.na(clinical$os_time_days)), "\n")
cat("stage missing:", sum(is.na(clinical$stage)), "MSI:", paste(names(table(clinical$msi_status)), table(clinical$msi_status), sep = "=", collapse = ", "), "\n")
cat("wrote:", outdir, "\n")
