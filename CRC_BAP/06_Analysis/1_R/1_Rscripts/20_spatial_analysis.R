#!/usr/bin/env Rscript
# 20_spatial_analysis.R
# Spatial test of the frozen K=2 ecotypes on public Visium CRC sections.
#
# Scores are expressed in the BULK reference space: every programme score is the mean, over the
# member genes, of gene-wise z-scores derived from the TCGA discovery log-expression matrix.
# The same statistic is computed for the 624 bulk tumours, so spot scores and bulk scores are on one
# scale and are NOT standardised within a section (that would force a ~50/50 call in every section).
#
# Datasets: GSE283052 (12 primary tumours), GSE267401 (2 patients, primary/metastasis pairs).
#
# Output: 06_Analysis/10_Spatial/*.tsv   09_Figures/plots/plot_25 .. plot_27

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(mclust)
})

WORK <- "04_Processed_Data/spatial"
OUT <- "06_Analysis/10_Spatial"
FIG <- "09_Figures/plots"
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source("06_Analysis/1_R/1_Rscripts/utils.R")

save_plot <- function(plot_fun, file_base, width = 7, height = 5) {
  # The composites are laid out at a fixed printed width, so a panel's authored size decides its printed
  # type size: a panel authored at 3.4 in prints its 11 pt text at >=7 pt in a two-column figure.
  # CRC_PANEL_MAX_W caps the authored width (inches; 0 = off) and CRC_PANEL_SCALE multiplies both
  # dimensions, so a re-authored panel set needs no edits to the individual plot calls.
  scale <- as.numeric(Sys.getenv("CRC_PANEL_SCALE", "1"))
  max_w <- as.numeric(Sys.getenv("CRC_PANEL_MAX_W", "0"))
  width <- width * scale
  height <- height * scale
  if (max_w > 0 && width > max_w) {
    shrink <- max_w / width
    width <- width * shrink
    height <- height * shrink
  }
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height)
  plot_fun()
  grDevices::dev.off()
  grDevices::jpeg(paste0(file_base, ".jpg"), width = width, height = height, units = "in", res = 300, quality = 95)
  plot_fun()
  grDevices::dev.off()
}

## 1. Bulk reference: gene-wise mean/SD and the same programme statistic -------
obj <- readRDS("04_Processed_Data/processed_bulk.rds")
bulk_expr <- obj$log_expr
gene_mu <- rowMeans(bulk_expr, na.rm = TRUE)
gene_sd <- apply(bulk_expr, 1, stats::sd, na.rm = TRUE)
gene_sd[!is.finite(gene_sd) | gene_sd < 1e-6] <- 1
cat("bulk reference:", nrow(bulk_expr), "genes x", ncol(bulk_expr), "samples\n")

msig <- msigdbr::msigdbr(db_species = "HS", collection = "H")
gene_sets <- split(msig$gene_symbol, msig$gs_name)
sym2ens <- unique(msig[, c("gene_symbol", "ensembl_gene")])
sym2ens <- sym2ens[!is.na(sym2ens$ensembl_gene) & sym2ens$ensembl_gene != "", , drop = FALSE]
sym2ens <- sym2ens[!duplicated(sym2ens$gene_symbol), , drop = FALSE]

bulk_z <- (bulk_expr - gene_mu) / gene_sd
ens_sets <- lapply(gene_sets, function(g) {
  e <- sym2ens$ensembl_gene[match(g, sym2ens$gene_symbol)]
  intersect(unique(e[!is.na(e)]), rownames(bulk_z))
})
ens_sets <- ens_sets[lengths(ens_sets) >= 5]
bulk_scores <- do.call(cbind, lapply(ens_sets, function(g) colMeans(bulk_z[g, , drop = FALSE], na.rm = TRUE)))
cat("bulk programme scores:", nrow(bulk_scores), "samples x", ncol(bulk_scores), "programmes\n")

labels_bulk <- fread("06_Analysis/03_Clustering/ecotype_labels.tsv", data.table = FALSE)
bulk_label <- labels_bulk$ecotype[match(rownames(bulk_scores), labels_bulk$sample_id)]
ok_bulk <- !is.na(bulk_label)
centroids <- t(sapply(c("E1", "E2"), function(ec) colMeans(bulk_scores[ok_bulk & bulk_label == ec, , drop = FALSE])))
rownames(centroids) <- c("E1", "E2")   # rows = ecotype, columns = programme

# honesty check: does this simplified statistic reproduce the GSVA-based labels?
call_bulk <- rownames(centroids)[apply(cor(t(bulk_scores[ok_bulk, , drop = FALSE]), t(centroids)), 1, which.max)]
agree_bulk <- mean(call_bulk == bulk_label[ok_bulk])
ari_bulk <- mclust::adjustedRandIndex(call_bulk, bulk_label[ok_bulk])
cat(sprintf("bulk self-check: agreement with GSVA labels %.3f, ARI %.3f\n", agree_bulk, ari_bulk))

## 2. Archive handling -------------------------------------------------------
dataset_root <- function(dataset) {
  for (cand in c(file.path("03_Raw_Data", dataset), file.path("03_Raw_Data", "spatial", dataset))) {
    if (file.exists(file.path(cand, paste0(dataset, "_RAW.tar")))) return(cand)
  }
  NULL
}

datasets <- c("GSE283052", "GSE267401")  # GSE284061 excluded: archive corrupt at GEO
roots <- setNames(lapply(datasets, dataset_root), datasets)
datasets <- datasets[!vapply(roots, is.null, logical(1))]
if (length(datasets) == 0) stop("no spatial archive found", call. = FALSE)

extract_members <- function(dataset) {
  root <- roots[[dataset]]
  tar_file <- file.path(root, paste0(dataset, "_RAW.tar"))
  work <- file.path(WORK, dataset)
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  members <- system2("tar", c("-tf", shQuote(tar_file)), stdout = TRUE, stderr = FALSE)
  want <- members[grepl("matrix\\.mtx\\.gz$|tissue_positions\\.csv\\.gz$|tissue_positions_list\\.csv\\.gz$|features\\.tsv\\.gz$|barcodes\\.tsv\\.gz$", members)]
  todo <- want[!file.exists(file.path(work, basename(want)))]
  if (length(todo) > 0) system2("tar", c("-xf", shQuote(tar_file), "-C", shQuote(work), shQuote(todo)), stdout = FALSE, stderr = FALSE)
  list(root = root, work = work, members = want)
}

read_positions <- function(path) {
  pos <- fread(path, header = FALSE, data.table = FALSE, fill = TRUE)
  if (is.character(pos[[1]][1]) && grepl("barcode", pos[[1]][1], ignore.case = TRUE)) {
    pos <- fread(path, header = TRUE, data.table = FALSE)
  } else {
    colnames(pos)[seq_len(4)] <- c("barcode", "in_tissue", "array_row", "array_col")
  }
  pos
}

## 3. Per-section scoring in bulk space --------------------------------------
score_section <- function(dataset, ex, section) {
  work <- ex$work
  mtx_path <- file.path(work, paste0(section, "matrix.mtx.gz"))
  pos_path <- file.path(work, paste0(section, "tissue_positions_list.csv.gz"))
  if (!file.exists(pos_path)) pos_path <- file.path(work, paste0(section, "tissue_positions.csv.gz"))
  feat_path <- file.path(work, paste0(section, "features.tsv.gz"))
  if (!file.exists(feat_path)) {
    sf <- list.files(ex$root, pattern = "features\\.tsv\\.gz$", full.names = TRUE)
    sf <- sf[!grepl("^GSM", basename(sf))]
    if (length(sf) == 0) sf <- list.files(ex$root, pattern = "features\\.tsv\\.gz$", full.names = TRUE)
    if (length(sf) > 0) feat_path <- sf[1]
  }
  bc_path <- file.path(work, paste0(section, "barcodes.tsv.gz"))
  if (!file.exists(bc_path)) {
    sb <- list.files(ex$root, pattern = paste0("^", dataset, "_barcodes\\.tsv\\.gz$"), full.names = TRUE)
    if (length(sb) > 0) bc_path <- sb[1]
  }
  if (!all(file.exists(c(mtx_path, pos_path, feat_path, bc_path)))) {
    cat("  skip", section, "(missing matrix/positions/features/barcodes)\n")
    return(NULL)
  }

  m <- Matrix::readMM(mtx_path)
  feat <- fread(feat_path, header = FALSE, data.table = FALSE)
  ids <- as.character(if (ncol(feat) >= 2) feat[[1]] else feat[[1]])
  symbols <- as.character(if (ncol(feat) >= 2) feat[[2]] else feat[[1]])
  if (nrow(m) != length(ids)) {
    cat("  skip", section, "(gene count mismatch)\n")
    return(NULL)
  }
  barcodes <- fread(bc_path, header = FALSE, data.table = FALSE)[[1]]
  if (length(barcodes) != ncol(m)) {
    cat("  skip", section, "(barcode/matrix mismatch:", length(barcodes), "vs", ncol(m), ")\n")
    return(NULL)
  }
  colnames(m) <- as.character(barcodes)

  # features may be Ensembl or symbols
  if (any(grepl("^ENSG", head(ids, 100)))) {
    gene_id <- sub("\\..*$", "", ids)
  } else {
    gene_id <- sym2ens$ensembl_gene[match(ids, sym2ens$gene_symbol)]
  }
  keep_gene <- !is.na(gene_id) & gene_id %in% rownames(bulk_z)
  if (sum(keep_gene) < 2000) {
    cat("  skip", section, "(only", sum(keep_gene), "genes mapped to the bulk reference)\n")
    return(NULL)
  }
  rownames(m) <- make.unique(ifelse(is.na(gene_id), ids, gene_id))

  pos <- read_positions(pos_path)
  pos <- pos[as.character(pos$in_tissue) %in% c("1"), , drop = FALSE]
  pos <- pos[pos$barcode %in% colnames(m), , drop = FALSE]
  if (nrow(pos) < 50) {
    cat("  skip", section, "(too few in-tissue spots:", nrow(pos), ")\n")
    return(NULL)
  }
  m <- m[, match(pos$barcode, colnames(m)), drop = FALSE]

  lib <- Matrix::colSums(m)
  lib[lib == 0] <- 1
  norm <- log1p(Matrix::t(Matrix::t(m) / lib * 1e4))

  scores <- sapply(ens_sets, function(g) {
    g <- intersect(g, rownames(norm))
    if (length(g) < 5) return(rep(NA_real_, ncol(norm)))
    sub <- as.matrix(norm[g, , drop = FALSE])
    z <- (sub - gene_mu[g]) / gene_sd[g]
    colMeans(z, na.rm = TRUE)
  })
  scores <- scores[, colSums(is.na(scores)) == 0, drop = FALSE]
  shared <- intersect(colnames(centroids), colnames(scores))
  if (length(shared) < 20) {
    cat("  skip", section, "(only", length(shared), "shared programmes)\n")
    return(NULL)
  }
  cen <- centroids[, shared, drop = FALSE]
  sim <- cor(t(scores[, shared, drop = FALSE]), t(cen), use = "pairwise.complete.obs")
  labels <- rownames(cen)[apply(sim, 1, which.max)]

  cat(sprintf("  scored %s: %d spots x %d genes, %.0f%% E2\n", section, nrow(pos), sum(keep_gene), 100 * mean(labels == "E2")))
  list(dataset = dataset, sample = section, labels = labels, sim = sim, barcodes = pos$barcode,
       row = as.numeric(pos$array_row), col = as.numeric(pos$array_col), scores = scores[, shared, drop = FALSE],
       stats = data.frame(dataset = dataset, sample = section, spots = nrow(pos),
                          genes = sum(keep_gene), fraction_E2 = mean(labels == "E2"),
                          mean_E2_similarity = mean(sim[, "E2"])))
}

res <- list()
for (dataset in datasets) {
  ex <- extract_members(dataset)
  sections <- sub("matrix\\.mtx\\.gz$", "", basename(ex$members[grepl("matrix\\.mtx\\.gz$", ex$members)]))
  cat(dataset, ": sections found =", length(sections), "\n")
  for (s in sections) {
    out <- try(score_section(dataset, ex, s), silent = TRUE)
    if (inherits(out, "try-error")) cat("  error in", s, ":", conditionMessage(attr(out, "condition")), "\n")
    else if (!is.null(out)) res[[length(res) + 1]] <- out
  }
}
if (length(res) == 0) stop("no spatial sections could be scored", call. = FALSE)
cat("scored sections:", length(res), "\n")

stats <- do.call(rbind, lapply(res, `[[`, "stats"))
saveRDS(res, file.path(OUT, "spatial_scores.rds"))

## Ecotype transfer diagnostics ---------------------------------------------
# The frozen bulk classifier is applied with two normalisations; both are reported so the
# failure mode is visible rather than hidden.
scores_filter <- function(m) {
  keep <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING")
  intersect(keep, colnames(m))
}
e2_progs <- scores_filter(res[[1]]$scores)

sec_scores <- do.call(rbind, lapply(res, function(r) colMeans(r$scores)))
rownames(sec_scores) <- stats$sample
bulk_e2 <- rowMeans(bulk_scores[, e2_progs, drop = FALSE])
sec_e2 <- rowMeans(sec_scores[, e2_progs, drop = FALSE])

transfer <- data.frame(
  analysis = c("bulk self-check (simplified statistic vs GSVA labels)",
               "spot calls: bulk-referenced normalisation",
               "spot calls: per-section normalisation"),
  statistic = c(sprintf("agreement %.3f; ARI %.3f", agree_bulk, ari_bulk),
                sprintf("%.0f%% of spots called E2", 100 * mean(stats$fraction_E2)),
                "not computed (forces ~50/50 by construction; invalid)"),
  verdict = c("statistic is not a faithful transfer of the classifier",
              "all spots fall on one side of the decision boundary (covariate shift)",
              "discarded"),
  stringsAsFactors = FALSE
)
fwrite(transfer, file.path(OUT, "spatial_transfer_diagnostics.tsv"), sep = "\t")
print(transfer)

spot_e2 <- do.call(c, lapply(res, function(r) rowMeans(r$scores[, e2_progs, drop = FALSE])))
summary_tab <- data.frame(
  dataset = stats$dataset,
  sample = stats$sample,
  spots = stats$spots,
  genes = stats$genes,
  section_E2_programme = round(as.numeric(sec_e2[stats$sample]), 3),
  stringsAsFactors = FALSE
)
fwrite(summary_tab, file.path(OUT, "spatial_section_programme_summary.tsv"), sep = "\t")

shift <- median(spot_e2) - median(bulk_e2[ok_bulk])
cat(sprintf("\nmedian E2 programme: bulk %.2f, spots %.2f, shift %.2f SD\n",
            median(bulk_e2[ok_bulk]), median(spot_e2), shift))
print(summary_tab)

## Figures: descriptive spatial mapping + explicit shift diagnostic ----------
top <- stats$sample[order(-stats$spots)][seq_len(min(8, nrow(stats)))]
map_df <- do.call(rbind, lapply(res, function(r) {
  if (!(r$sample %in% top)) return(NULL)
  data.frame(sample = r$sample, array_row = r$row, array_col = r$col,
             E2_programme = rowMeans(r$scores[, e2_progs, drop = FALSE]),
             stringsAsFactors = FALSE)
}))
map_df$sample <- factor(map_df$sample, levels = top)
save_plot(function() {
  print(ggplot(map_df, aes(array_col, -array_row, colour = E2_programme)) +
    geom_point(size = 0.35) +
    scale_colour_gradient2(low = "#2E86C1", mid = "grey90", high = "#C0392B",
                           midpoint = stats::median(map_df$E2_programme)) +
    facet_wrap(~sample, scales = "free", ncol = 4) + theme_bw(base_size = 8) +
    theme(aspect.ratio = 1, legend.position = "bottom") +
    labs(x = NULL, y = NULL, colour = "E2 programme",
         title = "Spatial map of the E2-defining programmes (EMT + angiogenesis + IL6-JAK-STAT3)"))
}, file.path(FIG, "plot_25"), 10, 7)

spot_df <- do.call(rbind, lapply(res, function(r) data.frame(
  sample = r$sample, value = rowMeans(r$scores[, e2_progs, drop = FALSE]), stringsAsFactors = FALSE)))
ord <- names(sort(tapply(spot_df$value, spot_df$sample, stats::median)))
spot_df$sample <- factor(spot_df$sample, levels = ord)
save_plot(function() {
  print(ggplot(spot_df, aes(sample, value, fill = sample)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.9) +
    theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = "E2 programme score per spot",
         title = "Within- and between-section variation of the E2 programme"))
}, file.path(FIG, "plot_26"), 9, 5)

dens <- rbind(
  data.frame(value = bulk_e2[ok_bulk], group = bulk_label[ok_bulk]),
  data.frame(value = spot_e2, group = "Visium spots")
)
save_plot(function() {
  print(ggplot(dens, aes(value, fill = group, colour = group)) +
    geom_density(alpha = 0.35, adjust = 1.2) + theme_bw() +
    scale_fill_manual(values = c(ECOTYPE_COL, `Visium spots` = "#7F7F7F")) +
    scale_colour_manual(values = c(ECOTYPE_COL, `Visium spots` = "#4D4D4D")) +
    labs(x = "E2 programme score (bulk z-scale)", y = "Density", fill = NULL, colour = NULL,
         title = sprintf("Covariate shift: spots sit %.1f SD below bulk samples", -shift)))
}, file.path(FIG, "plot_27"), 8, 5)

cat("\nsections:", nrow(stats), "| spots:", sum(stats$spots), "\n")
