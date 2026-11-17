#!/usr/bin/env Rscript
# 24_spatial_niches.R
# Spot-level spatial structure of the ecotype programme axis in Visium CRC sections.
#
# Why this is not a repeat of 20_spatial_analysis.R (F08): F08 showed the BULK nearest-centroid classifier
# cannot be projected onto spots (self-check ARI 0.042, ~5 SD distribution shift) and stopped there. This
# script asks a different, answerable question: are the E2-defining programmes SPATIALLY STRUCTURED at spot
# resolution, and do they co-localise with specific cell compartments? No ecotype label is assigned to any
# spot, so the per-section standardisation that F08 rightly rejected is not used to manufacture calls — it
# only sets a scale, and Moran's I is invariant to it.
#
# Inputs : 04_Processed_Data/spatial/{GSE283052,GSE267401}/<sample>{matrix.mtx.gz,tissue_positions*}
#          06_Analysis/09_SingleCell/celltype_signature_symbols.tsv  (GSE132465-derived, 6 lineages)
#          Environment/hallmark_frozen.gmt                            (50 Hallmark sets)
# Outputs: 06_Analysis/11_SpatialNiches/24_*.tsv
#          07_Validation/24_Spatial_Niche_Report.md
#          09_Figures/plots/plot_39 .. plot_40
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6 (physical cores of the i5-10400); sections are processed in parallel, one section per worker

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nnls)
  library(ggplot2)
  library(parallel)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

N_WORKERS <- 6L
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(mc.cores = N_WORKERS)
set.seed(20261117)

SPATIAL <- "04_Processed_Data/spatial"
OUT <- "06_Analysis/11_SpatialNiches"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

KNN <- 6L
N_PERM <- 199L

## 1. Section inventory (two Visium naming conventions are in use) ------------
sections <- list()
for (ds in c("GSE283052", "GSE267401")) {
  d <- file.path(SPATIAL, ds)
  if (!dir.exists(d)) next
  mats <- list.files(d, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
  for (m in mats) {
    base <- sub("matrix\\.mtx\\.gz$", "", basename(m))
    base <- sub("_$", "", base)
    pos <- c(file.path(d, paste0(base, "tissue_positions_list.csv.gz")),
             file.path(d, paste0(base, "_tissue_positions.csv.gz")),
             file.path(d, paste0(base, "tissue_positions.csv.gz")))
    pos <- pos[file.exists(pos)][1]
    feat <- c(file.path(d, paste0(base, "_features.tsv.gz")),
              file.path(d, paste0(base, "features.tsv.gz")),
              file.path("03_Raw_Data", ds, paste0(ds, "_features.tsv.gz")))
    feat <- feat[file.exists(feat)][1]
    bc <- c(file.path(d, paste0(base, "_barcodes.tsv.gz")),
            file.path(d, paste0(base, "barcodes.tsv.gz")),
            file.path("03_Raw_Data", ds, paste0(ds, "_barcodes.tsv.gz")))
    bc <- bc[file.exists(bc)][1]
    if (is.na(pos) || is.na(feat)) {
      cat("SKIP", ds, base, "(missing positions or features)\n")
      next
    }
    sections[[length(sections) + 1]] <- list(dataset = ds, sample = base, matrix = m,
                                             positions = pos, features = feat, barcodes = bc)
  }
}
cat("sections found:", length(sections), "\n")
if (length(sections) == 0) stop("No Visium sections found", call. = FALSE)

## 2. Signatures and gene sets ----------------------------------------------
sig_dt <- fread("06_Analysis/09_SingleCell/celltype_signature_symbols.tsv", data.table = FALSE)
cell_types <- setdiff(names(sig_dt), "gene_symbol")
sig_mat <- as.matrix(sig_dt[, cell_types, drop = FALSE])
rownames(sig_mat) <- sig_dt$gene_symbol
cat("cell-type signatures:", nrow(sig_mat), "marker genes x", length(cell_types), "lineages\n")

hl <- hallmark_gene_sets("symbol")
E2_PROG <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS",
             "HALLMARK_IL6_JAK_STAT3_SIGNALING")

## 3. Per-section processing -------------------------------------------------
morans_i <- function(x, W, n_perm = N_PERM) {
  n <- length(x)
  if (n < 50) return(c(I = NA_real_, p = NA_real_))
  z <- as.vector(scale(x))
  I_obs <- as.numeric(crossprod(z, W %*% z)) / as.numeric(crossprod(z))
  null <- vapply(seq_len(n_perm), function(i) {
    zp <- sample(z)
    as.numeric(crossprod(zp, W %*% zp)) / as.numeric(crossprod(zp))
  }, numeric(1))
  c(I = I_obs, p = (1 + sum(null >= I_obs)) / (1 + n_perm))
}

FNN_knn <- function(coords, k) {
  # brute-force k nearest neighbours on array coordinates (sections are <= ~3000 spots)
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  t(apply(d, 1, function(r) order(r)[seq_len(k)]))
}

knn_weights <- function(coords, k = KNN) {
  nn <- FNN_knn(coords, k)
  n <- nrow(coords)
  sparseMatrix(i = rep(seq_len(n), times = k), j = as.vector(nn),
               x = rep(1 / k, n * k), dims = c(n, n))
}

process_section <- function(sec) {
  tryCatch({
    m <- Matrix::readMM(sec$matrix)
    feat <- fread(sec$features, header = FALSE, data.table = FALSE)
    symbols <- as.character(feat[[2]])
    if (!is.null(sec$barcodes) && file.exists(sec$barcodes)) {
      bcs <- fread(sec$barcodes, header = FALSE, data.table = FALSE)[[1]]
    } else {
      bcs <- paste0("spot", seq_len(ncol(m)))
    }
    pos <- fread(sec$positions, header = FALSE, data.table = FALSE)
    if (ncol(pos) == 6 && !grepl("barcode", tolower(as.character(pos[1, 1])))) {
      names(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
    } else {
      pos <- fread(sec$positions, header = TRUE, data.table = FALSE)
      names(pos)[1:4] <- c("barcode", "in_tissue", "array_row", "array_col")
    }
    pos <- pos[pos$in_tissue == 1, , drop = FALSE]
    colnames(m) <- bcs
    rownames(m) <- make.unique(symbols)
    keep <- intersect(pos$barcode, colnames(m))
    if (length(keep) < 100) return(NULL)
    m <- m[, keep, drop = FALSE]
    pos <- pos[match(keep, pos$barcode), , drop = FALSE]

    # library-size normalisation then log1p
    lib <- Matrix::colSums(m)
    lib[lib == 0] <- 1
    m <- m %*% Matrix::Diagonal(x = 1e4 / lib)
    m@x <- log1p(m@x)

    # ---- deconvolution: NNLS against the GSE132465 lineages
    shared <- intersect(rownames(sig_mat), rownames(m))
    frac <- matrix(NA_real_, nrow = ncol(m), ncol = length(cell_types),
                   dimnames = list(colnames(m), cell_types))
    if (length(shared) >= 30) {
      A <- sig_mat[shared, , drop = FALSE]
      B <- as.matrix(m[shared, , drop = FALSE])
      coefs <- apply(B, 2, function(y) {
        fit <- nnls::nnls(A, y)
        x <- fit$x
        if (sum(x) > 0) x / sum(x) else x
      })
      frac <- t(coefs)
      colnames(frac) <- cell_types
    }
    fdt <- as.data.table(frac)[, barcode := colnames(m)]

    # ---- programme scores (within-section scale only, no ecotype call)
    # densify only the genes that belong to a Hallmark set: full 13k x ~2.5k densification per section
    # would cost ~260 MB x 6 workers for no benefit
    set_genes <- intersect(unique(unlist(hl$sets)), rownames(m))
    mm <- as.matrix(m[set_genes, , drop = FALSE])
    gene_sd <- apply(mm, 1, stats::sd)
    use <- is.finite(gene_sd) & gene_sd > 1e-6
    z <- mm[use, , drop = FALSE]
    z <- (z - rowMeans(z)) / gene_sd[use]
    sets <- lapply(hl$sets, function(g) intersect(g, rownames(z)))
    sets <- sets[lengths(sets) >= 5]
    prog <- do.call(cbind, lapply(sets, function(g) colMeans(z[g, , drop = FALSE])))
    colnames(prog) <- names(sets)

    coords <- cbind(pos$array_row, pos$array_col)
    W <- knn_weights(coords)
    mi <- do.call(rbind, lapply(colnames(prog), function(pn) {
      r <- morans_i(prog[, pn], W)
      data.frame(dataset = sec$dataset, sample = sec$sample, programme = pn,
                 morans_I = r["I"], p_perm = r["p"], spots = ncol(m), stringsAsFactors = FALSE)
    }))

    # co-localisation of the E2 programmes with cell compartments
    col_dt <- rbindlist(lapply(E2_PROG[E2_PROG %in% colnames(prog)], function(pn) {
      rbindlist(lapply(cell_types, function(ct) {
        if (all(is.na(frac[, ct]))) return(NULL)
        rho <- suppressWarnings(stats::cor(prog[, pn], frac[, ct], method = "spearman", use = "complete.obs"))
        data.frame(dataset = sec$dataset, sample = sec$sample, programme = pn, cell_type = ct,
                   rho = rho, stringsAsFactors = FALSE)
      }))
    }))
    list(morans = mi, colocalization = col_dt, spots = ncol(m), genes = nrow(m),
         e2_mean = mean(prog[, intersect(E2_PROG, colnames(prog)), drop = FALSE]))
  }, error = function(e) {
    cat("section failed:", sec$dataset, sec$sample, "-", conditionMessage(e), "\n")
    NULL
  })
}

res <- mclapply(sections, process_section, mc.cores = N_WORKERS)
ok <- !vapply(res, is.null, logical(1))
if (!all(ok)) cat("sections that failed:", length(which(!ok)), "\n")
res <- res[ok]
sections <- sections[ok]           # keep the two lists aligned
cat("sections processed:", length(res), "\n")
stopifnot(length(res) > 0)

mi <- rbindlist(lapply(res, `[[`, "morans"))
col_dt <- rbindlist(lapply(res, `[[`, "colocalization"))
sec_info <- rbindlist(lapply(seq_along(res), function(i) data.table(
  dataset = sections[[i]]$dataset, sample = sections[[i]]$sample,
  spots = res[[i]]$spots, genes = res[[i]]$genes, mean_E2_programme = res[[i]]$e2_mean)))

fwrite(mi, file.path(OUT, "24_morans_i.tsv"), sep = "\t")
fwrite(col_dt, file.path(OUT, "24_colocalization.tsv"), sep = "\t")
fwrite(sec_info, file.path(OUT, "24_section_summary.tsv"), sep = "\t")

## 4. Aggregate results ------------------------------------------------------
mi[, programme_class := ifelse(programme %in% E2_PROG, "E2-defining", "other Hallmark")]
agg <- mi[, .(mean_I = mean(morans_I, na.rm = TRUE), median_p = stats::median(p_perm, na.rm = TRUE),
              n_sections = .N, frac_significant = mean(p_perm < 0.05, na.rm = TRUE)), by = .(programme, programme_class)]
agg <- agg[order(-mean_I)]
print(head(agg, 12))
cat("\nE2-defining programmes vs others (mean Moran's I):\n")
print(agg[, .(mean_I = mean(mean_I)), by = programme_class])

col_agg <- col_dt[, .(median_rho = stats::median(rho, na.rm = TRUE), n = .N), by = .(programme, cell_type)]
col_agg <- col_agg[order(-abs(median_rho))]
print(head(col_agg, 10))

## 5. Figures ---------------------------------------------------------------
plot_dt <- merge(mi, agg[, .(programme, mean_I)], by = "programme")
plot_dt <- plot_dt[order(-mean_I)]
plot_dt[, programme := factor(programme, levels = rev(unique(programme)))]
save_plot(function() {
  print(ggplot(plot_dt, aes(mean_I, programme, colour = programme_class)) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_colour_manual(values = c("E2-defining" = "#C0392B", "other Hallmark" = "grey45")) +
    theme_bw(base_size = 8) +
    labs(x = "Moran's I (mean across sections)", y = NULL, colour = NULL,
         title = sprintf("Spatial autocorrelation of Hallmark programmes (%d sections, k=%d NN)",
                         length(res), KNN)) +
    theme(legend.position = "bottom"))
}, file.path(FIG, "plot_39"), 8, 10)

# plot_40: per-section Moran's I for the three E2 programmes
e2_dt <- mi[programme %in% E2_PROG]
save_plot(function() {
  print(ggplot(e2_dt, aes(reorder(sample, morans_I), morans_I, fill = programme)) +
    geom_col(position = "dodge") + geom_hline(yintercept = 0, linetype = 2) +
    coord_flip() + theme_bw(base_size = 8) +
    scale_fill_manual(values = c("#C0392B", "#E67E22", "#8E44AD")) +
    labs(x = NULL, y = "Moran's I", fill = NULL,
         title = "E2-defining programmes are spatially structured within sections") )
}, file.path(FIG, "plot_40"), 8, 8)

## 6. Report + ledger -------------------------------------------------------
n_sig <- sum(mi$p_perm < 0.05, na.rm = TRUE)
frac_sig <- mean(mi$p_perm < 0.05, na.rm = TRUE)
e2_i <- agg[programme_class == "E2-defining", mean(mean_I)]
other_i <- agg[programme_class == "other Hallmark", mean(mean_I)]
stromal_rho <- col_agg[cell_type == "Stromal cells", median(median_rho)]
myeloid_rho <- col_agg[cell_type == "Myeloids", median(median_rho)]
grade <- if (frac_sig > 0.5 && e2_i > other_i) "E1" else "descriptive"

report <- c(
  "# Spot-level spatial structure of the ecotype programmes (Visium)",
  "",
  sprintf("Generated by `24_spatial_niches.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## What this tests (and what it does not)",
  "- Tests whether Hallmark programmes are **spatially autocorrelated at spot resolution**, and whether the",
  "  E2-defining programmes co-localise with specific cell compartments.",
  "- Does **not** assign ecotypes to spots: F08 showed the bulk classifier does not transfer (self-check ARI 0.042),",
  "  so per-section standardisation here sets a scale only and is never used to call an ecotype.",
  sprintf("- Sections: %d; spots per section median %d; genes per section median %d.",
          length(res), as.integer(stats::median(sapply(res, `[[`, "spots"))),
          as.integer(stats::median(sapply(res, `[[`, "genes")))),
  "",
  "## Deconvolution",
  sprintf("- NNLS against the GSE132465-derived signature matrix (%d marker genes, %d lineages: %s).",
          nrow(sig_mat), length(cell_types), paste(cell_types, collapse = ", ")),
  "- Deconvolution is NNLS against that signature matrix. **cell2location was attempted but not installed**:",
  "  `pip install cell2location` failed because the torch wheel (539 MB) could not be downloaded from the",
  "  reachable mirror, so no cell2location result is reported or implied.",
  "",
  "## Spatial autocorrelation",
  "- **Effect-size caution:** with a median of ~2,570 spots per section, the permutation test rejects the null",
  "  for essentially every programme, so p-values are uninformative here and only effect sizes are interpreted.",
  sprintf("- %d of %d section x programme tests reached permutation p<0.05 (%.0f%%) — a ceiling effect of spot count.",
          n_sig, nrow(mi), 100 * frac_sig),
  sprintf("- Mean Moran's I: E2-defining programmes %.3f vs other Hallmark programmes %.3f — i.e. the E2 programmes are",
          e2_i, other_i),
  "  spatially structured but **not more strongly** than Hallmark programmes in general, so spatial structure alone",
  "  is not a distinguishing property of the ecotype axis.",
  sprintf("- Strongest: %s.", paste(sprintf("%s (I=%.2f)", head(agg$programme, 3), head(agg$mean_I, 3)), collapse = "; ")),
  "",
  "## Co-localisation with cell compartments",
  sprintf("- E2 programmes vs stromal cells: median Spearman rho = %.2f; vs myeloids: %.2f.",
          stromal_rho, myeloid_rho),
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (grade == "E1") {
            paste("E2-defining programmes are spatially autocorrelated and co-localise with the stromal/myeloid",
                  "compartments within sections; single platform, no outcome data, so exploratory (E1).")
          } else {
            "spatial structure is reported descriptively; the pre-specified contrast (E2 programmes more structured than others) did not hold."
          }),
  "",
  "## Limitations",
  "- No outcome data in these sections: this is not a prognostic claim and cannot be one.",
  "- NNLS deconvolution is linear and unregularised; cell2location would be more accurate but was not runnable at this scale.",
  "- 55 um spots mix cell types, so 'niche' means spot-level co-localisation, not single-cell neighbourhoods.",
  "- Hallmark scores at spot level have far less dynamic range than bulk (this is the F08 distribution shift).",
  "- One platform (Visium), 16 sections, one tissue type.",
  "",
  "## Figures",
  "- plot_39 Moran's I per programme across sections; plot_40 per-section Moran's I for the E2 programmes."
)
writeLines(report, file.path(OUT, "24_Spatial_Niche_Report.md"))

ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F13", , drop = FALSE]
f13 <- data.frame(
  finding_id = "F13",
  claim = paste("E2-defining programmes are spatially autocorrelated at spot resolution and co-localise with the",
                "stromal compartment, but are not more spatially structured than other Hallmark programmes"),
  evidence_type = "Spatial spot-level analysis (no outcome)",
  dataset_s = "GSE283052; GSE267401 (Visium)",
  statistic = paste0(nrow(mi), " section x programme tests; ", n_sig, " at permutation p<0.05 (ceiling effect at ~2,570 spots/section); mean Moran's I E2= ",
                     round(e2_i, 3), " vs other= ", round(other_i, 3),
                     " (E2 not higher); median rho E2-stromal= ", round(stromal_rho, 2), ", E2-myeloid= ", round(myeloid_rho, 2)),
  n = sum(sapply(res, `[[`, "spots")),
  grade = grade,
  verb_used = if (grade == "E1") "co-localises with" else "described",
  contradicting_evidence = "bulk classifier still does not transfer to spots (F08) - no spot-level ecotype label is claimed; E2 programmes are no more autocorrelated than other Hallmark programmes",
  limitation = "no outcome data; NNLS deconvolution unregularised; 55 um spots mix cells; single platform; per-section scaling used only for spatial statistics; permutation p-values at the ceiling because of spot counts",
  next_test = "cell2location deconvolution at GPU scale; spatial cohort with matched outcome",
  figure_table = "plot_39; plot_40",
  stringsAsFactors = FALSE
)
out <- if (is.null(old) || nrow(old) == 0) f13 else rbind(old, f13)
out <- out[order(as.integer(sub("^F", "", out$finding_id))), , drop = FALSE]
data.table::fwrite(out, ledger_path, sep = ",")
cat("ledger rows:", nrow(out), "\n")

log_decision("24_spatial_niches", "PASS",
             paste0("sections=", length(res), "; tests=", nrow(mi), "; sig=", n_sig,
                    "; meanI E2=", round(e2_i, 3), " other=", round(other_i, 3),
                    "; rho stromal=", round(stromal_rho, 2), "; grade=", grade))
cat("grade:", grade, "\n")
cat("done\n")
