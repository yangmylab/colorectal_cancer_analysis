# 02_signatures_deconvolution.R
# Score gene-set signatures (GSVA ssGSEA) and optionally deconvolve cell types (NNLS).
# Worker count: configured in run_pipeline.R.

run_signatures <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(GSVA)
    library(pheatmap)
    library(ggplot2)
  })

  obj <- readRDS(file.path(cfg$paths$processed, "processed_bulk.rds"))
  log_expr <- obj$log_expr
  clin <- obj$clinical

  read_gmt <- function(path) {
    lines <- readLines(path)
    sets <- list()
    for (l in lines) {
      f <- strsplit(l, "\t")[[1]]
      if (length(f) >= 3) sets[[f[1]]] <- unique(f[-c(1, 2)])
    }
    sets
  }

  id_type <- if (all(grepl("^ENSG", head(rownames(log_expr), 100)))) "ensembl" else "symbol"
  sig_path <- cfg$inputs$signatures
  if (!is.null(sig_path) && length(sig_path) > 0 && !is.na(sig_path) && file.exists(sig_path)) {
    gene_sets <- read_gmt(sig_path)
    first_gene <- unique(unlist(gene_sets))[1]
    custom_is_ensembl <- grepl("^ENSG", first_gene)
    if (id_type == "ensembl" && !custom_is_ensembl) {
      msig <- msigdbr::msigdbr(db_species = "HS", collection = "H")
      map <- unique(msig[, c("gene_symbol", "ensembl_gene")])
      map <- map[!is.na(map$ensembl_gene) & map$ensembl_gene != "", , drop = FALSE]
      mapv <- stats::setNames(map$ensembl_gene, map$gene_symbol)
      gene_sets <- lapply(gene_sets, function(g) unique(unname(mapv[g[g %in% names(mapv)]])))
    }
    sig_source <- sig_path
  } else {
    msig <- msigdbr::msigdbr(db_species = "HS", collection = "H")
    key <- if (id_type == "ensembl" && "ensembl_gene" %in% colnames(msig)) "ensembl_gene" else "gene_symbol"
    gene_sets <- split(msig[[key]], msig$gs_name)
    sig_source <- paste0("MSigDB Hallmark (msigdbr; ", key, ")")
  }

  gene_sets <- lapply(gene_sets, function(g) intersect(unique(g), rownames(log_expr)))
  sizes <- lengths(gene_sets)
  gene_sets <- gene_sets[sizes >= 5 & sizes <= 500]
  if (length(gene_sets) < 2) stop("Fewer than 2 usable gene sets after filtering", call. = FALSE)

  scores <- tryCatch(
    {
      param <- GSVA::gsvaParam(log_expr, gene_sets, kcdf = "Gaussian", minSize = 5, maxSize = 500, absRanking = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    },
    error = function(e) {
      log_decision("02_signatures", "FALLBACK", paste("gsvaParam failed:", conditionMessage(e), "; trying ssgseaParam"))
      param <- GSVA::ssgseaParam(log_expr, gene_sets, normalize = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    }
  )

  sig_dir <- file.path(cfg$paths$analysis, "05_Immune_Analysis")
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.frame(signature = rownames(scores), scores, check.names = FALSE), file.path(sig_dir, "signature_scores.tsv"), sep = "\t")

  deconv <- NULL
  ct_path <- cfg$inputs$celltype_signatures
  if (!is.null(ct_path) && length(ct_path) > 0 && !is.na(ct_path) && file.exists(ct_path)) {
    ct <- read_matrix_any(ct_path)
    genes <- intersect(rownames(ct), rownames(log_expr))
    ct <- ct[genes, , drop = FALSE]
    deconv <- matrix(NA_real_, nrow = ncol(ct), ncol = ncol(log_expr), dimnames = list(colnames(ct), colnames(log_expr)))
    for (j in seq_len(ncol(log_expr))) {
      b <- log_expr[genes, j]
      fit <- nnls::nnls(A = as.matrix(ct), b = as.numeric(b))
      x <- fit$x
      if (sum(x) > 0) x <- x / sum(x)
      deconv[, j] <- x
    }
    data.table::fwrite(data.frame(cell_type = rownames(deconv), deconv, check.names = FALSE), file.path(sig_dir, "deconvolution_nnls.tsv"), sep = "\t")
    deconv_note <- "NNLS fractions from user-supplied cell-type signatures"
  } else {
    deconv_note <- "No cell-type signature matrix supplied; signature scores are marker scores, not absolute fractions"
  }

  top_sigs <- names(sort(apply(scores, 1, stats::var), decreasing = TRUE))[seq_len(min(30, nrow(scores)))]
  ann_col <- data.frame(group = clin$group, row.names = clin$sample_id)
  save_plot(function() {
    pheatmap::pheatmap(scores[top_sigs, , drop = FALSE], show_colnames = FALSE, annotation_col = ann_col, main = "Top variable signatures")
  }, file.path(sig_dir, "02_signature_heatmap"), 9, 7)

  log_decision("02_signatures", "PASS", paste("source:", sig_source, "; sets:", length(gene_sets), "; deconvolution:", deconv_note))
  saveRDS(list(scores = scores, gene_sets = gene_sets, deconv = deconv, source = sig_source, note = deconv_note), file.path(cfg$paths$processed, "signature_scores.rds"))
  invisible(TRUE)
}
