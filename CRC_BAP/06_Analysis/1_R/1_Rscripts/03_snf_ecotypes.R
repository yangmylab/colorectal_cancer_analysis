# 03_snf_ecotypes.R
# Multi-view similarity fusion, consensus clustering, stability and nearest-centroid classifier.
# Worker count: configured in run_pipeline.R.

run_ecotypes <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(SNFtool)
    library(mclust)
    library(cluster)
    library(ggplot2)
    library(pheatmap)
  })

  bulk <- readRDS(file.path(cfg$paths$processed, "processed_bulk.rds"))
  sig <- readRDS(file.path(cfg$paths$processed, "signature_scores.rds"))
  scores <- sig$scores
  clin <- bulk$clinical
  samples <- intersect(colnames(scores), clin$sample_id)
  scores <- scores[, samples, drop = FALSE]
  clin <- clin[match(samples, clin$sample_id), , drop = FALSE]

  views <- list(RNA = t(scores))
  for (nm in c("methylation", "protein")) {
    path <- cfg$inputs[[nm]]
    if (!is.null(path) && length(path) > 0 && !is.na(path) && file.exists(path)) {
      m <- read_matrix_any(path)
      m <- m[, intersect(colnames(m), samples), drop = FALSE]
      if (ncol(m) >= 5 && nrow(m) >= 20) {
        m <- top_variable(m, 1000)
        m <- zscore_rows(m)
        views[[nm]] <- t(m)
      }
    }
  }

  # Single-cell deconvolution is kept as a descriptive layer, not a clustering view:
  # K=5 on RNA+celltype had ARI 0.885 but only 0.553 classifier CV, so it is not deployable.

  common <- Reduce(intersect, lapply(views, rownames))
  views <- lapply(views, function(v) v[common, , drop = FALSE])
  combined <- do.call(cbind, lapply(views, function(v) scale(v)))
  combined[!is.finite(combined)] <- 0

  make_W <- function(vlist) {
    Ws <- lapply(vlist, function(v) SNFtool::affinityMatrix(as.matrix(stats::dist(v)), K = 20, sigma = 0.5))
    if (length(Ws) == 1) Ws[[1]] else SNFtool::SNF(Ws, K = 20, t = 3)
  }

  evaluate_K <- function(W, combined, K) {
    cl <- tryCatch(SNFtool::spectralClustering(W, K, type = 3), error = function(e) NULL)
    if (is.null(cl)) return(NULL)
    sil <- tryCatch(mean(cluster::silhouette(cl, stats::dist(combined))[, 3]), error = function(e) NA_real_)
    n <- nrow(combined)
    B <- max(20, cfg$params$consensus_bootstraps)
    cons <- matrix(0, n, n)
    cnt <- 0
    for (b in seq_len(B)) {
      idx <- sample(n, max(2, floor(0.8 * n)))
      km <- tryCatch(kmeans(combined[idx, , drop = FALSE], centers = K, nstart = 10, iter.max = 100), error = function(e) NULL)
      if (is.null(km)) next
      lab <- rep(NA_integer_, n)
      lab[idx] <- km$cluster
      if (length(idx) < n) {
        d <- as.matrix(stats::dist(rbind(km$centers, combined[-idx, , drop = FALSE])))
        d <- d[seq_len(K), K + seq_len(n - length(idx)), drop = FALSE]
        lab[-idx] <- apply(d, 2, which.min)
      }
      cons <- cons + outer(lab, lab, "==")
      cnt <- cnt + 1
    }
    consensus <- if (cnt > 0) cons / cnt else matrix(NA_real_, n, n)
    dimnames(consensus) <- list(rownames(combined), rownames(combined))
    pac <- if (cnt > 0) mean(consensus > 0.1 & consensus < 0.9, na.rm = TRUE) else NA_real_
    aris <- numeric(20)
    for (b in seq_len(20)) {
      idx <- sample(n, max(2, floor(0.8 * n)))
      clb <- tryCatch(SNFtool::spectralClustering(W[idx, idx, drop = FALSE], K, type = 3), error = function(e) NULL)
      if (!is.null(clb)) aris[b] <- ari(clb, cl[idx])
    }
    list(
      stats = data.frame(K = K, silhouette = sil, min_group_fraction = min(table(cl)) / n, pac = pac, mean_ari = mean(aris, na.rm = TRUE), n_groups = length(unique(cl))),
      consensus = consensus,
      labels = cl
    )
  }

  Ks <- seq(cfg$params$k_min, cfg$params$k_max)
  W <- make_W(views)
  res <- lapply(Ks, function(k) evaluate_K(W, combined, k))
  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) == 0) stop("No K could be evaluated", call. = FALSE)
  stab <- do.call(rbind, lapply(res, `[[`, "stats"))

  ok <- stab[stab$min_group_fraction >= cfg$params$min_group_fraction & is.finite(stab$mean_ari), , drop = FALSE]
  if (nrow(ok) == 0) ok <- stab
  ok <- ok[order(-ok$mean_ari, -ok$silhouette, ok$K), , drop = FALSE]
  best <- ok[1, , drop = FALSE]
  best_res <- res[[match(best$K, stab$K)]]
  exploratory <- FALSE

  if (best$mean_ari < cfg$params$stability_ari_min) {
    log_decision("03_ecotypes", "UNSTABLE", paste("best K", best$K, "ARI", round(best$mean_ari, 3), "; trying RNA-only"))
    W_rna <- make_W(views["RNA"])
    res_rna <- lapply(Ks, function(k) evaluate_K(W_rna, combined, k))
    res_rna <- res_rna[!vapply(res_rna, is.null, logical(1))]
    stab_rna <- do.call(rbind, lapply(res_rna, `[[`, "stats"))
    ok_rna <- stab_rna[stab_rna$min_group_fraction >= cfg$params$min_group_fraction & is.finite(stab_rna$mean_ari), , drop = FALSE]
    if (nrow(ok_rna) == 0) ok_rna <- stab_rna
    ok_rna <- ok_rna[order(-ok_rna$mean_ari, -ok_rna$silhouette, ok_rna$K), , drop = FALSE]
    best_rna <- ok_rna[1, , drop = FALSE]
    if (best_rna$mean_ari > best$mean_ari) {
      W <- W_rna
      stab <- stab_rna
      res <- res_rna
      best <- best_rna
      best_res <- res[[match(best$K, stab$K)]]
    }
    if (best$mean_ari < cfg$params$stability_ari_min) {
      exploratory <- TRUE
      best_K <- 4
      best_res <- evaluate_K(W, combined, best_K)
      log_decision("03_ecotypes", "FACTOR_CONTINUUM_FALLBACK", paste("no K reached ARI", cfg$params$stability_ari_min, "; using K=4 as exploratory"))
    }
  }

  ecotype <- paste0("E", best_res$labels)
  names(ecotype) <- rownames(combined)

  X <- t(scores[, rownames(combined), drop = FALSE])
  levels_ec <- sort(unique(ecotype))
  centroids <- t(vapply(levels_ec, function(e) colMeans(X[ecotype == e, , drop = FALSE]), numeric(ncol(X))))
  rownames(centroids) <- levels_ec
  colnames(centroids) <- colnames(X)

  cv_acc <- numeric(10)
  n <- nrow(X)
  for (r in seq_len(10)) {
    idx <- sample(n, floor(0.8 * n))
    tr <- setdiff(seq_len(n), idx)
    cent_r <- t(vapply(levels_ec, function(e) {
      rows <- tr[ecotype[tr] == e]
      if (length(rows) == 0) rep(NA_real_, ncol(X)) else colMeans(X[rows, , drop = FALSE])
    }, numeric(ncol(X))))
    if (anyNA(cent_r)) next
    d <- as.matrix(stats::dist(rbind(cent_r, X[idx, , drop = FALSE])))
    d <- d[seq_along(levels_ec), length(levels_ec) + seq_along(idx), drop = FALSE]
    pred <- levels_ec[apply(d, 2, which.min)]
    cv_acc[r] <- mean(pred == ecotype[idx])
  }

  profile <- aggregate(X, by = list(ecotype = ecotype), mean)
  rownames(profile) <- profile$ecotype
  profile <- profile[, -1, drop = FALSE]
  top_sigs <- apply(profile, 1, function(x) paste(names(sort(x, decreasing = TRUE))[1:min(3, length(x))], collapse = ";"))

  out_dir <- file.path(cfg$paths$analysis, "03_Clustering")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  labels_df <- data.frame(sample_id = names(ecotype), ecotype = as.character(ecotype), stringsAsFactors = FALSE)
  data.table::fwrite(labels_df, file.path(out_dir, "ecotype_labels.tsv"), sep = "\t")
  data.table::fwrite(stab, file.path(out_dir, "stability.tsv"), sep = "\t")
  prof_df <- data.frame(ecotype = rownames(profile), top_signatures = top_sigs, profile, check.names = FALSE)
  data.table::fwrite(prof_df, file.path(out_dir, "ecotype_signature_profile.tsv"), sep = "\t")

  classifier <- list(centroids = centroids, features = colnames(X), levels = levels_ec, cv_accuracy = mean(cv_acc, na.rm = TRUE), method = "nearest_centroid", exploratory = exploratory)
  saveRDS(classifier, file.path(out_dir, "classifier.rds"))

  ann <- data.frame(ecotype = as.character(ecotype), stringsAsFactors = FALSE)
  rownames(ann) <- rownames(combined)
  ann_colors <- list(ecotype = stats::setNames(grDevices::rainbow(length(levels_ec)), levels_ec))
  save_plot(function() {
    pheatmap::pheatmap(best_res$consensus, show_rownames = FALSE, show_colnames = FALSE, annotation_col = ann, annotation_colors = ann_colors, main = paste("Consensus matrix, K =", best$K))
  }, file.path(out_dir, "03_consensus_heatmap"), 8, 7)

  pca <- stats::prcomp(combined, scale. = FALSE)
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], ecotype = ecotype)
  save_plot(function() {
    print(ggplot(pca_df, aes(PC1, PC2, colour = ecotype)) + geom_point(size = 2, alpha = 0.8) + theme_bw() + labs(title = "SNF view PCA"))
  }, file.path(out_dir, "03_ecotype_pca"), 7, 5)

  log_decision("03_ecotypes", if (exploratory) "EXPLORATORY" else "PASS", paste("K =", best$K, "; bootstrap ARI =", round(best$mean_ari, 3), "; CV accuracy =", round(mean(cv_acc, na.rm = TRUE), 3)))
  saveRDS(
    list(labels = ecotype, scores = scores, views = views, W = W, combined = combined, K = best$K, stability = stab, consensus = best_res$consensus, classifier = classifier, profile = profile, exploratory = exploratory, clinical = clin),
    file.path(cfg$paths$processed, "ecotypes.rds")
  )
  invisible(TRUE)
}
