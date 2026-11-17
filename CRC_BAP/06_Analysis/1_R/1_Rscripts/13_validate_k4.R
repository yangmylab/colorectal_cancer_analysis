#!/usr/bin/env Rscript
# 13_validate_k4.R
# Train the K=4 nearest-centroid classifier in TCGA and test it in GSE39582/GSE17536.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(GSVA)
  library(SNFtool)
  library(yaml)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")
cfg <- yaml::read_yaml("config/pipeline_config.yaml")
ec <- readRDS(file.path(cfg$paths$processed, "ecotypes.rds"))
surv <- readRDS(file.path(cfg$paths$processed, "survival_results.rds"))
sig <- readRDS(file.path(cfg$paths$processed, "signature_scores.rds"))

# K=4 labels in TCGA
lab4 <- SNFtool::spectralClustering(ec$W, 4, type = 3)
names(lab4) <- rownames(ec$combined)
lab4 <- stats::setNames(paste0("E", lab4), rownames(ec$combined))
X <- t(ec$scores[, names(lab4), drop = FALSE])
levels4 <- sort(unique(lab4))
cent4 <- t(vapply(levels4, function(e) colMeans(X[lab4 == e, , drop = FALSE]), numeric(ncol(X))))
rownames(cent4) <- levels4
colnames(cent4) <- colnames(X)

# TCGA K=4 Cox (complete cases used for the primary K=2 model)
dat <- surv$dat
dat$ecotype <- factor(lab4[as.character(dat$sample_id)], levels = levels4)
dat <- dat[!is.na(dat$ecotype), , drop = FALSE]
covs <- c("age", "sex", "stage_num", "msi_status", "sidedness")
covs <- covs[covs %in% names(dat)]
covs <- covs[vapply(dat[covs], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
for (cv in covs) if (is.character(dat[[cv]])) dat[[cv]] <- factor(dat[[cv]])
fb <- as.formula(paste("Surv(time, event) ~", if (length(covs)) paste(covs, collapse = " + ") else "1"))
ff <- as.formula(paste("Surv(time, event) ~ ecotype", if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
b <- coxph(fb, data = dat); f <- coxph(ff, data = dat)
a <- tryCatch(anova(b, f, test = "LRT"), error = function(e) NULL)
p_tcga <- if (!is.null(a)) as.data.frame(a)[2, "Pr(>|Chi|)"] else NA_real_

rows <- list(data.frame(cohort = "TCGA-COAD/READ", n = nrow(dat), events = sum(dat$event), c_index = as.numeric(concordance(f)$concordance), lrt_p = p_tcga, reference = "E1", stringsAsFactors = FALSE))

# Validation cohorts
val_dir <- cfg$inputs$validation_dir
for (cd in list.dirs(val_dir, recursive = FALSE, full.names = TRUE)) {
  ep <- file.path(cd, "expression.tsv"); cp <- file.path(cd, "clinical.tsv")
  if (!file.exists(ep) || !file.exists(cp)) next
  ve <- read_matrix_any(ep); vc <- fread(cp, data.table = FALSE)
  genes <- intersect(rownames(ve), unique(unlist(sig$gene_sets)))
  vg <- lapply(sig$gene_sets, function(g) intersect(g, genes)); vg <- vg[lengths(vg) >= 5]
  if (length(vg) < 2) next
  vs <- tryCatch(GSVA::gsva(GSVA::gsvaParam(ve, vg, kcdf = "Gaussian", minSize = 5, maxSize = 500, absRanking = TRUE), verbose = FALSE), error = function(e) NULL)
  if (is.null(vs)) next
  feats <- intersect(rownames(vs), colnames(cent4)); if (length(feats) < 2) next
  cent <- cent4[, feats, drop = FALSE]
  Xv <- t(vs[feats, , drop = FALSE])
  d <- as.matrix(dist(rbind(cent, Xv))); d <- d[seq_len(nrow(cent)), nrow(cent) + seq_len(nrow(Xv)), drop = FALSE]
  pred <- rownames(cent)[apply(d, 2, which.min)]; names(pred) <- colnames(vs)
  vc$ecotype <- pred[as.character(vc$sample_id)]
  if ("stage" %in% names(vc) && !"stage_num" %in% names(vc)) vc$stage_num <- as.numeric(factor(vc$stage, levels = c("I", "II", "III", "IV")))
  vc$time <- vc$os_time_days; vc$event <- vc$os_event
  vc <- vc[!is.na(vc$time) & !is.na(vc$event) & vc$time > 0 & !is.na(vc$ecotype), , drop = FALSE]
  if (nrow(vc) < 20 || sum(vc$event) < 5) next
  vc$ecotype <- factor(vc$ecotype, levels = levels4)
  vc_covs <- covs[covs %in% names(vc)]
  vc_covs <- vc_covs[vapply(vc[vc_covs], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
  for (cv in vc_covs) if (is.character(vc[[cv]])) vc[[cv]] <- factor(vc[[cv]])
  vfb <- as.formula(paste("Surv(time, event) ~", if (length(vc_covs)) paste(vc_covs, collapse = " + ") else "1"))
  vff <- as.formula(paste("Surv(time, event) ~ ecotype", if (length(vc_covs)) paste("+", paste(vc_covs, collapse = " + ")) else ""))
  vb <- tryCatch(coxph(vfb, data = vc), error = function(e) NULL)
  vf <- tryCatch(coxph(vff, data = vc), error = function(e) NULL)
  if (is.null(vf)) next
  va <- tryCatch(anova(vb, vf, test = "LRT"), error = function(e) NULL)
  vp <- if (!is.null(va)) as.data.frame(va)[2, "Pr(>|Chi|)"] else NA_real_
  rows[[length(rows) + 1]] <- data.frame(cohort = basename(cd), n = nrow(vc), events = sum(vc$event), c_index = as.numeric(concordance(vf)$concordance), lrt_p = vp, reference = "E1", stringsAsFactors = FALSE)
}

out <- do.call(rbind, rows)
dir.create("07_Validation", recursive = TRUE, showWarnings = FALSE)
fwrite(out, "07_Validation/K4_validation.tsv", sep = "\t")
cat("K=4 validation:\n"); print(out)
