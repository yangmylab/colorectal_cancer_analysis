# 04_survival_validation.R
# Cox modelling, calibration, decision-curve analysis and frozen-classifier external validation.

run_survival_validation <- function(cfg, synthetic = FALSE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(glmnet)
    library(ggplot2)
  })

  ec <- readRDS(file.path(cfg$paths$processed, "ecotypes.rds"))
  sig <- readRDS(file.path(cfg$paths$processed, "signature_scores.rds"))
  cl <- ec$clinical
  lab <- ec$labels
  cl$ecotype <- as.character(lab[cl$sample_id])
  cl <- cl[!is.na(cl$ecotype), , drop = FALSE]

  endpoint <- "css"
  if (!all(c("css_time_days", "css_event") %in% names(cl)) || sum(cl$css_event, na.rm = TRUE) < cfg$params$min_events) endpoint <- "os"
  if (!all(c(paste0(endpoint, "_time_days"), paste0(endpoint, "_event")) %in% names(cl))) stop("No usable survival endpoint (css/os)", call. = FALSE)
  time_col <- paste0(endpoint, "_time_days")
  event_col <- paste0(endpoint, "_event")
  dat <- cl
  dat$time <- as.numeric(dat[[time_col]])
  dat$event <- as.integer(dat[[event_col]])
  dat <- dat[!is.na(dat$time) & !is.na(dat$event) & dat$time > 0, , drop = FALSE]
  dat$ecotype <- factor(dat$ecotype)

  if ("stage" %in% names(dat)) {
    if (is.character(dat$stage) || is.factor(dat$stage)) {
      dat$stage_num <- as.numeric(factor(dat$stage, levels = c("I", "II", "III", "IV")))
    } else {
      dat$stage_num <- as.numeric(dat$stage)
    }
  }
  covs <- c("age", "sex", "stage_num", "msi_status", "sidedness")
  covs <- covs[covs %in% names(dat)]
  covs <- covs[vapply(dat[covs], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
  for (cv in covs) if (is.character(dat[[cv]])) dat[[cv]] <- factor(dat[[cv]])

  model_vars <- c("time", "event", "ecotype", covs)
  n_before <- nrow(dat)
  dat <- dat[stats::complete.cases(dat[, model_vars, drop = FALSE]), , drop = FALSE]
  log_decision("04_survival", "COMPLETE_CASE", paste("excluded", n_before - nrow(dat), "rows with missing model covariates"))
  if (nrow(dat) < 20) stop("Fewer than 20 complete cases for survival modelling", call. = FALSE)

  f_base <- as.formula(paste("Surv(time, event) ~", if (length(covs) > 0) paste(covs, collapse = " + ") else "1"))
  f_full <- as.formula(paste("Surv(time, event) ~ ecotype", if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""))
  base <- survival::coxph(f_base, data = dat, x = TRUE)
  full <- survival::coxph(f_full, data = dat, x = TRUE)
  lrt <- tryCatch(stats::anova(base, full, test = "LRT"), error = function(e) NULL)
  lrt_df <- if (!is.null(lrt)) as.data.frame(lrt) else NULL
  lrt_p <- if (!is.null(lrt_df) && nrow(lrt_df) >= 2 && "Pr(>|Chi|)" %in% names(lrt_df)) lrt_df[2, "Pr(>|Chi|)"] else NA_real_
  c_full <- as.numeric(survival::concordance(full)$concordance)
  c_base <- as.numeric(survival::concordance(base)$concordance)
  ph <- tryCatch(survival::cox.zph(full), error = function(e) NULL)
  ph_p <- if (!is.null(ph)) min(ph$table[, "p"], na.rm = TRUE) else NA_real_

  coef_tab <- data.frame(term = names(coef(full)), HR = exp(coef(full)), lower = exp(confint(full))[, 1], upper = exp(confint(full))[, 2], p = summary(full)$coefficients[, "Pr(>|z|)"], stringsAsFactors = FALSE)

  out_dir <- file.path(cfg$paths$analysis, "06_Other_Analysis")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(coef_tab, file.path(out_dir, "cox_model.tsv"), sep = "\t")
  perf <- data.frame(endpoint = endpoint, n = nrow(dat), events = sum(dat$event), c_index_base = c_base, c_index_full = c_full, lrt_p = lrt_p, ph_min_p = ph_p, covariates = paste(covs, collapse = ","))
  data.table::fwrite(perf, file.path(out_dir, "survival_performance.tsv"), sep = "\t")

  # Calibration at 5 years (approximate; censoring handled by observed-status subset)
  cal_tab <- NULL
  pred5 <- NULL
  obs5 <- NULL
  sf <- tryCatch(survival::survfit(full, newdata = dat), error = function(e) NULL)
  if (!is.null(sf)) {
    tm <- 1825
    sv <- sf$surv
    if (is.matrix(sv) && length(sf$time) > 0) {
      idx <- max(which(sf$time <= tm))
      if (length(idx) == 0 || !is.finite(idx)) {
        pred5 <- rep(0, nrow(dat))
      } else {
        pred5 <- as.numeric(1 - sv[idx, ])
      }
    } else if (length(sv) > 0) {
      pred5 <- rep(as.numeric(1 - sv[length(sv)]), nrow(dat))
    } else {
      pred5 <- rep(NA_real_, nrow(dat))
    }
    obs5 <- ifelse(dat$time <= tm & dat$event == 1, 1, ifelse(dat$time >= tm, 0, NA))
    cal <- data.frame(pred = pred5, obs = obs5)
    cal <- cal[is.finite(cal$pred) & !is.na(cal$obs), , drop = FALSE]
    if (nrow(cal) >= 20) {
      cal$decile <- cut(cal$pred, breaks = 10, include.lowest = TRUE)
      cal_tab <- aggregate(cbind(pred, obs) ~ decile, data = cal, mean)
      data.table::fwrite(cal_tab, file.path(out_dir, "calibration_5y.tsv"), sep = "\t")
      save_plot(function() {
        print(ggplot(cal_tab, aes(pred, obs)) + geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(size = 2) + geom_line() + theme_bw() + labs(x = "Predicted 5-year risk", y = "Observed 5-year risk", title = "Calibration (5-year)"))
      }, file.path(out_dir, "04_calibration_5y"), 6, 5)
    }
  }

  # Decision-curve analysis (approximate, complete 5-year status only)
  if (!is.null(pred5) && !is.null(obs5)) {
    ok <- !is.na(obs5)
    if (sum(ok) >= 20) {
      thresholds <- seq(0.10, 0.50, 0.05)
      nb <- vapply(thresholds, function(pt) {
        tp <- mean(pred5[ok] >= pt & obs5[ok] == 1)
        fp <- mean(pred5[ok] >= pt & obs5[ok] == 0)
        tp - fp * (pt / (1 - pt))
      }, numeric(1))
      dca <- data.frame(threshold = thresholds, net_benefit = nb)
      data.table::fwrite(dca, file.path(out_dir, "dca_5y.tsv"), sep = "\t")
      save_plot(function() {
        print(ggplot(dca, aes(threshold, net_benefit)) + geom_line() + geom_point() + theme_bw() + labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve (5-year, approximate)"))
      }, file.path(out_dir, "04_dca_5y"), 6, 5)
    }
  }

  # Kaplan-Meier by ecotype
  km <- survival::survfit(Surv(time, event) ~ ecotype, data = dat)
  save_plot(function() {
    plot(km, col = seq_along(levels(dat$ecotype)), lwd = 2, xlab = "Days", ylab = "Survival probability", main = paste("Kaplan-Meier by ecotype (", endpoint, ")", sep = ""))
    legend("topright", legend = levels(dat$ecotype), col = seq_along(levels(dat$ecotype)), lwd = 2, bty = "n")
  }, file.path(out_dir, "04_km_ecotype"), 7, 5)

  # External validation with frozen classifier
  val_rows <- list()
  val_dir <- cfg$inputs$validation_dir
  if (!is.null(val_dir) && length(val_dir) > 0 && !is.na(val_dir) && dir.exists(val_dir)) {
    cohorts <- list.dirs(val_dir, recursive = FALSE, full.names = TRUE)
    for (cd in cohorts) {
      ep <- file.path(cd, "expression.tsv")
      cp <- file.path(cd, "clinical.tsv")
      if (!file.exists(ep) || !file.exists(cp)) next
      ve <- read_matrix_any(ep)
      vc <- data.table::fread(cp, data.table = FALSE)
      genes <- intersect(rownames(ve), unique(unlist(sig$gene_sets)))
      vg <- lapply(sig$gene_sets, function(g) intersect(g, genes))
      vg <- vg[lengths(vg) >= 5]
      if (length(vg) < 2) next
      vs <- tryCatch({
        param <- GSVA::gsvaParam(ve, vg, kcdf = "Gaussian", minSize = 5, maxSize = 500, absRanking = TRUE)
        GSVA::gsva(param, verbose = FALSE)
      }, error = function(e) NULL)
      if (is.null(vs)) next
      feats <- intersect(rownames(vs), ec$classifier$features)
      if (length(feats) < 2) next
      cent <- ec$classifier$centroids[, feats, drop = FALSE]
      Xv <- t(vs[feats, , drop = FALSE])
      d <- as.matrix(stats::dist(rbind(cent, Xv)))
      d <- d[seq_len(nrow(cent)), nrow(cent) + seq_len(nrow(Xv)), drop = FALSE]
      pred_lab <- rownames(cent)[apply(d, 2, which.min)]
      names(pred_lab) <- colnames(vs)
      vc$ecotype <- pred_lab[as.character(vc$sample_id)]
      if ("stage" %in% names(vc) && !"stage_num" %in% names(vc)) vc$stage_num <- as.numeric(factor(vc$stage, levels = c("I", "II", "III", "IV")))
      for (cv in covs) if (cv %in% names(vc) && is.character(vc[[cv]])) vc[[cv]] <- factor(vc[[cv]])
      vc$time <- as.numeric(vc[[paste0(endpoint, "_time_days")]])
      vc$event <- as.integer(vc[[paste0(endpoint, "_event")]])
      vc <- vc[!is.na(vc$time) & !is.na(vc$event) & vc$time > 0 & !is.na(vc$ecotype), , drop = FALSE]
      if (nrow(vc) < 20 || sum(vc$event) < 5) next
      vc$ecotype <- factor(vc$ecotype, levels = levels(dat$ecotype))
      vc_covs <- covs[covs %in% names(vc)]
      vc_covs <- vc_covs[vapply(vc[vc_covs], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
      vfit <- tryCatch(survival::coxph(as.formula(paste("Surv(time, event) ~ ecotype", if (length(vc_covs) > 0) paste("+", paste(vc_covs, collapse = " + ")) else "")), data = vc), error = function(e) NULL)
      if (is.null(vfit)) next
      coefs <- summary(vfit)$coefficients
      ecoef <- grep("^ecotype", rownames(coefs), value = TRUE)[1]
      hr_val <- if (!is.na(ecoef)) exp(coefs[ecoef, "coef"]) else NA_real_
      p_val <- if (!is.na(ecoef)) coefs[ecoef, "Pr(>|z|)"] else NA_real_
      disc_coef <- coef(full)[grep("^ecotype", names(coef(full)))[1]]
      dir_ok <- !is.na(hr_val) && !is.na(disc_coef) && sign(log(hr_val)) == sign(disc_coef)
      val_rows[[length(val_rows) + 1]] <- data.frame(cohort = basename(cd), n = nrow(vc), events = sum(vc$event), c_index = as.numeric(survival::concordance(vfit)$concordance), hr = hr_val, p_hr = p_val, direction_consistent = dir_ok, endpoint = endpoint)
    }
  }
  val_tab <- if (length(val_rows) > 0) do.call(rbind, val_rows) else data.frame(cohort = "none", n = 0, events = 0, c_index = NA_real_, hr = NA_real_, p_hr = NA_real_, direction_consistent = NA, endpoint = endpoint)
  data.table::fwrite(val_tab, file.path(out_dir, "validation_performance.tsv"), sep = "\t")

  # Evidence ledger: never write synthetic results as biological findings
  if (!synthetic) {
    ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
    old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
    if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% c("F01", "F03"), , drop = FALSE]
    sig <- !is.na(lrt_p) && lrt_p < 0.05
    f01 <- data.frame(
      finding_id = "F01",
      claim = if (sig) paste("Ecotype is associated with", toupper(endpoint), "in discovery") else paste("No significant incremental association between ecotype and", toupper(endpoint), "after adjustment"),
      evidence_type = "Cox model",
      dataset_s = "TCGA-COAD/READ",
      statistic = paste0("C-index=", round(c_full, 3), "; base=", round(c_base, 3), "; LRT p=", signif(lrt_p, 3)),
      n = nrow(dat),
      grade = if (sig) "E1" else "not_supported",
      verb_used = if (sig) "associated with" else "not associated",
      contradicting_evidence = "none",
      limitation = "retrospective public data; TCGA-COAD/READ only; OS not CSS",
      next_test = "independent external validation",
      figure_table = "Figure_2_km; Table_survival_performance",
      stringsAsFactors = FALSE
    )
    rows <- f01
    replicated <- any(val_tab$c_index > 0.55 & !is.na(val_tab$p_hr) & val_tab$p_hr < 0.05 & val_tab$direction_consistent, na.rm = TRUE)
    if (replicated) {
      vbest <- val_tab[which.max(val_tab$c_index), , drop = FALSE]
      f03 <- data.frame(
        finding_id = "F03",
        claim = paste("Ecotype association is replicated in", vbest$cohort),
        evidence_type = "External validation",
        dataset_s = vbest$cohort,
        statistic = paste0("C-index=", round(vbest$c_index, 3), "; HR=", round(vbest$hr, 3), "; p=", signif(vbest$p_hr, 3)),
        n = vbest$n,
        grade = "E2",
        verb_used = "replicated in",
        contradicting_evidence = "none",
        limitation = "cohort/platform heterogeneity",
        next_test = "orthogonal protein/spatial validation",
        figure_table = "Table_validation",
        stringsAsFactors = FALSE
      )
      rows <- rbind(f01, f03)
    } else if (any(val_tab$cohort != "none")) {
      v <- val_tab[1, , drop = FALSE]
      f03 <- data.frame(
        finding_id = "F03",
        claim = paste("Ecotype association did not replicate significantly in", v$cohort),
        evidence_type = "External validation",
        dataset_s = v$cohort,
        statistic = paste0("C-index=", round(v$c_index, 3), "; HR=", round(v$hr, 3), "; p=", signif(v$p_hr, 3)),
        n = v$n,
        grade = "not_supported",
        verb_used = "not replicated",
        contradicting_evidence = "external cohort",
        limitation = "cohort/platform heterogeneity; ecotype HR not significant",
        next_test = "additional independent cohorts",
        figure_table = "Table_validation",
        stringsAsFactors = FALSE
      )
      rows <- rbind(f01, f03)
    }
    if (!is.null(old) && nrow(old) > 0) rows <- rbind(old, rows)
    data.table::fwrite(rows, ledger_path, sep = ",")
  } else {
    syn_dir <- file.path(cfg$paths$analysis, "06_Other_Analysis")
    data.table::fwrite(data.frame(note = "SYNTHETIC DRY RUN — NOT BIOLOGICAL EVIDENCE", c_index = c_full, lrt_p = lrt_p, n = nrow(dat), events = sum(dat$event)), file.path(syn_dir, "synthetic_survival_summary.tsv"), sep = "\t")
  }

  log_decision("04_survival", "PASS", paste("endpoint:", endpoint, "; n:", nrow(dat), "; events:", sum(dat$event), "; C-index:", round(c_full, 3), "; LRT p:", signif(lrt_p, 3), "; PH min p:", signif(ph_p, 3)))
  saveRDS(list(perf = perf, coef = coef_tab, cal = cal_tab, validation = val_tab, dat = dat, km = km, full = full, base = base), file.path(cfg$paths$processed, "survival_results.rds"))
  invisible(TRUE)
}
