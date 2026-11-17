#!/usr/bin/env Rscript
# 23_kinase_activity.R
# Kinase-substrate activity inference in CPTAC-COAD: does the E2-vs-E1 difference involve kinase activity
# beyond what kinase abundance explains?
#
# Why this script exists: 21_phospho_analysis.R reported differential phosphosites but explicitly withheld
# any kinase-activity claim, because no kinase-substrate library could be obtained offline. That gap is now
# closed with KSEAapp::KSData (PhosphoSitePlus-derived; 19,739 kinase-substrate pairs, 387 kinases) — an
# offline prior, no network call.
#
# Inputs : 03_Raw_Data/CPTAC_COAD/data_phosphoprotein_quantification.txt  (site matrix, re-filtered as in 21)
#          03_Raw_Data/CPTAC_COAD/data_protein_quantification.txt        (kinase abundance)
#          07_Validation/21_phospho_ecotype_labels.tsv                   (frozen phospho ecotype calls)
#          KSEAapp::KSData                                               (kinase-substrate prior)
# Outputs: 07_Validation/23_kinase_activity.tsv, 23_kinase_activity_report.md
#          09_Figures/plots/plot_37 .. plot_38
#
# Parallelism (one level only, BLAS pinned to 1 thread per worker):
#   workers = 6 (physical cores of the i5-10400); work split = per-kinase activity + activity/abundance r

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(mclust)
  library(parallel)
})

source("06_Analysis/1_R/1_Rscripts/utils.R")

N_WORKERS <- 6L
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(mc.cores = N_WORKERS)
set.seed(20261117)

IN <- "03_Raw_Data/CPTAC_COAD"
OUT <- "07_Validation"
FIG <- "09_Figures/plots"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## 1. Site matrix, filtered and imputed exactly as in 21_phospho_analysis.R ---------
ph <- fread(file.path(IN, "data_phosphoprotein_quantification.txt"), data.table = FALSE, check.names = FALSE)
site_ids <- ph$ENTITY_STABLE_ID
genes <- ph$GENE_SYMBOL
residue <- toupper(gsub("[^A-Za-z0-9]", "", sub("^p", "", ph$PHOSPHOSITE)))
site <- as.matrix(ph[, -(1:3), drop = FALSE])
storage.mode(site) <- "numeric"
rownames(site) <- make.unique(site_ids)

miss <- rowMeans(is.na(site))
keep <- which(miss <= 0.50)               # identical threshold to 21 (fixed a priori there)
site <- site[keep, , drop = FALSE]
genes <- genes[keep]
residue <- residue[keep]
cat("sites analysed:", nrow(site), "of", length(miss), "\n")

impute_knn <- function(m) {
  ok <- tryCatch(impute::impute.knn(m, k = 10, rowmax = 0.5, colmax = 0.9, rng.seed = 20261117),
                 error = function(e) NULL)
  if (!is.null(ok)) return(ok$data)
  for (i in seq_len(nrow(m))) {
    na <- is.na(m[i, ])
    if (any(na)) m[i, na] <- stats::median(m[i, ], na.rm = TRUE)
  }
  m
}
site_z <- t(scale(t(impute_knn(site))))
site_z <- site_z[is.finite(rowSums(site_z)), , drop = FALSE]
idx <- match(rownames(site_z), rownames(site))
sites <- data.table(site_id = rownames(site_z), gene = toupper(genes[idx]), res = residue[idx])

## 2. Kinase-substrate prior (offline) ---------------------------------------
data(KSData, package = "KSEAapp")
ks <- as.data.table(KSData)
ks <- unique(ks[, .(kinase = toupper(GENE),
                    substrate = toupper(SUB_GENE),
                    res = toupper(gsub("[^A-Za-z0-9]", "", SUB_MOD_RSD)))])
ks <- ks[nzchar(kinase) & nzchar(substrate) & nzchar(res)]
cat("prior:", nrow(ks), "kinase-substrate pairs |", length(unique(ks$kinase)), "kinases\n")

pairs <- merge(sites, ks, by.x = c("gene", "res"), by.y = c("substrate", "res"), allow.cartesian = TRUE)
cat("measured pairs:", nrow(pairs), "| sites:", length(unique(pairs$site_id)),
    "| kinases:", length(unique(pairs$kinase)), "\n")
sets <- split(pairs$site_id, pairs$kinase)
MIN_SUBSTRATES <- 5L
sets <- lapply(sets, unique)
sets <- sets[lengths(sets) >= MIN_SUBSTRATES]
cat("kinases with >=", MIN_SUBSTRATES, "measured substrates:", length(sets), "\n")

## 3. Per-sample kinase activity ---------------------------------------------
# vapply rather than cbind(mclapply(...)): a missing substrate set must fail loudly instead of being
# silently dropped from the matrix (that silently reduced 85 kinases to 37 during development).
row_of <- setNames(seq_len(nrow(site_z)), rownames(site_z))
act <- vapply(sets, function(ids) {
  rr <- row_of[ids]
  rr <- rr[!is.na(rr)]
  if (length(rr) < MIN_SUBSTRATES) return(rep(NA_real_, ncol(site_z)))
  colMeans(site_z[rr, , drop = FALSE], na.rm = TRUE)
}, numeric(ncol(site_z)))
act <- t(act)                                   # kinases x samples (limma layout)
act <- act[rowSums(is.na(act)) == 0, , drop = FALSE]
cat("activity matrix:", nrow(act), "kinases x", ncol(act), "tumours\n")

## 4. Ecotype labels (from 21) and activity test -----------------------------
lab <- fread(file.path(OUT, "21_phospho_ecotype_labels.tsv"), data.table = FALSE)
samples <- intersect(colnames(act), lab$sample_id)   # act is kinases x samples after the transpose
grp <- lab$ecotype_phospho[match(samples, lab$sample_id)]
ok <- !is.na(grp)
cat("tumours with an ecotype label:", sum(ok), "| E1", sum(grp[ok] == "E1"), "E2", sum(grp[ok] == "E2"), "\n")

design <- model.matrix(~ 0 + factor(grp[ok], levels = c("E1", "E2")))
colnames(design) <- c("E1", "E2")
fit <- eBayes(contrasts.fit(lmFit(act[, samples[ok], drop = FALSE], design),
                            makeContrasts(E2 - E1, levels = design)))
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
act_tab <- data.table(kinase = rownames(tt), logFC_activity = tt$logFC, t = tt$t,
                      p = tt$P.Value, fdr = tt$adj.P.Val, n_substrates = lengths(sets[rownames(tt)]))
zE1 <- rowMeans(act[, samples[ok][grp[ok] == "E1"], drop = FALSE])
zE2 <- rowMeans(act[, samples[ok][grp[ok] == "E2"], drop = FALSE])
sd_all <- apply(act[, samples[ok], drop = FALSE], 1, stats::sd)
act_tab[, cohens_d := (zE2[kinase] - zE1[kinase]) / sd_all[kinase]]

## 5. Kinase abundance (protein) ---------------------------------------------
prot <- fread(file.path(IN, "data_protein_quantification.txt"), data.table = FALSE, check.names = FALSE)
pm <- as.matrix(prot[, -1, drop = FALSE])
storage.mode(pm) <- "numeric"
rownames(pm) <- sub("\\|.*$", "", prot[[1]])
if (any(duplicated(rownames(pm)))) {
  v <- apply(pm, 1, function(x) stats::var(x, na.rm = TRUE))
  pm <- pm[order(rownames(pm), -v), , drop = FALSE]
  pm <- pm[!duplicated(rownames(pm)), , drop = FALSE]
}
kin_in_prot <- intersect(act_tab$kinase, rownames(pm))
cat("kinases also quantified at the protein level:", length(kin_in_prot), "\n")

prot_diff <- do.call(rbind, mclapply(kin_in_prot, function(k) {
  x1 <- pm[k, samples[ok][grp[ok] == "E1"]]
  x2 <- pm[k, samples[ok][grp[ok] == "E2"]]
  ok1 <- is.finite(x1); ok2 <- is.finite(x2)
  if (sum(ok1) < 10 || sum(ok2) < 10) return(NULL)
  wt <- stats::t.test(x2[ok2], x1[ok1])
  data.frame(kinase = k,
             prot_diff = mean(x2[ok2]) - mean(x1[ok1]),
             prot_p = wt$p.value,
             r_activity_abundance = suppressWarnings(stats::cor(act[k, samples[ok]], as.numeric(pm[k, samples[ok]]), use = "complete.obs")),
             stringsAsFactors = FALSE)
}, mc.cores = N_WORKERS))
act_tab <- merge(act_tab, prot_diff, by = "kinase", all.x = TRUE)
act_tab[, activity_not_explained := !is.na(prot_p) & prot_p > 0.05]

## 6. Results -----------------------------------------------------------------
hits <- act_tab[fdr < 0.05]
hits_unexplained <- hits[!is.na(prot_p) & prot_p > 0.05]
cat("\nkinases tested:", nrow(act_tab), "| activity FDR<0.05:", nrow(hits),
    "| of those with unchanged abundance:", nrow(hits_unexplained), "\n")
if (nrow(hits)) print(head(hits[order(-abs(cohens_d)), c("kinase", "cohens_d", "fdr", "n_substrates", "prot_diff", "prot_p")], 12))

cor_eff <- suppressWarnings(stats::cor(act_tab$logFC_activity, act_tab$prot_diff, use = "complete.obs"))
fwrite(act_tab[order(p)], file.path(OUT, "23_kinase_activity.tsv"), sep = "\t")

## 7. Figures -----------------------------------------------------------------
act_tab[, sig := ifelse(fdr < 0.05 & activity_not_explained, "activity only (FDR<0.05, abundance ns)",
                        ifelse(fdr < 0.05, "activity + abundance (FDR<0.05)", "ns"))]
top_lab <- head(act_tab[order(p)][fdr < 0.10], 12)
save_plot(function() {
  print(ggplot(act_tab, aes(cohens_d, -log10(p), colour = sig)) +
    geom_point(size = 1.4, alpha = 0.8) +
    scale_colour_manual(values = c("activity only (FDR<0.05, abundance ns)" = "#C0392B",
                                   "activity + abundance (FDR<0.05)" = "#E67E22", ns = "grey70")) +
    geom_text_repel(data = top_lab, aes(label = kinase), size = 2.4, max.overlaps = 20) +
    theme_bw(base_size = 9) +
    labs(x = "Kinase activity effect, Cohen's d (E2 vs E1)", y = expression(-log[10]~P),
         colour = NULL,
         title = sprintf("Kinase-substrate activity inference (%d kinases, >=%d substrates)",
                         nrow(act_tab), MIN_SUBSTRATES)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)))
}, file.path(FIG, "plot_37"), 8, 6)

both <- act_tab[!is.na(prot_diff)]
save_plot(function() {
  print(ggplot(both, aes(prot_diff, cohens_d, colour = sig)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
    scale_colour_manual(values = c("activity only (FDR<0.05, abundance ns)" = "#C0392B",
                                   "activity + abundance (FDR<0.05)" = "#E67E22", ns = "grey70")) +
    geom_text_repel(data = head(both[order(p)][fdr < 0.10], 10), aes(label = kinase), size = 2.3, max.overlaps = 15) +
    theme_bw(base_size = 9) +
    labs(x = "Kinase protein abundance difference (E2 - E1, log ratio)",
         y = "Kinase activity effect (Cohen's d)",
         colour = NULL,
         title = sprintf("Activity vs abundance for %d kinases (r = %.2f)", nrow(both), cor_eff)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7)))
}, file.path(FIG, "plot_38"), 8, 6)

## 8. Report + ledger ---------------------------------------------------------
grade <- if (nrow(hits_unexplained) > 0) "E1" else "not_supported"
report <- c(
  "# Kinase-substrate activity inference (CPTAC-COAD)",
  "",
  sprintf("Generated by `23_kinase_activity.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Prior and data",
  sprintf("- Kinase-substrate prior: KSEAapp::KSData (PhosphoSitePlus), %d pairs / %d kinases, **offline**.",
          nrow(ks), length(unique(ks$kinase))),
  sprintf("- Measured pairs after matching to the %d filtered phosphosites: %d pairs / %d sites.",
          nrow(site), nrow(pairs), length(unique(pairs$site_id))),
  sprintf("- Kinases with >=%d measured substrates: %d (of %d with any substrate). Median substrates per testable kinase: %d.",
          MIN_SUBSTRATES, nrow(act_tab), length(sets), as.integer(stats::median(act_tab$n_substrates))),
  "",
  "## Activity differs between ecotypes?",
  sprintf("- Tumours with labels: %d (E1 %d / E2 %d), activity scored as the mean z-score of a kinase's measured substrate sites.",
          sum(ok), sum(grp[ok] == "E1"), sum(grp[ok] == "E2")),
  sprintf("- Kinases at FDR<0.05: **%d of %d**; %d of those have unchanged kinase abundance (p>0.05).",
          nrow(hits), nrow(act_tab), nrow(hits_unexplained)),
  if (nrow(hits)) {
    paste0("- Largest effects: ",
           paste(sprintf("%s (d=%.2f, FDR=%.3g, %d substrates)", hits$kinase[order(-abs(hits$cohens_d))][1:min(6, nrow(hits))],
                         hits$cohens_d[order(-abs(hits$cohens_d))][1:min(6, nrow(hits))],
                         hits$fdr[order(-abs(hits$cohens_d))][1:min(6, nrow(hits))],
                         hits$n_substrates[order(-abs(hits$cohens_d))][1:min(6, nrow(hits))]),
                collapse = "; "), ".")
  } else "- No kinase reached FDR<0.05.",
  "",
  "## Activity versus abundance",
  sprintf("- Across the %d kinases quantified at both levels, activity effect vs abundance difference r = %.2f.",
          nrow(both), cor_eff),
  sprintf("- %d kinases change in activity without a detectable change in their own protein level.",
          nrow(hits_unexplained)),
  "",
  "## Grade",
  sprintf("- **%s** — %s", grade,
          if (grade == "E1") {
            paste("kinase activity (inferred from substrate phosphorylation) separates the ecotypes for a subset of kinases,",
                  "and some of those are not explained by kinase abundance; single cohort, inferred rather than measured activity,",
                  "so exploratory (E1) only.")
          } else {
            paste("no kinase activity signal survives FDR control, so the phosphoproteome difference (F09) is not attributable",
                  "to a specific kinase in this cohort at this depth.")
          }),
  "",
  "## Limitations",
  "- Activity is **inferred** from substrate phosphorylation, not measured (no kinase-activity assay).",
  "- PhosphoSitePlus substrate coverage is incomplete and biased toward well-studied kinases; median substrates per testable kinase is low.",
  "- Substrate sites are reused across kinases (hierarchical kinase families), so scores are not independent.",
  "- One cohort, 96 tumours with labels; no external phosphoproteome exists in this project to replicate.",
  "- Kinase abundance is mass-spectrometry protein level, which is itself noisy at this depth.",
  "",
  "## Figures",
  "- plot_37 kinase-activity volcano; plot_38 activity effect vs kinase abundance difference."
)
writeLines(report, file.path(OUT, "23_kinase_activity_report.md"))

ledger_path <- file.path("13_Reports", "Findings_Ledger.csv")
old <- if (file.exists(ledger_path)) data.table::fread(ledger_path, data.table = FALSE) else NULL
if (!is.null(old) && nrow(old) > 0) old <- old[!old$finding_id %in% "F12", , drop = FALSE]
f12 <- data.frame(
  finding_id = "F12",
  claim = if (nrow(hits_unexplained) > 0) {
    "A subset of kinases shows ecotype-associated activity that is not explained by kinase abundance (inferred)"
  } else {
    "No kinase-activity signal separates the ecotypes beyond kinase abundance at this depth"
  },
  evidence_type = "Kinase-substrate activity inference",
  dataset_s = "CPTAC-COAD (PDC000116) phosphoproteome",
  statistic = paste0(nrow(act_tab), " kinases tested (>=", MIN_SUBSTRATES, " measured substrates); ",
                     nrow(hits), " at FDR<0.05; ", nrow(hits_unexplained),
                     " with unchanged kinase abundance; activity-vs-abundance r=", round(cor_eff, 2),
                     if (nrow(hits)) paste0("; top: ", paste(head(hits$kinase[order(-abs(hits$cohens_d))], 5), collapse = ", ")) else ""),
  n = sum(ok),
  grade = grade,
  verb_used = if (grade == "E1") "associated with" else "not associated",
  contradicting_evidence = "not applicable",
  limitation = "activity inferred from substrate phosphorylation, not measured; PhosphoSitePlus coverage biased and substrates shared across kinases; single cohort; no external phosphoproteome",
  next_test = "kinase-activity assay or phosphoproteome cohort with matched abundance; replicate in an independent CPTAC-like dataset",
  figure_table = "plot_37; plot_38",
  stringsAsFactors = FALSE
)
out <- if (is.null(old) || nrow(old) == 0) f12 else rbind(old, f12)
out <- out[order(as.integer(sub("^F", "", out$finding_id))), , drop = FALSE]
data.table::fwrite(out, ledger_path, sep = ",")
cat("ledger rows:", nrow(out), "\n")

log_decision("23_kinase", "PASS",
             paste0("kinases tested=", nrow(act_tab), "; FDR<0.05=", nrow(hits),
                    "; activity-only=", nrow(hits_unexplained), "; activity-vs-abundance r=", round(cor_eff, 2),
                    "; grade=", grade))
cat("grade:", grade, "\n")
cat("done\n")
