# Overnight analysis status — CRC Track A

**Date:** 2026-11-17
**Project:** `/run/media/meng/Seqdata/colorectal_cancer/CRC_BAP`

## Bottom line

A **real discovery analysis on TCGA-COAD/READ is complete**.

- **624 patients** with primary CRC; 15,000 genes; 128 OS events in the clinical table, 79 in the complete-case survival model.
- Stable **K=2 ecotype** solution (bootstrap ARI 0.975; classifier CV 0.979).
- **Primary ecotype–OS association is not significant after adjustment**: C-index 0.763 → 0.772; LRT p=0.107.
- **K=4 exploratory signal**: C-index 0.785; LRT p=0.024; **post hoc**, E1 only, requires external validation.
- **External validation completed in GSE39582 and GSE17536**: GSE39582 C-index 0.609, HR 1.16, p=0.333; GSE17536 C-index 0.765, HR 1.53, p=0.0785; two-cohort meta HR 1.25 (95% CI 0.97–1.61), p=0.081. Direction consistent but **not significant**; **not replicated**; no E2 claim.
- **K=4 exploratory partition not replicated**: GSE39582 p=0.090; GSE17536 p=0.200.
- No wet-lab validation.

## Data acquired

- TCGA-COAD/READ RNA-seq: **701 files, 2.96 GB**, GDC Data Release 46.0.
- 647 primary-tumour files collapsed to **624 patients** (highest total TPM per patient).
- 51 solid-tissue normal, 1 metastatic, 2 recurrent files downloaded but not used in discovery.
- External validation: **GSE39582** series matrix complete (165 MB); 585 tumours; 23,963 Ensembl-mapped genes; 194 OS events.
- External validation: **GSE17536** complete (50 MB); 177 tumours; 23,963 Ensembl-mapped genes; 73 OS events.
- **GSE14333** complete (81 MB); DFS-only endpoint, not used for the OS validation.
- Clinical: 458 COAD / 166 READ; stage I 96, II 196, III 150, IV 64, missing 118; sex M 330/F 291/missing 3; OS time missing 66.

## Results

| Item | Result |
|---|---|
| Samples / patients | 624 / 624 |
| Genes | 15,000 |
| QC outliers | 14 / 624 (2.2%) |
| Signatures | 50 Hallmark sets |
| Ecotypes | K=2: E1 355, E2 269 |
| E1 top programmes | Hedgehog, pancreas beta cells, MYC targets V2 |
| E2 top programmes | EMT, angiogenesis, IL6–JAK–STAT3 |
| Bootstrap ARI | 0.975 |
| Classifier CV accuracy | 0.979 |
| Survival n / events | 476 / 79 |
| C-index base → full | 0.763 → 0.772 |
| LRT p (K=2) | 0.107 |
| K=4 C-index / p | 0.785 / 0.024 (post hoc, E1) |
| GSE39582 validation C-index / HR / p | 0.609 / 1.16 / 0.333 (not significant) |
| GSE17536 validation C-index / HR / p | 0.765 / 1.53 / 0.0785 (not significant) |
| Two-cohort meta-analysis | HR 1.25 (0.97–1.61); p=0.081 |
| K=4 validation (GSE39582 / GSE17536) | p=0.090 / 0.200 (not replicated) |
| Outlier-excluded C-index / p | 0.773 / 0.151 |

## Cross-platform layers (new, 2026-11-17)

| Layer | Dataset | Result |
|---|---|---|
| Single cell (descriptive) | GSE132465 (63,689 cells, 6 lineages) | used for NNLS deconvolution of TCGA only |
| **Protein** | CPTAC-COAD via PDC000116 / cBioPortal (8,067 proteins, 97 samples) | mRNA-vs-protein **ecotype agreement 85.4 %**, ARI 0.496, median programme r=0.41 (IFN-α 0.86, IFN-γ 0.80, EMT 0.71); OS underpowered (8 events): HR 2.44, p=0.286 |
| **Methylation** | TCGA-COAD/READ 450K (458 files, 6.0 GB → 19,269 promoter genes, 393 tumours) | ecotype **does not reproduce**: ARI vs methylation cluster 0.006; 5/19,269 genes FDR<0.05 and \|logFC\|>0.5; global beta p=0.113; CIMP-high 74.5 % E1 vs 84.3 % E2; cluster–OS HR 0.74, p=0.398 |
| **Spatial (descriptive)** | GSE283052 (12 sections) + GSE267401 (4 sections) = 16 sections / 43,818 spots | bulk ecotype classifier **does not transfer** (self-check ARI 0.042; all spots on one side of the boundary; scores ≈5 SD below bulk). E2 programmes vary as intra-tumour niches; no ecotype claim |

Reading of the layers: the K=2 axis is reproducible at the **protein** level and is **not** a
methylation-defined class. The methylome therefore does not explain the EMT/angiogenic programme of E2.

## Packages installed (with user approval)

`minfi`, `IlluminaHumanMethylation450kanno.ilmn12.hg19` (+ dependencies `illuminaio`, etc.).
No Python packages were installed; the existing `sc-analysis` conda environment (numpy/pandas/scanpy) was used.

## Gate status

| Gate | Status |
|---|---|
| G0 intake | PASS |
| G1 annotation | PASS (624 samples, 624 donors, 0 ambiguities) |
| G2 design | PASS |
| G3 final | WARN (null findings F01/F03/F04 ungraded; E1 for K=4 exploratory) |

## What remains

1. **Additional external validation**: GSE39582 and GSE17536 tested (not significant); GSE14333 is complete but DFS-only, so it cannot serve the OS validation — use it only for a secondary DFS analysis.
2. **K=4 hypothesis**: already tested in GSE39582 (p=0.090) and GSE17536 (p=0.200); not replicated. Do not present it as confirmatory.
3. **Protein/methylation layers**: done (CPTAC-COAD protein, TCGA 450K methylation) — see the cross-platform table above.
4. **Spatial layer**: done descriptively (16 sections, 43,818 spots) — see `CRC_BAP/07_Validation/Spatial_Report.md`. A genuine spatial ecotype claim needs a spot-level classifier trained on single-cell/spatial data (cell2location/RCTD or CNSA CNP0003321).
5. **Phosphoproteome**: CPTAC phosphoproteome table downloaded (`data_phosphoprotein_quantification.txt`), kinase-activity analysis not yet run.
6. **Manuscript**: current result is a stable K=2 axis with a null primary survival association and a negative external validation; the framing must reflect that.
7. **Optional packages**: MOFA2, ConsensusClusterPlus, decoupleR, minfi, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC, survminer are still not installed; the pipeline used available fallbacks.

## Figures

The final figure set is numbered `plot_1` … `plot_15` (PDF + JPEG, 300 dpi) in `CRC_BAP/09_Figures/plots/`, with
`plot_manifest.tsv` mapping each number to its source file and description.

| Plot | Content |
|---|---|
| plot_1 | GSE132465 single-cell marker signatures |
| plot_2–plot_5 | QC: PCA uncorrected / batch-corrected, library size, top-500 variable features |
| plot_6 | Hallmark signature heatmap |
| plot_7 | Consensus matrix (selected K) |
| plot_8 | SNF ecotype PCA |
| plot_9 | Ecotype stability across K |
| plot_10 | Kaplan–Meier by ecotype |
| plot_11 | 5-year calibration |
| plot_12 | 5-year decision curve |
| plot_13 | Mean TCGA cell-type fractions (NNLS) |
| plot_14 | K sensitivity of the survival model |
| plot_15 | External-validation forest plot |
| plot_16–plot_19 | CPTAC-COAD protein layer: programme heatmap, mRNA-vs-protein concordance, KM, ecotype agreement |
| plot_20–plot_24 | TCGA 450K methylation layer: CIMP, global methylation, top DMPs, KM, cluster concordance |
| plot_25–plot_27 | Visium spatial layer: E2-programme maps, within- vs between-section variation, covariate-shift diagnostic |

## Files

- `CRC_BAP/09_Figures/plots/plot_1` … `plot_15` (+ `plot_manifest.tsv`)
- `CRC_BAP/13_Reports/Final_Report.md`
- `CRC_BAP/13_Reports/Findings_Ledger.csv`
- `CRC_BAP/06_Analysis/03_Clustering/stability.tsv`
- `CRC_BAP/06_Analysis/03_Clustering/ecotype_labels.tsv`
- `CRC_BAP/06_Analysis/06_Other_Analysis/survival_performance.tsv`
- `CRC_BAP/08_Robustness/Robustness_Report.md`
- `CRC_BAP/Logs/DECISION_LOG.tsv`
- `CRC_BAP/99_Synthetic_DryRun/` — synthetic code validation only, not biological evidence

## Notes

- No software was installed during this run.
- The synthetic dry run remains isolated and is not part of the biological results.
- The READY watcher was stopped after the real run.
- GSE39582 is complete and stored at `03_Raw_Data/GSE39582/GSE39582_series_matrix.txt.gz`; its validation inputs are under `03_Raw_Data/validation/GSE39582/`.
