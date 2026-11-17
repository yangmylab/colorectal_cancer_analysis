# CRC Track A — governed analysis project

This directory is the `/bioanalyst` BAP project for the cross-platform CRC ecotype study.

## Status

- **Real discovery analysis completed on TCGA-COAD/READ.**
- G0 PASS, **G1 PASS** (624 tumours / 624 donors), G2 PASS, G3 WARN (null F01/F03 ungraded; E1 for K=4 exploratory).
- Stable K=2 ecotype (ARI 0.975); primary ecotype–OS association **not significant** (C-index 0.763 → 0.772; LRT p=0.107).
- K=4 exploratory signal (C-index 0.785; LRT p=0.024; post hoc).
- **GSE39582 external validation completed**: C-index 0.609; HR 1.16; p=0.333; not significant.
- **GSE17536 external validation completed**: C-index 0.765; HR 1.53; p=0.0785; not significant.
- Two-cohort meta-analysis: HR 1.25 (95% CI 0.97–1.61); p=0.081.
- **K=4 exploratory partition did not replicate**: GSE39582 p=0.090; GSE17536 p=0.200.
- No wet-lab validation.

## Key results

- E1 (355 tumours): Hedgehog, pancreas beta cells, MYC targets V2.
- E2 (269 tumours): EMT, angiogenesis, IL6–JAK–STAT3.
- 476 complete-case patients, 79 OS events.
- 14/624 QC outliers.
- Synthetic dry run is separate under `99_Synthetic_DryRun/` and is not biological evidence.

## Run the pipeline

```bash
cd CRC_BAP
bash run_all.sh config/pipeline_config.yaml
```

Current inputs:

- `03_Raw_Data/expression.tsv` — 15,000 genes × 624 samples
- `03_Raw_Data/clinical.tsv` — 624 rows

## Add external validation

GSE39582 is already built under `03_Raw_Data/validation/GSE39582/` and was analysed
(C-index 0.609; HR 1.16; p=0.333; not significant).

For additional cohorts (GSE17536, GSE14333), adapt the GEO converter:

```bash
Rscript 06_Analysis/1_R/1_Rscripts/10_geo_series_to_inputs.R \
  --series-matrix 03_Raw_Data/<cohort>/<cohort>_series_matrix.txt.gz \
  --outdir 03_Raw_Data/validation/<cohort>
# fill survival fields, then rerun 04_survival_validation.R
```

## Pipeline

| Step | Script | Function |
|---|---|---|
| 00 | `06_Analysis/1_R/1_Rscripts/00_validate_inputs.R` | schema + G1 gate |
| 01 | `01_bulk_qc.R` | normalization, QC, PCA, batch |
| 02 | `02_signatures_deconvolution.R` | GSVA ssGSEA, optional NNLS |
| 03 | `03_snf_ecotypes.R` | SNF, consensus clustering, classifier |
| 04 | `04_survival_validation.R` | Cox, calibration, DCA, validation |
| 05 | `05_figures_report.R` | figures, decision log, report |
| 07 | `07_robustness.R` | K and outlier sensitivity |
| 99 | `99_synthetic_dryrun.R` | synthetic code test only |

## Key files

- `13_Reports/Final_Report.md`
- `13_Reports/Findings_Ledger.csv`
- `06_Analysis/03_Clustering/`
- `06_Analysis/06_Other_Analysis/`
- `08_Robustness/Robustness_Report.md`
- `Logs/DECISION_LOG.tsv`

## Data provenance

- GDC Data Release 46.0; 701 RNA-seq files; size/MD5 verified.
- 647 primary-tumour files → 624 patients.
- Stage missing 118; OS time missing 66; MSI/purity unavailable.
- GSE39582 series matrix complete; validation inputs under `03_Raw_Data/validation/GSE39582/`.
- GSE17536 series matrix complete; validation inputs under `03_Raw_Data/validation/GSE17536/`.
- GSE14333 complete; DFS-only, not used for the OS validation.

## Notes

- No software was installed during this run.
- Optional packages (MOFA2, ConsensusClusterPlus, decoupleR, minfi, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC, survminer) are not installed; the pipeline uses available fallbacks.
