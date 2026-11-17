# QC Report — CRC Track A

QC date: 2026-11-17 01:05:20
Script: `06_Analysis/1_R/1_Rscripts/01_bulk_qc.R`
Samples in: 120 → after QC: 120 (no samples dropped at this step)
Features: 2000 → 2000
Normalization: quantile-normalized log2
Batch: corrected (limma removeBatchEffect; 3 levels)

## Per-sample summary
Table: `05_QC/Sample_QC.tsv`.
Outliers flagged: 1 ( 0.8 % )

## Decision
QC checkpoint: PASS
Outliers were flagged, not silently removed. Sensitivity analysis will exclude them.

## Limitations
Purity and tissue composition were not available at QC time; they are handled in the model/sensitivity analysis.
