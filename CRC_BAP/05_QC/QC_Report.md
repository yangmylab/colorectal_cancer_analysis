# QC Report — CRC Track A

QC date: 2026-11-17 02:14:55
Script: `06_Analysis/1_R/1_Rscripts/01_bulk_qc.R`
Samples in: 624 → after QC: 624 (no samples dropped at this step)
Features: 15000 → 15000
Normalization: quantile-normalized log2
Batch: corrected (limma removeBatchEffect; 2 levels)

## Per-sample summary
Table: `05_QC/Sample_QC.tsv`.
Outliers flagged: 14 ( 2.2 % )

## Decision
QC checkpoint: PASS
Outliers were flagged, not silently removed. Sensitivity analysis will exclude them.

## Limitations
Purity and tissue composition were not available at QC time; they are handled in the model/sensitivity analysis.
