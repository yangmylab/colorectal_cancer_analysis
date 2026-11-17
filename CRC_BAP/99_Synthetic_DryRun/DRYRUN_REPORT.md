# Synthetic dry-run report

**SYNTHETIC PIPELINE VALIDATION — NOT BIOLOGICAL EVIDENCE.**

This run used simulated data solely to verify that every module executes and that the decision logic works.

## Results
- Samples: 120
- Features: 2000
- Signatures scored: 12
- Selected K: 4
- Best bootstrap ARI: 1
- Classifier CV accuracy: 1
- Survival endpoint: css
- Events: 115
- C-index (full): 0.674
- C-index (base): 0.636
- LRT p: 0.000178
- Validation cohorts: 1

## Decision
The pipeline ran end-to-end. Synthetic outputs are isolated under `99_Synthetic_DryRun/` and must never be reported as findings.

## What remains blocked
The real analysis is blocked at G1 because no CRC datasets are present in `03_Raw_Data/`.
