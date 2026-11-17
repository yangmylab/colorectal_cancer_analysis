# Robustness and sensitivity report — CRC Track A

Generated: 2026-11-17 08:43:28

## Primary result
- K = 2 , n = 476 , events = 79
- C-index base/full: 0.763 / 0.772
- LRT p: 0.107

## K sensitivity (same outcome, same covariates)

| K | n | events | C-index base | C-index full | LRT p | min group fraction |
|---|---|---|---|---|---|---|
| 2 | 476 | 79 | 0.763 | 0.772 | 0.107 | 0.414 |
| 3 | 476 | 79 | 0.763 | 0.784 | 0.055 | 0.252 |
| 4 | 476 | 79 | 0.763 | 0.785 | 0.0236 | 0.166 |
| 5 | 476 | 79 | 0.763 | 0.78 | 0.0618 | 0.105 |

## Outlier-exclusion sensitivity
- Excluding 14 QC outliers: n = 463 , events = 77
- C-index base/full: 0.767 / 0.773 ; LRT p: 0.151

## Decision
A finer K (4) improves C-index by 0.013 and reaches LRT p < 0.05; report as exploratory and consider a pre-specified follow-up.

## Interpretation
The primary K=2 solution is retained because it is the selected stability solution in the discovery data. Finer partitions are reported as exploratory sensitivity analyses.

## Limitations
- TCGA-COAD/READ only; OS not CSS; MSI and purity unavailable.
- No independent external validation yet.
- Hallmark programmes only; no protein, methylation, single-cell or spatial layers in this run.
