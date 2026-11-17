# K=4 validation report — CRC Track A

**Status:** K=4 exploratory signal was tested in two independent cohorts and did **not** replicate significantly.

## Method

- K=4 labels were derived in TCGA-COAD/READ with SNF spectral clustering (`K=4`).
- A nearest-centroid classifier was trained on TCGA Hallmark scores.
- The classifier was frozen and applied to GSE39582 and GSE17536.
- Cox model: ecotype + age + sex + stage (+ MSI/sidedness where available).
- LRT compares the ecotype model with the covariate-only base model.

## Results

| Cohort | n | events | C-index | LRT p | Interpretation |
|---|---|---|---|---|---|
| TCGA-COAD/READ (discovery) | 476 | 79 | 0.785 | 0.024 | discovery signal |
| GSE39582 | 573 | 190 | 0.614 | 0.090 | not significant |
| GSE17536 | 177 | 73 | 0.772 | 0.200 | not significant |

## Conclusion

The K=4 partition does not replicate at p<0.05 in either independent cohort. It remains an **E1 exploratory** discovery finding and must not be presented as a validated taxonomy. The primary K=2 solution also did not replicate significantly (see `Validation_Report.md`).
