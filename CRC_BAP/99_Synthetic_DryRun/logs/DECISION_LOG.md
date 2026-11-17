# Decision Log

Generated: 2026-11-17 01:05:36

| time | step | decision | rationale |
|---|---|---|---|
| 2026-11-17 01:05:18 | 00_validate | SYNTHETIC_G1 | synthetic annotation checked; project G1 remains BLOCKED (no real data) |
| 2026-11-16 17:05:18 | 00_validate | PASS | samples: 120 genes: 2000 |
| 2026-11-16 17:05:20 | 01_bulk_qc | PASS | normalization: quantile-normalized log2 ; outliers: 1 of 120 ; batch: corrected (limma removeBatchEffect; 3 levels) |
| 2026-11-16 17:05:26 | 02_signatures | PASS | source: 99_Synthetic_DryRun/data/signatures.gmt ; sets: 12 ; deconvolution: No cell-type signature matrix supplied; signature scores are marker scores, not absolute fractions |
| 2026-11-16 17:05:34 | 03_ecotypes | PASS | K = 4 ; bootstrap ARI = 1 ; CV accuracy = 1 |
| 2026-11-16 17:05:36 | 04_survival | PASS | endpoint: css ; n: 120 ; events: 115 ; C-index: 0.674 ; LRT p: 0.000178 ; PH min p: 0.404 |
