# Decision Log

Generated: 2026-11-17 08:44:13

| time | step | decision | rationale |
|---|---|---|---|
| 2026-11-17 02:14:45 | 00_validate | PASS | samples: 624 genes: 15000 |
| 2026-11-16 18:14:55 | 01_bulk_qc | PASS | normalization: quantile-normalized log2 ; outliers: 14 of 624 ; batch: corrected (limma removeBatchEffect; 2 levels) |
| 2026-11-16 18:15:30 | 02_signatures | PASS | source: MSigDB Hallmark (msigdbr; ensembl_gene) ; sets: 50 ; deconvolution: No cell-type signature matrix supplied; signature scores are marker scores, not absolute fractions |
| 2026-11-16 18:15:56 | 03_ecotypes | PASS | K = 2 ; bootstrap ARI = 0.975 ; CV accuracy = 0.979 |
| 2026-11-16 18:15:56 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-16 18:15:57 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-16 18:15:58 | 05_report | PASS | final report written |
| 2026-11-16 18:22:54 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-16 18:22:54 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-16 18:26:52 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-16 18:27:27 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-16 18:29:55 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-16 18:30:29 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-17 00:24:20 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-17 00:24:59 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-17 00:25:16 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-17 00:25:56 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
| 2026-11-17 00:35:47 | 02_signatures | PASS | source: MSigDB Hallmark (msigdbr; ensembl_gene) ; sets: 50 ; deconvolution: NNLS fractions from user-supplied cell-type signatures |
| 2026-11-17 00:38:13 | 02_signatures | PASS | source: MSigDB Hallmark (msigdbr; ensembl_gene) ; sets: 50 ; deconvolution: NNLS fractions from user-supplied cell-type signatures |
| 2026-11-17 00:38:53 | 03_ecotypes | PASS | K = 5 ; bootstrap ARI = 0.885 ; CV accuracy = 0.553 |
| 2026-11-17 00:39:06 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-17 00:39:46 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.769 ; LRT p: 0.198 ; PH min p: 0.247 |
| 2026-11-17 00:41:24 | 05_report | PASS | final report written |
| 2026-11-17 00:42:36 | 03_ecotypes | PASS | K = 2 ; bootstrap ARI = 0.976 ; CV accuracy = 0.986 |
| 2026-11-17 00:42:37 | 04_survival | COMPLETE_CASE | excluded 53 rows with missing model covariates |
| 2026-11-17 00:43:18 | 04_survival | PASS | endpoint: os ; n: 476 ; events: 79 ; C-index: 0.772 ; LRT p: 0.107 ; PH min p: 0.192 |
