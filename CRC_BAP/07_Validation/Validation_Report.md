# Validation report — CRC Track A

**Status:** two independent cohorts analysed (GSE39582, GSE17536). Neither the primary K=2 association nor the exploratory K=4 partition replicated at p<0.05.

## 1. Discovery result

- K=2 ecotype, n=476 complete cases, 79 OS events.
- Base C-index 0.763; full C-index 0.772; LRT p=0.107.
- Discovery ecotype E2 vs E1 HR 1.45 (95% CI 0.92–2.28; p=0.106).
- Primary K=2 ecotype–OS association: **not supported**.

## 2. External validation — K=2 classifier

| Cohort | n | events | C-index | ecotype HR (E2 vs E1) | p | direction |
|---|---|---|---|---|---|---|
| GSE39582 | 573 | 190 | 0.609 | 1.16 | 0.333 | consistent |
| GSE17536 | 177 | 73 | 0.765 | 1.53 | 0.0785 | consistent |

**Two-cohort fixed-effect meta-analysis:** HR **1.25** (95% CI 0.97–1.61), **p=0.081**; heterogeneity Q=0.98, tau²=0.

- Direction is consistent in both cohorts, but neither is significant at p<0.05.
- **Conclusion:** the primary K=2 association is **not replicated** at the pre-specified significance level. No E2 claim.

## 3. External validation — K=4 exploratory partition

| Cohort | n | events | C-index | LRT p |
|---|---|---|---|---|
| TCGA-COAD/READ (discovery) | 476 | 79 | 0.785 | 0.024 |
| GSE39582 | 573 | 190 | 0.614 | 0.090 |
| GSE17536 | 177 | 73 | 0.772 | 0.200 |

**Conclusion:** K=4 does **not** replicate significantly. It remains E1 exploratory.

## 4. Cohort characteristics

- GSE39582: 585 tumours, 23,963 Ensembl-mapped genes, 194 OS events (190 in model), 77 MSI / 459 MSS, stage missing 8.
- GSE17536: 177 tumours, 23,963 Ensembl-mapped genes, 73 OS events, 55 DSS events, 36 DFS events, stage complete.

## 5. Pre-specified validation rules

- Freeze the classifier and cutoffs before touching validation data.
- Score Hallmark programmes with the same Ensembl mapping.
- Fit Cox with ecotype + available covariates; report C-index, HR, calibration, DCA.
- E2 requires same direction **and** p<0.05 in an independent cohort.
- Do not refit the classifier on validation data.

## 6. GSE14333 (not used for OS validation)

- GSE14333 series matrix complete (81 MB; GPL570).
- It provides **DFS only** (`DFS_Time`, `DFS_Cens`), not OS.
- It was therefore **not** used for the primary OS validation. It can support a secondary DFS analysis in a follow-up.

## 7. Next action

- Optional: build GSE14333 for a secondary DFS analysis (do not relabel DFS as OS).
- If no cohort replicates the OS association, report the stable K=2 axis as a descriptive taxonomy without a survival claim.
- Do not present K=4 as confirmatory.

Signed: Reasonix agent · 2026-11-17
