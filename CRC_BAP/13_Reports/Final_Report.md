# Final Report — CRC Track A

Generated: 2026-11-17 (updated after the proteome, methylation, spatial, phosphoproteome and DFS layers)

## 1. Research question and hypothesis
See `00_Project_Design/Analysis_Plan.md` §1–2.

## 2. Datasets and provenance
See `01_Dataset_Discovery/Dataset_Inventory.csv` and `Dataset_Selection_Report.md`.

## 3. Sample annotation
Samples: 624 ; donors: 624
G1 annotation gate: see `02_Metadata/Metadata_Audit_Report.md`.

## 4. Design and statistical plan
See `00_Project_Design/Analysis_Plan.md` §4.

## 5. Quality control
Outliers flagged: 14 of 624
See `05_QC/QC_Report.md`.

## 6. Methods
Scripts: `06_Analysis/1_R/1_Rscripts/` (`21_phospho_analysis.R` = phosphoproteome, `22a_/22b_` = GSE14333 DFS); logs: `Logs/`; environment: `Environment/`.

## 7. Main findings
The graded ledger is authoritative: `13_Reports/Findings_Ledger.csv` (F01–F10).

| ID | Claim (abbreviated) | Grade |
|---|---|---|
| F01 | No incremental association of the K=2 ecotype with OS after adjustment (C-index 0.772, LRT p = 0.107) | not_supported |
| F02 | K=4 ecotype associated with OS in discovery — post hoc | E1 |
| F03 | K=2 association did not replicate in GSE39582 / GSE17536 (meta HR 1.25, p = 0.081) | not_supported |
| F04 | K=4 association did not replicate | not_supported |
| F05 | Single-cell reference resolves six lineages; TCGA deconvolution is descriptive only | descriptive |
| F06 | Ecotypes largely reproducible at the protein level (83.3 % agreement, ARI 0.439) | E3 |
| F07 | Ecotype does not reproduce at the DNA-methylation level (ARI 0.006) | not_supported |
| F08 | Bulk classifier does not transfer to Visium spots; spatial data descriptive only | not_supported |
| F09 | Ecotypes reproduce at the phosphoproteome level (87.5 % agreement, ARI 0.558, permutation p < 0.001) | E3 |
| F10 | Ecotypes associated with DFS in GSE14333 (adjusted HR 2.37, p = 0.004) — secondary endpoint | E1 |
| F11 | The disease-free association replicates in independent GSE39582 (RFS stage I-III, adjusted HR 1.50, p = 0.018) | E2 |

## 8. Biological interpretation
The stable axis is a two-way split between a Hedgehog / pancreas-beta / MYC-V2 programme (E1) and an
EMT / angiogenesis / IL6–JAK–STAT3 programme (E2). It survives a change of assay twice — at the protein and
at the phosphoproteome level — and fails to survive one: promoter methylation carries essentially no
information about the partition. That pattern (transcriptome + proteome + phosphoproteome agree, methylation
and bulk-to-spot projection do not) is the result; it is explicitly **not** a claim that ecotypes stratify
overall survival, which is null in discovery and in external validation.

## 9. Independent validation
- **OS (primary endpoint):** GSE39582 C-index 0.609 (p = 0.333), GSE17536 C-index 0.765 (p = 0.079);
  fixed-effect meta HR 1.25 (0.97–1.61), p = 0.081. Direction consistent, significance not reached.
- **DFS (secondary endpoint, GSE14333, n = 226, 50 events):** adjusted HR 2.37 (1.27–4.40), LRT p = 0.0043;
  stage II–III subset HR 2.15; 500-resample bootstrap 100 % of HRs above 1. Graded E1, kept strictly
  separate from the OS analysis — DFS and OS are never pooled.
- **RFS replication (GSE39582, independent cohort, stage I–III, n = 483, 137 events):** adjusted HR
  1.50 (1.07–2.10), p = 0.018; all stages HR 1.55 (1.12–2.15), p = 0.008; bootstrap 99 % of resamples
  above 1. GSE17536 RFS stage I–III HR 2.51 (1.14–5.53), p = 0.022. Fixed-effect summary with GSE14333
  HR 1.67 (1.24–2.24), p = 0.0007, no heterogeneity (p = 0.21) — graded **E2** (F11).
- **K=4 partition** did not replicate in either OS cohort (p = 0.090 / 0.200) and remains post hoc (E1).

## 10. Robustness
Outlier-excluded discovery model C-index 0.773 (LRT p = 0.151); K-sensitivity table in `08_Robustness/`;
classifier self-checks (bulk ARI 0.042 on spatial transfer) are reported rather than hidden; the phospho
concordance is tested against a 1000-permutation null (95th percentile ARI 0.033 vs observed 0.558).

## 11. Limitations
Retrospective public data; computational-only; no wet-lab validation.
Cross-platform layers come from small cohorts (CPTAC-COAD: 96 tumours, 8 deaths; GSE14333: 226 tumours,
50 DFS events, stage I–III only). The phospho matrix is 59.8 % missing before filtering. Kinase-activity
inference was not performed because no substrate library is available offline. The ecotype classifier
cannot be projected onto spatial spots, so the spatial layer stays descriptive.

## 12. Reproducibility package
Scripts, configs, decision log and figure manifest are in this project, now under git version control.
A Zenodo DOI is pending. Reproduce the new layers in order: `21_phospho_analysis.R` → `22a_build_gse14333_dfs.R`
→ `22b_dfs_validation.R` → `14_collect_plots.R`.

## 13. Proposed follow-up experiments
IHC/CODEX validation of ecotype markers; organoid perturbation of prioritized targets; kinase-substrate
activity scoring once an annotation library is installed; a spot-level (deconvolution-trained) spatial
classifier; an independent DFS cohort (e.g. GSE38832) to test F10.

## 14. Reviewer scorecard
| Class | Issue | Evidence | Corrective action | Status |
|---|---|---|---|---|
| Critical | None yet | — | — | — |
| Major | No wet-lab validation | design | state limitation | open |
| Major | Prognostic claim unsupported | F01, F03 | report null OS honestly; DFS only E1 | addressed |
| Minor | Public-data heterogeneity | design | sensitivity analyses | open |
| Minor | Kinase activity not inferable offline | F09 | documented, not claimed | open |

**Verdict:** findings graded and reported within their evidence; manuscript drafting is the next gate (G3).

## 15. Evidence-grade distribution
- E3 (cross-platform convergence): F06 (protein), F09 (phosphoproteome)
- **E2 (replicated association): F11 (disease-free/RFS, GSE39582)**
- E1 (exploratory / post hoc): F02 (K=4 discovery), F10 (DFS discovery, secondary endpoint)
- descriptive: F05 (single-cell deconvolution)
- not_supported: F01, F03, F04 (OS association), F07 (methylation), F08 (spatial transfer)
- The dissociation is the message: the ecotype axis reproduces across assays and predicts disease-free
  survival, but no OS association replicates in any cohort tested.
