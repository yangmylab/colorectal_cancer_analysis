# Methylation layer — TCGA-COAD/READ 450K

**Date:** 2026-11-17
**Data source:** GDC open-access methylation beta values (SeSAMe "Methylation Beta Value",
Illumina Human Methylation 450), projects TCGA-COAD and TCGA-READ, Data Release 46.0.
**Question:** does the frozen K=2 transcriptomic ecotype reproduce in DNA methylation?

## 1. Inputs

| Item | Value |
|---|---|
| Beta-value files downloaded | 458 (410 primary tumours, 6.0 GB) |
| Probes with ≥90 % completeness | 402,779 |
| Primary tumours in the matrix | 393 |
| Promoter probes mapped to a gene | 176,236 → **19,269 genes** |
| Samples with an RNA ecotype label | 388 |

Probes were mapped with the 450K annotation (TSS1500/TSS200/1stExon/5'UTR); a probe annotated
to several genes is assigned to the first named gene. Gene selection is therefore
**annotation-based, not variance-based**, so the differential-methylation test is not circular.

## 2. Results

| Metric | Value |
|---|---|
| CIMP panel genes available | 8 / 8 |
| CIMP-high / CIMP-low samples | 313 / 80 |
| CIMP-high fraction, E1 | 74.5 % |
| CIMP-high fraction, E2 | 84.3 % |
| Mean global beta, E1 vs E2 | 0.4756 vs 0.4808 (Wilcoxon p = 0.113) |
| Differentially methylated genes (FDR<0.05, \|logFC\|>0.5) | **5 of 19,269** |
| ARI, RNA ecotype vs unsupervised methylation cluster (K=2) | **0.006** |
| Methylation cluster and OS (n = 336, 84 events) | HR 0.74 (0.38–1.48), p = 0.398; C-index 0.511 → 0.738 with covariates |

Top differential genes (E2 vs E1, promoter): LMO7, C1orf226, TMIGD1, IQCK, TM4SF4 — all with
\|logFC\| < 0.4 on the M-value scale, i.e. small effects.

## 3. Interpretation

- **The K=2 ecotype is a transcriptional axis that does not replicate at the DNA-methylation level.**
  Unsupervised methylation clustering gives a partition essentially unrelated to the RNA ecotypes
  (ARI 0.006), and only 5 of 19,269 promoter genes differ between ecotypes after FDR control.
- Global methylation and CIMP status differ only modestly (CIMP-high 74.5 % in E1 vs 84.3 % in E2;
  global beta p = 0.113), so the EMT/angiogenic programme of E2 is not driven by a distinct
  methylation class.
- The methylation clusters themselves carry no prognostic signal (HR 0.74, p = 0.398).

## 4. Limitations

- 450K array covers promoters and CpGs, not the full methylome; EMT/angiogenesis programmes can be
  regulated without promoter methylation changes.
- Only 393 of 624 discovery tumours had usable methylation data; the smaller sample limits power for
  small effects.
- Sex-chromosome and polymorphic probes were not filtered beyond the completeness threshold.

## 5. Figures

- `plot_20` — CIMP score by ecotype
- `plot_21` — global methylation by ecotype
- `plot_22` — top differentially methylated genes (heatmap)
- `plot_23` — KM by methylation cluster
- `plot_24` — RNA ecotype vs methylation cluster

## 6. Files

- `03_Raw_Data/TCGA/methylation/` (raw beta files, matrix, gene-level matrix)
- `06_Analysis/08_Methylation/methylation_summary.tsv`
- `06_Analysis/08_Methylation/differential_methylation_ecotype.tsv`
- `06_Analysis/08_Methylation/methylation_survival.tsv`
- `06_Analysis/2_Python/1_PythonScripts/16_fetch_tcga_methylation.py`, `19_build_gene_methylation.py`
- `06_Analysis/1_R/1_Rscripts/19_methylation_analysis.R`

Signed: Reasonix agent · 2026-11-17
