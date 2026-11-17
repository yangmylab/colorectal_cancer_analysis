# Proteomics validation — CPTAC-COAD

**Date:** 2026-11-17
**Data source:** CPTAC prospective colon cancer (Vasaikar et al., *Cell* 2019), study PDC000116,
downloaded as processed tables from the cBioPortal data hub (`coad_cptac_2019`).
**Question:** do the frozen K=2 transcriptomic ecotypes reproduce at the protein level?

## 1. Inputs

| Layer | Genes / features | Samples |
|---|---|---|
| Protein (log2 TMT, cBioPortal `data_protein_quantification.txt`) | 8,067 | 97 |
| mRNA (RSEM, `data_mrna_seq_v2_rsem.txt`) | 13,481 | 107 |
| Both layers | — | 96 |

## 2. Method (pre-specified, classifier frozen)

1. Hallmark gene sets (MSigDB H, 50 sets) were scored per sample on mRNA and, independently, on protein:
   per-gene z-score across samples, then the mean over the measured members of the set (≥5 members).
2. The **frozen** K=2 nearest-centroid classifier from the TCGA discovery analysis was applied to both
   score matrices (no refitting, no sample-level labels reused).
3. mRNA-assigned and protein-assigned ecotypes were compared (agreement, adjusted Rand index).
4. OS was taken from the CPTAC clinical table (`OS_MONTHS`, `OS_STATUS`).

## 3. Results

| Metric | Value |
|---|---|
| Samples with both layers | 96 |
| Protein-assigned E1 / E2 | 49 / 47 |
| mRNA vs protein ecotype agreement | **85.4 %** |
| Adjusted Rand index (mRNA vs protein) | **0.496** |
| Median programme correlation (mRNA vs protein) | **r = 0.41** |
| Top concordant programmes | IFN-α response r=0.86; IFN-γ response r=0.80; EMT r=0.71; angiogenesis r=0.68 |

Survival in CPTAC-COAD is **not informative**: only **8 deaths** among 96 patients with protein data.

| Analysis | n | events | C-index | HR (E2 vs E1) | 95 % CI | p |
|---|---|---|---|---|---|---|
| CPTAC-COAD protein ecotype, OS | 96 | 8 | 0.596 | 2.44 | 0.47–12.59 | 0.286 |

## 4. Interpretation

- The **immune (IFN) and stromal (EMT/angiogenesis)** programmes that define E2 are the most
  reproducible across RNA and protein, and most samples receive the same ecotype from both layers.
- The concordance is high but imperfect (ARI 0.50): post-transcriptional regulation and the much smaller
  protein gene space (8,067 vs 13,481) both attenuate agreement.
- The protein-level survival test is **underpowered** (8 events) and is reported as descriptive only.
  It is directionally consistent with discovery (HR > 1 for E2) but cannot confirm or refute the
  transcriptomic null result.

## 5. Limitations

- CPTAC and TCGA are different patient populations with different treatment eras.
- Protein programme scores are means over measured members only; coverage differs per set.
- No phosphosite-level inference was attempted here (the phosphoproteome table is downloaded and held
  for a follow-up kinase-activity analysis).

## 6. Figures

- `plot_16` — protein Hallmark programme heatmap by protein-assigned ecotype
- `plot_17` — programme-level mRNA vs protein concordance
- `plot_18` — KM by protein-assigned ecotype (descriptive; 8 events)
- `plot_19` — ecotype assignment agreement, mRNA vs protein

## 7. Files

- `07_Validation/proteomics_validation.tsv`
- `07_Validation/proteomics_validation_metrics.tsv`
- `07_Validation/proteomics_programme_correlations.tsv`
- `03_Raw_Data/CPTAC_COAD/`
- `06_Analysis/1_R/1_Rscripts/17_proteomics_validation.R`

Signed: Reasonix agent · 2026-11-17
