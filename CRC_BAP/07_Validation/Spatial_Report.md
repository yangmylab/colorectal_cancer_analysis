# Spatial layer — Visium CRC sections

**Date:** 2026-11-17
**Datasets:** archives supplied by Prof. Meng, both verified (`tar -tf` member count + `gzip -t` on members):

| Dataset | Sections | In-tissue spots | Genes per matrix |
|---|---|---|---|
| GSE283052 (12 primary CRC tumours) | 12 | 26,766 | 36,601 |
| GSE267401 (2 patients, primary/metastasis pairs) | 4 | 17,052 | 37,082 |
| **total** | **16** | **43,818** | 13,562 mapped to the bulk reference |

`GSE284061` remains unusable: its GEO archive is corrupt at the source — members fail `gzip -t` and the
md5 is identical across four independent downloads (curl, wget, aria2c, and the user's own copy).

## 1. Method

Spots were scored in the **bulk reference space**: each Hallmark programme score is the mean, over its
member genes, of gene-wise z-scores whose mean/SD come from the TCGA discovery log-expression matrix
(15,000 genes × 624 tumours). The same statistic was computed for the 624 bulk tumours, so spot and
bulk scores share a scale and are **not** standardised within a section.

## 2. The bulk ecotype classifier does not transfer to spots

Two normalisations were tried before any spatial claim was made; both fail, and both failures are
recorded (`06_Analysis/10_Spatial/spatial_transfer_diagnostics.tsv`) rather than hidden:

| Test | Result |
|---|---|
| Self-check: does the simplified statistic reproduce the GSVA-based bulk labels (624 tumours)? | agreement 0.604, **ARI 0.042** — only the majority class is reproduced, so the statistic is not a faithful transfer of the classifier |
| Spot calls with per-section normalisation | forces ≈50/50 E1/E2 in every section by construction — **invalid**, discarded |
| Spot calls with bulk-referenced normalisation | 100 % of spots, and all 16 sections, land on one side of the decision boundary |
| Distributional position | spot E2-programme scores sit **≈5 SD below** the bulk distribution (bulk median ≈ 0.0, spots ≈ −5.0) |

This is the expected consequence of resolution: a 55 µm spot holds a few thousand transcripts from a
mixture of cell types, with no compartment comparable to a whole tumour. A spot-level ecotype
classifier must be **trained** on single-cell or spatial data; bulk centroids cannot be projected onto
spots.

## 3. What the spatial data do show (descriptive only)

- The E2-defining programmes (EMT + angiogenesis + IL6–JAK–STAT3) vary **spatially within every
  section**, in coherent patches rather than uniform fields (`plot_25`). The metastasis section CTC21M
  has an extended E2-high region plus a distinct E2-low region; primary sections P02 and P11 are
  comparatively uniform.
- Between-section variation is smaller than within-section variation (`plot_26`): at spot resolution the
  programme axis is organised as **intra-tumour niches**, not as whole-tumour classes.
- No spatial claim about ecotypes is made. Per the pre-specified rule (PROTOCOL_TrackA §7), spatial data
  are used descriptively when they do not validate the bulk taxonomy — which is this case.

## 4. Figures

- `plot_25` — spatial map of the E2-defining programmes per spot (8 largest sections)
- `plot_26` — within- vs between-section variation of the E2 programme
- `plot_27` — covariate-shift diagnostic: spot vs bulk E1/E2 distributions

## 5. Files

- `06_Analysis/10_Spatial/spatial_section_programme_summary.tsv`
- `06_Analysis/10_Spatial/spatial_transfer_diagnostics.tsv`
- `06_Analysis/10_Spatial/spatial_scoring_selfcheck.tsv`
- `06_Analysis/10_Spatial/spatial_scores.rds` (per-spot programme scores, centroid correlations)
- `06_Analysis/1_R/1_Rscripts/20_spatial_analysis.R`
- `03_Raw_Data/GSE283052/`, `03_Raw_Data/GSE267401/`

## 6. Next step for a genuine spatial claim

Train the spot-level classifier on the single-cell reference (GSE132465) — cell2location or RCTD
deconvolution followed by programme scoring — or use a spatial cohort with matched single-cell data
(CNSA `CNP0003321`, the spatial arm of GSE225857) and test whether the spatial ecotype predicts outcome.
That is a new analysis, not a re-use of the bulk centroids.

Signed: Reasonix agent · 2026-11-17
