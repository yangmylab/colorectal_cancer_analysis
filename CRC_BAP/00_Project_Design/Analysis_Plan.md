# Analysis Plan — CRC cross-platform ecotype taxonomy (Track A)

> Governed by the `/bioanalyst` BAP. Sections 1–2 gate **G0**; sections 3–6 gate **G2**.
> G1 is a hard sample-annotation gate. This project is currently **G1 BLOCKED** because no CRC
> data are present on disk; no main analysis has been run.

---

## 1. Research question

- **Biological question:** Can colorectal cancer (CRC) be decomposed into a small number of reproducible, spatially organised tumour–immune ecotypes that refine prognosis and therapy-relevant associations beyond stage, MSI and existing bulk classifications (CMS/CRIS/PDS)?
- **System / model:** human colorectal adenocarcinoma (colon and rectum), primary tumours, with available public transcriptome, proteome, methylation, single-cell and spatial data.
- **Exposure / intervention / independent variable:** ecotype assignment derived from single-cell-anchored transcriptional, epigenetic and proteomic programmes.
- **Primary outcome:** cancer-specific survival (CSS); overall survival (OS) where CSS is unavailable.
- **Secondary outcomes (exploratory):** OS, RFS, progression/metastasis association, anti-EGFR/chemotherapy/ICI response associations, morphology–ecotype prediction, circulating ecotype detection.
- **Mode:** hypothesis-driven for the ecotype–prognosis question; exploratory for therapy, morphology, circulating and microbiome extensions.
- **Deliverable:** reproducible classifier, benchmark, figures/tables, manuscript draft and portal specification.
- **What would falsify the hypothesis:** ecotypes fail to reach bootstrap stability (ARI ≤ 0.70) or fail to replicate in ≥ 2 independent cohorts; or the ecotype term adds no discrimination/calibration/decision-curve benefit over the base model in any external cohort.

## 2. Hypothesis

- **Primary hypothesis:** adding ecotype assignment to a base model (age, sex, stage, MSI, sidedness, CMS/CRIS where available) improves CSS discrimination, calibration or decision-curve benefit in an independent cohort.
- **Alternative hypotheses worth testing:** (i) ecotypes are protein-level states that replicate in PXD079707/CPTAC; (ii) ecotypes map to recurrent spatial niches; (iii) ecotypes are detectable in cfDNA/EV signals; (iv) ecotypes associate with therapy response; (v) a simpler factor-continuum explains the data better than discrete ecotypes.
- **Prior evidence:** CMS/CRIS/PDS coexist and remain platform-dependent (PMID 38319359; PMID 41607760); a 4.27M-cell CRC immune atlas defines immune phenotypes (PMID 41529665); the 1,000-CRC proteogenomic resource reports that current classifications are incomplete (PRIDE PXD079707); progression and CRLM proteogenomic cohorts exist (PMID 38861363; PMID 41195591).
- **Expected effect size (for planning):** HR ≈ 1.4–1.6 per adverse ecotype vs reference; unknown until data are analysed. Planning assumes ≥150 CSS events in discovery.

### Datasets

| Role | Accession | Repository | Publication (DOI/PMID) | n | Groups | Notes |
|---|---|---|---|---|---|---|
| Discovery (bulk multi-omic) | TCGA-COAD / TCGA-READ | GDC | PMID 22810696 | 461 + 172 cases; 701 RNA-seq files; 2,202 methylation files | Tumour vs normal; stage/MSI | Primary discovery; open tiers only |
| Proteomic anchor | PXD079707 | PRIDE | not indexed as of 2026-11-17 | 1,000 treatment-naive CRCs | Tumour; clinical follow-up | DIA proteomics; WSIs; CC0; verify metadata |
| Proteomic validation | CPTAC colon (PDC000116, confirm) | PDC | PMID 31031003 | 110 tumours | Tumour vs normal | Proteome/phosphoproteome |
| Progression | PXD052391 / IPX0002188000 | ProteomeXchange/iProx | PMID 38861363 | 148 patients; 435 samples | Nontumour → IEN → infiltration → advanced | WES restricted (HRA007452) |
| Metastasis | IPX0007391001 / PRJCA020890 | iProx/GSA | PMID 41195591 | 34 patients; 102 samples | Primary/normal/liver metastasis | Matched CRLM |
| Validation (bulk) | GSE39582 | GEO | NA (verify) | 585 | Tumour; survival | Core external validation |
| Validation (bulk) | GSE17536 / GSE14333 / GSE38832 | GEO | NA (verify) | 177 / 290 / 122 | Tumour; survival | Independent cohorts; check overlap |
| Reference (single-cell) | GSE132465 / GSE144735 / GSE178341 / GSE201348 / GSE312983 | GEO | PMID 32451460 (GSE132465); others verify | 23 / 6 / 181 / 72 / 189 | Tumour, normal, polyp, age | Cell-state and chromatin reference |
| Spatial validation | GSE225857 / GSE283052 / GSE320042 / GSE284989 | GEO | NA (verify) | 8 / 12 / 38 / 22 | Primary/metastasis | Multi-platform validation |

Full provenance: `01_Dataset_Discovery/Dataset_Inventory.csv`.
Rationale for selection and rejections: `01_Dataset_Discovery/Dataset_Selection_Report.md`.

---

## 3. Sample annotation status

- Annotation table: `02_Metadata/Sample_Annotation.csv` — **624** samples, **624** donors.
- Ambiguities: `02_Metadata/Ambiguous_Samples.csv` — **0 open** rows.
- Audit report: `02_Metadata/Metadata_Audit_Report.md`.
- **G1 status:** **PASS (2026-11-17)** — TCGA-COAD/READ downloaded from GDC; 647 primary-tumour files collapsed to 624 patients.

| Group | n samples | n donors | biological replicates | technical replicates |
|---|---|---|---|---|
| Primary tumour | 624 | 624 | 624 | 0 |
| Solid tissue normal | 51 files | — | — | — (not used) |
| Metastasis / recurrent | 3 files | — | — | — (not used) |

**G1 resolution:** G1 passed. Stage is missing for 118 patients, OS time for 66, and MSI/purity are unavailable. Normal samples are not used as controls. The 23 duplicate/alternate primary-tumour samples were collapsed by highest total TPM per patient.

---

## 4. Statistical analysis plan

### 4.1 Design

- **Design type:** cross-sectional multi-cohort integrative analysis (retrospective public data).
- **Unit of inference:** biological sample = patient × tumour/normal; cells and spatial spots are never treated as independent patients.
- **Independence:** one primary tumour per patient in discovery; paired normal/tumour and matched primary/metastasis are encoded explicitly; technical replicates collapsed by pre-specified rules.
- **Discovery vs validation split:** TCGA-COAD/READ (+ PXD079707 for protein) is discovery; GSE39582, GSE17536, GSE14333, GSE38832, GSE87211, GSE103479 and CPTAC are external validation; the classifier and cutoffs are frozen before validation is touched.
- **K selection:** among K with bootstrap ARI ≥ 0.70 and minimum group size ≥ 5%, prefer the highest mean ARI; ties are broken by silhouette, then by the smaller K.

### 4.2 Model specification

```
Response:            CSS time-to-event (primary); OS/RFS (secondary)
Model:               ~ ecotype + age + sex + stage + msi_status + sidedness (+ CMS where available)
Random effects:      none for discovery; dataset as random effect in cross-cohort meta-analysis
Primary contrast:    adverse ecotype(s) vs reference ecotype
Test:                Wald/LRT for Cox; C-index difference with bootstrap CI
Multiplicity:        BH FDR 0.05 within families; primary endpoint alpha 0.05
Assumptions:         proportional hazards; independence; adequate events per variable
Known violations:    stage/MSI may be partially confounded with ecotype; purity and batch may differ by cohort
```

### 4.3 Confounders

| Variable | Available? | Complete? | Included in model? | Rationale |
|---|---|---|---|---|
| age | yes | to verify | yes | prognostic |
| sex | yes | to verify | yes | prognostic |
| stage | yes | to verify | yes | dominant prognostic factor |
| MSI/MMR | yes | to verify | yes | immune biology and prognosis |
| sidedness | yes | to verify | yes | biology and outcome |
| purity | estimated | to verify | yes (sensitivity) | stromal/immune confounding |
| batch/cohort | yes | to verify | yes (correction/stratum) | technical confounding |
| treatment | partial | often missing | no in primary; yes in therapy extension | indication confounding |

- **Fully confounded with group (effect not estimable):** to be determined at G2 after metadata audit.
- **Deliberately excluded and why:** post-treatment variables and mediators are excluded from the primary prognostic model.

### 4.4 Feasibility verdict

With ~600 discovery tumours and ≥150 CSS events, the design can estimate a 3–5 group ecotype model with 5–7 covariates and moderate effects (HR ≈ 1.4–1.6). It cannot support causal claims, treatment-effect claims, or cell-level inference. Proteomic and spatial extensions are validation layers, not discovery. If events are fewer than 50, the primary analysis will be reduced to Kaplan–Meier/log-rank plus univariable Cox and explicitly labelled underpowered.

### 4.5 Analyses explicitly NOT run

- Cell-level differential expression or cell-level survival (pseudo-replication).
- “Optimal” survival cutoffs fitted on the evaluation cohort.
- Causal claims from observational data without MR/colocalisation.
- Drug efficacy claims from CMap/GDSC associations.
- Treating inferred ligand–receptor signalling as demonstrated communication.

### 4.6 Pre-specified exclusions and QC thresholds

- Exclude non-primary tumours from discovery (keep metastasis for the extension).
- Exclude samples with missing stage or outcome from the primary survival analysis; report counts.
- Exclude samples with extreme purity or expression outliers under pre-specified rules (Mahalanobis distance, PCA distance, detection rate); every exclusion gets an ID and reason.
- Collapse technical replicates by mean (log-scale) or sum (counts), documented.
- Do not exclude samples based on ecotype or outcome.

---

## 5. Analysis pipeline (planned steps)

| # | Step | Script (planned) | Inputs | Outputs | Method / key parameters |
|---|---|---|---|---|---|
| 00 | ingest, checksums, G1 annotation gate | `06_Analysis/1_R/1_Rscripts/00_validate_inputs.R` | `03_Raw_Data/` | `02_Metadata/Sample_Annotation.csv`, `Logs/` | schema checks; `bap.py gate g1` |
| 01 | bulk QC and normalisation | `01_bulk_qc.R` | expression + clinical | `05_QC/`, `04_Processed_Data/` | voom/limma or log2; quantile; ComBat; PCA |
| 02 | signature scoring and deconvolution | `02_signatures_deconvolution.R` | processed expression + GMT | `06_Analysis/05_Immune_Analysis/` | GSVA ssGSEA; NNLS where a signature matrix exists |
| 03 | multi-view integration and ecotypes | `03_snf_ecotypes.R` | RNA/methylation/protein views | `06_Analysis/03_Clustering/` | SNF; consensus k-means; bootstrap ARI; nearest-centroid classifier |
| 04 | survival and external validation | `04_survival_validation.R` | labels + clinical + validation cohorts | `06_Analysis/06_Other_Analysis/`, `07_Validation/` | Cox/glmnet; C-index; calibration; DCA; frozen classifier |
| 05 | figures and report | `05_figures_report.R` | all results | `09_Figures/`, `10_Tables/`, `13_Reports/` | PDF + JPEG; tables; decision log |
| 99 | synthetic dry run | `99_synthetic_dryrun.R` | synthetic data | `99_Synthetic_DryRun/` | validates code only; not biological evidence |

Parallelisation: 8 workers by default (16 physical cores; one level only; BLAS pinned to 1 per worker), configurable in `config/pipeline_config.yaml`.

---

## 6. Reporting and validation plan

- **Reporting guideline:** TRIPOD+AI for prediction models; REMARK for prognostic markers; STROBE-like reporting for observational data.
- **Validation strategy:** external cohorts with pre-specified direction, effect correlation and C-index retention; classifier frozen before validation; report calibration and DCA, not only AUC.
- **Robustness menu:** purity adjustment, batch-correction method, platform restriction, cohort exclusion, MSI/stage subgroup analysis, alternative K, alternative integration method.
- **Figure plan:** Figure 1 study design/data matrix; Figure 2 single-cell reference; Figure 3 integration/ecotypes/stability; Figure 4 biology; Figure 5 spatial; Figure 6 clinical/therapy; Figure 7 drug/resource; Figures 8–10 extensions.
- **Evidence-grade expectation:** E2 for replicated ecotype–prognosis associations; E3 only if transcriptome + proteome + spatial/clinical converge; E4 is not achievable without experimental validation.

---

## Sign-off

| Gate | Status | Date | Note |
|---|---|---|---|
| G0 intake | PASS | 2026-11-17 | Question, hypothesis and outcomes specified |
| G1 annotation | PASS | 2026-11-17 | 624 tumours / 624 donors; 0 ambiguities |
| G2 design & statistics | PASS (design) | 2026-11-17 | Model, confounders, feasibility and validation specified; execution awaits G1 |
| G3 final review | PENDING | — | Requires completed analysis |
