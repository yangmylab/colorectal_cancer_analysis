# Track A protocol: cross-platform, single-cell-anchored CRC ecotype taxonomy

**Version:** 2.0 (2026-11-17; expanded after the database/literature scan)
**Status:** protocol for approval; no data downloaded, no software installed
**Chosen constraints:** public data only; computational validation only; target journals *Journal of Translational Medicine*, *Computational and Structural Biotechnology Journal*, *Briefings in Bioinformatics*
**Companion files:** `CRC_multiomics_strategy.md`, `DATASETS.md`, `DATASETS_EXTENDED.md`, `LITERATURE_GAPS_2024_2026.md`, `METHODS_TOOLBOX_2026.md`, `SOFTWARE.md`, `PILOT_PLAN.md`

---

## 1. Decision summary

| Item | Decision |
|---|---|
| Study type | Retrospective integrative analysis of public human CRC data |
| Central question | Are there reproducible, spatially organized tumor-immune ecotypes that refine CRC prognosis beyond stage/MSI/CMS/CRIS? |
| Discovery | TCGA-COAD/READ (multi-omic) + CPTAC colon (proteome/phosphoproteome) |
| Proteomic anchor (1,000) | PXD079707 (DIA proteomics, targeted genomics, WSIs, clinical) |
| Progression | PXD052391 / IPX0002188000 (148 patients, 435 samples) |
| Metastasis | IPX0007391001 / PRJCA020890 (34 patients, 102 samples) |
| Morphology | TCGA diagnostic slides + PXD079707 WSIs; CHIEF/MUSK features |
| Circulating | cfDNA methylation/5hmC and EV/RNA cohorts |
| Benchmark | CMS/CRIS/PDS stability + integration/deconvolution/niche/survival methods |
| Reference | Public scRNA/snRNA/scATAC atlases |
| Spatial validation | Public Visium/Stereo-seq/high-definition spatial cohorts |
| External validation | GEO bulk cohorts + CPTAC + cBioPortal (mutation/CNV context) |
| Primary endpoint | Cancer-specific survival (CSS); OS where CSS is unavailable |
| Primary comparison | Ecotype + base model vs. base model alone (stage, age, sex, MSI, sidedness, CMS/CRIS) |
| Validation capacity | Computational only; no wet-lab claims |
| Journal targets | *J Transl Med* (clinical translation), *CSBJ* (computational framework), *Briefings in Bioinformatics* (reproducible framework/resource) |

**Framing sentence for the paper.**
> We present a reproducible, cross-platform framework that converts single-cell programs into bulk/spatial ecotypes, quantifies the instability of existing CRC classifications, and identifies ecotypes with external prognostic and therapeutic relevance.

This framing is honest about public-data limitations and positions the work as a reproducibility/translation study rather than a mechanistic discovery.

---

## 2. Study question, population and endpoints

### 2.1 PICO
- **Population:** adults with primary colorectal adenocarcinoma (colon and rectum), treatment-naive, with bulk transcriptome and clinical follow-up.
- **Index:** ecotype assignment derived from multi-view programs.
- **Comparator:** existing classifications (CMS, CRIS, PDS) and clinical variables.
- **Outcomes:** CSS (primary), OS/RFS (secondary), therapy-response associations (exploratory).

### 2.2 Primary hypothesis
> Adding ecotype assignment to a base model (age, sex, stage, MSI, sidedness, CMS/CRIS) improves discrimination and calibration for CSS in an independent cohort.

### 2.3 Secondary hypotheses
1. Ecotypes are reproducible across platforms (RNA, methylation, proteome, phosphoproteome) and cohorts.
2. Ecotypes map to recurrent spatial niches (TLS, invasive front, hypoxic/immune-excluded, CAF-rich).
3. Ecotypes differ in candidate therapeutic vulnerabilities.
4. Ecotypes refine risk within stage II/III and within MSS tumors.

### 2.4 Inclusion / exclusion
**Include:**
- Primary tumor samples with RNA-seq (or array) and survival/time-to-event data.
- Stage I-IV if stage is recorded; treatment-naive primary tumors.
- Consistent sample-type coding (TCGA 01; GEO primary tumor annotations).

**Exclude:**
- Metastatic samples from the discovery cohort (analyzed separately for metastasis extension).
- Normals, cell lines, PDX, organoids.
- Duplicated patients/aliquots; samples with missing stage or follow-up.
- Low-quality samples passing pre-specified QC failures.

**Handling missing modalities:**
- Multi-view integration uses samples with ≥2 views; complete-case sensitivity analysis.
- Methylation/proteome subsets analyzed separately and integrated at the factor level.
- Do not impute entire omic views.

---

## 3. Project structure and data flow

```
colorectal_cancer/
  CRC_multiomics_strategy.md
  PROTOCOL_TrackA.md
  PILOT_PLAN.md
  DATASETS.md
  SOFTWARE.md
  analysis/
    README.md
    config/study_config.yaml
    metadata/sample_sheet_template.tsv
    scripts/                 # planned module entry points (see scripts/README.md)
    data/                    # user-downloaded data (not created by assistant)
    processed/
    reference/
    results/
    figures/
    manifests/
```

**Data flow.**
`raw -> QC -> harmonized matrices -> reference signatures -> deconvolution/program scores -> multi-view factors -> ecotypes -> spatial validation -> clinical models -> drug mapping -> figures/resource`

**Rules.**
- Never modify raw data; all derived objects are versioned.
- Every table/figure must be regenerable from scripts + manifest.
- Lock the analysis cohort before outcome modeling.
- Record all parameters and seeds in `study_config.yaml`.

---

## 4. Phase-by-phase plan

### P0. Setup, manifest and pre-registration (weeks 1-2)

**Tasks.**
1. Create the project skeleton and `study_config.yaml`.
2. Build the sample sheet: one row per sample, with patient ID, cohort, tissue, platform, modality, stage, MSI, sidedness, age, sex, treatment, outcomes, purity, batch.
3. Download Tier 1 datasets manually (see `DATASETS.md` and `DATASETS_EXTENDED.md`), including a PXD079707 subset, PXD052391 and IPX0007391001; record checksums.
4. Check patient overlap across cohorts.
5. Write and register the statistical analysis plan (OSF) before outcome analysis.
6. Freeze the cohort definition and primary endpoint.

**Outputs.**
- `metadata/manifest.tsv`, `metadata/sample_sheet.tsv`, `metadata/overlap_report.tsv`
- Registered SAP (PDF/DOI)
- `config/study_config.yaml` (locked)

**QC gate G0.**
- All discovery samples traceable to source files.
- No patient overlap between discovery and validation.
- Primary endpoint and model covariates frozen.

---

### P1. Bulk preprocessing and harmonization (weeks 2-6)

**P1.1 TCGA-COAD/READ.**
- RNA-seq: STAR counts from GDC; gene-level TPM; log2(TPM+1); filter low-expression genes; PCA/outlier audit.
- Methylation: IDAT -> beta (minfi); BMIQ; M-values; DMP/DMR (limma/bumphunter); CIMP classification.
- CNV: ASCAT/CNVkit; gene-level calls; burden scores.
- Mutation: MAF aggregation (maftools); driver annotation; MSI status.
- RPPA: normalize; map antibodies to genes; use as orthogonal protein view.
- Clinical: GDC clinical XML/API; harmonize stage (AJCC 7/8), vital status, follow-up.
- Slides: optional; keep metadata for future pathology correlation.

**P1.2 Proteogenomics.**
- **PXD079707:** process DIA proteomics (Spectronaut/DIA-NN), protein-level normalization, gene mapping; 1,000 treatment-naive CRCs; match to clinical/WSI metadata.
- **CPTAC colon:** proteome/phosphoproteome normalization; independent validation; kinase activity.
- **PXD052391 / IPX0002188000:** progression cohort (nontumor, intraepithelial neoplasia, infiltration, advanced); stage-aware normalization.
- **IPX0007391001 / PRJCA020890:** matched primary/normal/liver metastasis; paired analysis.
- **PXD019103:** orthogonal colon proteome validation.
- Confirm accessions, data levels, processed tables and licenses; document.

**P1.3 GEO validation cohorts.**
- Arrays: RMA/quantile normalization; re-annotate probes to current gene symbols; collapse duplicates.
- Counts: `ComBat-seq` for cohort/batch; `limma::removeBatchEffect` for arrays as sensitivity.
- Extract stage, MSI, sidedness, treatment, survival from series metadata and publications.
- Flag cohorts with missing/ambiguous covariates.

**P1.4 Purity and confounding.**
- ESTIMATE/ABSOLUTE purity; include purity as covariate.
- Register batch, cohort, platform, scan date (if available).
- Do not correct the outcome; correct only technical variation.

**Outputs.**
- `processed/rna/`, `processed/methylation/`, `processed/proteome/`, `processed/cnv/`, `processed/mutation/`
- `metadata/harmonized_clinical.tsv`
- QC report: sample counts, missingness, PCA, cross-platform correlations.

**QC gate G1.**
- ≥500 discovery tumors with RNA + clinical; ≥300 with proteome (PXD079707 subset + CPTAC).
- Protein-level ecotype program scores show pre-specified positive concordance with RNA-based scores.
- Progression and metastasis cohorts are processed and stage-annotated.
- Cross-platform RNA-protein correlation for shared genes reported (expect positive but modest).
- No residual batch separation on PCA after correction (or documented residual confounding).
- Outlier rate <5%, with reasons recorded.

---

### P2. Single-cell reference atlas (weeks 4-8)

**P2.1 Data and QC.**
- Ingest GSE132465, GSE144735, GSE178341, GSE166555, GSE201348, GSE231559, GSE312983, GSE285117.
- Standard QC: doublets, ambient RNA, mitochondrial/ribosomal thresholds; retain metadata (patient, tissue, platform).
- Harmonize cell-type labels to a common ontology.

**P2.2 Integration and annotation.**
- `scVI`/`scANVI` integration; Harmony as sensitivity.
- Canonical markers + `CellTypist`; manual curation.
- Malignant-cell inference: `inferCNV`, `CopyKAT`; separate malignant from normal epithelium.
- Compartment-specific subclustering: epithelial, fibroblast, myeloid, T/NK, B/plasma, endothelial, glial/neuronal.

**P2.3 Program discovery.**
- NMF/cNMF per compartment; select robust programs (bootstrap stability, ≥20 genes).
- Epithelial programs: stem/proliferative, regenerative, EMT, metabolic, secretory, antigen-presentation.
- Fibroblast programs: myofibroblast, inflammatory, antigen-presenting, matrix-remodeling.
- Myeloid programs: monocyte-derived, tissue-resident, immunosuppressive, antigen-presenting.
- Lymphoid programs: cytotoxic, exhausted, memory, Treg, B/plasma, TLS-associated.
- scATAC (GSE201348): ArchR/Signac; gene activity; regulons; link chromatin programs to RNA programs.

**P2.4 Signature construction.**
- Program gene lists (top weighted genes) and scores.
- TLS signature, invasive-front signature, hypoxia, immune exclusion.
- Export `reference/signatures/*.gmt` and `reference/programs/*.tsv`.

**QC gate G2.**
- Cell-type proportions consistent with known CRC biology.
- Marker expression validates annotations.
- Programs stable across atlases (Jaccard overlap of top genes >0.3 between independent atlases, or documented differences).
- No single dataset dominates integration (report dataset mixing).

---

### P3. Bulk deconvolution and program scoring (weeks 6-9)

**Tasks.**
- Bulk deconvolution: CIBERSORTx, EPIC, quanTIseq, MCP-counter, xCell; consensus rank scores.
- Purity adjustment: report raw and adjusted fractions.
- Program scoring: ssGSEA/GSVA, AUCell/UCell; z-score within cohort; batch-aware.
- Compare deconvolution estimates with pathology/purity where available.
- Build a "program matrix" for integration: samples x programs.

**Outputs.**
- `results/deconvolution/` (per method + consensus)
- `results/programs/` (scores per sample/cohort)
- Benchmark report: method concordance, correlation with purity.

**QC gate G3a.**
- Major compartment fractions concordant across methods (Spearman >0.7).
- Program scores reproducible across technical replicates/cohorts.
- No program score explained almost entirely by purity (>0.9 correlation flagged).

---

### P4. Multi-view integration and ecotype discovery (weeks 8-12)

**P4.1 Views.**
- RNA: top robust programs + selected variable genes.
- Methylation: CIMP status, DMP scores, selected probes.
- Proteome: protein modules (WGCNA/PCA).
- Phosphoproteome: kinase/pathway activity scores.
- CNV: burden and arm-level events.

**P4.2 Integration.**
- Primary: MOFA+ (15 factors; sample-level; missing views allowed).
- Sensitivity: iClusterBayes, SNF, DIABLO.
- Report variance explained per view/factor; factor stability (bootstrap).

**P4.3 Ecotype clustering.**
- Consensus clustering on selected factors (K=2-8; 1,000 bootstraps).
- Metrics: PAC, silhouette, bootstrap ARI, group size.
- Select K by stability + interpretability + minimum 5% per group + clinical separation.
- Label ecotypes by top programs (do not name before data).

**P4.4 Comparison with existing classifications.**
- Compute CMS (CMSclassifier/NeoCMS), CRIS, PDS where possible.
- Report ARI/NMI, reclassification tables, and survival of discordant cases.
- Quantify instability of CMS/CRIS across cohorts/platforms (a key result).

**P4.5 Locked classifier.**
- Train a parsimonious classifier on program scores (nearest shrunken centroid or multinomial elastic net).
- Nested cross-validation; export weights and thresholds.
- Freeze before external validation.

**Outputs.**
- `results/integration/mofa/`, `results/ecotypes/`
- Ecotype labels, factor loadings, stability plots, CMS/CRIS comparison.
- Locked classifier (`models/ecotype_classifier.*`).

**Decision gate G3.**
- Ecotypes stable (bootstrap ARI >0.70).
- ≥3 ecotypes, each ≥5% of samples, with distinct biology.
- Replication in ≥2 independent cohorts (ARI >0.5 with discovery).
- If G3 fails: adjust views/K; if still failing, pivot to Track D (benchmark/resource).

---

### P5. Biological characterization (weeks 10-14)

**Tasks.**
- GSEA/GSVA (Hallmark, Reactome, KEGG).
- TF activity: decoupleR (DoRothEA, PROGENy).
- Kinase activity: KSEA/PhosphoSitePlus; phosphosite enrichment.
- Metabolism: scFEA/Compass proxies; expression-based flux.
- Subgroup biology: MMRd vs MMRp, EOCRC (<50), right vs left, stage.
- Marker panel selection for future IHC (3-5 per ecotype).
- Optional: neoantigen/immunopeptidome overlay; microbial peptide signals.

**Outputs.**
- `results/biology/` (pathway, TF, kinase, metabolism tables)
- Ecotype marker panels and signatures.
- Figure 4 data.

**QC gate G4a.**
- Enrichment results survive FDR <0.05 and are not driven by a single cohort.
- Marker panels validated in single-cell reference (specificity).

---

### P6. Spatial validation (weeks 12-16)

**P6.1 Processing.**
- Per dataset: platform-specific QC, normalization, spot/bin filtering.
- Deconvolution with cell2location/RCTD using reference signatures.
- Ecotype scoring per spot (program scores; classifier transfer).

**P6.2 Niche analysis.**
- Niche detection: CellCharter; neighborhood enrichment: squidpy.
- Ligand-receptor: CellChat/NicheNet.
- TLS score, invasive-front score, hypoxia, immune exclusion.
- Co-localization of ecotype programs with spatial domains.

**P6.3 Replication.**
- Test ecotype-niche associations across ≥2 datasets/platforms.
- Meta-analyze effect sizes; report heterogeneity.
- Compare bulk-predicted ecotype with spatial ecotype score.

**Outputs.**
- `results/spatial/` (per dataset + meta-analysis)
- Spatial niche map; ecotype-niche coupling.
- Figure 5 data.

**QC gate G4.**
- Ecotype-specific niches replicate in ≥2 datasets.
- Spatial marker patterns consistent with bulk predictions.
- If G4 fails: use spatial data descriptively; do not claim spatial validation.

---

### P7. Clinical validation (weeks 14-20)

**P7.1 Endpoints.**
- Primary: CSS (TCGA); OS where CSS is unavailable.
- Secondary: RFS, OS.
- Exploratory: response to ICI, anti-EGFR, adjuvant chemotherapy.

**P7.2 Models.**
- Base model: age + sex + stage + MSI + sidedness + CMS/CRIS (if available).
- Full model: base + ecotype.
- Algorithms: penalized Cox (elastic net), random survival forest, gradient boosting.
- Nested cross-validation for tuning; bootstrap for optimism.

**P7.3 Metrics.**
- C-index; 5-year time-dependent AUC.
- Calibration slope/intercept and calibration plots.
- Decision-curve analysis (net benefit across thresholds).
- NRI/IDI versus base model.
- Kaplan-Meier by ecotype and risk group.

**P7.4 External validation.**
- Apply locked classifier to GSE39582, GSE17536, GSE14333, GSE38832, GSE87211, GSE103479, CPTAC.
- Report discrimination, calibration, and DCA per cohort.
- Meta-analyze across cohorts; report heterogeneity.

**P7.5 Subgroups and sensitivity.**
- Stage II/III; MSS; EOCRC; sidedness; purity; batch-correction method; cohort exclusion.

**Outputs.**
- `results/clinical/` (model objects, metrics, calibration, DCA)
- Figure 6 data; Tables 1-3.

**Decision gate G5.**
- Incremental value in ≥1 independent cohort (likelihood-ratio or DCA).
- Calibration acceptable (slope 0.8-1.2; no systematic over/under-estimation).
- If G5 fails: reframe as a biological/resource paper; target *CSBJ*/*Briefings in Bioinformatics*; report null clinical result honestly.

---

### P8. Therapeutic vulnerability and drug repurposing (weeks 18-22)

**Tasks.**
- Build ecotype expression signatures (up/down).
- CMap/LINCS connectivity scoring.
- GDSC/CTRP/PRISM cell-line sensitivity; map cell lines to ecotype-like states using CCLE expression.
- DepMap CRISPR dependencies for ecotype markers/pathways.
- Filter by clinical stage, Blood-Brain Barrier, known toxicity.
- Prioritize 3-5 candidates per ecotype with mechanism.

**Outputs.**
- `results/drug/` (ranked compounds, evidence table)
- Figure 7 data; supplementary drug table.

**Framing.**
- Hypothesis-generating only; no claims of efficacy.
- State that experimental validation is required.

---

### P9. Reproducibility, resource and portal (weeks 20-26)

**Tasks.**
- Nextflow/Snakemake workflow for all modules.
- Containers (Docker/Singularity); `renv.lock`, `environment.yml`.
- Seeds, parameters, checksums, session info.
- Shiny/Cellxgene portal: upload expression -> ecotype score; explore atlas.
- GitHub repository + Zenodo DOI.
- Reproducibility audit: fresh clone, rerun a subset, compare outputs.

**Outputs.**
- `workflow/`, `containers/`, `portal/`, `manifests/`
- Zenodo release; portal URL.

**QC gate G6.**
- A fresh environment reproduces key figures/tables from raw/processed inputs.
- All code and metadata public (or available on request for any restricted data).

---

### P10. Manuscript preparation (weeks 24-36)

**Tasks.**
- Finalize figures/tables.
- Write abstract, introduction, results, discussion, methods.
- Prepare cover letter and journal-specific highlights.
- Post preprint (bioRxiv/medRxiv).
- Submit; prepare point-by-point response.

**Outputs.**
- Manuscript draft, figures, supplementary tables, cover letter.
- Preprint DOI.

---

### P11. Extension phases (E1-E7; optional, weeks 40-70)

**E1. Proteogenomic anchoring and progression/metastasis (weeks 40-46).**
- PXD079707: protein-level ecotype scores; concordance with RNA-based scores; survival associations.
- CPTAC: independent proteome/phosphoproteome validation; kinase activity.
- PXD052391: trajectory across nontumor -> intraepithelial neoplasia -> infiltration -> advanced; stage-specific ecotype shifts.
- IPX0007391001: matched primary/normal/liver metastasis; ecotype retention vs switch; comparison with C1/C2 proteomic subtypes.
- **Outputs:** Figure 8; Tables 4-5.

**E2. Pathology foundation-model morphology (weeks 44-50).**
- WSI QC/tiling; CHIEF/MUSK/UNI/CONCH frozen features; attention MIL; nested CV; external validation (TCGA <-> PXD079707).
- Predict ecotypes and prognosis; compare with transcriptomic/spatial ecotypes.
- **Outputs:** Figure 9a.

**E3. Circulating ecotype (weeks 46-52).**
- cfDNA methylation/5hmC deconvolution; EV/RNA signatures; association with tissue ecotype, stage, MRD and outcome.
- **Outputs:** Figure 9b.

**E4. Therapy-response mapping (weeks 46-54).**
- Anti-EGFR (GSE140973), chemotherapy (GSE103479/GSE106584/GSE156915), TNT (GSE339362), METIMMOX (GSE313035), MMRd ICI (GSE333300/GSE330748).
- Pre-specified response definitions; adjusted analyses; exploratory only.
- **Outputs:** Figure 6b; Table 3.

**E5. Cross-platform and method benchmark (weeks 48-56).**
- CMS/CRIS/PDS stability; integration/deconvolution/niche/survival method comparisons; runtime.
- **Outputs:** Figure 10; supplementary benchmark.

**E6. Causal triangulation (optional; weeks 52-58).**
- MR/colocalization for ecotype driver genes; eQTL/pQTL resources; no causal claims without robust instruments.

**E7. Microbiome/microbial peptides (optional; weeks 52-60).**
- curatedMetagenomicData; MetaPhlAn 4/HUMAnN 3; MMUPHin; ANCOM-BC2; microbial peptides from CPTAC/PXD079707; plasma microbial cfRNA.
- Companion or supplementary; no causal claims.

**Extension decision gate G7.**
- Each extension must meet its pre-specified criterion; otherwise report null or move to supplement.
- The core paper (Aims 1-4) is not delayed by extensions.

---

## 5. Statistical analysis plan (SAP)

### 5.1 Primary analysis
- **Hypothesis:** full model (base + ecotype) improves CSS discrimination/calibration vs base model.
- **Test:** likelihood-ratio test for nested Cox models; report C-index difference with bootstrap CI.
- **Alpha:** 0.05 (two-sided).
- **Effect measure:** HR per ecotype vs reference; adjusted survival curves.

### 5.2 Secondary analyses
- OS and RFS (same framework).
- Ecotype association with stage, MSI, sidedness, EOCRC (chi-square/Fisher; multinomial regression).
- Spatial niche association (mixed-effects models with dataset as random effect).
- Drug connectivity (rank-based; FDR).

### 5.3 Multiplicity
- Primary endpoint: no adjustment.
- Secondary endpoints: Benjamini-Hochberg FDR <0.05 within families.
- Exploratory: report nominal p-values with explicit "exploratory" label.

### 5.4 Missing data
- Clinical covariates: multiple imputation (`mice`, 20 imputations) if missing <20%; complete-case sensitivity.
- Outcomes: exclude samples with missing time/event; report.
- Omic views: no imputation of entire views; use available-view integration and sensitivity.

### 5.5 Model validation
- Internal: bootstrap optimism correction (1,000 resamples); nested CV.
- External: pre-specified cohorts; no recalibration unless explicitly reported.
- Calibration: slope/intercept; calibration plots by decile.
- DCA: net benefit at threshold probabilities 0.1-0.5.

### 5.6 Power
- Assume 150 CSS events in discovery; detect HR 1.5 per ecotype vs rest at 80% power (two-sided 0.05) with ~600 samples.
- For external validation, report precision (CI width) rather than post-hoc power.
- Spatial: replicate in ≥2 datasets; report effect sizes and heterogeneity.

### 5.7 Reporting
- TRIPOD+AI for prediction models.
- REMARK for prognostic markers.
- Full cohort flow diagram; STROBE-like reporting for observational data.
- Report all pre-specified analyses and negative results.

---

## 6. Operational ecotype definition (prevents post-hoc storytelling)

1. **Input:** MOFA+ factors with ≥5% variance explained in ≥2 views and interpretable loadings; program scores.
2. **Clustering:** consensus clustering, K=2-8, 1,000 bootstraps; compute PAC, silhouette, ARI.
3. **Selection:** choose the smallest K with bootstrap ARI >0.70, all groups ≥5%, and distinct top programs; if multiple K qualify, prefer the one with better clinical separation after adjustment.
4. **Labeling:** name each ecotype by its top three programs and dominant TME state; do not use mechanism names unless supported.
5. **Classifier:** nearest shrunken centroid or multinomial elastic net on program scores; nested CV; export weights.
6. **Replication:** apply locked classifier; require ARI >0.5 with discovery in ≥2 cohorts.
7. **Comparison:** CMS/CRIS/PDS concordance and reclassification; survival of discordant cases.
8. **Stopping rule:** if no stable solution, do not force ecotypes; pivot to a continuum/factor-based paper or Track D.

---

## 7. Spatial analysis plan (detail)

**Per dataset.**
1. QC: spot/bin counts, mitochondrial fraction, tissue coverage; remove low-quality spots.
2. Normalization: platform-appropriate (SCTransform/sctransform or scanpy normalization).
3. Deconvolution: cell2location (reference signatures), RCTD/SpatialDecon as sensitivity.
4. Ecotype scoring: program scores + locked classifier; smooth across neighbors.
5. Niche detection: CellCharter; cluster spots into niches; annotate by cell composition.
6. Neighborhood enrichment: squidpy; test ecotype-score enrichment in niches.
7. Ligand-receptor: CellChat/NicheNet; focus on ecotype-specific interactions.
8. Domain scores: TLS, invasive front, hypoxia, immune exclusion, CAF-rich.
9. Replication: meta-analysis across datasets; report effect sizes/heterogeneity.
10. Visualization: spatial maps, niche composition, representative sections.

**Interpretation rules.**
- Spatial data validate bulk ecotypes; they do not define them.
- Require replication in ≥2 datasets/platforms for any spatial claim.
- Report platform and spot-size differences.

---

## 8. Drug repurposing plan (detail)

- **Signature:** ecotype up/down genes (top 150 each; FDR <0.05).
- **CMap/LINCS:** query; connectivity score; rank compounds.
- **Cell-line sensitivity:** GDSC/CTRP/PRISM; map CCLE cell lines to ecotype-like states; test association.
- **Dependencies:** DepMap CRISPR; ecotype markers/pathways.
- **Filtering:** clinical stage, known CRC activity, toxicity, BBB (if relevant).
- **Prioritization:** mechanism plausibility + existing clinical agent + ecotype specificity.
- **Output:** table with compound, target, direction, evidence, caveats.
- **Claim:** hypothesis-generating; no efficacy claim.

---

## 9. Reproducibility and resource

- **Workflow:** Nextflow DSL2 or Snakemake; one module per phase; parameterized by `study_config.yaml`.
- **Containers:** Docker/Singularity images per module; pin versions.
- **R:** `renv.lock`; **Python:** `environment.yml`/`requirements.txt`.
- **Seeds:** fixed for integration, clustering, CV, bootstrap.
- **Checksums:** manifest for raw and processed files.
- **Documentation:** README, methods, parameter table, troubleshooting.
- **Portal:** Shiny app for ecotype scoring; Cellxgene for atlas exploration.
- **Release:** GitHub + Zenodo DOI; code license (MIT/BSD-3); data according to source licenses.

---

## 10. Figures, tables and manuscript

### Main figures (10)
1. Study design, cohort flow and data matrix (including proteogenomic/progression/metastasis resources).
2. Single-cell reference atlas and programs.
3. Multi-view integration, ecotype discovery and stability.
4. Ecotype biology (pathways, TFs, kinases, metabolism; subgroups).
5. Spatial validation and niches.
6. Clinical validation (survival, calibration, DCA, external cohorts) and therapy associations.
7. Therapeutic vulnerabilities and resource/portal.
8. Progression and metastasis axis (PXD052391, IPX0007391001, GSE201348, GSE225857).
9. Morphology-molecular and circulating ecotypes (WSI features, cfDNA/EV).
10. Benchmark and reproducibility (CMS/CRIS/PDS instability, method comparison, portal).

### Tables
- Table 1: cohort characteristics (discovery + validation + proteogenomic cohorts).
- Table 2: ecotype marker programs and comparison with CMS/CRIS/PDS.
- Table 3: prioritized drug candidates, therapy associations and evidence.
- Table 4: progression/metastasis cohort characteristics and ecotype transitions.
- Table 5: morphology and circulating ecotype performance.
- Table 6: benchmark results (stability, runtime, calibration).

### Supplementary
- QC reports; deconvolution benchmarks; sensitivity analyses; all signatures; negative results; code/manifest; cohort overlap.

### Manuscript structure
- Abstract: background (subtype instability), approach (cross-platform ecotypes), results (stability, spatial niches, clinical utility), conclusion (reproducible framework), data/code availability.
- Introduction: CRC burden; heterogeneity; limits of CMS/CRIS/PDS; why single-cell-anchored, cross-platform ecotypes.
- Results: 5 subsections (resource, atlas, ecotypes, spatial, clinical/therapeutic).
- Discussion: comparison with prior atlases; limitations (public data, retrospective, computational only); future validation.
- Methods: detailed, reproducible.
- Data/code availability.

---

## 11. Journal tailoring

| Journal | What to emphasize | What to avoid | Suggested framing |
|---|---|---|---|
| *Journal of Translational Medicine* | Clinical translation: risk stratification, stage II/III decisions, MSS/immune-hot, therapy selection, drug hypotheses | Mechanistic claims without wet lab | "A reproducible ecotype framework with translational utility" |
| *Computational and Structural Biotechnology Journal* | Computational framework, benchmarks, scalability, runtime, reproducibility, resource | Clinical overclaiming | "A cross-platform multi-omic integration framework for CRC ecotypes" |
| *Briefings in Bioinformatics* | Method/resource contribution, benchmarking, usability, code/portal, comparison with existing tools | Purely descriptive biology | "A reproducible framework and benchmark for single-cell-anchored cancer ecotypes" |

**Common requirements.**
- Data/code availability statement.
- Limitations paragraph that explicitly addresses public-data and computational-only constraints.
- Clear distinction between association and causation.
- Calibration/DCA, not just AUC.
- Pre-registration/SAP mention (strengthens rigor).

---

## 12. Reviewer-proofing checklist

- [ ] Patient-level splits; test sets locked before modeling.
- [ ] No feature selection or tuning on the full dataset.
- [ ] Nested cross-validation; bootstrap optimism correction.
- [ ] External validation in ≥3 independent cohorts.
- [ ] Calibration and decision-curve analysis reported.
- [ ] Batch correction evaluated with negative controls; sensitivity analyses.
- [ ] Purity adjustment and sensitivity.
- [ ] Multiple-testing correction; all pre-specified analyses reported.
- [ ] Class imbalance and missing data addressed.
- [ ] CMS/CRIS/PDS comparison with ARI/NMI and reclassification.
- [ ] Spatial claims replicated in ≥2 datasets.
- [ ] Drug claims labeled hypothesis-generating.
- [ ] Code, manifests, versions, seeds available.
- [ ] Limitations explicitly stated (retrospective, public data, no wet lab).

---

## 13. Risk register and decision points

| Risk | Trigger | Response |
|---|---|---|
| Ecotypes unstable | ARI <0.70 or <3 groups | Adjust views/K; use factor continuum; pivot to Track D |
| No clinical increment | G5 fails | Reframe as biological/resource paper; report null result |
| Batch confounding | PCA separation persists | ComBat/MMUPHin; include batch; restrict to one platform for primary analysis |
| CPTAC mismatch | Few overlapping genes/samples | Use CPTAC as orthogonal validation only; do not integrate at sample level |
| Spatial heterogeneity | No replication | Use spatial descriptively; remove spatial validation claim |
| Reviewer demands wet lab | Top-tier review | Target computational journals; prepare collaboration plan |
| Scooped | Preprint conflict | Emphasize reproducibility/clinical utility; release resource early |
| Data overlap | Same patients across cohorts | Deduplicate; exclude overlapping cohort from validation |
| PXD079707 is new/unpublished | File completeness/metadata unclear | Verify; use CPTAC as fallback; do not rely solely on it |
| PXD052391 WES restricted | Cannot access WES | Use open proteome/phosphoproteome; request WES only if essential |
| Pathology FM compute/domain shift | GPU/time or poor external performance | Frozen features; CHIEF/MUSK/UNI/CONCH benchmark; external slide sets |
| Liquid biopsy heterogeneity | Analyte/platform differences | Pre-specify analyte-specific features; no pooling; external validation |
| Scope creep | Extensions delay core paper | Submit core at week 36; extensions as follow-up/revision |

---

## 14. Timeline (weeks)

| Weeks | Phase | Milestone |
|---|---|---|
| 1-2 | P0 | Manifest, SAP, pre-registration |
| 2-6 | P1 | Bulk harmonization; QC |
| 4-8 | P2 | Single-cell reference atlas |
| 6-9 | P3 | Deconvolution/program scores |
| 8-12 | P4 | Ecotypes and locked classifier |
| 10-14 | P5 | Biological characterization |
| 12-16 | P6 | Spatial validation |
| 14-20 | P7 | Clinical validation |
| 18-22 | P8 | Drug repurposing |
| 20-26 | P9 | Workflow/portal/release |
| 24-36 | P10 | Manuscript and submission |

Overlap is intentional; the pilot (weeks 1-3) runs in parallel with P0/P1.

**Total effort:** the core analysis is ~36 weeks; add 3-6 months for manuscript revision and resubmission (~9-12 months overall).

### Extension timeline (weeks 40-70; optional)
| Weeks | Extension | Milestone |
|---|---|---|
| 40-46 | E1 | Proteogenomics; progression/metastasis |
| 44-50 | E2 | Pathology foundation model |
| 46-52 | E3 | Circulating ecotype |
| 46-54 | E4 | Therapy-response mapping |
| 48-56 | E5 | Cross-platform/method benchmark |
| 52-58 | E6 | Causal triangulation (optional) |
| 52-60 | E7 | Microbiome/microbial peptides (optional) |
| 60-70 | — | Extended manuscript/revision |

**Rule:** submit the core paper at week 36 if the core criteria pass. Extensions are a follow-up or a revision strength; they must not delay the core paper.

---

## 15. Immediate next actions

1. Approve this protocol and the analysis cohort definition.
2. Download Tier 1 datasets manually (TCGA, PXD079707 subset, PXD052391, IPX0007391001, GSE39582, GSE132465, GSE225857); record checksums.
3. Approve the software list; install only after confirmation.
4. Run the 3-week pilot (`PILOT_PLAN.md`) before full data acquisition, including a PXD079707 subset.
5. Decide whether extensions (E1-E7) are in the same paper or a follow-up; freeze the core scope.
6. Register the SAP on OSF after the pilot confirms feasibility.
7. Freeze `study_config.yaml` and begin P1.

---

## Appendix A. Minimum viable analysis (if time is limited)

If the full protocol is too ambitious, the minimum publishable core is:
1. TCGA + CPTAC + PXD079707 subset + GSE39582 + GSE132465.
2. Single-cell programs -> bulk deconvolution -> MOFA+ -> 3-4 ecotypes.
3. Stability + comparison with CMS/CRIS.
4. External prognostic validation in GSE39582/GSE17536.
5. One spatial dataset for descriptive validation.
6. Drug repurposing table.
This supports a *CSBJ*/*J Transl Med* paper; the full protocol supports *Briefings in Bioinformatics* or a higher-tier translational journal. If only one extension is possible, prioritize PXD079707 proteomic anchoring or the progression cohort (PXD052391).

---

## Appendix B. What would raise the ceiling

- Wet-lab validation of 2-3 ecotype markers (IHC/CODEX) or drug hits.
- A prospective clinical cohort.
- A unique resource (portal + benchmark) with broad usability.
- Mechanistic follow-up of one ecotype (e.g., CAF-myeloid circuit).
- Full integration of the 1,000-CRC proteogenome (PXD079707) with single-cell/spatial data.
- A pathology foundation-model ecotype classifier validated on external WSIs.
- A circulating ecotype signature validated for MRD/outcome.
Without these, keep claims calibrated and target the chosen journals.
