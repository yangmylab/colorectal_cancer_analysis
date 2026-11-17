# Methods and tools toolbox (2023-2026)

**Prepared:** 2026-11-17
**Purpose:** select the analysis stack for Track A, with fit, caveats, validation requirements and compute.
**Rule:** every tool class gets one primary method + at least one sensitivity method; no single-method claims.

---

## 1. Selection principles

1. **Benchmark, don't adopt blindly.** New methods (foundation models, spatial integration) must be compared against established baselines on CRC data.
2. **Match the data structure.** Bulk cohorts, single-cell atlases, spatial sections and proteomics each need different normalization and validation.
3. **External validation beats internal metrics.** A method that improves cross-cohort replication is preferred over one that maximizes within-cohort fit.
4. **Interpretability matters.** Ecotypes must be explainable by programs/markers, not only by latent factors.
5. **Reproducibility is part of the method.** Containers, seeds, lockfiles and manifests.
6. **Pre-specify.** Choose primary tools and sensitivity analyses before outcome modeling.

---

## 2. Data access and harmonization

| Tool | Role | Notes |
|---|---|---|
| `gdc-client` | TCGA/GDC download | Manifest + checksums; open tiers only |
| `sra-tools` (`prefetch`, `fasterq-dump`) | SRA/ENA raw data | Only if reprocessing |
| `pdc-client` / PDC portal | CPTAC proteomics | Prefer processed tables |
| PRIDE API / `pridepy` | ProteomeXchange datasets | PXD079707, PXD052391, PXD019103 |
| `GEOparse` / GEO FTP | GEO series matrices | Record platform and probe annotation |
| `curatedMetagenomicData` | Uniform microbiome profiles | Batch-aware meta-analysis |
| `cBioPortal` API | Mutation/CNV/clinical | No expression for most MSK/GENIE samples |
| `ExperimentHub` / `AnnotationHub` | Reference annotations | Version-pin |

**Caveat:** sample overlap across GEO series is common; check patient IDs before using a cohort as validation.

---

## 3. Bulk transcriptome, epigenome and proteome

| Task | Primary | Sensitivity / alternatives | Caveats |
|---|---|---|---|
| Quantification | STAR + featureCounts; Salmon | Kallisto, RSEM | Gene-level vs isoform-level |
| Normalization | TPM/log2; RMA for arrays | FPKM-UQ, VST | Do not compare raw TPM across pipelines |
| Batch correction | `ComBat-seq` (counts), `limma::removeBatchEffect` (arrays) | RUV, SVA | Report uncorrected results; batch can be confounded with outcome |
| Differential expression | `limma`/`DESeq2`/`edgeR` | Wilcoxon, permutation | FDR; adjust purity/stage |
| Methylation | `minfi` + BMIQ; `limma`/`bumphunter`/`DMRcate` | `sesame`, `wateRmelon` | Beta vs M-values; cell-type heterogeneity |
| CNV | ASCAT, CNVkit, PureCN | FACETS | Purity/ploidy |
| Mutation/MSI | `maftools`, MSIsensor | MMR gene expression/IHC | MSI calling differs by tool |
| Proteomics | DIA-NN, Spectronaut, FragPipe; `MSstats` | MaxQuant | Normalization, missingness (MNAR), batch |
| Phosphoproteomics | Site-level normalization; PTM-SEA, KSEA, RoKAI | GSEA on phosphosites | Site localization, protein abundance confounding |
| Pathway/TF activity | `GSVA`/`fgsea`; `decoupleR` (DoRothEA, PROGENy) | GSEA desktop | Redundancy; directionality |
| Purity | ESTIMATE, ABSOLUTE | pathology estimate | Include as covariate |

**Validation:** cross-platform correlation (RNA-protein), negative-control probes, batch-PCA, outlier audit.

---

## 4. Single-cell and single-nucleus

| Task | Primary | Sensitivity | Caveats |
|---|---|---|---|
| QC | `scanpy`/`Seurat`; doublet removal (`scDblFinder`/`Solo`); ambient correction (`CellBender`) | `EmptyDrops` | Over-filtering removes malignant/neutrophil cells |
| Integration | `scVI`/`scANVI` | Harmony, scPoli, scArches, BBKNN | Batch vs biology; report mixing metrics |
| Annotation | markers + `CellTypist` | manual, scGPT/Geneformer (benchmark) | Foundation models can mislabel rare states |
| Malignant inference | `inferCNV`, `CopyKAT` | `Numbat`, `SCEVAN` | Normal-cell contamination; diploid tumors |
| Programs | NMF/cNMF | PCA, ICA, meta-program | Stability across atlases |
| Chromatin | `ArchR`/`Signac`; gene activity; SCENIC/pySCENIC | Cicero | Sparse data; peak-to-gene linking |
| Multi-modal | `MultiVI` | `totalVI`, `scMoMaT` | Missing modalities |
| Reference mapping | `scArches`, `scPoli` | `scANVI` | Query batch effects |

**Foundation models (scGPT, Geneformer, scFoundation, CellPLM):** use as **feature extractors/benchmarks**, not as unquestioned annotators. Validate against marker-based labels and external atlases. Recent work shows they capture organized biological knowledge but limited regulatory logic; benchmark before use.

---

## 5. Spatial transcriptomics and spatial multi-omics

| Task | Primary | Sensitivity | Caveats |
|---|---|---|---|
| Processing | Space Ranger / platform pipeline | manual | Platform version matters |
| Deconvolution | `cell2location` | RCTD, SpatialDecon, BayesPrism, InstaPrism, spOT-NMF | Reference dependence; spot size |
| Niche detection | `CellCharter` | BANKSY, CellNiche, GraphST, SpaGCN, UTAG | Number of niches; stability |
| Neighborhood analysis | `squidpy` | `sosta` | Radius/graph choice |
| Cell-cell communication | `CellChat` v2 | LIANA+, NicheNet, Connectome | Ligand-receptor database bias; spatial distance |
| Multi-omics alignment | `SpatialGlue`, `DePass`, `MISO` | Riemannian alignment, HESTIA | Platform pairing; missing modalities |
| Image-guided enhancement | `CoxFormer`, image-guided spatial omics | HESTIA | Generative methods need validation |
| Whole-genome spatial | RAEFISH (PMID 41038164) | Stereo-seq, Visium HD, Xenium, CosMx | Not yet widely available |
| Spatial CRISPR | SPAC-seq/TARDIS (PMID 42190664) | — | Perturbation, not observational |

**Benchmark evidence:** four subcellular platforms differ in sensitivity, diffusion, segmentation and clustering (PMID 41107232). Therefore:
- Do not pool platforms without platform-aware normalization.
- Replicate niche claims in ≥2 datasets.
- Report platform, spot/bin size and segmentation method.

**CRC-specific spatial resources:** GSE225857, GSE283052, GSE320042 (+ GSE320038/GSE320040), GSE284989, GSE280311/GSE280312, GSE270745, GSE217414; plus HTAN/10x/CROST/SpatialDB for additional sections.

---

## 6. Multi-omics integration and subtyping

| Method | Type | Strengths | Caveats |
|---|---|---|---|
| `MOFA+` | Factor analysis | Missing views, interpretable factors, variance explained | Linear; factor number; scale |
| `MEFISTO` | Spatial/temporal MOFA | Smooth factors over space/time | Requires coordinates |
| `iClusterBayes` | Bayesian latent factors | Joint clustering | Compute; model assumptions |
| `SNF`/`SNFtool` | Similarity fusion | Simple, robust | Graph construction |
| `DIABLO`/`mixOmics` | Supervised multi-block | Integrates with outcome | Supervised leakage risk |
| `MIDAS` | Multi-omics autoencoder | Scalable | Interpretability |
| `MultiVI` | Single-cell multi-modal | Missing modalities | Single-cell focus |
| `GLUE` | Graph-linked embedding | Cross-modality | Requires prior graph |
| `PINSPlus` | Perturbation clustering | Stability | Bulk focus |
| `CIMLR` | Multi-kernel clustering | Multi-omic | Compute |
| `NEMO` | Network-based | Similarity fusion | Graph dependence |

**Recommended:** MOFA+ primary; iClusterBayes/SNF/DIABLO sensitivity; consensus clustering with bootstrap ARI; compare to CMS/CRIS/PDS.

---

## 7. Clinical prediction, survival and causal inference

| Task | Primary | Sensitivity | Caveats |
|---|---|---|---|
| Survival model | Penalized Cox (`glmnet`) | Random Survival Forest, XGBoost/LightGBM, DeepSurv | Nonlinear models need calibration |
| Discrimination | C-index, time-dependent AUC | Uno's C, IBS | Report CIs |
| Calibration | Calibration slope/intercept, plots | Brier score | Recalibration must be reported |
| Clinical utility | Decision-curve analysis | Net reclassification improvement | Thresholds pre-specified |
| Explainability | SHAP, SurvSHAP(t) | permutation importance | Correlated features |
| Validation | Nested CV + external cohorts | bootstrap optimism | No leakage |
| Causal | MR (`TwoSampleMR`), colocalization (`coloc`) | SMR/HEIDI, mediation | Instrument validity, pleiotropy |
| QTL resources | GTEx, eQTL Catalogue, UKB-PPP, deCODE | — | Tissue/cell-type mismatch |

**Reporting:** TRIPOD+AI, REMARK, STROBE.

---

## 8. Pathology foundation models and imaging

| Tool | Role | Notes |
|---|---|---|
| CHIEF (PMID 39232164) | General pathology foundation model | Weakly supervised; strong domain-shift performance |
| MUSK (PMID 39779851) | Vision-language pathology model | Image-text; biomarker/outcome prediction |
| UNI / CONCH | Pathology foundation models | Feature extraction; check license |
| Prov-GigaPath | Whole-slide foundation model | Large-scale; compute |
| HIPT / CLAM | WSI multiple-instance learning | Baselines |
| OpenSlide, TIAToolbox, PathML, MONAI | WSI I/O and preprocessing | Tissue detection, tiling |
| TCGA diagnostic slides (625) | H&E resource | Slide quality varies |
| PXD079707 WSIs (1,000) | Proteogenomic-matched H&E | Direct morphology-proteome link |

**Use:** extract tile/slide features -> attention MIL or linear probe -> predict ecotypes/prognosis. **Validate:** external slide sets, calibration, no patient overlap. **Caveat:** foundation-model features are not biologically interpretable; use them as one layer, not the sole evidence.

---

## 9. Microbiome and metabolome

| Task | Primary | Sensitivity | Caveats |
|---|---|---|---|
| Profiling | MetaPhlAn 4; HUMAnN 3 | mOTUs, Kraken2 | Database/version effects |
| Batch correction | MMUPHin | ComBat | Study confounding |
| Differential abundance | ANCOM-BC2 | LEfSe, MaAsLin2, Wilcoxon | Compositionality; zero inflation |
| Metabolite prediction | MelonnPan | — | Prediction, not measurement |
| Classification | SIAMCAT | random forest | Overfitting; confounders |
| Metabolomics | Metabolomics Workbench / MetaboLights; METASPACE, SpaceM | — | Sparse human CRC tissue data |

**Use:** optional companion; pre-specify confounders (BMI, diet, antibiotics, geography); no causal claims.

---

## 10. Drug repurposing and perturbation

| Resource | Role | Caveats |
|---|---|---|
| CMap / LINCS | Connectivity signatures | Cell-line context |
| GDSC / CTRP / PRISM | Cell-line sensitivity | 2D vs 3D; tissue mismatch |
| DepMap | CRISPR dependencies | Essentiality vs selectivity |
| DREIMT | Immune drug repurposing | Database bias |
| Open Targets / DrugBank | Target evidence | Annotation lag |
| scPerturb / GEARS / CPA | Perturbation prediction | Benchmark before use |

**Use:** hypothesis-generating; require orthogonal evidence; no efficacy claims.

---

## 11. Workflow and reproducibility

| Tool | Role |
|---|---|
| Nextflow / Snakemake | Pipeline orchestration |
| Docker / Singularity | Containers |
| `renv` / conda / pip-tools | Environment lockfiles |
| `targets` (R) | Dependency-aware pipelines |
| Quarto / R Markdown | Reproducible reports |
| Zenodo | DOI release |
| Cellxgene / Shiny | Atlas/portal |
| `sessioninfo` / `reticulate` | Version capture |

---

## 12. Recommended Track A stack

### Primary
- Bulk: STAR/Salmon -> `limma`/`DESeq2`; `ComBat-seq`; `GSVA`; `decoupleR`.
- Single-cell: `scVI`/`scANVI`; `CellTypist`; `inferCNV`/`CopyKAT`; NMF; `SCENIC`.
- Deconvolution: `cell2location` (spatial) + consensus bulk deconvolution (CIBERSORTx/EPIC/quanTIseq/MCP-counter/xCell).
- Integration: `MOFA+`; `ConsensusClusterPlus`; bootstrap ARI.
- Spatial: `CellCharter` + `squidpy` + `CellChat` v2.
- Clinical: penalized Cox + RSF; calibration/DCA; nested CV.
- Pathology: CHIEF/MUSK features + attention MIL.
- Drug: CMap/LINCS + GDSC/CTRP/PRISM + DepMap.

### Sensitivity
- Harmony instead of scVI; RCTD/BayesPrism instead of cell2location; iClusterBayes/SNF/DIABLO instead of MOFA+; BANKSY/GraphST instead of CellCharter; XGBoost/DeepSurv instead of Cox; UNI/CONCH instead of CHIEF/MUSK.

### Benchmark
- Integration: MOFA+ vs iClusterBayes vs SNF vs DIABLO.
- Deconvolution: cell2location vs RCTD vs BayesPrism/InstaPrism; synthetic spots.
- Niche: CellCharter vs BANKSY vs CellNiche.
- Survival: Cox vs RSF vs XGBoost/DeepSurv; calibration and DCA.
- Pathology: CHIEF vs MUSK vs UNI/CONCH.

---

## 13. Compute and data-size estimates

| Task | Resource | Time |
|---|---|---|
| Bulk harmonization (20 cohorts) | 16 cores, 64 GB | 1-3 days |
| Single-cell integration (1-4M cells) | 32 cores, 128 GB; GPU 24 GB recommended | 1-2 weeks |
| NMF/programs | 16-32 cores | 2-5 days |
| MOFA+ | 16-32 cores, 64-128 GB | 1-3 days |
| Consensus clustering (1,000 bootstraps) | 32 cores | 3-7 days |
| Spatial deconvolution (40 sections) | 16-32 cores, 64 GB; GPU optional | 1-2 weeks |
| Pathology feature extraction (1,000 WSIs) | GPU 24 GB | 3-10 days |
| Survival/modeling | 16 cores | 2-5 days |
| Microbiome (optional) | 16 cores | 3-7 days |

**Storage:** 2-4 TB (raw + processed + containers); processed matrices only for most atlases.

---

## 14. Validation rules by tool class

| Class | Minimum validation |
|---|---|
| Deconvolution | Method concordance (Spearman >0.7), synthetic-spot benchmark, purity sensitivity |
| Integration | Bootstrap ARI >0.70, external replication ARI >0.50, batch-mixing metrics |
| Clustering | Stability selection, ≥5% per group, interpretable markers |
| Spatial niches | Replication in ≥2 datasets/platforms; report platform effects |
| Cell-cell communication | Spatial distance-aware, ligand-receptor database sensitivity |
| Survival | Nested CV, calibration, DCA, external cohorts |
| Pathology FM | External slide set, calibration, no patient overlap |
| Microbiome | Confounder adjustment, batch correction, negative controls |
| Drug | Orthogonal evidence, mechanism, no efficacy claim |
| Causal | Instrument strength, pleiotropy, colocalization |

---

## 15. What not to do

- Do not use a single deconvolution method as ground truth.
- Do not pool spatial platforms without platform-aware normalization.
- Do not use CMS/CRIS labels as ground truth; treat them as comparators.
- Do not use foundation-model embeddings as the only evidence for a biological claim.
- Do not select features or tune on the full dataset.
- Do not report drug candidates as validated treatments.
- Do not make causal claims from observational public data without MR/colocalization support.
- Do not ignore batch, purity, stage, MSI and ancestry confounding.

---

## 16. Key method references (verified)

- PMID 41107232 - spatial platform benchmark, *Nat Commun* 2025
- PMID 39232164 - CHIEF, *Nature* 2024
- PMID 39779851 - MUSK, *Nature* 2025
- PMID 41038164 - RAEFISH, *Cell* 2025
- PMID 42190664 - SPAC-seq/TARDIS, *Cell* 2026
- PMID 40488620 - ProCanFDL, *Cancer Discov* 2025
- PMID 42020427 - CellNiche, *Nat Commun* 2026
- PMID 42568193 - image-guided spatial omics, *Bioinformatics* 2026
- PMID 42754563 - CoxFormer, *Nat Commun* 2026
- PMID 42786276 - DePass, *Nat Cell Biol* 2026
- PMID 42986684 - spOT-NMF, *Brief Bioinform* 2026
- PMID 42412825 - Riemannian spatial multi-omics alignment, *Bioinformatics* 2026
- PMID 42721220 - sosta, *PLoS Comput Biol* 2026
- PMID 41802282 - imaging-to-omics prediction review, *Brief Bioinform* 2026
