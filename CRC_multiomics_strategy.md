# Colorectal cancer multi-omics analysis strategy (public-data study)

**Prepared:** 2026-11-17
**Version:** 2.0 (expanded after the 2026 database and literature scan)
**Scope:** discovery/validation using public human CRC data only (no new wet-lab generation assumed)
**Companion files:** `DATASETS.md` (manual download checklist), `DATASETS_EXTENDED.md` (additional verified resources), `LITERATURE_GAPS_2024_2026.md` (gap matrix), `METHODS_TOOLBOX_2026.md` (new methods), `SOFTWARE.md` (environment checklist), `PROTOCOL_TrackA.md` (step-by-step protocol), `PILOT_PLAN.md` (3-week de-risking pilot), `analysis/` (project skeleton)

---

## 1. Executive summary

**Selected study (Track A; confirmed 2026-11-17).** A cross-platform, single-cell-anchored **"ecotype" taxonomy of colorectal cancer (CRC)**, defined by tumor-intrinsic programs plus spatially organized immune/stromal states, validated across transcriptome, proteome, phosphoproteome, epigenome, single-cell and spatial data, and tested for prognostic and therapeutic utility against existing classifications (CMS, CRIS, PDS).

**2026 update: why the design changed.** Three new resources make a substantially stronger study possible:
- **PXD079707** (PRIDE, CC0): 1,000 treatment-naive CRCs with DIA proteomics, targeted genomic profiling, whole-slide histology and clinical follow-up, explicitly designed to refine existing classifications and enable AI-based histology.
- **PXD052391 / IPX0002188000**: 148 patients and 435 samples spanning nontumor, intraepithelial neoplasia, infiltration and advanced CRC, with proteome and phosphoproteome.
- **IPX0007391001 / PRJCA020890**: 34 patients and 102 samples of matched primary CRC, normal colon and liver metastasis.

The study now has **proteomic depth (1,000 cases)**, a **progression axis (435 samples)** and a **metastasis axis (102 samples)**, plus 625 TCGA diagnostic slides, 1,000 matched WSIs, and large cfDNA/EV cohorts for a circulating extension. The novelty shifts from "another proteomic subtype" to a **cross-platform ecotype framework with single-cell, spatial, morphology, circulating and therapy-response validation**.

**Working title.**
> *Resolving colorectal cancer heterogeneity across transcriptome, proteome, epigenome and space: a reproducible ecotype taxonomy with clinical utility*

**One-sentence aim.**
> Determine whether CRC can be decomposed into a small number of reproducible, spatially organized tumor-immune ecotypes that are conserved across platforms and cohorts and that outperform existing bulk subtypes for prognosis and therapy selection.

**Why this is publishable now.**
1. CRC has excellent public multi-omic coverage (TCGA, CPTAC, >4 million single cells, dozens of spatial sections, large bulk validation cohorts).
2. Existing taxonomies (CMS, CRIS, PDS) are unstable across platforms/cohorts and are not clinically deployed; a rigorous, cross-platform, single-cell-anchored alternative addresses a recognized problem rather than adding another descriptive atlas.
3. Spatial and proteomic validation are feasible with public data, which is what distinguishes a "reanalysis" from a resource-grade study.
4. The analysis is fully reproducible and can include a deployable classifier and portal, which strengthens the clinical-translational case.
5. The 1,000-CRC proteogenomic resource (PXD079707) removes the previous sample-size ceiling; CPTAC, the 148-patient progression cohort and the 34-patient CRLM cohort provide independent protein-level validation.
6. Pathology foundation models (CHIEF, MUSK) plus 1,000 matched WSIs make a morphology-molecular layer feasible without new tissue collection.
7. Large cfDNA/5hmC and EV/RNA cohorts make a circulating-ecotype extension feasible, linking tissue taxonomy to non-invasive monitoring.
8. The 2026 literature explicitly asks for reproducibility and clinical utility rather than more atlases; the project is framed around that question.

**Primary deliverables.**
- Harmonized, versioned multi-omic resource and sample manifest.
- Proteogenomic ecotype anchor (PXD079707 + CPTAC + progression + CRLM cohorts).
- Reproducible ecotype classifier (R/Python package + locked model weights).
- Spatial niche map linking ecotypes to TLS, invasive front, hypoxia and immune exclusion, replicated across platforms.
- Progression and metastasis ecotype axis (nontumor -> adenoma -> carcinoma -> liver metastasis).
- Morphology-molecular ecotype model using pathology foundation models.
- Circulating ecotype signature from cfDNA/EV data.
- External validation of prognosis and therapy-response associations.
- Ecotype-to-drug vulnerability map (hypothesis-generating).
- Benchmark of CMS/CRIS/PDS stability and of integration/deconvolution/niche methods.
- Open portal, code, and intermediate results.

**Target journals (selected).** Because the study is computational-only, it is optimized for *Journal of Translational Medicine* (clinical translation), *Computational and Structural Biotechnology Journal* (framework/benchmarks), and *Briefings in Bioinformatics* (reproducible framework/resource). *Genome Medicine* / *eBioMedicine* remains a stretch only if external validation is exceptionally strong; *Gut* / *Nature Communications* / *Cancer Cell* would require experimental validation.

---

## 2. What is already established (2023-2026 landscape)

The field is no longer short of descriptive atlases. The strategic question is what remains unresolved.

| Area | Landmark public work | What it established | What remains open |
|---|---|---|---|
| Single-cell atlas | 4.27M-cell CRC immune atlas, 1,670 samples, spatial and protein profiling (*Cancer Cell* 2026, PMID 41529665) | Immune phenotypes (desert/B/T/myeloid-enriched), myeloid programs, antigen-presenting neutrophils, IL-1 niches | Cross-platform reproducibility; prospective clinical deployment; therapy-response validation |
| Age continuum | 189-tumor snRNA-seq atlas across eight decades (GSE312983) | Early-onset CRC (EOCRC) has distinct cell states | Bulk/proteome/spatial validation; whether EOCRC programs are prognostic or actionable |
| MMRd vs MMRp | 181-sample single-cell atlas (GSE178341) | Immune microenvironment differs by MMR status | Why a subset of MSS tumors is immune-hot; predictive biomarkers beyond MSI |
| Polyp-to-cancer | scRNA + scATAC continuum (GSE201348) | Cell-state and chromatin changes across adenoma-carcinoma sequence | Spatial/proteomic confirmation; targets for interception |
| Isoform/neoantigen | Long-read isoform atlas (*Cell Genomics* 2024, PMID 39216476); immunopeptidomics of MSI/MSS (2026, PMID 42106149) | Tumor-specific isoforms and noncanonical antigens; MSS is not antigen-poor | Integration with subtype/ecotype and spatial context |
| Proteogenomics | CPTAC colon cancer, 110 tumors (*Cell* 2019, PMID 31031003) | Proteomic/phosphoproteomic layers; Rb phosphorylation; glycolysis linked to immune evasion in MSI-H | No modern multi-cohort proteome-anchored taxonomy; limited spatial linkage |
| Population-scale proteogenomics | PXD079707 (1,000 treatment-naive CRCs, 2026); PXD052391 (148 patients/435 samples); IPX0007391001 (34 patients/102 samples) | Proteomic depth, progression axis and metastasis axis | Not yet integrated with single-cell, spatial, morphology, circulating or therapy-response data |
| Spatial biology | Visium/Stereo-seq/Xenium studies: liver metastasis (GSE225857), ecotypes (GSE320042), cholesterol/ICI resistance (GSE284989), TLS/CCL19 (*Cancer Cell* 2024, PMID 39137726), neutrophil niches | Niche-level organization matters; TLS and CAF states are therapy-relevant | Small samples, single-platform; no cross-modality ecotype validation |
| CAF/immune mechanisms | COL10A1+ CAF (2025), CTHRC1+ CAF-WNT5A (2025), Th17-TWEAK (2024), NNMT CAF (*Nature* 2025, PMID 40702186) | Specific CAF/immune circuits drive progression and immunosuppression | Which circuits are ecotype-defining and reproducible across cohorts |
| Molecular taxonomy | CMS1-4, CRIS A-D, PDS1-3; review PMID 38319359 | Multiple transcriptomic classifications | Instability, TME confounding, lack of cross-platform validation, no clinical adoption |
| Microbiome | Multi-cohort metagenomics (curatedMetagenomicData), microbial peptides reanalyzed from CPTAC (PMID 40263502) | Species and functions associated with CRC | Causal integration with host multi-omics; confounding; paired metabolomics |

**Take-home:** the novel space is not "another atlas" but **integration, reproducibility and clinical utility** across modalities and cohorts.

---

## 3. Gap analysis (the paper's justification)

**G1. Subtype instability.** CMS is dominated by TME signal; CRIS removes stromal contamination but is less validated; PDS describes a continuum. Assignments vary by platform, normalization and cohort. No classification has demonstrated stable, cross-platform performance with a locked, deployable model.

**G2. Missing single-cell-to-bulk-to-space bridge.** Single-cell atlases and spatial studies are usually interpreted separately from large bulk/proteomic cohorts. There is no systematic mapping of single-cell-derived programs onto proteome, methylation and spatial data in one framework.

**G3. Proteome and phosphoproteome are under-integrated.** CPTAC colon cancer is a rich anchor but has not been integrated with modern single-cell/spatial references and large transcriptomic validation cohorts.

**G4. Clinical utility is asserted, not demonstrated.** Most subtype papers report associations; few provide external calibration, decision-curve benefit, or comparison with stage/MSI in independent cohorts.

**G5. Microbiome-host integration is largely correlative.** Public metagenomes are abundant, but batch-aware, confounder-controlled integration with host transcriptome/proteome and metabolome is rare.

**Positioning.** The study should explicitly frame itself as a **reproducibility and translation study**, not as a discovery atlas. Negative/robustness results (e.g., CMS instability) are part of the contribution.

---

## 4. Central hypothesis and specific aims

**Central hypothesis.**
> CRC contains a limited set of reproducible tumor-immune ecotypes, encoded by coordinated transcriptional, epigenetic and proteomic programs and organized into recurrent spatial niches. These ecotypes are conserved across platforms and cohorts, refine prognosis beyond stage/MSI/CMS/CRIS, and map to actionable therapeutic vulnerabilities.

**Aim 1. Build the harmonized multi-omic resource.**
Integrate TCGA-COAD/READ (RNA, methylation, CNV, mutation, RPPA, miRNA, slides), CPTAC colon (proteome, phosphoproteome), large GEO bulk cohorts, single-cell/single-nucleus atlases, and spatial transcriptomics into a versioned, patient-level resource with a confounder/batch registry.

**Aim 2. Define and stress-test ecotypes.**
Derive single-cell-anchored programs and multi-view factors; cluster into ecotypes; quantify stability (bootstrap ARI), cross-platform concordance, and reclassification relative to CMS/CRIS/PDS.

**Aim 3. Validate ecotypes in space.**
Map ecotypes to spatial spots/sections; identify recurrent niches (TLS, invasive front, hypoxic/immune-excluded, CAF-rich); test replication across Visium, Stereo-seq and Xenium-like data.

**Aim 4. Demonstrate clinical and therapeutic relevance.**
Test prognostic and therapy-response associations in independent cohorts with calibration and decision-curve analysis; build an ecotype-drug vulnerability map; release a locked classifier and portal.

**Aim 5. Anchor ecotypes in proteogenomics and map progression/metastasis.**
Use PXD079707 (1,000 treatment-naive CRCs), CPTAC colon, PXD052391 (148 patients, 435 progression samples) and IPX0007391001 (34 patients, 102 primary/normal/liver-metastasis samples) to test whether ecotypes are protein-level states, whether they change along the adenoma-carcinoma-metastasis axis, and whether primary ecotypes predict metastasis niches.

**Aim 6. Link morphology to molecular ecotypes.**
Extract pathology foundation-model features (CHIEF/MUSK/UNI/CONCH) from TCGA diagnostic slides and PXD079707 WSIs; train attention-based multiple-instance models to predict ecotypes and prognosis; validate on external slide sets.

**Aim 7. Derive and validate circulating ecotypes.**
Use cfDNA methylation/5hmC and EV/RNA cohorts (e.g., GSE263030, GSE97923, GSE188627, GSE188433, GSE315715) to test whether tissue ecotypes leave detectable circulating signatures and whether these track stage, MRD and outcome.

**Aim 8. Map ecotype-therapy associations.**
Test ecotype associations with anti-EGFR resistance (GSE140973), chemotherapy (GSE103479, GSE106584, GSE156915), total neoadjuvant therapy (GSE339362), METIMMOX (GSE313035) and MMRd ICI contexts (GSE333300, GSE330748). Label all as exploratory.

**Aim 9. Benchmark and release a reproducible resource.**
Quantify CMS/CRIS/PDS instability across ~20 cohorts and 3-4 platforms; benchmark integration, deconvolution and niche methods; release the locked classifier, benchmark, portal, code and manifests.

**Go/no-go criteria.**
Core:
- Ecotypes must be stable (bootstrap ARI > 0.70) and reproducible in at least two independent cohorts/platforms.
- At least one ecotype must add prognostic information beyond stage + MSI + CMS/CRIS (likelihood-ratio or DCA benefit).
- Spatial ecotype-niche coupling must replicate in at least two spatial datasets.
Extensions (independent; failure does not kill the core):
- Proteome: ecotype program scores correlate with protein/phosphoprotein modules in PXD079707/CPTAC (pre-specified positive concordance).
- Progression/metastasis: at least one ecotype shows a stage/metastasis trend replicated in >=2 cohorts.
- Morphology: ecotype prediction AUC >0.65 in an external slide set.
- Circulating: at least one circulating feature associates with tissue ecotype and outcome in >=2 cohorts.
- Therapy: at least one pre-specified ecotype-therapy association replicates; otherwise report null.
If core criteria fail, pivot to Track D (benchmark/resource paper) or Track C (microbiome companion).

---

## 5. Study design and dataset matrix

### 5.1 Design
- **Discovery:** TCGA-COAD/READ + CPTAC colon (multi-omic, survival).
- **Reference:** single-cell/single-nucleus atlases (cell states and programs).
- **Spatial validation:** Visium/Stereo-seq/high-definition spatial series.
- **External validation:** large GEO bulk cohorts + cBioPortal/MSK-IMPACT/GENIE + CPTAC proteome.
- **Orthogonal/optional:** microbiome, metabolomics, imaging.

### 5.2 Dataset matrix (verified accessions and counts as of 2026-11-17)

| Role | Dataset | Content | n | Notes |
|---|---|---|---|---|
| Discovery | TCGA-COAD / TCGA-READ (GDC) | RNA-seq, 450K methylation, CNV, WXS/WGS, miRNA, RPPA, H&E slides, clinical | 461 COAD + 172 READ cases; 701 gene-expression files, 2,202 methylation files, 495 RPPA, 41 ATAC-seq, 1,972 slide-image files (625 diagnostic-slide files) | Primary multi-omic discovery; survival |
| Discovery (proteome) | CPTAC colon cancer (Vasaikar 2019; PDC000116 - confirm on PDC) | Proteome + phosphoproteome, matched normal | 110 prospectively collected colon tumors | Proteome/kinase anchor |
| Discovery (proteome, 1,000) | PXD079707 (PRIDE; CC0) | DIA proteomics, targeted genomics, whole-slide histology, clinical follow-up | 1,000 treatment-naive CRCs | Proteomic anchor; AI histology; confirm processed tables |
| Progression (proteome) | PXD052391 / IPX0002188000 | Proteome + phosphoproteome, nontumor -> intraepithelial neoplasia -> infiltration -> advanced | 148 patients, 435 samples | WES restricted (HRA007452) |
| Metastasis (proteome) | IPX0007391001 / PRJCA020890 | Genomics, transcriptome, proteome, phosphoproteome; primary/normal/liver metastasis | 34 patients, 102 samples | CRLM axis |
| Reference atlas | GSE132465 | scRNA-seq, 23 Korean patients | 23 patients, 63,689 cells | Core cell-state reference |
| Reference atlas | GSE144735 | scRNA-seq, 6 Belgian patients, core/border/normal | 6 patients, 27,414 cells | Regional heterogeneity |
| Reference atlas | GSE178341 | scRNA-seq MMRd/MMRp atlas | 181 samples | Immune phenotype reference |
| Reference atlas | GSE166555 | scRNA-seq normal/adjacent tumor | 25 samples | Normal/adjacent baseline |
| Reference atlas | GSE188711 | scRNA-seq left- vs right-sided CRC | 6 samples | Sidedness |
| Reference atlas | GSE201348 | scRNA + scATAC, polyp-adenoma-CRC continuum | 72 samples | Progression and chromatin |
| Reference atlas | GSE231559 | scRNA-seq, MD Anderson CRC Moon Shot | 26 samples | Independent reference |
| Reference atlas | GSE312983 | snRNA-seq across age continuum | 189 tumors | EOCRC |
| Reference atlas | GSE285117 | Neutrophil-focused series from the 4.27M-cell study | series-level | Myeloid states |
| Spatial | GSE225857 | scRNA + spatial, primary CRC and liver metastases | 8 samples | Metastasis niches |
| Spatial | GSE283052 | Spatial transcriptomics of CRC tissues | 12 samples | Primary tumors |
| Spatial | GSE320042 | Spatial ecotypes, high-resolution | 38 samples | Direct ecotype comparison |
| Spatial | GSE284989 | Spatial genomics, cholesterol/ICI resistance | 22 samples | Therapy resistance |
| Spatial | GSE280311 / GSE280312 | High-definition spatial immune profiling | 12 / 16 samples | Immune niches |
| Spatial | GSE270745 | Spatial atlas, DMBT1/dysplasia | 6 samples | Early lesions |
| Spatial | GSE217414 | Liver metastases from CRC | 4 samples | Metastasis |
| Validation (bulk) | GSE39582 | Affymetrix colon cancer cohort | 585 samples | Large survival validation |
| Validation (bulk) | GSE17536 | Moffitt colon cancer | 177 samples | Survival |
| Validation (bulk) | GSE14333 | 290 primary CRCs | 290 samples | Survival |
| Validation (bulk) | GSE38832 | Colon cancer | 122 samples | Survival |
| Validation (bulk) | GSE87211 | Rectal cancer +/- controls | 363 samples | Rectal validation |
| Validation (bulk) | GSE103479 | Stage II/III treated vs untreated | 156 samples | Adjuvant chemotherapy |
| Validation (bulk, new) | GSE241101 | Multiregional primary CRC | 286 samples | Intratumor heterogeneity |
| Validation (bulk, new) | GSE159216 | CRC liver metastases | 283 samples | Metastasis |
| Validation (bulk, new) | GSE248381 | Slow-cycling phenotype | 258 samples | Subtype validation |
| Validation (bulk, new) | GSE161158 | SMAD4-modulated profile | 250 samples | Stage II/III DFS |
| Validation (bulk, new) | GSE233616 | Array + methylation subtypes | 188 samples | Epigenome/subtyping |
| Validation (bulk, new) | GSE140973 | Anti-EGFR secondary resistance | 114 samples | Therapy resistance |
| Validation (bulk, new) | GSE333907 | Early-onset CRC cohort | 83 samples | EOCRC |
| Validation (bulk, new) | GSE313035 | METIMMOX | 67 samples | Immunotherapy |
| Validation (bulk, new) | GSE339362 | Total neoadjuvant therapy response | 56 samples | Therapy |
| Circulating extension | GSE263030 | cfDNA 5hmC, occult CRC | 602 samples | Early detection |
| Circulating extension | GSE97923 | cfDNA methylation, COAD | 773 samples | Circulating methylation |
| Circulating extension | GSE188627 | EV-miRNA, metastatic CRC | 290 samples | EV signature |
| Circulating extension | GSE188433 | Plasma EV long RNA, CRC/adenoma | 123 samples | Early detection |
| Circulating extension | GSE315715 | Circulating mRNA | 56 samples | Early detection |
| Circulating extension | GSE216561 | Plasma cfRNA human/microbial | 113 samples | Multi-omic plasma |
| Pathology extension | PXD079707 WSIs + TCGA diagnostic slides | H&E matched to proteome/omics | 1,000 + 625 | Morphology-molecular |
| Validation (clinical-genomic) | cBioPortal: MSK-IMPACT / AACR GENIE / TCGA PanCancer | Mutations, CNV, limited outcomes | cohort-dependent | Metastatic/real-world context |
| Optional (microbiome) | curatedMetagenomicData (Bioconductor) | Uniform taxonomic/functional profiles from published CRC case-control cohorts | dozens of cohorts | Batch-aware meta-analysis only |
| Optional (metabolome) | Metabolomics Workbench / MetaboLights | Sparse human CRC tissue/serum metabolomics | study-dependent | Targeted validation, not discovery |
| Optional (imaging) | TCGA diagnostic slides; TCIA/IDC; Human Protein Atlas | H&E, IHC | 625 diagnostic slides | Morphology/IHC orthogonal check |
| Normal reference | HCA Gut Cell Atlas | Normal colon/rectum cell atlas | public | Baseline signatures |

**Overlap strategy.** Patient-level overlap across GEO series must be checked (same institution/cohort re-used). Deduplicate before validation; never let a validation cohort share patients with discovery. Full inventory, additional cohorts and access notes are in `DATASETS_EXTENDED.md`.

### 5.3 Sample-size logic
- Discovery: ~600 tumors with ≥150 death events (TCGA-COAD/READ is adequate for ecotype discovery; PXD079707 provides proteome in 1,000 and CPTAC in 110).
- Validation: ≥1,000 tumors across ≥4 independent cohorts; ≥300 for proteome validation; ≥148 patients/435 samples for progression; ≥34 patients/102 samples for metastasis.
- Spatial: ≥40 sections across ≥3 platforms for niche replication.
- Power for ecotype-prognosis: detect HR ~1.4 per ecotype at 80% power with ~150 events; pre-specify.

---

## 6. Analysis plan (modules M0-M18)

### M0. Data governance and harmonization
- Version-pinned manifest: accession, download date, checksum, license/DUA.
- Sample sheet with ontology: tissue (tumor/normal/metastasis), stage, MSI/MMR, sidedness, age, sex, treatment, batch, purity.
- Harmonize gene identifiers (GENCODE), re-annotate probes, parse TCGA barcodes (01 tumor, 11 normal; exclude 06 metastatic unless analyzing metastasis).
- Register confounders and batch variables before analysis; lock test sets.
- **QC gate:** no patient overlap between discovery and validation; all sample-level metadata traceable.

### M1. Bulk preprocessing and QC
- **RNA-seq:** STAR/Salmon quantification, TPM/log2, gene filtering, PCA/outlier detection, `ComBat-seq` or `limma::removeBatchEffect` for cohort effects; purity estimates (ABSOLUTE/ESTIMATE).
- **Methylation:** IDAT -> beta (minfi), BMIQ normalization; DMP/DMR (limma/bumphunter); CIMP classification; M-values for integration.
- **Proteome/phosphoproteome:** log2, median/quantile normalization, missing-value filtering (MNAR-aware), phosphosite normalization to protein abundance, batch correction. Sources: PXD079707 (1,000), CPTAC (110), PXD052391 (148 patients/435 samples), IPX0007391001 (34 patients/102 samples), PXD019103.
- **CNV/mutation:** ASCAT/CNVkit; MAF aggregation (maftools); MSI calling; driver annotation.
- **QC gates:** cross-platform correlation (RNA-protein for shared genes), negative-control probes, batch-PCA, outlier audit.

### M2. Single-cell reference construction
- QC: doublet removal, ambient RNA correction, mitochondrial/ribosomal filters.
- Integration: `scVI`/`scANVI` (batch-aware), Harmony as sensitivity.
- Annotation: canonical markers + `CellTypist`; manual curation; malignant-cell inference (`inferCNV`, `CopyKAT`).
- Program discovery: NMF/cNMF on epithelial, fibroblast, myeloid, T/NK, B/plasma, endothelial, glial/neuronal cells.
- Chromatin: scATAC (GSE201348) with ArchR/Signac; gene activity; regulons (SCENIC/pySCENIC).
- **Outputs:** cell-type signatures, malignant meta-programs, CAF/immune states, TLS signature.

### M3. Deconvolution and program scoring
- Bulk: CIBERSORTx, EPIC, quanTIseq, MCP-counter, xCell; consensus across methods (rank-based).
- Spatial: cell2location, RCTD, SpatialDecon; benchmark on synthetic mixtures/simulated spots.
- Program scores: ssGSEA/GSVA, AUCell/UCell; purity-adjusted where appropriate.
- **QC gate:** deconvolution concordance across methods (Spearman > 0.7 for major compartments) and against pathology estimates where available.

### M4. Multi-view integration and ecotype discovery
- Views: RNA (program scores + top variable genes), methylation (DMP/CIMP), proteome, phosphoproteome, CNV, RPPA.
- Methods: MOFA+ (primary), iClusterBayes, SNF, DIABLO (sensitivity).
- Clustering: consensus clustering over K=2-8; stability selection; bootstrap ARI; silhouette.
- Comparison: NMI/ARI with CMS, CRIS, PDS; reclassification tables; survival of discordant cases.
- Ecotype definition: tumor-intrinsic program (e.g., stem/regenerative, EMT, metabolic, immune-secretory) x TME state (immune-hot/cold, CAF-rich, TLS-high).
- **Go/no-go:** stable ecotypes; replication in >=2 cohorts; interpretable biology.

### M5. Biological interpretation
- GSEA/GSVA (Hallmark, Reactome, KEGG).
- Transcription-factor activity: decoupleR (DoRothEA, PROGENy).
- Kinase activity: KSEA/PhosphoSitePlus, phosphosite enrichment.
- Metabolism: scFEA/Compass; expression-based flux proxies.
- Subgroup biology: MMRd vs MMRp, EOCRC vs late-onset, right vs left, primary vs metastasis.
- Optional: neoantigen/immunopeptidome overlay; microbial peptide signals.
- **Output:** ecotype biology atlas, marker panels, pathway/kinase hypotheses.

### M6. Spatial validation
- Map ecotypes to spots/sections (cell2location/RCTD), then:
  - neighborhood enrichment and niche detection (`squidpy`, `CellCharter`);
  - ligand-receptor and signaling (`CellChat`, `NicheNet`);
  - TLS score, invasive-front score, hypoxia, immune exclusion;
  - colocalization of ecotype-defining programs with spatial domains.
- Replicate across Visium, Stereo-seq and high-definition spatial data.
- **Go/no-go:** ecotype-specific niches replicate in >=2 datasets; marker spatial patterns match bulk predictions.

### M7. Clinical utility
- Endpoints: OS, cancer-specific survival, RFS; response to ICI, anti-EGFR (RAS/BRAF WT), adjuvant chemotherapy.
- Models: penalized Cox, random survival forests, gradient boosting; nested cross-validation.
- Metrics: time-dependent AUC, C-index, calibration, decision-curve analysis, NRI/IDI versus stage + MSI + CMS/CRIS.
- External validation in GEO/CPTAC/cBioPortal; subgroup analyses (stage II/III, MSI, EOCRC, sidedness).
- Report per TRIPOD+AI; avoid overclaiming from retrospective public data.
- **Go/no-go:** incremental prognostic value in at least one independent cohort; calibrated risk groups.

### M8. Therapeutic vulnerability and drug repurposing
- Phospho-kinase dependencies; pathway activation scores.
- CMap/LINCS signatures; GDSC/CTRP/PRISM sensitivity; DepMap CRISPR dependencies.
- Build ecotype-drug association matrix; prioritize 3-5 testable hypotheses.
- **Framing:** hypothesis-generating; explicitly state that experimental validation is required.

### M9. Microbiome/metabolome integration (optional companion)
- Profiling: MetaPhlAn 4, HUMAnN 3; batch-aware CLR; MMUPHin; ANCOM-BC2.
- Multi-cohort meta-analysis with random effects; confounder adjustment (BMI, diet, antibiotics, geography, batch).
- Host-microbe integration: matched cohorts where available; otherwise ecological correlation with explicit limitations.
- Metabolite prediction (MelonnPan) and validation in public metabolomics.
- **Decision:** if signals are robust and survive confounder sensitivity, develop as a companion paper; do not dilute Track A.

### M10. Reproducibility and resource release
- Workflow: Nextflow/Snakemake; containers (Docker/Singularity); `renv`/conda lockfiles.
- Code and intermediate results on GitHub + Zenodo DOI.
- Interactive portal (Shiny/Cellxgene) for ecotype scoring and exploration.
- Document versions, parameters, seeds, and negative results.

### M11. Proteogenomic anchoring and progression/metastasis
- PXD079707: process DIA proteomics (Spectronaut/DIA-NN), protein-level normalization, gene mapping; score ecotype programs at the protein level; test concordance with RNA-based ecotypes.
- CPTAC: independent proteome/phosphoproteome validation; kinase activity.
- PXD052391: place samples along nontumor -> intraepithelial neoplasia -> infiltration -> advanced; test ecotype program trajectories and stage-specific vulnerabilities.
- IPX0007391001: matched primary/normal/liver metastasis; test ecotype retention vs switch, niche remodeling and the C1/C2 proteomic subtype relationship.
- **Outputs:** protein-level ecotype scores, progression/metastasis maps; Figure 8.

### M12. Pathology foundation-model morphology
- WSI QC/tiling (OpenSlide, TIAToolbox, CLAM); tissue detection.
- Feature extraction with CHIEF, MUSK, UNI and CONCH (frozen features; benchmark).
- Models: attention multiple-instance learning, linear probe, penalized Cox; nested CV; external slide sets (TCGA vs PXD079707).
- Test ecotype prediction, prognostic value and correlation with transcriptomic/spatial ecotypes.
- **Outputs:** morphology-ecotype model; Figure 9a.

### M13. Circulating ecotype
- cfDNA methylation/5hmC: normalize; deconvolve tissue-of-origin using reference methylomes; compute ecotype scores.
- EV/RNA: normalize; derive and score signatures.
- Test association with tissue ecotypes (matched/semi-matched), stage, MRD and outcome; pre-specify features.
- **Outputs:** circulating ecotype signature; Figure 9b.

### M14. Therapy-response mapping
- Cohorts: GSE140973 (anti-EGFR), GSE103479/GSE106584/GSE156915 (chemo), GSE339362 (TNT), GSE313035 (METIMMOX), GSE333300/GSE330748 (MMRd ICI).
- Pre-specify response definitions; test ecotype association; adjust for stage, MSI, RAS/BRAF; label exploratory.
- **Outputs:** therapy-ecotype map; Figure 6b/Table 3.

### M15. Cross-platform and method benchmark
- Subtyping: assign CMS/CRIS/PDS across cohorts; ARI/NMI; survival of discordant cases.
- Integration: MOFA+ vs iClusterBayes vs SNF vs DIABLO; stability, runtime, external replication.
- Deconvolution: cell2location vs RCTD vs BayesPrism/InstaPrism; synthetic spots.
- Niche: CellCharter vs BANKSY vs CellNiche.
- Survival: Cox vs RSF vs XGBoost/DeepSurv; calibration/DCA.
- **Outputs:** benchmark tables/figures; Figure 10.

### M16. Causal triangulation (optional)
- MR/colocalization for ecotype driver genes using eQTL/pQTL resources (GTEx, eQTL Catalogue, UKB-PPP, deCODE).
- Mediation analyses only where assumptions are defensible.
- **Outputs:** exploratory causal evidence table; no causal claims without robust instruments.

### M17. Microbiome and microbial peptides (optional)
- curatedMetagenomicData; MetaPhlAn 4/HUMAnN 3; MMUPHin; ANCOM-BC2; MelonnPan.
- Microbial peptides reanalyzed from CPTAC/PXD079707 proteomes; plasma microbial cfRNA (GSE174302, GSE216561, GSE264208).
- Confounder adjustment (BMI, diet, antibiotics, geography); no causal claims.
- **Outputs:** companion analysis or supplementary module.

### M18. Extended resource and portal
- Locked classifier, benchmark results, portal, code and manifests.
- Interactive scoring for user expression/proteome data.
- **Outputs:** Zenodo release, portal, reproducible workflow.

---

## 7. Statistical rigor and reporting

**Pre-registration (OSF).** Primary endpoint, ecotype K range, validation cohorts, metrics, and subgroup analyses. This is a strong differentiator for a public-data paper.

**Data hygiene.**
- Patient-level splits only; lock the test set before modeling.
- Deduplicate cohorts and check for sample overlap.
- No feature selection on the full dataset; nest all tuning.

**Statistics.**
- FDR (Benjamini-Hochberg) for omics; permutation tests for clustering stability.
- Bootstrap confidence intervals; report effect sizes with CIs.
- Sensitivity analyses: purity adjustment, batch correction method, platform, MSI status, stage, age, sidedness.
- Missing data: multiple imputation for clinical variables; report missingness.
- Calibration and decision curves, not just AUC.

**Reporting standards.**
- TRIPOD+AI (prediction models), REMARK (prognostic markers), MINSEQE/MIAME-like metadata for omics.
- Full code, manifests, and session info.

**Claim discipline.**
- Observational public data support association, not causality.
- Avoid "drives", "causes", "validated clinically" unless prospectively tested.
- Use "reproducible", "externally validated association", "hypothesis-generating vulnerability".

---

## 8. Figures and manuscript outline

### Main figures (10)
1. **Study design and data landscape:** dataset matrix, sample overlap, modality coverage, analysis workflow, proteogenomic/progression/metastasis resources.
2. **Single-cell reference atlas:** cell states, malignant programs, CAF/immune states, chromatin programs.
3. **Multi-view integration and ecotype discovery:** MOFA+ factors, consensus clustering, stability, comparison/reclassification vs CMS/CRIS/PDS.
4. **Ecotype biology:** pathways, TFs, kinases, metabolism; MMRd/MSS, EOCRC, sidedness; marker panels.
5. **Spatial validation:** ecotype mapping, niches, TLS/invasive front, ligand-receptor, cross-platform replication.
6. **Clinical utility and therapy:** survival, external validation, calibration, decision curves, therapy-response associations.
7. **Therapeutic vulnerabilities and resource:** ecotype-drug map, kinase dependencies, classifier/portal schematic.
8. **Progression and metastasis:** nontumor -> adenoma -> carcinoma -> liver metastasis (PXD052391, IPX0007391001, GSE201348, GSE225857).
9. **Morphology and circulating ecotypes:** WSI foundation-model features, morphology-ecotype prediction, cfDNA/EV circulating signatures.
10. **Benchmark and reproducibility:** CMS/CRIS/PDS instability, integration/deconvolution/niche/survival method comparison, portal/resource.

**Supplementary:** QC, deconvolution benchmarks, cohort table, sensitivity analyses, all marker/signature lists, negative results, data/code availability.

### Manuscript outline
- Abstract (250 words): problem, approach, key findings, clinical implication, availability.
- Introduction (~800 words): CRC burden, heterogeneity, limits of CMS/CRIS/PDS, why cross-platform ecotypes.
- Results (5 subsections mirroring Aims 1-4 + resource).
- Discussion (~1,000 words): interpretation, comparison with prior atlases, limitations (public data, retrospective, no wet lab), future validation.
- Methods (~3,000 words): data sources, QC, integration, statistics, reproducibility.
- Data/code availability; author contributions; competing interests.

### Title options
1. *A cross-platform ecotype taxonomy of colorectal cancer resolves subtype instability and maps therapeutic vulnerabilities.*
2. *Single-cell-anchored, spatially validated ecotypes of colorectal cancer outperform bulk molecular subtypes for prognostication.*
3. *Integrative multi-omic and spatial atlas of colorectal cancer defines reproducible tumor-immune ecotypes with clinical utility.*

---

## 9. Novelty and positioning

**What is new (v2.0).**
1. **Proteogenomic anchor at population scale:** first ecotype taxonomy anchored in the 1,000-CRC PXD079707 resource, CPTAC, the 148-patient progression cohort and the 34-patient CRLM cohort.
2. **Cross-platform reproducibility:** the same ecotypes are tested in RNA, protein, phosphoprotein, methylation and spatial data, not just one modality.
3. **Single-cell anchoring:** bulk subtypes are re-derived from cell-state programs, reducing TME confounding.
4. **Explicit reproducibility benchmark:** CMS/CRIS/PDS instability is quantified across ~20 cohorts and 3-4 platforms, with a locked classifier.
5. **Multi-platform spatial meta-atlas:** ecotype-associated niches are replicated across Visium, Stereo-seq, Xenium and CosMx-like data.
6. **Morphology-molecular ecotypes:** pathology foundation models link H&E to ecotypes and prognosis.
7. **Circulating ecotypes:** cfDNA/EV signatures link tissue ecotypes to non-invasive monitoring.
8. **Therapy-response map:** ecotype associations across anti-EGFR, chemotherapy, TNT and ICI contexts.
9. **Progression/metastasis axis:** adenoma-carcinoma-metastasis ecotype trajectory.
10. **Deployable resource:** locked classifier, benchmark, portal, code and manifests.

**What not to claim.**
- Do not claim a new "atlas" as the main novelty; atlases are abundant.
- Do not claim causality or clinical readiness from retrospective public data.
- Do not claim pan-cancer generality without pan-cancer validation.
- Do not compete with PXD079707 on "proteomic subtypes"; build on it.

**Differentiation from closest work.**
- vs. PXD079707: that is a proteomic resource with AI histology; this study adds single-cell anchoring, cross-platform reproducibility, spatial niche validation, circulating signatures, therapy mapping and a locked classifier.
- vs. 4.27M-cell atlas: that study is immune-centric and descriptive; this study is cross-modality, taxonomy-focused and clinically validated.
- vs. CPTAC: CPTAC has 110 tumors; this study uses 1,000 primary tumors plus progression and CRLM cohorts.
- vs. CMS/CRIS/PDS: this study explicitly tests and addresses their instability and clinical utility gap.
- vs. single spatial studies: those are single-platform and descriptive; this study is a multi-platform meta-analysis with replication.
- vs. pathology foundation-model papers: those focus on pan-cancer diagnosis/prognosis; this study is ecotype-aware CRC morphology.

---

## 10. Target journals and publication strategy

| Tier | Journals | Requirements / fit |
|---|---|---|
| Selected target | *Journal of Translational Medicine*, *Computational and Structural Biotechnology Journal*, *Briefings in Bioinformatics* | Computational-only, public-data study; each venue emphasizes a different strength (clinical translation / framework / resource-benchmark) |
| Stretch (no new data) | *Genome Medicine*, *eBioMedicine*, *BMC Medicine* | Only with exceptionally strong external validation and a unique resource |
| Requires experimental validation | *Gut*, *Nature Communications*, *Cancer Cell*, *Nature Cancer* | Not feasible without wet-lab collaboration |
| Alternative framing | *Genomics Proteomics Bioinformatics*, *Database*, *Molecular Oncology* | If the paper becomes a benchmark/resource or a focused translational study |

**Strategy.**
- Post a bioRxiv/medRxiv preprint at submission to establish priority.
- Tailor the title/abstract to the venue: clinical translation (*J Transl Med*), computational framework/benchmarks (*CSBJ*), or reproducible resource (*Briefings in Bioinformatics*).
- Prepare a cover letter emphasizing reproducibility, external validation and clinical decision impact, not "another atlas".
- Include a limitations paragraph that pre-empts "public data only" criticism and states that no wet-lab validation was performed.
- If reviewers request wet-lab validation, have a collaboration plan ready (organoids/CODEX/IHC on a small cohort).

---

## 11. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Batch/cohort confounding dominates | High | Harmonization, ComBat/MMUPHin, negative controls, sensitivity analyses, report uncorrected results |
| Ecotypes unstable/overfit | Medium | Bootstrap stability, locked test sets, external replication, pre-registered K range |
| Weak incremental clinical value | Medium | Focus on stage II/III and MSI subgroups; use DCA/NRI; report null results honestly |
| Proteome sample size limited | High | Use CPTAC as anchor/validation, not discovery; emphasize cross-platform concordance |
| Spatial data small/heterogeneous | High | Use spatial for validation and niche discovery; require replication in >=2 datasets |
| Microbiome confounding | High | Keep as optional companion; confounder-adjusted meta-analysis; no causal claims |
| Scooped by a larger atlas | Medium | Preprint early; differentiate on reproducibility/clinical utility, not cell count |
| Data access/DUA delays | Medium | Start with open data; list gated sources separately; no controlled data required for the core |
| Reviewer demands wet lab | High for top tier | Target computational-friendly journals; prepare collaboration plan |
| PXD079707 is a new/unpublished resource | Medium | Verify file completeness, metadata and license; use CPTAC as fallback; do not rely solely on it |
| PXD052391 WES is restricted | Low | Use open proteome/phosphoproteome; request WES only if essential |
| Pathology FM compute/domain shift | Medium | Use frozen features; benchmark CHIEF/MUSK/UNI/CONCH; external slide sets |
| Liquid biopsy analyte heterogeneity | High | Pre-specify analyte-specific features; no pooling; external validation |
| New-cohort overlap | Medium | Patient-level deduplication; check institution/cohort |
| Scope creep from extensions | High | Complete core Aims 1-4 first; treat extensions as pre-specified optional modules |

---

## 12. Timeline, compute and team

### 12-month timeline
| Month | Milestone |
|---|---|
| 1 | Data acquisition, manifest, QC, pre-registration |
| 2 | Bulk harmonization; single-cell reference atlas |
| 3 | Deconvolution; multi-view integration |
| 4 | Ecotype discovery and stability |
| 5 | Biological interpretation |
| 6 | Spatial validation |
| 7 | Clinical validation |
| 8 | Drug repurposing and sensitivity analyses |
| 9 | Figures; classifier packaging |
| 10 | Portal; reproducibility audit |
| 11 | Manuscript drafting |
| 12 | Submission and revision |

### Extension timeline (months 13-18; optional, separate paper or revision)
| Month | Milestone |
|---|---|
| 13-14 | Proteogenomic anchoring; progression/metastasis axis |
| 14-15 | Pathology foundation-model morphology |
| 15-16 | Circulating ecotype |
| 15-17 | Therapy-response mapping; benchmark |
| 17-18 | Extended manuscript/revision |

**Strategy note:** submit the core paper at month 12 if the core criteria pass. Do not delay the core paper for extensions; use them as a follow-up or a revision strength.

### Compute
- **Disk:** 1.6 TB free locally; plan for 2-4 TB (raw + intermediate + containers). Processed matrices only for most atlases.
- **RAM:** 64 GB minimum; 128 GB recommended for MOFA+/single-cell integration.
- **CPU:** 16-32 cores; GPU optional (scVI, deep learning).
- **Runtime:** bulk integration days; single-cell integration 1-2 weeks; spatial deconvolution days per section.

### Team/skills
- Bioinformatics (multi-omics integration, single-cell, spatial).
- Biostatistics (survival, prediction models, calibration).
- Pathology (spatial annotation, IHC validation).
- Clinical oncology (endpoint definitions, interpretation).
- Data engineering (workflows, portal).

### Ethics and data use
- Public de-identified data; comply with each repository's terms.
- No attempt at re-identification.
- For gated data (dbGaP/EGA/Broad Single Cell Portal), obtain DUAs; do not include controlled data in the core plan.

---

## 13. Alternative tracks (not selected)

**Track B. Adenoma-carcinoma-metastasis axis atlas.**
Focus on progression and EOCRC, using GSE201348, GSE312983, GSE225857, spatial data and TCGA/CPTAC. Higher biological novelty; heavier harmonization; smaller clinical validation.

**Track C. Microbiome-metabolite-host multi-omic integration.**
Use curatedMetagenomicData + public metabolomics + host transcriptome/proteome. High novelty but paired data are scarce and confounding is severe; best as a companion paper.

**Track D. Benchmark/resource paper.**
Benchmark integration/deconvolution/subtyping methods on CRC, release a portal and best-practice guidelines. Lower biological novelty; strong fit for bioinformatics journals; fast and low-risk.

---

## 14. Decisions and immediate next actions

**Decisions (2026-11-17; v2.0).**
- **Track A selected:** cross-platform, single-cell-anchored ecotype taxonomy.
- **Validation:** computational only (public data); no wet-lab or clinical-efficacy claims.
- **Target tier:** *Journal of Translational Medicine*, *Computational and Structural Biotechnology Journal*, *Briefings in Bioinformatics*.
- **v2.0 expansion:** PXD079707 (1,000 CRCs), PXD052391 (progression), IPX0007391001 (CRLM), pathology foundation models, circulating ecotypes, therapy mapping and a methods benchmark.
- **Core-first rule:** complete Aims 1-4 and submit; extensions are optional and must not delay the core paper.

**Next actions.**
1. **Approve the v2.0 protocol** in `PROTOCOL_TrackA.md` and the cohort definition in `analysis/config/study_config.yaml`.
2. **Download Tier 1 manually** from `DATASETS.md`/`DATASETS_EXTENDED.md`, including a PXD079707 subset, PXD052391 and IPX0007391001; record checksums (no automatic downloads).
3. **Run the 3-week pilot** in `PILOT_PLAN.md` on TCGA + GSE39582 + GSE132465 + one GSE225857 section + a PXD079707 subset.
4. **Approve the software list** in `SOFTWARE.md`; install only after confirmation.
5. **Decide whether extensions are in the same paper or a follow-up**; freeze the core scope.
6. **Freeze the config and register the SAP** before outcome modeling.

---

## Appendix. Key references (PubMed IDs)

- 41529665 - 4.27M-cell CRC immune atlas, *Cancer Cell* 2026
- 39216476 - isoform-resolution CRC atlas, *Cell Genomics* 2024
- 42106149 - MSI/MSS immunopeptidomics, *Mol Cell Proteomics* 2026
- 31031003 - CPTAC colon proteogenomics, *Cell* 2019
- 40263502 - microbial peptides from CPTAC CRC proteomics, *Sci Rep* 2025
- 38319359 - molecular pathological classification update, *Virchows Arch* 2024
- 39137726 - CCL19+ fibroblasts/TLS in CRLM, *Cancer Cell* 2024
- 38335276 - Th17/TWEAK EMT in CRLM, *Cancer Res* 2024
- 40702186 - NNMT CAF immunosuppression, *Nature* 2025
- 40826474 - COL10A1+ CAF, *J Exp Clin Cancer Res* 2025
- 40050872 - CTHRC1+ CAF/WNT5A, *J Transl Med* 2025
- 40511612 - ANGPTL4 proteogenomic target, *JNCI* 2025
- 39644773 / 40409349 / 41320324 - microbiome-metabolite examples (context only)
