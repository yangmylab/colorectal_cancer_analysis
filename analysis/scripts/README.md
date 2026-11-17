# Planned analysis modules

**Status:** interfaces only. No script has been written or executed yet. Each module will be implemented after protocol approval and the pilot.

---

## Naming and interface rules

```
NN_<module>.<ext>          # entry point
NN_<module>_utils.<ext>    # helper functions
```

Every entry point:
1. reads `../config/study_config.yaml`;
2. takes explicit input/output paths as arguments (or derives them from config);
3. writes logs to `../logs/`;
4. writes a version/session file;
5. is idempotent;
6. never writes to `../data/raw/`.

Languages: R for bulk/methylation/survival/MOFA/ecotypes; Python for single-cell, spatial, deep learning, and portal back end. Mixed modules exchange `.h5ad`, `.rds`, `.tsv`, `.parquet`.

---

## Module list

| Script | Language | Input | Output | Purpose |
|---|---|---|---|---|
| `00_manifest.R` | R | raw files + metadata | `manifests/manifest.tsv`, checksums | Inventory, checksums, sample sheet validation |
| `01_bulk_rnaseq_qc.R` | R | TCGA/GEO RNA | `processed/rna/`, QC report | Normalization, filtering, batch, purity |
| `02_methylation_qc.R` | R | TCGA IDAT/beta | `processed/methylation/` | BMIQ, DMP/DMR, CIMP, M-values |
| `03_proteome_qc.R` | R | CPTAC tables | `processed/proteome/`, `processed/phospho/` | Normalization, site-level processing |
| `04_cnv_mutation.R` | R | TCGA/GEO CNV/MAF | `processed/cnv/`, `processed/mutation/` | CNV burden, drivers, MSI |
| `05_reference_atlas.py` | Python | scRNA/snRNA/scATAC | `reference/atlas/` | QC, integration, annotation, malignant inference |
| `06_programs.R` | R | reference atlas | `reference/programs/`, `reference/signatures/` | NMF programs, signatures, marker panels |
| `07_deconvolution.R` | R | bulk + signatures | `results/deconvolution/` | CIBERSORTx/EPIC/quanTIseq/MCP-counter/xCell consensus |
| `08_program_scoring.R` | R | bulk + signatures | `results/programs/` | ssGSEA/GSVA/AUCell scores |
| `09_mofa_integration.R` | R | views | `results/integration/mofa/` | MOFA+, factor stability, variance explained |
| `10_ecotypes.R` | R | factors/programs | `results/ecotypes/`, `models/` | Consensus clustering, stability, classifier |
| `11_interpretation.R` | R | ecotypes + omics | `results/biology/` | GSEA, TF/kinase activity, metabolism |
| `12_spatial_validation.py` | Python | spatial data + signatures | `results/spatial/` | Deconvolution, niches, neighborhoods, LR |
| `13_clinical_models.R` | R | ecotypes + clinical | `results/clinical/` | Cox/RSF, calibration, DCA, external validation |
| `14_drug_repurposing.R` | R | ecotype signatures | `results/drug/` | CMap/LINCS, GDSC/CTRP/PRISM, DepMap |
| `15_portal/` | R/Python | results/models | `portal/` | Shiny/Cellxgene app |
| `16_figures.R` | R | results | `figures/` | Publication figures |
| `17_proteogenomics.R` | R | PXD079707/CPTAC/PXD052391/IPX0007391001 | `processed/proteome/`, `results/proteogenomics/` | Protein-level ecotype scores; progression/metastasis |
| `18_progression_metastasis.R` | R | progression/CRLM cohorts | `results/progression/` | Stage trajectory; primary-metastasis mapping |
| `19_pathology_fm.py` | Python | WSIs + slide metadata | `results/pathology/` | CHIEF/MUSK/UNI/CONCH features; attention MIL |
| `20_circulating.py` | Python | cfDNA/EV cohorts | `results/circulating/` | Circulating ecotype scores |
| `21_therapy_response.R` | R | therapy cohorts | `results/therapy/` | Ecotype-therapy associations |
| `22_benchmark.R` | R | all results | `results/benchmark/` | CMS/CRIS/PDS stability; method comparisons |
| `23_causal.R` | R | QTL/GWAS summary stats | `results/causal/` | Optional MR/colocalization |
| `24_microbiome.R` | R | metagenomes/cfRNA | `results/microbiome/` | Optional microbiome/microbial peptides |
| `workflow/main.nf` | Nextflow | all | pipeline | Orchestration and containers |

---

## Module contracts

### `00_manifest.R`
- Validates required columns, unique IDs, event coding, date formats.
- Computes MD5 for every raw file.
- Fails loudly on duplicate sample IDs or missing patient IDs.

### `01_bulk_rnaseq_qc.R`
- Input: counts/TPM + sample sheet.
- Output: normalized matrix, PCA, outlier table, purity table.
- Must report pre/post-batch PCA and cross-platform correlations.

### `05_reference_atlas.py`
- Input: processed matrices (`.h5ad`) + metadata.
- Output: integrated `.h5ad`, cell-type labels, malignant calls, dataset-mixing metrics.
- Must save the exact package versions and seeds.

### `09_mofa_integration.R`
- Input: one matrix per view, samples in rows.
- Output: MOFA object, factor scores, variance explained, bootstrap stability.
- Must handle missing views explicitly and report which samples use which views.

### `10_ecotypes.R`
- Input: factor scores/program scores.
- Output: labels, consensus matrices, ARI/PAC, classifier object.
- Must not use outcome data for clustering or feature selection.

### `13_clinical_models.R`
- Input: locked ecotype labels + clinical + cohort indicator.
- Output: model objects, C-index, time-AUC, calibration, DCA, NRI/IDI.
- Must use nested CV and a locked external test set.

### `12_spatial_validation.py`
- Input: spatial data + reference signatures.
- Output: deconvolution, ecotype scores, niches, neighborhood enrichment, LR interactions.
- Must record platform, spot size, and normalization.

---

## Minimal pilot implementation

For the 3-week pilot, only these are needed:
`00`, `01` (TCGA + GSE39582), `05`, `06`, `07`, `08`, `09`, `10`, a reduced `12` on one section, and a small `17` (PXD079707 subset, optional).

Everything else is deferred until the pilot passes.

---

## Definition of done (per module)

- [ ] Runs from a clean environment using the config.
- [ ] Produces the documented outputs.
- [ ] Writes a log and version file.
- [ ] Has a small synthetic test or a documented sanity check.
- [ ] Is reviewed by a second person.
- [ ] Is committed with a clear message.
