# Software and environment checklist

**Project:** Colorectal cancer multi-omics strategy
**Prepared:** 2026-11-17
**Status:** NOTHING HAS BEEN INSTALLED. Per project standing instruction, this is the proposed environment only. Confirm before any installation.

---

## 0. Recommended environment

- **OS:** Linux (current machine: linux/amd64)
- **Hardware:** 16-32 CPU cores; 64 GB RAM minimum (128 GB recommended); 2-4 TB free disk; optional NVIDIA GPU for scVI/deep learning.
- **Languages:** R >= 4.4, Python >= 3.11
- **Managers:** `renv` (R), `conda`/`mamba` (Python + CLI), `Nextflow` + `Singularity`/`Docker` for workflows.
- **Version control:** git + GitHub; Zenodo DOI for releases.
- **Reproducibility:** lockfiles (`renv.lock`, `environment.yml`), container images, session info, fixed seeds.

---

## 1. R / Bioconductor

### Core data handling and statistics
- `tidyverse`, `data.table`, `readr`, `janitor`
- `Biobase`, `SummarizedExperiment`, `MultiAssayExperiment`
- `survival`, `survminer`, `glmnet`, `rms`, `timeROC`, `survivalROC`
- `lme4`, `emmeans`, `broom`, `broom.mixed`
- `missForest` / `mice` (missing data)

### Bulk transcriptomics and differential analysis
- `limma`, `edgeR`, `DESeq2`
- `sva` (ComBat/ComBat-seq)
- `GSVA`, `GSEABase`, `fgsea`
- `decoupleR` (TF/pathway activity)
- `ComplexHeatmap`, `EnhancedVolcano`, `ggplot2`, `patchwork`

### Methylation
- `minfi`, `IlluminaHumanMethylation450kanno.ilmn12.hg19` (and EPIC annotation if needed)
- `bumphunter`, `DMRcate`, `wateRmelon` (BMIQ)
- `sesame` (alternative)

### Copy number, mutation, purity
- `maftools`
- `ESTIMATE`
- `ABSOLUTE` / `PureCN` (if WES/WGS available)
- `CNVkit` is Python; R-side plotting via `ComplexHeatmap`

### Deconvolution and immune scoring
- `immunedeconv` (wraps CIBERSORTx/EPIC/quanTIseq/MCP-counter/xCell)
- `EPIC`
- `xCell`
- `IOBR` (optional, convenience wrapper)

### Single-cell
- `Seurat` (v5)
- `Signac` (scATAC)
- `ArchR` (scATAC, alternative)
- `inferCNV`
- `CopyKAT`
- `CellChat`
- `SCENIC` / `SCENIC+` (regulons)
- `monocle3` / `slingshot` (trajectories, optional)

### Multi-omics integration
- `MOFA2`
- `mixOmics` (DIABLO, sPLS)
- `iClusterBayes`
- `SNFtool`
- `ConsensusClusterPlus`
- `coca` (optional)

### Spatial (R-side)
- `Space Ranger` outputs are handled in Python; R-side visualization via `SpatialExperiment`, `spotlight`, `SPOTlight` (deconvolution)
- `SpatialDecon`

### Microbiome (optional)
- `curatedMetagenomicData`
- `MMUPHin`
- `ANCOMBC` (ANCOM-BC2)
- `phyloseq`

### Utilities
- `renv`, `sessioninfo`, `targets` (pipeline), `quarto`

---

## 2. Python

### Core
- `numpy`, `pandas`, `scipy`, `scikit-learn`, `statsmodels`
- `matplotlib`, `seaborn`, `plotly`
- `joblib`, `tqdm`

### Single-cell
- `scanpy`, `anndata`
- `scvi-tools` (scVI, scANVI, totalVI, MultiVI)
- `harmonypy`
- `scarches` (reference mapping, optional)
- `decoupler` (Python version)
- `gseapy`

### Spatial
- `squidpy`
- `cell2location`
- `cellcharter`
- `tangram-sc`
- `spatialdata`, `spatialdata-io`

### Survival / prediction
- `lifelines`
- `scikit-survival`
- `shap`
- `optuna` (hyperparameter tuning)

### Statistics / multiple testing
- `statsmodels` (FDR), `scikit-posthocs`

### Workflow / IO
- `snakemake` or Nextflow (Java/Groovy)
- `pyyaml`, `click`, `papermill` (optional)

---

## 3. Command-line tools

### Data access
- `gdc-client` (TCGA/GDC)
- `sra-tools` (`prefetch`, `fasterq-dump`) for raw SRA/ENA
- `pdc-client` (if available/needed for PDC)
- `awscli` / `gsutil` (for cloud-hosted atlases, optional)

### Bulk RNA-seq (only if reprocessing)
- `fastp` / `trim_galore`
- `STAR` or `salmon`
- `subread` (`featureCounts`)
- `qualimap` / `multiqc`

### Single-cell / spatial (only if reprocessing)
- `cellranger` / `spaceranger` (10x; license required)
- `STARsolo` (alternative)

### Microbiome (optional)
- `MetaPhlAn` 4
- `HUMAnN` 3
- `KneadData` (host depletion)

### Workflow and containers
- `nextflow` or `snakemake`
- `singularity` / `apptainer` or `docker`
- `git`, `quarto`, `pandoc`

---

## 4. Installation plan (do not run until approved)

### 4.1 R
```r
# Preferred: renv per project
install.packages("renv")
renv::init()
# Bioconductor
install.packages("BiocManager")
BiocManager::install(c(
  "SummarizedExperiment", "MultiAssayExperiment",
  "limma", "edgeR", "DESeq2", "sva",
  "minfi", "bumphunter", "DMRcate", "wateRmelon",
  "GSVA", "fgsea", "decoupleR",
  "ComplexHeatmap", "maftools", "ESTIMATE",
  "MOFA2", "mixOmics", "iClusterBayes", "SNFtool",
  "ConsensusClusterPlus", "Seurat", "Signac", "inferCNV",
  "CellChat", "SCENIC", "curatedMetagenomicData", "MMUPHin", "ANCOMBC"
))
# CRAN
install.packages(c(
  "tidyverse", "data.table", "survival", "survminer", "glmnet", "rms",
  "timeROC", "survivalROC", "missForest", "mice", "immunedeconv",
  "EPIC", "xCell", "SPOTlight", "SpatialDecon", "renv", "sessioninfo",
  "targets", "quarto"
))
```

### 4.2 Python (conda/mamba)
```bash
mamba create -n crc_multiomics python=3.11 -y
mamba activate crc_multiomics
pip install numpy pandas scipy scikit-learn statsmodels matplotlib seaborn \
  scanpy anndata scvi-tools harmonypy decoupler gseapy squidpy cell2location \
  cellcharter tangram-sc spatialdata spatialdata-io lifelines scikit-survival \
  shap optuna pyyaml
```

### 4.3 CLI
```bash
# GDC
# Download the official gdc-client binary for linux/amd64 from
# https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
# SRA
mamba install -c bioconda sra-tools -y
# Workflow
mamba install -c bioconda nextflow -y
# Containers (if not already present)
# Install singularity/apptainer per site policy
```

### 4.4 Optional
```bash
# Cellranger/spaceranger require 10x registration and are only needed for raw reprocessing.
# MetaPhlAn/HUMAnN only if the microbiome companion is pursued.
```

---

## 5. Environment validation

After installation, run a minimal smoke test:
```bash
Rscript -e 'library(MOFA2); library(Seurat); library(limma); sessionInfo()'
python -c "import scanpy, scvi, squidpy, cell2location; print('ok')"
gdc-client --version
nextflow -version
```
Record outputs in `metadata/environment_versions.txt`.

---

## 6. Reproducibility rules

- Pin versions in `renv.lock` and `environment.yml`.
- Use container images for each major module (bulk, single-cell, spatial).
- Fix seeds for clustering, dimensionality reduction and cross-validation.
- Store intermediate objects (`.rds`, `.h5ad`, `.zarr`) with checksums.
- Do not upgrade packages mid-analysis; use separate environments per module if needed.
- Archive the final environment with the Zenodo release.

---

## 7. Notes and constraints

- No installation has been performed. Confirm the study track and this list before installing.
- Some tools have license constraints: `cellranger`/`spaceranger` (10x), `CIBERSORTx` (academic registration), `SCENIC`/`pySCENIC` (check license), `GSEA` (registration for desktop version; `fgsea` is open).
- Prefer open-source equivalents when possible: `fgsea` instead of GSEA desktop, `decoupleR` instead of proprietary pathway tools.
- If GPU is unavailable, `scVI` still runs on CPU, but integration of >1M cells may be slow; subsample or use `scArches` reference mapping.
- `MOFA+` memory scales with samples x features; use program scores or top-variable features rather than full transcriptomes.
