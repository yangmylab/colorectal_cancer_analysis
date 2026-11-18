# Environment

## Verified environment (2026-11-18)

- **R:** 4.6.1 (2026-06-24)
- **R packages installed:** 1,754
- **Core R packages available:** Seurat, limma, edgeR, DESeq2, survival, glmnet, rms, mice, NMF, SNFtool, mixOmics, GSVA, fgsea, msigdbr, ComplexHeatmap, pheatmap, data.table, ggplot2, preprocessCore, sva, GEOquery, CellChat, mclust, nnls, riskRegression, prodlim, future, yaml, jsonlite, **minfi, bumphunter, IlluminaHumanMethylation450kanno.ilmn12.hg19, hgu133plus2cdf, decoupleR, KSEAapp, survminer, pdftools, magick**.
- **Python:** 3.11.18 in conda env `sc-analysis` (`/home/meng/miniconda3/envs/sc-analysis/bin/python`): scanpy 1.11.5, anndata 0.12.19, cell2location 0.1.5, scvi-tools 1.4.2, pyro-ppl 1.8.6, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.1, torch 2.15.0 (CPU only). Base `python3` lacks numpy.
- **Hardware:** **6 physical cores / 12 SMT threads** (`nproc` reports 12 — do not pass it to heavy work), 61 GiB RAM, 1.6 TB free on `/run/media/meng/Seqdata`.
- **Figure toolchain:** `pdflatex` with TikZ, `lualatex`, `ghostscript`, `poppler-utils` (`pdftoppm`, `pdfinfo`, `pdftotext`), `pandoc`.

## Optional packages not installed

MOFA2, ConsensusClusterPlus, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC. (`decoupleR`, `minfi` and `survminer` from an earlier version of this list are installed — see above.)

**Fallbacks implemented:** SNF/NMF for multi-view integration, GSVA ssGSEA for signatures, manual consensus clustering, Cox/C-index/calibration/DCA for survival.

## Installation history

Packages were installed during the 2026-11-17 sessions (with the user's approval), so this is no longer a
no-install environment: `minfi` 1.58.0, `bumphunter` 1.54.0, the 450K annotation 0.6.1, `SNFtool` 2.3.1,
`hgu133plus2cdf` 2.18.0, `decoupleR` 2.17.0, `KSEAapp` 2.0; `magick` 2.9.1 and `xml2` 1.6.0 were rebuilt against
the system ImageMagick 7 / libxml2; the conda env `sc-analysis` was recreated and `cell2location` 0.1.5 installed
(it requires `pyro-ppl==1.8.6` — 1.9.2 breaks its import). Any *new* package still needs explicit approval first.

## Capturing the environment

```bash
Rscript -e 'writeLines(capture.output(sessionInfo()), "Environment/env_R.txt")'
Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "Environment/packages_R.tsv", sep="\t", row.names=FALSE, quote=FALSE)'
```
