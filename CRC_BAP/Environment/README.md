# Environment

## Verified environment (2026-11-17)

- **R:** 4.6.1 (2026-06-24)
- **R packages installed:** 1,696
- **Core R packages available:** Seurat, limma, edgeR, DESeq2, survival, glmnet, rms, mice, NMF, SNFtool, mixOmics, GSVA, fgsea, msigdbr, ComplexHeatmap, pheatmap, data.table, ggplot2, preprocessCore, sva, GEOquery, CellChat, mclust, nnls, riskRegression, prodlim, future, yaml, jsonlite.
- **Python:** 3.11 in conda env `sc-analysis` (numpy, pandas, scipy, scikit-learn confirmed). Base `python3` is 3.14 without numpy.
- **Hardware:** 16 cores, 61 GiB RAM, 1.6 TB free on `/run/media/meng/Seqdata`.

## Optional packages not installed

MOFA2, ConsensusClusterPlus, decoupleR, minfi, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC, survminer.

**Fallbacks implemented:** SNF/NMF for multi-view integration, GSVA ssGSEA for signatures, manual consensus clustering, Cox/C-index/calibration/DCA for survival.

## No installations

No R or Python packages were installed during this run. If the user wants the optional packages, list them and install only after explicit approval (project standing instruction).

## Capturing the environment

```bash
Rscript -e 'writeLines(capture.output(sessionInfo()), "Environment/env_R.txt")'
Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "Environment/packages_R.tsv", sep="\t", row.names=FALSE, quote=FALSE)'
```
