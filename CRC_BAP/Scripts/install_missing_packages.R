#!/usr/bin/env Rscript
# install_missing_packages.R
# Rebuilds the R side of this project's environment after the drift discovered 2026-11-17:
#   - `magick` was compiled against ImageMagick 6 and could not load, which broke GSVA (and therefore
#     `02_signatures_deconvolution.R` and the validation scoring in `04_survival_validation.R`);
#   - minfi, SNFtool, GEOquery, msigdbdf, hgu133plus2cdf and the 450K annotation were absent.
#
# Usage:  Rscript Scripts/install_missing_packages.R
# Sources: CRAN + Bioconductor via the Tsinghua mirrors (bioconductor.org and NCBI are unreachable here).
# Parallelism: Ncpus = 6 (physical cores of the i5-10400; nproc reports 12 SMT threads).

options(timeout = 900)
options(Ncpus = 6L)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
cat("BiocManager version:", as.character(BiocManager::version()), "\n")
print(BiocManager::repositories())

install_cran <- function(pkgs) {
  for (p in pkgs) {
    cat("\n=== CRAN:", p, "===\n")
    if (requireNamespace(p, quietly = TRUE)) cat("already available\n") else install.packages(p, Ncpus = 6L)
  }
}

install_bioc <- function(pkgs) {
  cat("\n=== Bioconductor:", paste(pkgs, collapse = ", "), "===\n")
  BiocManager::install(pkgs, update = FALSE, ask = FALSE, Ncpus = 6L)
}

# 1. magick must be rebuilt first: GSVA's namespace chain fails on the ImageMagick 6 ABI.
install.packages("magick", Ncpus = 6L)

# 2. CRAN packages
install_cran(c("SNFtool", "msigdbdf"))

# 3. Bioconductor packages (CEL processing + GEO metadata + methylation layer)
install_bioc(c("GEOquery", "hgu133plus2cdf", "minfi",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19"))

cat("\n=== verification ===\n")
for (p in c("magick", "GSVA", "SNFtool", "msigdbdf", "GEOquery", "affy", "hgu133plus2cdf",
            "minfi", "bumphunter", "IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  ok <- tryCatch({suppressPackageStartupMessages(loadNamespace(p)); TRUE}, error = function(e) {
    cat("  FAILED:", p, "-", conditionMessage(e), "\n"); FALSE
  })
  cat(sprintf("%-45s %s\n", p, ifelse(ok, "LOADS", "BROKEN")))
}
