# Handoff — Project Handoff Document

**Project:** CRC Track A — cross-platform ecotype taxonomy (colorectal cancer)
**Root:** `/run/media/meng/Seqdata/colorectal_cancer` (pipeline in `CRC_BAP/`)
**Updated:** 2026-11-17 (session that repaired the environment and added the RFS replication, F11)

---

## 1. Objective (What we are doing)

Build a **reproducible, cross-platform taxonomy of colorectal cancer** from public multi-omic data,
starting from bulk transcriptomes (TCGA-COAD/READ discovery + GEO external validation) and testing the
resulting "ecotypes" at protein, phosphoprotein, DNA-methylation, single-cell and spatial levels.
The scientific question is whether tumour-immune ecotypes refine prognosis beyond stage/MSI/CMS/CRIS,
and whether they are reproducible across platforms rather than artefacts of one assay.
Deliverables: a graded findings ledger, a numbered figure set, validation/robustness reports, and a
manuscript whose every claim is licensed by the evidence grade (E1–E3; E4 unreachable without wet lab).
Success criteria: stable taxonomy in discovery **and** independent validation, honest reporting of
negative results, and no claim beyond the data (no causal or clinical-utility claims).

**Current headline:** the axis reproduces across assays (protein E3, phosphoproteome E3) and **predicts
disease-free survival in three cohorts (E2, F11)**, while **no OS association replicates (F01/F03)**.
That dissociation — recurrence yes, death no — is the story the manuscript has to tell.

---

## 2. Completed Items (What has been done)

### 2.1 Data acquisition (all on disk, verified)

| Dataset | Content | Location | Size |
|---|---|---|---|
| TCGA-COAD/READ RNA-seq | 701 files → 624 patients × 15,000 genes | `CRC_BAP/03_Raw_Data/TCGA/` | 9.1 GB |
| GSE39582 | 585 tumours, 23,963 Ensembl genes, 194 OS events, **536 RFS / 147 events** | `03_Raw_Data/validation/GSE39582/` | 165 MB |
| GSE17536 | 177 tumours, 73 OS events, **145 RFS / 36 events** | `03_Raw_Data/validation/GSE17536/` | 208 MB |
| GSE14333 | 290 tumours (GPL570), **DFS only** (226 with DFS, 50 events) | `03_Raw_Data/GSE14333/` | 81 MB |
| GSE38832 | **122 CEL files (GPL570) — expression built, NO metadata yet** | `03_Raw_Data/GSE38832/` | 620 MB |
| GSE132465 single cell | 63,689 cells, 6 lineages | `03_Raw_Data/GSE132465/` | 1.2 GB |
| CPTAC-COAD (PDC000116) | proteome 8,067 × 97; phosphoproteome 31,262 sites × 97; mRNA 13,481 × 107 | `03_Raw_Data/CPTAC_COAD/` | 29 MB |
| TCGA 450K methylation | 402,779 probes × 393 samples → 19,269 gene promoters | `03_Raw_Data/TCGA/methylation/` | 6.3 GB |
| GSE283052 / GSE267401 Visium | 16 sections / 43,818 spots | `03_Raw_Data/GSE283052|GSE267401/` | 1.0 GB |

### 2.2 Discovery (TCGA-COAD/READ)

- QC 14/624 outliers (2.2 %); 50 Hallmark programmes; SNF integration.
- **K=2 selected**: E1 355 / E2 269, bootstrap ARI 0.975, nearest-centroid CV 0.979.
  E1 = Hedgehog / pancreas-beta / MYC-V2; E2 = EMT / angiogenesis / IL6–JAK–STAT3.
- OS (n=476, 79 events, adjusted): C-index 0.763 → 0.772, LRT **p = 0.107** → F01 `not_supported`.
- K=4 post hoc: C-index 0.785, p = 0.024 → F02 **E1**.

### 2.3 External validation

- **OS (F03/F04, `not_supported`)**: GSE39582 C-index 0.609/p 0.333; GSE17536 0.765/p 0.0785;
  meta HR 1.25 (0.97–1.61), p 0.081. K=4 did not replicate (p 0.090/0.200).
- **DFS discovery (F10, E1)**: GSE14333, n=226, 50 events, adjusted HR 2.37 (1.27–4.40), LRT p 0.0043.
- **RFS replication (F11, E2) — NEW**: `22d_dfs_replication.R`. GSE39582 stage I–III adjusted
  **HR 1.50 (1.07–2.10), p 0.018** (n=483, 137 events; C-index 0.545→0.570); all stages 1.55 (1.12–2.15),
  p 0.008; unadjusted stage I–III 1.48 (1.06–2.07), p 0.023. Bootstrap 500 × 6 workers: median 1.48
  (1.06–2.13), **99 % above 1**. GSE17536 RFS stage I–III HR 2.51 (1.14–5.53), p 0.022 (all stages 2.83,
  p 0.0056). Fixed-effect with GSE14333: HR 1.67 (1.24–2.24), p 0.0007, heterogeneity p 0.21.
  Primary model is stage I–III by user decision (curative-intent, matches GSE14333 which has no stage IV).

### 2.4 Cross-platform layers

- **Single cell (F05, descriptive):** NNLS deconvolution; RNA+cell-type K=5 rejected (CV 0.553).
- **Protein (F06, E3) → RE-GRADED NUMBERS:** 83.3 % mRNA-vs-protein agreement, ARI 0.439, median
  programme r 0.40 (IFN-α 0.87, IFN-γ 0.81, EMT 0.71). Cox identical to before (C-index 0.596,
  HR 2.44, p 0.286). The change from 85.4 %/0.496 is a *gene-set mapping* artefact — see §5.13.
- **Methylation (F07, `not_supported`):** ARI vs methylation cluster 0.006; 5/19,269 genes FDR<0.05.
- **Spatial (F08, `not_supported`):** transfer self-check ARI 0.042; spots ≈5 SD below bulk; descriptive.
- **Phosphoproteome (F09, E3):** 11,334 of 31,262 sites (≤50 % missing); 87.5 % agreement with the
  protein call, ARI 0.558, 1000-permutation null 95th pct 0.033, p 0.000999; 3,095 sites FDR<0.05
  (538 with |log2FC|>0.5); 2,554-gene phospho/protein discordance screen (276 FDR<0.05).
  Survival reported as **not estimable** (Cox separation with 8 deaths).
  Kinase-activity inference **deliberately not claimed** (no substrate library offline).

### 2.5 Environment repaired this session (was blocking everything)

| Symptom | Root cause | Fix |
|---|---|---|
| `GSVA` would not load | `magick.so` built against ImageMagick **6**; system ships **7** | `install.packages("magick")` → 2.9.1 |
| `GEOquery` lazy-load failed | `xml2.so` linked `libxml2.so.2`; system has SONAME `.so.16`; **conda's `xml2-config` 2.13.9 shadowed the system 2.15.2** | rebuild with `PATH=/usr/local/bin:/usr/bin:/bin` |
| `minfi`, 450K anno, `SNFtool`, `hgu133plus2cdf` absent | never installed on this box | installed via tuna Bioconductor/CRAN |
| `msigdbr()` could not fetch | `msigdbdf` not on tuna CRAN, 404 in CRAN Archive, cran.r-project.org + r-universe unreachable | frozen GMT `Environment/hallmark_frozen.gmt` + `hallmark_gene_sets()` in `utils.R` |
| conda env `sc-analysis` gone | deleted after the previous session | recreated with `--override-channels` (py3.11.18, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, sklearn 1.9.1, scanpy 1.11.5) |

Install scripts kept: `CRC_BAP/Scripts/install_missing_packages.R`,
`CRC_BAP/Scripts/export_frozen_hallmark_gmt.R`.

### 2.6 Outputs

- **Figures:** `09_Figures/plots/plot_1` … `plot_36` (PDF + JPEG 300 dpi), `plot_manifest.tsv` (36 entries).
  Numbering assigned in `14_collect_plots.R`; append at `plot_37+`.
  plot_16–19 protein; 20–24 methylation; 25–27 spatial; **28–32 phospho; 33–34 GSE14333 DFS; 35–36 GSE39582 RFS**.
- **Ledger:** `13_Reports/Findings_Ledger.csv` — F01–F11 (12 fields). F09/F10/F11 are written by their
  scripts and rewritten in place; do not hand-edit them.
- **Reports:** `13_Reports/Final_Report.md`; `07_Validation/{Validation,K4_validation,Proteomics,
  Methylation,Spatial}_Report.md`; `21_Phospho_Report.md`; `22b_DFS_Report.md`; `22d_DFS_Replication_Report.md`.
- **Decision log:** `Logs/DECISION_LOG.tsv` — now also records 17_proteomics, 21_phospho, 22a/22b_gse14333,
  22c_gse38832 (BLOCKED), 22d_replication.
- **Git:** commits `15231d6` → `e64e9c8` → `5c1d04b` → (this session's final commit).

---

## 3. Current Blockers (Where we are stuck)

1. **GSE38832 has no metadata.** The RAW archive is CEL-only; `GSE38832_series_matrix.txt.gz` is absent and
   NCBI is unreachable from this host (IPv4/IPv6, http/https, eutils — all fail; Chinese mirrors work).
   Expression is already built, so the cohort runs the moment the file lands in
   `03_Raw_Data/GSE38832/`. **Needs a manual download by the user.**
2. **`msigdbdf` cannot be installed** (not on any reachable mirror). `msigdbr()` therefore cannot fetch;
   everything now goes through the frozen GMT. If a future step needs a *different* MSigDB collection
   (not Hallmark), it must be added to the GMT manually.
3. **GSE284061 archive is corrupt at the GEO source** (md5 `e3ee0da0cd0a77f2973f3c8fe75aa8a9` reproduced by
   three downloaders). Do not retry; the authors must re-upload.
4. **The bulk classifier cannot be projected onto Visium spots** (self-check ARI 0.042) — spatial stays
   descriptive (F08).
5. **CPTAC-COAD survival is uninformative** (8 deaths) — protein/phospho layers test reproducibility only.
6. **450K methylation covers promoters only** — EPIC/WGBS needed to exclude non-promoter regulation.
7. **Kinase-activity inference needs a substrate library** (`decoupleR`/KSEAapp/PhosR or a KEA table).
   Not installed, not downloadable → not claimed in F09.
8. **Manuscript not started** (`12_Manuscript/` empty). `Authors.xlsx` is still an unfilled template.
9. **Optional packages still absent** (would need approval): MOFA2, ConsensusClusterPlus, decoupleR,
   iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC, survminer.
10. **`README.md`/`Analysis_Plan.md`/`pipeline_config.yaml` are stale** about the machine: they say
    16 physical cores / `n_workers: 8`; the box has **6 physical / 12 SMT**. Change the config before
    claiming the pipeline's parallelism in a methods section.

---

## 4. Next Steps (What to do next)

**P0 — manuscript (owner: Prof. Meng + agent; trigger: now)**
`CRC_BAP/12_Manuscript/` around: stable K=2 axis; protein + phospho reproducibility (E3); methylation and
spatial non-reproduction; **null OS (F01/F03) vs replicated disease-free association (F11 E2)** — the
dissociation is the paper. Resume from the ledger F01–F11 and `plot_manifest.tsv` (36 figures). Fill
`Authors.xlsx`.

**P0 — GSE38832 DFS (owner: user file + agent; blocked by #1)**
Drop `GSE38832_series_matrix.txt.gz` into `03_Raw_Data/GSE38832/`, then
`Rscript 06_Analysis/1_R/1_Rscripts/22c_build_gse38832_dfs.R` (auto-detects endpoint keys/units, asserts no
OS columns) and add a `22e_*` analysis mirroring `22d_dfs_replication.R` with the GSE38832 files.

**P1 — kinase-substrate activity (needs approval)**
Install `decoupleR` or place a KEA/PhosphoSitePlus table in `03_Raw_Data/CPTAC_COAD/`; re-run
`21_phospho_analysis.R` and add the activity matrix as plot_37. Do not claim activity until then.

**P1 — genuine spatial claim** (blocked by #4): deconvolve spots with GSE132465 (needs the Python env, now
restored), train a spot-level classifier, test outcome.

**P2 — more OS/DFS cohorts:** GSE87211, GSE103479 (protocol) — NCBI unreachable, so manual download.
**P2 — environment housekeeping:** `msigdbdf` via another route; optional packages above; refresh
`Environment/packages_R.tsv` and add a conda YAML export now that `sc-analysis` exists again.

---

## 5. Pitfalls & Anti-patterns (Never repeat these)

1. **Never standardise scores within a spatial section** — it manufactures ~50/50 splits.
2. **Self-check a transferred classifier on its own domain first** (ARI 0.042 = stop).
3. **Never select features by variance then test them** (produced 10,644 fake DMPs).
4. **Never `topTable(coef=<group>)` on `~0 + group`** — use `makeContrasts()` + `contrasts.fit()`.
5. **Check `dim()` before `dist()`** — cluster samples with `dist(scale(t(mat)))`.
6. **Verify archives** (`gzip -t`, member count); identical md5 from three downloaders = broken source.
7. **Do not project bulk centroids onto Visium spots.**
8. **`pkill -f` kills your own shell** — kill by PID; and no `rm -rf "$VAR"` (the host refuses
   non-literal deletes; use a fresh directory instead).
9. **Never pipe an inline interpreter** (`Rscript -e ... | head`) — write `$TMPDIR/*.R` and run the file.
10. **Do not renumber figures** — `plot_N` is assigned in `14_collect_plots.R`; append at `plot_37+`.
11. **Never present K=4 or the K=5 single-cell view as confirmatory.**
12. **Never let DFS/RFS enter the OS path** — builders assert no `os_*`/`css_*`/`rfs_*` columns and
    `04_survival_validation.R` only scans `03_Raw_Data/validation/`.
13. **Never quote a separated Cox coefficient as an effect size** (8 deaths → HR 5e8, CI 0–Inf).
14. **A gene-set *id mapping* can move results.** The protein layer's agreement changed 85.4 %→83.3 %
    purely because symbol sets came from `org.Hs.eg.db` instead of msigdbr's column; 2 of 96 borderline
    argmax calls flipped. When you change how sets are built, re-run and re-quote the numbers.
15. **`msigdbdf`, `GEOquery`, `minfi`, `GSVA` may be missing again** — this environment drifts. Check with
    `Scripts/install_missing_packages.R` before assuming a pipeline step can run.
16. **conda needs `--override-channels` with the tuna mirrors**; otherwise it contacts
    `repo.anaconda.com` for ToS notices and dies (network unreachable), even when `-c` mirrors are given.
17. **GEO RAW tars may not be gzipped despite the name** — `tar -tf` works, `tar -xzf` fails with
    "gzip: stdin: not in gzip format". Use `tar -xf`.
18. **ABI drift breaks R packages silently.** ImageMagick 6→7 and libxml2 `.so.2`→`.so.16` both produced
    "cannot open shared object file" at *load* time; and conda's `xml2-config` in PATH was silently used
    for a system build. Rebuild the affected package with a minimal `PATH`.
19. **R matrix division is not column-broadcast** — `sum_z / cnt[rows, ]` fails; use `as.numeric(cnt[rows, 1])`.
20. **NA in a matrix makes `sd()` NA and `if()` throw inside every forked worker** ("all scheduled cores
    encountered errors"). Guard with `is.finite()`.
21. **Never hand-edit script-managed ledger rows** (F09/F10/F11) — they are delete-and-append on each run.

---

## 6. Key Resources & Environment Info

**R (verified 2026-11-17 after repair)**
- R 4.6.1; `.libPaths()` = `/usr/local/lib/R/site-library`, `/usr/lib/R/site-library`, `/usr/lib/R/library`.
- Working and version-checked: GSVA 2.6.6, magick 2.9.1, xml2 1.6.0, GEOquery 2.80.1, affy 1.90.0,
  hgu133plus2cdf 2.18.0, hgu133plus2.db 3.13.0, minfi 1.58.0, bumphunter 1.54.0,
  IlluminaHumanMethylation450kanno.ilmn12.hg19 0.6.1, SNFtool 2.3.1, org.Hs.eg.db 3.23.1, limma 3.68.5,
  survival 3.8.12, Seurat 6.0.0, NMF 0.28, mixOmics 6.36.0, ComplexHeatmap 2.28.0, mclust 6.1.3,
  nnls 1.6, sva 3.60.0, preprocessCore 1.74.0, impute 1.86.0, msigdbr 26.1.1 (loads, but cannot fetch).
- **Absent:** `msigdbdf`, `decoupleR`, `SPOTlight`, `immunedeconv`, `EPIC`, `xCell`, `estimate`,
  `timeROC`, `survminer`, `ReactomePA`, `MOFA2`, `ConsensusClusterPlus`, `iClusterBayes`.
- CRAN mirror `https://mirrors.tuna.tsinghua.edu.cn/CRAN`; BioC mirror
  `https://mirrors.tuna.tsinghua.edu.cn/bioconductor` (bioconductor.org itself is **blocked**).

**Python**
- `sc-analysis` at `/home/meng/miniconda3/envs/sc-analysis/bin/python` — 3.11.18, numpy 2.4.6,
  pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.1, scanpy 1.11.5, anndata 0.12.19.
- Create/repair with **`--override-channels`** and the tuna conda-forge + pkgs/main mirrors.

**Hardware / parallelism**
- **6 physical cores / 12 SMT threads** (`nproc` = 12; never pass `detectCores()` to heavy work),
  61 GiB RAM, 1.6 TB free on `/run/media/meng/Seqdata`.
- One level of parallelism, **6 workers**, `OMP/OPENBLAS/MKL_NUM_THREADS=1` set at the top of each script.
  `config/pipeline_config.yaml` still says `n_workers: 8` — **stale, fix before quoting methods**.

**Network reality on this host (2026-11-17)**
- Reachable: `mirrors.tuna.tsinghua.edu.cn`, `mirrors.aliyun.com`, `baidu.com`.
- Blocked: `ftp.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`, `eutils.ncbi.nlm.nih.gov`, `www.ebi.ac.uk`,
  `bioconductor.org`, `cran.r-project.org`, `zenodo.org`, `ghproxy.net`, `r-universe.dev`.
- Consequence: **all GEO-sequence downloads are manual**; plan around cohorts already on disk.

**Entry points (relative to `CRC_BAP/`)**
- Pipeline: `./run_all.sh` → `run_pipeline.R` (now runnable again: GSVA restored).
- Discovery/validation: `06_Analysis/1_R/1_Rscripts/{00,01,02,03,04,05,07,10,11,12,13,14}_*.R`.
- Protein: `17_proteomics_validation.R`; methylation: `19_methylation_analysis.R`;
  spatial: `20_spatial_analysis.R`; phospho: `21_phospho_analysis.R`;
  GSE14333 DFS: `22a_build_gse14333_dfs.R` + `22b_dfs_validation.R`;
  **GSE38832 build: `22c_build_gse38832_dfs.R`**; **RFS replication: `22d_dfs_replication.R`**.
- Environment: `Scripts/install_missing_packages.R`, `Scripts/export_frozen_hallmark_gmt.R`,
  `Environment/hallmark_frozen.gmt` (50 Hallmark sets, Ensembl).
- Figure manifest: `14_collect_plots.R` (run last; idempotent; safe self-copy).

**Data provenance**
- TCGA RNA + 450K: GDC API `https://api.gdc.cancer.gov/` (Data Release 46.0).
- GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSE…>/suppl/` (**manual download only**).
- CPTAC-COAD: cBioPortal datahub `coad_cptac_2019` via
  `https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/coad_cptac_2019/<file>`
  (the `media` host serves LFS content). Study PDC000116.
- Hallmark sets: frozen GMT in `Environment/` (no network).

**Cross-machine resume**
1. `Rscript Scripts/install_missing_packages.R` (needs <=6 cores, ~25 min).
2. Recreate `sc-analysis` with `--override-channels` + tuna mirrors.
3. Copy `03_Raw_Data/` (≈12 GB) or re-download manually (NCBI is blocked here).
4. `./run_all.sh`; then layer scripts 17 → 19 → 20 → 21 → 22a → 22b → 22d → 14.
5. Expect `plot_manifest.tsv` = 36 entries and `Findings_Ledger.csv` = F01–F11.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

Restate the state before doing anything (stable K=2 axis; protein + phospho reproducible E3; methylation and
spatial negative; **OS null everywhere; disease-free survival replicates, F11 E2**), confirm the blockers —
especially whether `GSE38832_series_matrix.txt.gz` has been supplied — and resume from the P0 items in §4.
Standing rules: list any missing dataset for manual download and wait (NCBI is blocked, so this is normal
here); do not install software without explicit approval; never mix DFS/RFS with OS.
