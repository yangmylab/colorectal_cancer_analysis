# Handoff — Project Handoff Document

**Project:** CRC Track A — cross-platform ecotype taxonomy (colorectal cancer)
**Root:** `/run/media/meng/Seqdata/colorectal_cancer` (pipeline in `CRC_BAP/`)
**Updated:** 2026-11-17 (session that added the phosphoproteome layer and the GSE14333 DFS analysis)

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

---

## 2. Completed Items (What has been done)

### 2.1 Data acquisition (all on disk, verified)

| Dataset | Content | Location | Size |
|---|---|---|---|
| TCGA-COAD/READ RNA-seq | 701 files → 624 patients × 15,000 genes | `CRC_BAP/03_Raw_Data/TCGA/` | 9.1 GB |
| GSE39582 | 585 tumours, 23,963 Ensembl genes, 194 OS events | `CRC_BAP/03_Raw_Data/validation/GSE39582/` | 165 MB |
| GSE17536 | 177 tumours, 73 OS events | `CRC_BAP/03_Raw_Data/validation/GSE17536/` | 208 MB |
| GSE14333 | 290 tumours (GPL570), **DFS only** — now built | `CRC_BAP/03_Raw_Data/GSE14333/` | 81 MB |
| GSE132465 single cell | 63,689 cells, 6 lineages | `CRC_BAP/03_Raw_Data/GSE132465/` | 1.2 GB |
| CPTAC-COAD (PDC000116) | proteome 8,067 × 97; **phosphoproteome 31,262 sites × 97**; mRNA 13,481 × 107 | `CRC_BAP/03_Raw_Data/CPTAC_COAD/` | 29 MB |
| TCGA 450K methylation | 402,779 probes × 393 samples → 19,269 gene promoters | `CRC_BAP/03_Raw_Data/TCGA/methylation/` | 6.3 GB |
| GSE283052 Visium | 12 tumours, 26,766 spots | `CRC_BAP/03_Raw_Data/GSE283052/` | 669 MB |
| GSE267401 Visium | 4 sections, 17,052 spots | `CRC_BAP/03_Raw_Data/GSE267401/` | 359 MB |

### 2.2 Discovery analysis (TCGA-COAD/READ, real data)

- QC 14/624 outliers (2.2 %); 50 Hallmark programmes scored; SNF multi-view integration.
- **K=2 selected**: E1 355 / E2 269, bootstrap ARI 0.975, nearest-centroid CV 0.979.
  E1 = Hedgehog / pancreas-beta / MYC-V2; E2 = EMT / angiogenesis / IL6–JAK–STAT3.
- Survival (n=476, 79 OS events, adjusted): C-index 0.763 → 0.772, LRT **p = 0.107** → F01 `not_supported`.
- K=4 exploratory: C-index 0.785, p = 0.024 → F02 **E1**.

### 2.3 External validation, OS (frozen classifier, no refitting)

- GSE39582 C-index 0.609 / p 0.333; GSE17536 C-index 0.765 / p 0.0785; meta **HR 1.25 (0.97–1.61), p 0.081**.
- K=4 did not replicate (p 0.090 / 0.200). F03/F04 `not_supported`; no E2 claim anywhere.

### 2.4 Cross-platform layers

- **Single cell (descriptive, F05):** NNLS deconvolution (epithelial 0.546, stromal 0.267, myeloid 0.123,
  B 0.044, T 0.015, mast 0.006). RNA+cell-type K=5 reached ARI 0.885 but CV 0.553 → rejected, not deployable.
- **Protein (F06, E3):** 85.4 % mRNA-vs-protein ecotype agreement, ARI 0.496, median programme r 0.41
  (IFN-α 0.86, IFN-γ 0.80, EMT 0.71). OS uninformative (8 deaths).
- **Methylation (F07, `not_supported`):** ARI vs unsupervised methylation cluster 0.006; 5/19,269 genes
  FDR<0.05 with |logFC|>0.5; global beta p 0.113; CIMP-high 74.5 % (E1) vs 84.3 % (E2).
- **Spatial (F08, `not_supported`):** bulk classifier does not transfer to spots (self-check ARI 0.042;
  100 % of spots on one side; ≈5 SD shift). E2 programmes form intra-tumour niches (descriptive only).
- **Phosphoproteome (F09, E3) — NEW this session:** `21_phospho_analysis.R`.
  31,262 sites → 11,334 kept at ≤50 % missing (threshold fixed a priori; full sensitivity ladder written),
  kNN (k=10) imputation, sites averaged per gene (3,414 genes) so hyper-phosphorylated proteins do not
  dominate. 49/50 Hallmark programmes scorable in all three CPTAC layers; frozen classifier recovers 49/50
  features. **Phospho-assigned ecotypes 47 E1 / 50 E2, 87.5 % agreement with the protein call, ARI 0.558,
  1000-permutation null 95th percentile 0.033, p = 0.000999.** Differential phosphorylation: 3,095 sites
  FDR<0.05 (538 with |log2FC|>0.5); 296 with mRNA-assigned labels. Phospho/protein discordance: 2,554 genes,
  276 FDR<0.05. Median programme r: phospho–protein 0.22, phospho–mRNA 0.19, protein–mRNA 0.42.
  Survival **not estimable** (Cox separated with 8 deaths) — reported as such, not as HR 5e8.
- **DFS (F10, E1) — NEW this session:** `22a_build_gse14333_dfs.R` + `22b_dfs_validation.R`.
  GSE14333 rebuilt from the series matrix (GEOquery is absent; 54,674 GPL570 probes → 22,167 symbols by
  highest-expressed probe). Values already log2 — no transform. DFS recorded for 226 of 290 tumours
  (all 61 stage IV patients are `NA`), 50 events, median follow-up 38.5 months. Endpoint parses correctly
  (events rise by stage: 2/41 I, 14/94 II, 34/91 III). Adjusted (age, sex, stage) **HR 2.37 (1.27–4.40),
  LRT p = 0.0043**, C-index 0.686 → 0.711; unadjusted 2.82 (1.52–5.24); stage II–III 2.15 (1.15–4.01);
  500-resample bootstrap median 2.42 (1.27–5.25), **100 % of resamples > 1**; cox.zph min p 0.08.
  Kept strictly out of the OS analysis by construction (builder refuses to write OS/CSS/RFS columns).

### 2.5 Outputs

- **Figures:** `CRC_BAP/09_Figures/plots/plot_1` … `plot_34` (PDF + JPEG 300 dpi) with
  `plot_manifest.tsv` (34 entries). Numbering is assigned in `14_collect_plots.R`; do not renumber.
  plot_16–19 protein; plot_20–24 methylation; plot_25–27 spatial; **plot_28–32 phospho; plot_33–34 DFS**.
- **Ledger:** `13_Reports/Findings_Ledger.csv` — F01–F10, 12 fields. F09/F10 are written by their scripts
  (rewritten in place), so the ledger is reproducible, not hand-edited.
- **Reports:** `13_Reports/Final_Report.md` (rewritten to the F01–F10 state);
  `07_Validation/{Validation,K4_validation,Proteomics,Methylation,Spatial}_Report.md`,
  `07_Validation/21_Phospho_Report.md`, `07_Validation/22b_DFS_Report.md`;
  `08_Robustness/Robustness_Report.md`.
- **Decision log:** `CRC_BAP/Logs/DECISION_LOG.tsv` (added rows for 21_phospho, 22a_gse14333, 22b_gse14333).
- **Version control: DONE** — `git init` + commit `15231d6` (177 files; raw/intermediate data ignored by
  `.gitignore`), plus this session's commit.

### 2.6 Bugs found and fixed (do not reintroduce)

1. `limma` with `~0 + group`: `topTable(coef="E2")` tests the group **mean** → use `makeContrasts()` + `contrasts.fit()`.
2. Hierarchical clustering: `dist()` over genes instead of samples (`hclust(dist(scale(t(mat))))`).
3. Variance-based gene selection followed by differential testing = circular (produced 10,644 fake DMPs).
4. CPTAC protein IDs are `SYMBOL|SYMBOL`; split on `|`.
5. `minfi::getAnnotation()` needs the annotation package attached and uses `Relation_to_Island`.
6. Visium GEO matrices are the full 4,992-spot array; join positions by barcode.
7. **New:** R matrix division is not column-broadcast — `sum_z / cnt[rows, ]` fails; use `as.numeric(cnt[rows, 1])`.
8. **New:** 3,691 CPTAC protein rows contain NA, so `sd()` returns NA and `if (sd(x) == 0)` errors *inside
   every forked worker* ("all scheduled cores encountered errors"). Guard with `is.finite()`/`complete.cases`.
9. **New:** with 8 deaths a Cox fit can separate; `exp(coef)` becomes 5e8 with CI 0–Inf. Detect
   `|coef| > 4` or `SE > 5` and report "not estimable".
10. **New:** `msigdbr` 26.1.1 needs the `msigdbdf` data package or an online release download; it is
    offline here. Take the Hallmark sets from `04_Processed_Data/signature_scores.rds` (frozen discovery
    sets, Ensembl) and map to symbols with `org.Hs.eg.db`.
11. **New:** `GEOquery` and `minfi` are no longer installed — GSE14333 is parsed with `readLines` +
    `fread(skip=, nrows=)` instead.
12. **New:** `cut -d,` on the findings ledger looks corrupt because quoted fields contain commas —
    validate with a real CSV parser.

---

## 3. Current Blockers (Where we are stuck)

1. **GSE284061 archive is corrupt at the GEO source** (md5 `e3ee0da0cd0a77f2973f3c8fe75aa8a9` reproduced by
   curl/wget/aria2c/manual download; control files verify). Do not retry; the authors must re-upload.
2. **The bulk ecotype classifier cannot be projected onto Visium spots** (self-check ARI 0.042). No spatial
   ecotype claim is possible with the current classifier.
3. **CPTAC-COAD survival is uninformative** (8 deaths in 96–97 profiled tumours) — protein and phospho
   layers test reproducibility, not prognosis.
4. **K=4 remains post hoc** and does not replicate.
5. **450K methylation covers promoters/CpGs only** — cannot exclude non-promoter regulation; EPIC/WGBS needed.
6. **Kinase-activity inference is not possible offline**: no `decoupleR`/`KSEAapp`/`PhosR`, no local
   kinase–substrate library, and `msigdbr` cannot fetch anything. F09 therefore reports phosphosite and
   phospho/protein-discordance evidence only. Needs either an approved install or an approved download.
7. **GSE14333 has DFS only, and no DFS for stage IV** — F10 is a secondary-endpoint, stage I–III result and
   must never be pooled with OS.
8. **Environment drift since the previous session (not yet repaired — user deferred):**
   - conda env **`sc-analysis` no longer exists** (only `bioinfo`, `docking`, `pymol_env`); base
     `/home/meng/miniconda3/bin/python` has numpy 2.4.4 / pandas 3.0.2 / sklearn / scipy but **no
     scanpy/anndata**, so the Python steps (01/15/16/19) cannot be re-run as written.
   - **`GSVA` fails to load**: `magick.so` needs `libMagick++-6.Q16.so.9`, the system now ships
     ImageMagick 7. This breaks `02_signatures_deconvolution.R` and the validation scoring in
     `04_survival_validation.R` if re-run. New layer scripts deliberately avoid GSVA.
   - **`minfi`, `SNFtool` and `GEOquery` are not installed** although the previous handoff said they were.
9. **Manuscript not started** (`12_Manuscript/` empty). An empty author template `Authors.xlsx` (13 columns:
   姓名/Email/Mobile/Affiliations/ORCID) appeared on 2026-11-17 11:56, untracked.
10. **Optional packages not installed** (would need explicit approval): MOFA2, ConsensusClusterPlus,
    decoupleR, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate, timeROC, survminer, msigdbdf.

---

## 4. Next Steps (What to do next)

**P0 — manuscript framing (owner: Prof. Meng + agent; trigger: now)**
Write `CRC_BAP/12_Manuscript/` around: stable K=2 axis (ARI 0.975, CV 0.979), protein reproducibility
(85.4 %, F06 E3), **phospho reproducibility (87.5 %, ARI 0.558, F09 E3)**, methylation non-reproduction
(ARI 0.006, F07), spatial non-transfer (F08), null OS (LRT p 0.107; meta HR 1.25, p 0.081) and the
exploratory DFS signal (F10, E1). Resume from `13_Reports/Findings_Ledger.csv` (F01–F10) and
`09_Figures/plots/plot_manifest.tsv` (34 figures). `Authors.xlsx` needs filling.

**P1 — kinase-substrate activity scoring (owner: agent; needs an approved install/download)**
Either install `decoupleR` (then use `decoupleR::run_ksa`/KSEA-style scoring) or drop a KEA/PhosphoSitePlus
substrate table into `03_Raw_Data/CPTAC_COAD/`. Re-run `21_phospho_analysis.R` and add the activity matrix
as plot_35. Until then, do **not** claim kinase activity.

**P1 — replicate F10 in an independent DFS cohort (owner: agent)**
GSE38832 (protocol §, 122 samples) via `10_geo_series_to_inputs.R` + the `11_/12_` builder pattern, then a
`22c_*` DFS analysis. Keep DFS/OS separation absolute.

**P1 — a genuine spatial ecotype claim (blocked by #2)**
Train a spot-level classifier: deconvolve spots against GSE132465 (cell2location/RCTD), score programmes per
spot, then test outcome. Alternative: CNSA `CNP0003321`. Reuse `score_section()` in
`06_Analysis/1_R/1_Rscripts/20_spatial_analysis.R`.

**P2 — additional OS cohorts:** GSE38832, GSE87211, GSE103479 for the frozen classifier.

**P2 — environment repair (owner: Prof. Meng decision):** rebuild `magick` (fixes GSVA), recreate
`sc-analysis`, install `minfi` + `IlluminaHumanMethylation450kanno.ilmn12.hg19` + `SNFtool` + `GEOquery`
+ `msigdbdf`, then re-run `run_all.sh` to confirm end-to-end reproducibility.

---

## 5. Pitfalls & Anti-patterns (Never repeat these)

1. **Never standardise scores within a spatial section and then call ecotypes** — per-section z-scoring
   forces ~50/50 by construction.
2. **Always self-check a transferred classifier on its own training domain first** (ARI 0.042 = stop).
3. **Never select features by variance and then test them** — select by annotation, then test.
4. **Never use `topTable(coef=<group>)` on a `~0 + group` design** — use `makeContrasts()`.
5. **Check `dim()` before `dist()`** — clustering samples needs `dist(scale(t(mat)))`.
6. **Verify downloaded archives** (`gzip -t`, member count); identical md5 from three downloaders means the
   source is broken.
7. **Do not project bulk centroids onto Visium spots.**
8. **`pkill -f <pattern>` kills your own shell** — kill by PID or use a bracket pattern.
9. **Never pipe an inline interpreter** (`Rscript -e ... | head`) — the host rejects it. Write a probe to
   `$TMPDIR/*.R` and run the script file.
10. **Do not renumber figures** — `plot_N` is assigned in `14_collect_plots.R`; append at `plot_35+`.
11. **Do not present K=4 or the K=5 single-cell view as confirmatory.**
12. **Never let a DFS/other endpoint enter the OS path** — the GSE14333 builder asserts that no
    `os_*`/`css_*`/`rfs_*` columns exist, and `04_survival_validation.R` only scans
    `03_Raw_Data/validation/`.
13. **Never quote a separated Cox coefficient as an effect size** (8 deaths → HR 5e8, CI 0–Inf).
14. **Never assume `msigdbr`, `GEOquery`, `minfi` or `GSVA` work on this box** — check first; the
    environment has drifted and no installs are authorised.
15. **Never hand-edit `Findings_Ledger.csv`** — F09/F10 are written by their scripts (delete-and-append).

---

## 6. Key Resources & Environment Info

**Languages / environments (verified 2026-11-17, this session)**
- R **4.6.1**, `.libPaths()` = `/usr/local/lib/R/site-library`, `/usr/lib/R/site-library`, `/usr/lib/R/library`;
  **1,740** packages. Working: Seurat, limma, edgeR, DESeq2, survival, glmnet, rms, mice, NMF, mixOmics,
  fgsea, msigdbr (set table only, **cannot fetch**), ComplexHeatmap, pheatmap, data.table, ggplot2,
  ggrepel, **org.Hs.eg.db**, **hgu133plus2.db**, mclust, nnls, riskRegression, prodlim, impute, matrixStats.
  **Broken/absent:** `GSVA` (magick → ImageMagick 6 ABI), `minfi`, `SNFtool`, `GEOquery`, `decoupleR`,
  `msigdbdf`, `ReactomePA`, `survminer`, `timeROC`, `SPOTlight`.
- Python: **no `sc-analysis` env**. `python3` = miniconda base (numpy 2.4.4, pandas 3.0.2, scipy, sklearn,
  **no scanpy/anndata**). Conda envs present: `bioinfo`, `docking`, `pymol_env`.
- Hardware: **6 physical cores / 12 SMT threads** (`nproc` = 12; `detectCores()` must not be trusted),
  ~61 GiB RAM, 1.6 TB free on `/run/media/meng/Seqdata`.
- Parallel convention: one level only, **6 workers** for heavy work (`mclapply(mc.cores = 6)`), with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` set at the top of each script. Note
  `config/pipeline_config.yaml` still says `n_workers: 8` — stale, supersede by editing the config
  (which also documents "16 physical cores", likewise stale).
- No `renv.lock` / `environment.yml`; capture with `installed.packages()` → `Environment/packages_R.tsv`.

**Entry points (all paths relative to `CRC_BAP/`)**
- Whole pipeline: `./run_all.sh` (wraps `run_pipeline.R`). **Currently cannot complete** until GSVA is fixed.
- Discovery/validation: `06_Analysis/1_R/1_Rscripts/{00,01,02,03,04,05,07,10,11,12,13,14}_*.R`.
- Methylation: `2_Python/1_PythonScripts/{16,19}_*.py`; layer script `1_R/1_Rscripts/19_methylation_analysis.R`.
- Protein: `1_R/1_Rscripts/17_proteomics_validation.R`.
- Spatial: `1_R/1_Rscripts/20_spatial_analysis.R`.
- **Phospho (new):** `1_R/1_Rscripts/21_phospho_analysis.R` → `07_Validation/21_*`, plot_28–32.
- **DFS (new):** `1_R/1_Rscripts/22a_build_gse14333_dfs.R` then `22b_dfs_validation.R` →
  `03_Raw_Data/GSE14333/22a_*`, `07_Validation/22b_*`, plot_33–34.
- Figure manifest: `1_R/1_Rscripts/14_collect_plots.R` (run last; assigns `plot_N`, idempotent).

**Data provenance (for re-download)**
- TCGA RNA + 450K: GDC API `https://api.gdc.cancer.gov/` (Data Release 46.0, open access).
- GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSE…>/suppl/` (GSE39582, GSE17536, GSE14333, GSE132465,
  GSE283052, GSE267401).
- CPTAC-COAD: cBioPortal datahub `coad_cptac_2019` via
  `https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/coad_cptac_2019/<file>`
  (the `media` host serves LFS content; raw.githubusercontent serves pointers). Study PDC000116.
- Hallmark gene sets: reuse `04_Processed_Data/signature_scores.rds$gene_sets` (no network).

**Cross-machine resume**
1. Install R 4.6.x + packages from `Environment/packages_R.tsv`; recreate a Python env with
   numpy/pandas/scipy/scikit-learn/scanpy for the `*_*.py` steps.
2. Copy `CRC_BAP/03_Raw_Data/` (≈11 GB) or re-download from the URLs above.
3. Run `./run_all.sh`; then layer scripts in order 17 → 19 → 20 → 21 → 22a → 22b → 14.
4. Check `09_Figures/plots/plot_manifest.tsv` has 34 entries and `13_Reports/Findings_Ledger.csv` has F01–F10.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

The new session should restate the current state from this document before doing anything (stable K=2 axis;
protein **and** phosphoproteome reproducible; methylation not; spatial descriptive only; OS null in discovery
and external validation; DFS significant but E1 and secondary), confirm which blockers are still true —
especially the **environment drift in §3.8**, which was deferred by the user — and resume from the P0 items in
§4. If the previous session was interrupted, first record what actually completed before the interruption,
then continue. Standing rules: list any missing dataset for manual download and wait; do not install software
without explicit approval; never mix DFS with OS.
