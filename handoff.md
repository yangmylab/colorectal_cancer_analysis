# Handoff — Project Handoff Document

**Project:** CRC Track A — cross-platform ecotype taxonomy (colorectal cancer)
**Root:** `/run/media/meng/Seqdata/colorectal_cancer` (pipeline in `CRC_BAP/`)
**Updated:** 2026-11-17 (end of the session that added the protein, methylation and spatial layers)

---

## 1. Objective (What we are doing)

Build a **reproducible, cross-platform taxonomy of colorectal cancer** from public multi-omic data,
starting from bulk transcriptomes (TCGA-COAD/READ discovery + GEO external validation) and testing the
resulting "ecotypes" at protein, DNA-methylation, single-cell and spatial levels.
The scientific question is whether tumour-immune ecotypes refine prognosis beyond stage/MSI/CMS/CRIS,
and whether they are reproducible across platforms rather than artefacts of one assay.
Deliverables: a graded findings ledger, a numbered figure set, validation/robustness reports, and a
manuscript whose every claim is licensed by the evidence grade (E1–E4) in the ledger.
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
| GSE14333 | 81 MB, **DFS only** (no OS) | `CRC_BAP/03_Raw_Data/GSE14333/` | 81 MB |
| GSE132465 single cell | 63,689 cells, 6 lineages | `CRC_BAP/03_Raw_Data/GSE132465/` | 1.2 GB |
| CPTAC-COAD (PDC000116, via cBioPortal datahub) | proteome 8,067 × 97; phosphoproteome; mRNA 13,481 × 107 | `CRC_BAP/03_Raw_Data/CPTAC_COAD/` | 29 MB |
| TCGA 450K methylation | 458 beta files (410 primary tumours) → 402,779 probes × 393 samples | `CRC_BAP/03_Raw_Data/TCGA/methylation/` | 6.3 GB |
| GSE283052 Visium | 12 primary tumours, 26,766 in-tissue spots | `CRC_BAP/03_Raw_Data/GSE283052/` | 669 MB |
| GSE267401 Visium | 4 sections (2 patients, primary/metastasis), 17,052 spots | `CRC_BAP/03_Raw_Data/GSE267401/` | 359 MB |

### 2.2 Discovery analysis (TCGA-COAD/READ, real data)

- QC: 14/624 outliers (2.2 %); 50 Hallmark programmes scored by GSVA; SNF multi-view integration.
- **K=2 ecotype is the selected solution**: E1 355 / E2 269, bootstrap ARI **0.975**, nearest-centroid
  classifier CV accuracy **0.979**. E1 top programmes = Hedgehog / pancreas-beta / MYC-V2; E2 = EMT /
  angiogenesis / IL6–JAK–STAT3.
- Survival (n=476, 79 OS events, adjusted for age, sex, stage, sidedness): C-index 0.763 → **0.772**,
  LRT **p = 0.107** → the primary ecotype–OS association is **not supported** (F01, `not_supported`).
- **K=4 exploratory** (post hoc): C-index 0.785, LRT p = 0.024 → graded **E1 only** (F02).
- Robustness: outlier-excluded C-index 0.773 / p 0.151; K-sensitivity table in `08_Robustness/`.

### 2.3 External validation (frozen classifier, no refitting)

- GSE39582: C-index 0.609, HR 1.16, p 0.333. GSE17536: C-index 0.765, HR 1.53, p 0.0785.
- Two-cohort fixed-effect meta: **HR 1.25 (0.97–1.61), p = 0.081** — direction consistent, **not significant**.
- K=4 partition did not replicate: GSE39582 p 0.090, GSE17536 p 0.200.
- Verdict recorded as F03 / F04 `not_supported`; no E2 claim anywhere.

### 2.4 Cross-platform layers (this session)

- **Single cell (descriptive):** GSE132465 → NNLS deconvolution of TCGA; mean fractions epithelial
  0.546, stromal 0.267, myeloid 0.123, B 0.044, T 0.015, mast 0.006 (F05).
  A RNA+cell-type K=5 solution reached ARI 0.885 but classifier CV was only 0.553 → **rejected** as
  not deployable; the K=2 transcriptome-only taxonomy remains primary.
- **Protein (CPTAC-COAD):** frozen classifier applied to protein programme scores — mRNA-vs-protein
  ecotype agreement **85.4 %**, ARI **0.496**, median programme r **0.41** (IFN-α 0.86, IFN-γ 0.80,
  EMT 0.71, angiogenesis 0.68). OS is underpowered (only **8 deaths**): HR 2.44, p 0.286 — descriptive
  (F06, E3).
- **Methylation (TCGA 450K):** gene-level promoter beta for **19,269 genes × 393 tumours**
  (annotation-based selection, **not** variance-based). Result: the ecotype **does not reproduce** —
  ARI vs unsupervised methylation cluster **0.006**, only **5/19,269** genes at FDR<0.05 with
  |logFC|>0.5, global beta p 0.113, CIMP-high 74.5 % (E1) vs 84.3 % (E2), methylation-cluster OS
  HR 0.74, p 0.398 (F07, `not_supported`).
- **Spatial (GSE283052 + GSE267401, 16 sections / 43,818 spots):** the bulk classifier **does not
  transfer to spots**. Documented, not hidden: simplified-statistic self-check on the 624 bulk samples
  gives ARI **0.042**; per-section normalisation forces ~50/50 (invalid, discarded); bulk-referenced
  scoring puts **100 %** of spots on one side of the boundary and spot scores sit **≈5 SD below** the
  bulk distribution. Descriptive result: E2 programmes form **intra-tumour niches** (within-section
  variation > between-section variation); no spatial ecotype claim (F08, `not_supported`).

### 2.5 Outputs

- **Figures:** `CRC_BAP/09_Figures/plots/plot_1` … `plot_27` (PDF + JPEG 300 dpi) with
  `plot_manifest.tsv` mapping number → source → description. Numbering is assigned in
  `14_collect_plots.R`; do not renumber.
  - plot_1 single-cell signatures; plot_2–5 QC; plot_6–8 signatures/consensus/PCA; plot_9–15 stability,
    KM, calibration, DCA, deconvolution, K-sensitivity, validation forest;
    plot_16–19 protein layer; plot_20–24 methylation layer; plot_25–27 spatial layer.
- **Ledger:** `CRC_BAP/13_Reports/Findings_Ledger.csv` — F01–F08, 12 fields each, grades:
  F01/F03/F04/F07/F08 `not_supported`, F02 E1, F05 descriptive, F06 E3.
- **Reports:** `13_Reports/Final_Report.md`; `07_Validation/{Validation,K4_validation,Proteomics,
  Methylation,Spatial}_Report.md`; `08_Robustness/Robustness_Report.md`; `STATUS_OVERNIGHT.md`.
- **Decision log:** `CRC_BAP/Logs/DECISION_LOG.tsv` (every pipeline decision with timestamp).
- **Synthetic dry run** kept isolated in `99_Synthetic_DryRun/` — code validation only, never evidence.

### 2.6 Bugs found and fixed this session (do not reintroduce)

1. `limma` with a cell-means design (`~0 + group`): `topTable(coef="E2")` tests the group **mean**, not
   E2−E1 → had inflated "DMPs" to 10,644. Fixed with `makeContrasts(E2 - E1)` + `contrasts.fit`.
2. Hierarchical clustering: `dist()` was computed over **genes** instead of samples → `cutree` returned
   5,000 clusters. Correct form is `hclust(dist(scale(t(mat))))` for a genes × samples matrix.
3. Methylation gene selection was variance-based (top-60k probes) and then differentially tested —
   circular. Replaced with annotation-based promoter mapping (`probe2gene_promoter.tsv.gz`, 200,334
   probes → 20,164 genes → 19,269 usable genes), aggregated in Python from the full 485k-probe matrix.
4. CPTAC protein IDs are `SYMBOL|SYMBOL` — must be split on `|` before gene-set matching.
5. `minfi::getAnnotation()` needs the annotation package **attached** (`library()`), and the column is
   `Relation_to_Island` (not `Relation_to_UCSC_CG_Island`).
6. Visium GEO matrices are the **full 4,992-spot array**, not the filtered set; column names must come
   from `barcodes.tsv`, and positions must be joined **by barcode** (row orders differ).

---

## 3. Current Blockers (Where we are stuck)

1. **GSE284061 archive is corrupt at the GEO source.** Symptom: `gzip -t` fails on
   `GSM8676861_matrix.mtx.gz` and `GSM8676861_features.tsv.gz`; extracted from a 71,403,520-byte tar that
   matches the server `Content-Length`. Impact: that one spatial section is unusable. Already tried:
   `curl`, `wget`, `aria2c` (fresh and resumed), the per-sample URL
   (`/geo/samples/GSM8676nnn/GSM8676861/suppl/`), and the copy the user downloaded manually — all
   byte-identical (md5 `e3ee0da0cd0a77f2973f3c8fe75aa8a9`). Control files (`GSE39582`, `GSE17536`,
   `GSE283052_features.tsv.gz`) verify fine, so disk and network are sound. Do **not** retry the same
   downloads; the authors must re-upload.
2. **The bulk ecotype classifier cannot be projected onto Visium spots.** Symptom: ARI 0.042 self-check
   on bulk data; all spots on one side of the boundary; ≈5 SD distribution shift. Impact: no spatial
   ecotype claim is possible with the current classifier — the spatial layer stays descriptive.
3. **CPTAC-COAD survival is uninformative** (8 deaths in 97 profiled tumours) — the protein layer can
   test reproducibility of programmes, not prognosis.
4. **GSE14333 has DFS only** — cannot serve the OS validation; would need a separate DFS analysis and
   must not be relabelled as OS.
5. **450K methylation covers promoters/CpGs only** — cannot exclude regulation of the EMT/angiogenic
   programme outside promoter methylation; EPIC/WGBS would be needed.
6. **Manuscript not started** (`12_Manuscript/` empty); the framing must centre the stable K=2 axis and
   the negative survival/validation results.
7. **No version control**: the project is not a git repository, so `handoff.md` and code have no history
   (recommend `git init` + first commit).
8. **Optional packages not installed** (would need explicit user approval): MOFA2,
   ConsensusClusterPlus, decoupleR, iClusterBayes, SPOTlight, immunedeconv, EPIC, xCell, estimate,
   timeROC, survminer. Current pipeline uses SNF/NMF, GSVA ssGSEA, manual consensus clustering,
   Cox/C-index/calibration/DCA as fallbacks.

---

## 4. Next Steps (What to do next)

**P0 — manuscript framing (owner: Prof. Meng + agent; trigger: now)**
Write `CRC_BAP/12_Manuscript/` around: stable K=2 axis (ARI 0.975, CV 0.979), protein-level
reproducibility (85.4 % agreement), methylation non-reproduction (ARI 0.006), and the null prognostic
association (LRT p 0.107; meta HR 1.25, p 0.081). Resume from
`13_Reports/Findings_Ledger.csv` (F01–F08) and `09_Figures/plots/plot_manifest.tsv`.

**P0 — phosphoproteome layer (owner: agent; trigger: next session)**
`03_Raw_Data/CPTAC_COAD/data_phosphoprotein_quantification.txt` is already downloaded but unanalysed.
Plan: phosphosite → protein mapping, programme-level phospho scores, kinase-activity inference
(decoupleR/KSEA if installed), concordance with the protein layer. New script: `21_phospho_analysis.R`;
continue plot numbering at `plot_28`.

**P1 — a genuine spatial ecotype claim (owner: agent; blocked by #2)**
Train the spot-level classifier instead of projecting bulk centroids: deconvolve spots with the
GSE132465 single-cell reference (cell2location/RCTD), score programmes per spot, then test whether the
spatial ecotype predicts outcome. Alternative cohort with matched single-cell data: CNSA `CNP0003321`
(spatial arm of GSE225857). Resume from `06_Analysis/1_R/1_Rscripts/20_spatial_analysis.R` (scoring
functions `score_section()` and the bulk-reference construction at the top are reusable).

**P1 — GSE14333 secondary DFS analysis (owner: agent; blocked by #4)**
Build inputs from `03_Raw_Data/GSE14333/GSE14333_series_matrix.txt.gz`, label endpoint `DFS`
explicitly (`DFS_Time`, `DFS_Cens`), report separately from OS; never mix endpoints.

**P2 — additional external cohorts (owner: agent)**
Protocol lists GSE38832, GSE87211, GSE103479 for the frozen classifier. Add via
`10_geo_series_to_inputs.R` + the `11_/12_build_*_validation.R` pattern.

**P2 — environment hardening (owner: Prof. Meng decision)**
Approve installation of the optional packages listed in §3.8, then re-run the affected steps
(MOFA2 → multi-view factors; SPOTlight → spatial deconvolution; survminer → publication KM).

---

## 5. Pitfalls & Anti-patterns (Never repeat these)

1. **Never standardise scores within a spatial section and then call ecotypes.** Per-section z-scoring
   forces a ~50/50 split by construction — it manufactures "heterogeneity". Always score in a fixed
   reference space and report the distribution shift.
2. **Always self-check a transferred classifier on its own training domain first.** Compute ARI against
   the original labels on the bulk data *before* interpreting any new-domain output. A 0.042 ARI means
   stop, not "tune".
3. **Never select features by variance and then test them for differential signal.** That is circular
   (it produced 10,644 "DMPs"). Select by annotation or independence, then test.
4. **Never use `topTable(coef=<group>)` on a `~0 + group` design** — it tests the group mean. Use
   `makeContrasts()` + `contrasts.fit()`.
5. **Check `dim()` before `dist()`**: for a genes × samples matrix, clustering samples needs
   `dist(scale(t(mat)))`; getting it wrong silently produces nonsense cluster counts.
6. **Verify downloaded archives before analysis**: `tar -tf | wc -l` (member count) and `gzip -t` on a
   couple of members. A file can match `Content-Length` and still be corrupt; three downloaders
   reproducing identical md5 means the *source* is broken (GSE284061).
7. **Do not project bulk centroids onto Visium spots.** A 55 µm spot has no compartment structure
   comparable to a bulk tumour; train on single-cell/spatial data instead.
8. **`pkill -f <pattern>` kills your own shell** when the pattern appears in the command line (it
   terminated the tool call twice). Kill by PID, or use a bracket pattern, and start processes in a
   separate call.
9. **Never pipe an inline interpreter** (`Rscript -e ... | tail`): the host cannot audit it and rejects
   the call. Write the probe to `$TMPDIR/*.R` and run the script file.
10. **Do not renumber figures.** `plot_N` is assigned in `14_collect_plots.R`; new layers append at
    `plot_28+`.
11. **Do not present K=4 or the K=5 single-cell view as confirmatory** — both are post hoc / not
    deployable (classifier CV 0.553) and are graded E1 or rejected.

---

## 6. Key Resources & Environment Info

**Languages / environments**
- R **4.6.1** (2026-06-24), **1,705** packages installed. Key: Seurat, limma, edgeR, DESeq2, survival,
  glmnet, rms, mice, NMF, SNFtool, mixOmics, GSVA, fgsea, msigdbr, ComplexHeatmap, pheatmap, data.table,
  ggplot2, preprocessCore, sva, GEOquery, CellChat, mclust, nnls, riskRegression, prodlim, yaml,
  jsonlite, **minfi** and **IlluminaHumanMethylation450kanno.ilmn12.hg19** (installed this session, with
  user approval).
- Python **3.11** in conda env `sc-analysis` (`/home/meng/miniconda3/envs/sc-analysis/bin/python`;
  numpy 2.4.6, pandas 2.3.3, scipy, scikit-learn, scanpy). `base python3` is 3.14 **without** numpy —
  always use the conda interpreter.
- `aria2c` in conda env `aria2` (fast parallel fetch; note it produced a corrupt GEO tar once — verify
  archives after use).
- No `renv.lock` / `environment.yml`. Capture with:
  `Rscript -e 'ip <- installed.packages()[, c("Package","Version")]; write.table(ip, "Environment/packages_R.tsv", sep="\t", row.names=FALSE, quote=FALSE)'`
  and `conda env export -n sc-analysis > Environment/environment_sc-analysis.yml`.
- Hardware: 16 cores, ~61 GiB RAM, 1.6 TB free on `/run/media/meng/Seqdata`.

**Entry points**
- Whole pipeline: `cd CRC_BAP && ./run_all.sh` (wraps `run_pipeline.R`, logs to `Logs/run_all_*.log`).
- Step-wise (all under `CRC_BAP/06_Analysis/`):
  - `2_Python/1_PythonScripts/01_fetch_tcga_gdc.py` — TCGA RNA-seq fetch/build.
  - `1_R/1_Rscripts/00_validate_inputs.R`, `02_signatures_deconvolution.R`, `03_snf_ecotypes.R`,
    `04_survival_validation.R`, `05_figures_report.R`, `07_robustness.R`, `13_validate_k4.R`,
    `14_collect_plots.R` (assigns `plot_N`).
  - `2_Python/1_PythonScripts/15_build_celltype_signatures.py` — single-cell signatures.
  - `2_Python/1_PythonScripts/16_fetch_tcga_methylation.py` — methylation fetch + probe matrix.
  - `2_Python/1_PythonScripts/19_build_gene_methylation.py` — gene-level beta (needs
    `probe2gene_promoter.tsv.gz`, exported from the 450K annotation).
  - `1_R/1_Rscripts/17_proteomics_validation.R` — CPTAC protein layer (plot_16–19).
  - `1_R/1_Rscripts/19_methylation_analysis.R` — methylation layer (plot_20–24).
  - `1_R/1_Rscripts/20_spatial_analysis.R` — Visium layer + transfer diagnostics (plot_25–27).
  - `1_R/1_Rscripts/10_/11_/12_*.R` — GEO series → validation inputs.

**Data provenance (for re-download)**
- TCGA RNA: GDC API `https://api.gdc.cancer.gov/` (Data Release 46.0), open access.
  Methylation: same API, `data_type = "Methylation Beta Value"`, platform Illumina Human Methylation 450.
- GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSE…>/suppl/` (GSE39582, GSE17536, GSE14333,
  GSE132465, GSE283052, GSE267401).
- CPTAC-COAD: cBioPortal datahub `coad_cptac_2019` via
  `https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/coad_cptac_2019/<file>`
  (GitHub raw serves LFS pointers, the `media` host serves content). Provenance study PDC000116.

**Cross-machine resume**
1. Install R 4.6.x + the package list from `Environment/packages_R.tsv` (Bioconductor 3.23 for
   minfi/annotation); recreate conda env `sc-analysis` from the exported YAML.
2. Copy `CRC_BAP/03_Raw_Data/` (≈11 GB) or re-download from the provenance URLs above.
3. Run `./run_all.sh`; then the layer scripts in the order 17 → 19 → 20 (they read cached RDS/TSV).
4. Regenerate figures with `14_collect_plots.R` and check `09_Figures/plots/plot_manifest.tsv`.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

The new session should restate the current state from this document before doing anything (discovery K=2
+ nulls; protein reproducible; methylation not; spatial descriptive only), confirm which blockers are
still true, and resume from the P0 items in §4. If the previous session was interrupted, first record
what actually completed before the interruption, then continue. Note the standing rules: list any
missing dataset for manual download and wait; do not install software without explicit approval.
