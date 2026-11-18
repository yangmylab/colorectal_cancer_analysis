# Handoff — Project Handoff Document

**Project:** CRC Track A — cross-platform ecotype taxonomy (colorectal cancer)
**Root:** `/run/media/meng/Seqdata/colorectal_cancer` (pipeline and manuscript in `CRC_BAP/`)
**Updated:** 2026-11-18, third checkpoint: the manuscript revision driven by the pre-submission review, plus the
author's own manual edit of the docx. Since the second checkpoint: 30 of 44 panels re-authored at print width (four
figures now clear the 7 pt floor), one shared ecotype palette, QC panel corrections, the review package and its
13 closed findings, the author's structural edits ported into the markdown source, supplementary-panel
interpretation blended into the main-figure passages, an AI-accent removal pass, a Discussion comparison with
published classifications, and a Methods expansion from ~330 to 1,413 words.
**Git:** history from `15231d6`. Key commits: `f271129` figure typography, `ac1df66` metadata, `d9cefab` F13,
`c8c80b0` review, `cb62d25`/`60c43df` figure infrastructure, `c9d9462` wording + statements, `fb173bf` group sizes,
`ea159c4` re-authored panels, `77dfcc9` blended supplements, `cd210e7` AI-accent removal, `781ebd0` Discussion
comparison, `5d3cc27` expanded Methods. Working tree clean except the user's untracked `Authors.xlsx`.

---

## 1. Objective (What we are doing)

Build a **reproducible, cross-platform taxonomy of colorectal cancer** from public multi-omic data: a K=2
transcriptomic "ecotype" axis derived in TCGA-COAD/READ, frozen as a classifier, then tested without refitting at
protein, phosphoprotein, methylation, single-cell and spatial levels, and against disease-free and mortality
endpoints in independent cohorts.

The scientific question is whether the axis (i) survives a change of assay and (ii) predicts outcome **beyond
stage/MSI/CMS**, the published taxonomy. Deliverables: a graded findings ledger, a numbered figure set composed
into 6 main + 2 supplementary figures, validation/robustness reports, and a manuscript in which every claim is
licensed by the evidence grade (E1–E3) of the finding behind it.

**Success criteria:** stable taxonomy in discovery **and** independent validation; honest reporting of negative
results; no claim beyond the data (no causal or clinical-utility claims).

**Current headline (verified):** the axis reproduces across assays (protein F06 E3, phosphoproteome F09 E3) and is
**associated with disease-free survival in three cohorts (F11, E2)** — but it is **not incremental over CMS**
(F14), and **no mortality association replicates** (F01, F03, F15). Reproducibility is established; clinical
added value is not. That dissociation is the paper.

---

## 2. Completed Items (What has been done)

### 2.1 Data acquisition (all on disk, verified)

| Dataset | Content | Location |
|---|---|---|
| TCGA-COAD/READ RNA-seq | 624 patients × 15,000 genes (GDC DR 46.0) | `03_Raw_Data/TCGA/` |
| GSE39582 | 585 tumours; **536 RFS / 147 events**; **CMS for all 585**; MSI for 536 | `03_Raw_Data/validation/GSE39582/` |
| GSE17536 | 177 tumours; 73 OS events; 145 RFS / 36 events | `03_Raw_Data/validation/GSE17536/` |
| GSE14333 | 290 tumours (GPL570), DFS only: 226 with DFS, 50 events | `03_Raw_Data/GSE14333/` |
| GSE38832 | 122 CEL files → RMA matrix; **metadata fetched via GEO SOFT**; 92 with DFS (9 recurrences), 122 with DSS (28 deaths); **no age/sex** | `03_Raw_Data/GSE38832/` |
| CPTAC-COAD (PDC000116) | protein 8,067 × 97; phospho 31,262 sites × 97; mRNA 13,481 × 107 | `03_Raw_Data/CPTAC_COAD/` |
| TCGA 450K | 402,779 probes × 393 samples → 19,269 gene promoters | `03_Raw_Data/TCGA/methylation/` |
| GSE283052 / GSE267401 Visium | 16 sections / 43,818 spots | `03_Raw_Data/GSE283052|GSE267401/` |
| GSE132465 single cell | 63,689 cells, 6 lineages; annotation + TPM matrix | `03_Raw_Data/GSE132465/` |

### 2.2 Analysis results (ledger is authoritative: `13_Reports/Findings_Ledger.csv`, F01–F15)

| ID | Grade | Claim (abbreviated) |
|---|---|---|
| F01 | not_supported | No incremental association with OS after adjustment (C-index 0.772, LRT p=0.107) |
| F02 | E1 | K=4 associated with OS in discovery — **post hoc** |
| F03 | not_supported | OS association did not replicate (meta HR 1.25, p=0.081) |
| F04 | not_supported | K=4 did not replicate (p=0.090 / 0.200) |
| F05 | descriptive | Single-cell reference resolves six lineages; TCGA deconvolution descriptive only |
| F06 | E3 | Protein-level reproduction: **83.3% agreement, ARI 0.439**, median programme r 0.40 |
| F07 | not_supported | Methylation non-reproduction: **ARI 0.006**, 5/19,269 genes at FDR<0.05 |
| F08 | not_supported | Bulk classifier does not transfer to spots (self-check ARI 0.042; ~5 SD shift) |
| F09 | E3 | Phosphoproteome reproduction: **87.5%, ARI 0.558, permutation p=0.000999** |
| F11 | **E2** | Disease-free survival: GSE14333 HR 2.37; GSE39582 1.50 (p=0.018); GSE17536 2.51; meta **1.67 (1.24–2.24), p=0.0007** |
| F12 | E1 | Kinase activity: 9/78 at FDR<0.05, 5 abundance-independent; **inferred, not measured** |
| F13 | descriptive | E2 programmes spatially autocorrelated (I=0.528) but **not more than others (0.587)**; stromal co-localisation ρ=0.49 — **pattern reproduced** by a trained-reference cell2location check (matrix r 0.93), **absolute fractions method-dependent** |
| F14 | not_supported | **Not incremental over stage+MSI+CMS**: HR 1.42 (0.93–2.17), p=0.107; Cramér's V 0.52 |
| F15 | not_supported | No association with cancer-specific death in GSE38832 (HR 1.00, p=0.997) |

Ledger: 14 rows, 12 fields, ids in order, no duplicates; grades 7 not_supported / 2 E3 / 1 E2 / 2 E1 / 2 descriptive.

### 2.3 Key parameters and decisions (why the numbers are what they are)

- **Discovery:** 50 Hallmark programmes; SNF + consensus k-means; K=2 selected on stability (**bootstrap ARI 0.975**,
  classifier CV **0.979**); K=4 kept only as post hoc (F02) and never presented as confirmatory.
- **Phosphoproteome:** sites kept if quantified in ≥50% of tumours (**11,334 of 31,262**; threshold fixed before
  labels were read, full sensitivity ladder reported); kNN imputation (k=10); **sites averaged per gene** so
  hyper-phosphorylated proteins do not dominate → 3,414 genes.
- **Kinase activity:** KSEAapp `KSData` (PhosphoSitePlus, 19,739 pairs / 387 kinases); kinases require ≥5 measured
  substrates (78 testable, median 10).
- **Spatial:** NNLS against the GSE132465 signature (83 markers, 6 lineages); Moran's I on a k=6 neighbour graph,
  199 permutations; **no spot-level ecotype is ever assigned** (F08 constraint).
- **Outcomes:** disease-free endpoints (DFS/RFS) and mortality endpoints (OS/CSS/DSS) are **never pooled**. Primary
  disease-free model is **stage I–III adjusted for age, sex, stage**; all-stage is sensitivity. Bootstrap = 500
  resamples. F14 uses nested models: stage/age/sex → +MSI → +CMS → +MSI+CMS.
- **Parallelism everywhere:** one level, **6 workers** (6 physical cores; `nproc` reports 12 SMT), BLAS pinned to
  1 thread per worker.
- **Grading rule:** E3 = cross-platform convergence; E2 = replicated association; E1 = exploratory/post hoc;
  descriptive; not_supported. No claim exceeds its grade.

### 2.4 Environment repaired and extended (this and the previous session)

| Symptom | Root cause | Fix |
|---|---|---|
| `GSVA` would not load | `magick` built against ImageMagick 6, system has 7 | rebuilt magick 2.9.1 |
| `GEOquery` lazy-load failed | `xml2` linked `libxml2.so.2`; system has `.so.16`; **conda's `xml2-config` shadowed the system one** | rebuilt with `PATH=/usr/local/bin:/usr/bin:/bin` |
| minfi, 450K anno, SNFtool, hgu133plus2cdf absent | never installed | installed via tuna Bioconductor |
| `msigdbr()` could not fetch | `msigdbdf` on no reachable mirror | frozen GMT + `hallmark_gene_sets()` (offline) |
| conda env `sc-analysis` gone | deleted | recreated with `--override-channels`, py3.11.18 + scanpy 1.11.5 |
| `cell2location` absent | 539 MB torch wheel failed on aliyun | installed from tuna; needs **`pyro-ppl==1.8.6`** (1.9.2 breaks on `deep_getattr`) |

### 2.5 Deliverables (all on disk and committed)

- **Figures — 44 panels:** `09_Figures/plots/plot_1` … `plot_44` (PDF + JPEG 300 dpi) with
  `plot_manifest.tsv`.
- **Composed figures — 8:** `12_Manuscript/Figures/` → `Figure_1` … `Figure_6`,
  `Supplementary_Figure_1`, `Supplementary_Figure_2` (vector PDF + 300 dpi JPEG), panel letters A/B/C…, built from
  `12_Manuscript/figures_plan.tsv` by `Scripts/compose_figures.R`. Since 2026-11-18 the composer places the panel
  **PDF pages** with `pdflatex` + TikZ (true vector, selectable text; verified with `pdftotext`/`pdfimages`), uses
  **one scale factor per figure** so panel type size no longer varies within a composite, and sets panel letters at
  a fixed **8 pt Helvetica bold** in a reserved 4 mm strip above each panel. Options: `--mode=uniform|equal`,
  `--width_mm` (180), `--height_mm` (240), `--floor_pt` (7), `--gutter_mm`, `--letter_strip_mm`, `--letter_pt`.
  Requires `pdflatex` (TikZ) + `pdftoppm`.
- **Typography audit:** `12_Manuscript/compose_figures_typography.{tsv,md,docx}` — per panel: authored vs printed
  size, scale factor, printed body/title pt, floor pass/fail. Current state: **all 44 panels print at 2.8–5.4 pt,
  i.e. below the 7 pt floor** (full breakdown under §3 blocker 10).

| Figure | Panels | Plots | Role |
|---|---|---|---|
| 1 | 6 | 6, 7, 8, 9, 1, 13 | the axis in discovery |
| 2 | 6 | 16, 17, 19, 28, 29, 31 | reproduction at protein + phosphoproteome |
| 3 | 7 | 20, 21, 22, 24, 25, 26, 27 | failure: methylation, spatial transfer |
| 4 | 6 | 30, 32, 37, 38, 39, 40 | mechanism: phosphosites, kinase, spatial structure |
| 5 | 6 | 33, 34, 35, 36, 43, 44 | recurrence replicates, death does not |
| 6 | 5 | 10, 14, 15, 41, 42 | OS null, K=4 post hoc, no value over CMS |
| Supp 1 | 4 | 2, 3, 4, 5 | quality control (supports Figure 1) |
| Supp 2 | 4 | 11, 12, 18, 23 | diagnostics + negative controls (supports Figures 2, 3, 6) |

- **Manuscript:** `12_Manuscript/manuscript.md` (v0.2) + `manuscript.docx` (built with the user's reference
  template `MN/Paper_1/reference/tnr_fangsong_redcode.docx`; verified TNR + FangSong) + `claims_to_evidence.tsv`
  (14 claims with grade and citation caveat) + `figures_plan.tsv` (with `supports_main` column) + `README.md`.
  Results cites all 44 panels with **first citations strictly in order** (1A→…→6E→S1A→…→S2D).
- **Reports:** `13_Reports/Final_Report.md`; `07_Validation/{Validation,K4_validation,Proteomics,Methylation,
  Spatial}_Report.md`; `21_Phospho_Report.md`; `22b_DFS_Report.md`; `22d_DFS_Replication_Report.md`;
  `22e_GSE38832_Report.md`; `25_Incremental_Value_Report.md`; `06_Analysis/11_SpatialNiches/24_Spatial_Niche_Report.md`,
  `25_cell2location_check.md` (the 100-epoch prior-only run, kept as a superseded record) and
  `26_cell2location_reference.{md,tsv,_overall.tsv}` (the trained-reference verdict).
- **cell2location reference:** `Scripts/cell2location_reference.py` → `06_Analysis/09_SingleCell/
  cell2location_reference_signature.tsv` (2,007 genes x 6 lineages) + `cell2location_reference_summary.tsv`;
  cache `04_Processed_Data/singlecell/c2l_reference.h5ad` (10,187 cells, gitignored). Re-run the comparison with
  `Scripts/cell2location_check.py --prior reference --epochs 3000 --spots 1000`.
- **Verification tooling:** `Scripts/audit_analysis.py` → `13_Reports/Consistency_Audit.md` (**22/22 pass**);
  `Scripts/check_figure_citations.py` (**44/44 panels, order enforced**).
- **Pre-submission review (2026-11-18):** `13_Reports/Review_2026-11-18/Review_Report.{md,docx}` +
  `review_findings.tsv` — 26 severity-graded findings (5 P0, 13 P1, 8 P2) from the `peer-review`,
  `nature-statistics` and `paper-claim-audit` playbooks. Every number was checked against its source file: no
  fabricated value and no Results sentence exceeding its grade. The P0s are infrastructural — no references
  (R01), no availability/ethics/funding/contribution statements (R02), all 44 panels below the 7 pt floor (R03),
  no group sizes or model denominators (R04), no statistics in the legends (R05).
- **Revision pass 1 (2026-11-18, commits `cb62d25` + next):** figures infrastructure done — `save_plot`
  honours `CRC_PANEL_MAX_W`/`CRC_PANEL_SCALE` (tested: 9x7 in -> 3.4x2.64 in), one Okabe-Ito E1/E2 palette
  across 7 scripts, S1C re-coloured by project with a Wilcoxon subtitle, S1D annotation renamed `project`.
  Text: R06-R10, R22, R26 fixed; R09 reworded without inventing a test; `## Statements` added with TODO
  placeholders for funding, contributions, approvals and the code DOI.
  **Still open (see `13_Reports/Review_2026-11-18/fix_status.md`):** R01 bibliography, R03 re-render
  (run the figure scripts with `CRC_PANEL_MAX_W=3.4`, then `git checkout --` the data dirs 05_QC
  06_Analysis 07_Validation 08_Robustness 10_Tables 13_Reports and `Logs/DECISION_LOG.tsv` so only figures
  change, then recompose), R04 group sizes, R05 legends, R11 Methods, R15-R19 and R23-R25.

- **Next session, R01 (bibliography) — do this first, in one pass.** The manuscript still has no citations. Plan:
  (1) inventory every citable claim (CMS1-4, CRIS, PDS1-3, TCGA-COAD/READ, CPTAC-COAD, the six GEO cohorts,
  MSI/CIMP, and the methods: SNF, GSVA, limma, nnls, cell2location, KSEA/PhosphoSitePlus, MSigDB); (2) resolve each
  through OpenAlex/PubMed and keep only entries whose metadata returns — never write a reference from memory;
  (3) attach each to the single sentence it supports; (4) run `citation-audit`. Expect 45-60 entries.
- **Figure leftovers.** `20_spatial_analysis.R` aborts on a curl download of gene sets (the frozen GMT path is not
  reached in that script), so plot_25-27 keep their old 8-10 in width; fix the loader, re-run 20 and 14 with
  `CRC_PANEL_MAX_W=3.4`, recompose. A VISUAL pass over the regenerated composites is still owed — the re-render was
  verified by measurement and audits only, so overlapping labels in the 50-row heatmaps are not yet ruled out.

- **Decision log:** `Logs/DECISION_LOG.tsv`, 61 rows covering every step from `00_validate` to `25_incremental`.
- **Revision artefacts (2026-11-18):** `13_Reports/Review_2026-11-18/` holds `Review_Report.{md,docx}` (26 graded
  findings), `review_findings.tsv`, `fix_status.md` (per-finding state) and `humaniser_changelog.md` (AI-accent
  removal, then the Discussion comparison). `Scripts/rerender_panels.sh` re-authors the panels and restores the
  data artefacts; `Scripts/check_figure_citations.py` now enforces "main panels in order, each supplementary panel
  after its `supports_main` anchor".

### 2.6 Bugs found and fixed (do not reintroduce)

1. `limma` with `~0 + group`: `topTable(coef="E2")` tests the group mean → use `makeContrasts()` + `contrasts.fit()`.
2. `dist()` over genes instead of samples → `hclust(dist(scale(t(mat))))`.
3. Variance-based gene selection followed by differential testing = circular (10,644 fake DMPs).
4. CPTAC protein IDs are `SYMBOL|SYMBOL`; split on `|`.
5. Visium GEO matrices are the full 4,992-spot array; join positions by barcode.
6. R matrix division is not column-broadcast — use `as.numeric(cnt[rows, 1])`.
7. 3,691 CPTAC protein rows contain NA → `sd()` NA → `if()` throws **inside every forked worker**.
8. With 8 deaths a Cox fit can separate → detect `|coef|>4` / `SE>5` and report "not estimable".
9. `msigdbr` 26.1.1 needs `msigdbdf` or a download; neither available → frozen GMT.
10. `GEOquery`/`minfi` absent → GSE14333 parsed with `readLines` + `fread(skip=, nrows=)`.
11. `cbind(mclapply(...))` silently drops failing workers (85 kinases became 37) → use `vapply`.
12. Permutation p-values hit a ceiling with ~2,500 spots → interpret effect sizes only.
13. Adding a model spec broke a `factor()` level list → rows became `NA` and report bullets duplicated;
    `stopifnot(nrow(...) == 1)` now guards selections.
14. **"Batch-corrected" QC panel was a near-duplicate of the uncorrected one**: the "batch" was TCGA project
    (COAD vs READ), explaining ≤2.8% of any PC, and the corrected matrix was never used downstream. Now the
    project effect is *quantified* (panel 1B) and no correction is claimed.
15. Renaming a QC source file silently leaves a stale copy unless `14_collect_plots.R`'s mapping is updated too.

---

## 3. Current Blockers (Where we are stuck)

1. **`msigdbdf` cannot be installed** (tuna CRAN 404, CRAN/r-universe/zenodo unreachable). `msigdbr()` cannot
   fetch, so only the frozen Hallmark GMT is available. A *different* MSigDB collection would have to be added to
   `Environment/hallmark_frozen.gmt` manually.
2. **GSE284061 archive is corrupt at the GEO source** (md5 `e3ee0da0cd0a77f2973f3c8fe75aa8a9` from three
   downloaders). Do not retry; the authors must re-upload. One spatial section is unusable.
3. **The bulk classifier cannot be projected onto Visium spots** (self-check ARI 0.042). This is a design
   constraint, not a bug: spot-level labels would require a spot-trained classifier (see §4).
4. **CPTAC-COAD survival is uninformative** (8 deaths) — protein/phospho layers test reproducibility only.
5. ~~**F13's co-localisation statistic is method-dependent.**~~ **SETTLED 2026-11-18:** a reference trained on
   GSE132465 (2,007 genes, 10,187 cells, 250 epochs in 3.3 min) was run at 3,000 epochs x 1,000 spots over three
   sections. It reproduces the *spatial pattern* of the NNLS composition (matrix r 0.94/0.70/0.93, top-lineage
   agreement 0.97/0.74/0.99) but not the *absolute composition* (epithelial 0.39 vs 0.97; stromal 0.17 vs 0.01).
   So the earlier "did NOT corroborate (r 0.11)" reading is withdrawn as an under-training artefact, and F13's
   caveat narrows from "direction uncertain" to "absolute magnitude method- and gene-panel-dependent". Note that
   unregularised NNLS is itself degenerate here (~0.97 on one lineage, exactly 0 on several), which is why a
   per-lineage correlation is undefined for 9 of 18 lineage-section pairs. Evidence:
   `06_Analysis/11_SpatialNiches/26_cell2location_reference.{md,tsv,_overall.tsv}`.
6. **GSE38832's GEO record is thin**: 9 recurrences (cannot test F11), no age/sex, DSS only for mortality. A
   fuller follow-up would require the genuine series-matrix download (blocked FTP; manual download only).
7. **Only one cohort carries CMS/MSI** (GSE39582), so F14 (non-incrementality over CMS) is single-cohort by
   necessity — the main limitation on the clinical-value conclusion.
8. **`Authors.xlsx` is an empty template** (untracked) and no target journal has been chosen; both block
   submission.
9. ~~**Stale project metadata:** `config/pipeline_config.yaml` says `n_workers: 8` and `Analysis_Plan.md` says
   "16 physical cores".~~ **FIXED 2026-11-18:** both configs (real + synthetic) now set `n_workers: 6`,
   `Analysis_Plan.md` states 6 physical / 12 SMT, and `Environment/README.md` was corrected from "16 cores" /
   1,696 packages / "no installations" to the verified 6 physical cores, 1,754 packages and the real install
   history. `PILOT_PLAN.md` and `METHODS_TOOLBOX_2026.md` keep their 16-32 core *budgets* (they describe what each
   analysis needs, not this box) with a one-line caveat added.
10. **Panel legibility is part-fixed; four figures still miss the 7 pt floor.** After re-authoring 30 of 44 panels
    at `CRC_PANEL_MAX_W=3.4`, the printed body type is Fig 2 7.33, Fig 4 7.33, Fig 5 11.25, Fig 6 7.34 (all
    compliant), S2 6.35, Fig 3 4.01, Fig 1 3.10, S1 5.08. **21 of 44 panels meet the floor (was 0 of 44).**
    Measured: `12_Manuscript/compose_figures_typography.tsv`. Two causes remain: (a) `20_spatial_analysis.R`
    aborts on a curl download of gene sets, so plot_25–27 kept their 8–10 in width; (b) the discovery chain
    (02–05) was deliberately not re-run because consensus clustering can move the K=2 labels, so Fig 1, S1 and the
    04-derived panels keep their old size — see §5 pitfall 35.
11. **The manuscript still has no citations.** `13_Reports/Review_2026-11-18/` records R01 as open; the Discussion
    comparison added on 2026-11-18 names CMS/CRIS/PDS without bracketed references, so those paragraphs are where
    a bibliography is most needed. Plan in §4 P0-1.
12. **No visual pass over the regenerated panels.** The re-render was verified by measurement and by the two
    audits only; overlapping labels in the 50-row heatmaps are not yet ruled out.
13. **Legends still carry no n, test or interval definition**, and per-model n is missing beyond the group sizes now
    in Methods (review R05/R04-remainder).

---

## 4. Next Steps (What to do next)

**P0 — finish the submission package**

**P0-1 (agent, next session, do it in one pass).** Build the bibliography (review R01). Inventory every citable
claim, resolve each through OpenAlex/PubMed and keep only entries whose metadata returns, attach each to the single
sentence it supports, then run `citation-audit`. Expect 45–60 entries. Start with the three Discussion paragraphs
added on 2026-11-18 (`**Comparison with published classifications**`), which name CMS1–4, CRIS (A–D) and PDS
without citations.

**P0-2 (agent).** Finish the figure floor: fix the curl/gene-set loader in
`06_Analysis/1_R/1_Rscripts/20_spatial_analysis.R` (the frozen `Environment/hallmark_frozen.gmt` path is not
reached there), re-run 20 and 14 with `CRC_PANEL_MAX_W=3.4`, recompose, then look at all eight composites. Then
decide with the user whether to re-run the discovery chain for Figure 1 and S1, which risks moving the K=2 labels.

**P0-3 (agent).** Rewrite the twelve figure legends with n, test, interval and correction (review R05), and add
per-model n (R04 remainder). Resume from `12_Manuscript/figures_plan.tsv` and the per-layer reports.

**P0-4 (user).** Fill `Authors.xlsx` and choose the target journal; both gate formatting.

**P1 — settle the one soft statement (F13 co-localisation)** — **DONE 2026-11-18**
`Scripts/cell2location_reference.py` trains the reference on GSE132465 raw counts (2,007 genes, 10,187 cells, 250
epochs, 3.3 min); `Scripts/cell2location_check.py --prior reference --epochs 3000 --spots 1000` then runs three
sections. Outcome: the spatial pattern is reproduced (matrix r 0.93 median) but the absolute composition is not
(epithelial 0.39 vs 0.97), so F13's caveat was narrowed rather than removed, and the F13 ledger `limitation` field
now derives its cell2location sentence from `26_cell2location_reference_overall.tsv` instead of hardcoding a result.
What would remove the caveat entirely: agreement in absolute fractions, or a GPU-scale run over all 16 sections to
show the gap is a sampling/training artefact rather than a property of unregularised NNLS.

**P2 — optional strengthening**
**P2-1.** Fuller GSE38832 metadata (manual series-matrix download) to test F11 in a fourth cohort.
**P2-2.** A second cohort with CMS/MSI + a disease-free endpoint, to test whether F14 is cohort-specific.
**P2-3.** Refresh environment records: `Environment/packages_R.tsv` and a `conda env export` for `sc-analysis`.

**P2 — housekeeping**
**P2-4.** ~~Correct `n_workers` in `config/pipeline_config.yaml` and the "16 physical cores" line in `Analysis_Plan.md`.~~
   **DONE 2026-11-18** (both configs, Analysis_Plan, Environment/README plus caveats in PILOT_PLAN and
   METHODS_TOOLBOX_2026).

---

## 5. Pitfalls & Anti-patterns (Never repeat these)

1. **Never standardise scores within a spatial section** and then call ecotypes — it manufactures a ~50/50 split.
2. **Self-check a transferred classifier on its own training domain first** (ARI 0.042 = stop, not "tune").
3. **Never select features by variance then test them** — select by annotation, then test.
4. **Never `topTable(coef=<group>)` on `~0 + group`** — use `makeContrasts()` + `contrasts.fit()`.
5. **Check `dim()` before `dist()`** — cluster samples with `dist(scale(t(mat)))`.
6. **Verify archives** (`gzip -t`, member count); identical md5 from three downloaders means the source is broken.
7. **Do not project bulk centroids onto Visium spots.**
8. **`pkill -f` kills the calling shell** when the pattern is in its own command line (happened twice). Use
   `pgrep -f "[p]attern"` then `kill <pid>`. The host also refuses `rm -rf "$VAR"` — use fresh directories.
9. **Never pipe an inline interpreter** (`Rscript -e … | head`); write `$TMPDIR/*.R` and run the file.
10. **Do not renumber figures** — `plot_N` is assigned in `14_collect_plots.R`; append at `plot_45+`.
11. **Never present K=4 or the K=5 single-cell view as confirmatory.**
12. **Never let a disease-free endpoint enter the OS path** — builders assert no `os_*`/`css_*`/`rfs_*` columns;
    `04_survival_validation.R` only scans `03_Raw_Data/validation/`.
13. **Never quote a separated Cox coefficient as an effect size** (8 deaths → HR 5e8, CI 0–Inf).
14. **A gene-set id mapping can move results.** Changing symbol mapping moved protein agreement 85.4% → 83.3%
    (2 of 96 borderline argmax calls). Re-run and re-quote after any mapping change.
15. **Never state that a tool was installed without importing it.** The spatial report once claimed cell2location
    was installed when it was not.
16. **Never hand-edit script-managed ledger rows** (F09, F11, F12, F13, F14, F15) — each writer deletes and
    re-appends, and all writers now sort by finding id.
17. **`cbind(mclapply(...))` silently drops failures** — use `vapply`.
18. **Permutation p-values can hit a ceiling with large n** (2,500 spots → all 800 tests p<0.05). Report effect
    sizes, and say so.
19. **Do not "batch-correct" a biological covariate.** "COAD vs READ" is anatomy, not a technical batch
    (≤2.8% of any PC). Quantify such structure instead of removing it.
20. **If a before/after panel looks identical, the effect size is the finding** — plot the number, don't ship two
    look-alike scatterplots.
21. **After renaming a source panel, update `14_collect_plots.R`'s mapping** or a stale copy is silently reused.
22. **Adding a model specification breaks any `factor()` built from the old level list** — guard selections with
    `stopifnot(nrow(...) == 1)`.
23. **Check cross-artefact consistency by script, not by eye**: run
    `python3 Scripts/audit_analysis.py` and `python3 Scripts/check_figure_citations.py` before any commit that
    touches figures or the ledger.
24. **Never rescale panels independently inside a composite.** The old composer scaled every panel to a fixed
    900 px width, which made one figure's type size vary by 1.14×–1.67×; use one scale factor per figure
    (`--mode=uniform`).
25. **Printed type size is set by the panel's authored width, not by the montage.** A 9-in panel printed into a
    60 mm column is ~3 pt no matter how the composite is assembled. Check
    `12_Manuscript/compose_figures_typography.tsv` before calling any figure submission-ready — do not assume
    "it was drawn at 300 dpi, so it is fine".
26. **Reserve space for panel letters instead of overlaying them.** The old letter boxes sat on top of panel
    content; a 4 mm strip above each panel keeps them legible and identical across figures.
27. **`pdflatex` twice for the composite.** TikZ `remember picture`/`overlay` placement needs two passes; one pass
    silently produces mispositioned nodes.
28. **`groupby(...).groups` returns index labels, not positions** — mapping them through `ann.iloc[...]` throws on
    string ids. Build a `pd.Series(np.arange(n), index=...)` and map explicitly.
29. **Never assume a matrix's cell order matches the annotation.** Read the matrix header and reindex it to the
    annotation (`pd.Index(header).get_indexer(ann.index)`), asserting no `-1`, before attaching labels.
30. **cell2location's export keys vary by version.** In 0.1.5 + scvi-tools 1.4.2 the per-lineage profile is
    `varm['means_per_cluster_mu_fg']` (plus q05/q95), not `means_cell_abundance_w_sf`. Look the key up with a
    fallback list and error out loudly if none is present.
31. **A near-uniform cell2location fraction is a training-budget symptom, not a negative result.** 100 epochs with
    a marker-matrix prior returns ~1/6 per lineage; the same data at 3,000 epochs with a trained reference moves
    to r 0.93 composition agreement. Never quote the under-trained run as evidence against the NNLS fractions.
32. **Per-lineage correlations are undefined when the comparator is exactly zero.** Unregularised NNLS collapses
    onto one lineage (~0.97) and returns exactly 0 for several others, so `corr` on those is NaN, not 0. Report
    `n/a`, exclude them from the median, and add a composition-level statistic (matrix r over spots x lineages)
    that survives the degeneracy.
33. **`mclapply` forks inherit one RNG stream**, so a `set.seed()`-ed permutation test still moves individual
    p-values between runs (observed 0.02 → 0.005 for one section x programme) even though aggregates are stable.
    Re-running to change *text* is safe; re-running to regenerate *statistics* is not, so check which one you need.
34. **Generated reports inherit their text from the writer script, not from reality.** `24_Spatial_Niche_Report.md`
    kept saying cell2location "was attempted but not installed" after it had been installed. Fix the template in
    the script (and re-run); hand-editing the generated file only lasts until the next run.
35. **Do not re-run the discovery chain to change a figure.** `03_snf_ecotypes.R` derives the K=2 labels under
    `mclapply`, whose forks share one RNG stream; a re-run can move the labels and invalidate every downstream
    number. Re-author only the layer panels, or verify the labels are byte-identical afterwards.
36. **Treat the figure-citation sequence as an invariant when editing prose.** A single added "Figure 1C" shifts
    the first-citation order and fails `check_figure_citations.py`. Compare the extracted `(kind, number, panel)`
    sequence before and after any manuscript edit; the checker now also requires each supplementary panel to be
    cited after the main panel in its `supports_main` column.
37. **Write Methods from the code, never from memory.** Drafting the expanded Methods, two remembered values were
    wrong: the CIMP panel is **eight** genes (CACNA1G, CDKN2A, CRABP1, IGF2, MLH1, NEUROG1, RUNX3, SOCS1), not
    five; and the 0.70 stability floor gates *confirmatory* status rather than filtering candidate K (K is
    filtered by the 5% minimum group, then ordered by mean bootstrap ARI).
38. **The docx is a build artefact; the markdown is the source.** The author's manual docx edit had to be ported
    into `manuscript.md` before anything else, or the next `pandoc` rebuild would have discarded it. When the user
    says "work in the revised docx", extract it, diff against the markdown, and port the differences first.

---

## 6. Key Resources & Environment Info

**R (verified 2026-11-17)**
- R 4.6.1; `.libPaths()` = `/usr/local/lib/R/site-library`, `/usr/lib/R/site-library`, `/usr/lib/R/library`;
  1,740+ packages.
- Working: GSVA 2.6.6, magick 2.9.1, xml2 1.6.0, GEOquery 2.80.1, affy 1.90.0, hgu133plus2cdf 2.18.0,
  hgu133plus2.db 3.13.0, minfi 1.58.0, bumphunter 1.54.0, 450K annotation 0.6.1, SNFtool 2.3.1,
  org.Hs.eg.db 3.23.1, limma 3.68.5, survival 3.8.12, Seurat 6.0.0, NMF 0.28, mixOmics 6.36.0,
  ComplexHeatmap 2.28.0, mclust 6.1.3, nnls 1.6, sva 3.60.0, preprocessCore 1.74.0, impute 1.86.0,
  msigdbr 26.1.1 (loads, cannot fetch), **decoupleR 2.17.0**, **KSEAapp 2.0** (offline kinase–substrate prior).
- **Absent:** `msigdbdf`, `decoupleRdata`, `SPOTlight`, `immunedeconv`, `EPIC`, `xCell`, `estimate`, `timeROC`,
  `survminer`, `ReactomePA`, `MOFA2`, `ConsensusClusterPlus`, `iClusterBayes`.
- CRAN `https://mirrors.tuna.tsinghua.edu.cn/CRAN`; BioC
  `https://mirrors.tuna.tsinghua.edu.cn/bioconductor` (bioconductor.org itself is blocked).

**Python**
- `sc-analysis` at `/home/meng/miniconda3/envs/sc-analysis/bin/python` — 3.11.18, numpy 2.4.6, pandas 2.3.3,
  scipy 1.17.1, scikit-learn 1.9.1, scanpy 1.11.5, anndata 0.12.19, **cell2location 0.1.5**, scvi-tools 1.4.2,
  **pyro-ppl 1.8.6** (do not upgrade: 1.9.2 breaks cell2location's import), torch 2.15.0+cu132 (CPU-only box,
  `cuda available: False`).
- Recreate with `--override-channels` and the tuna conda-forge + pkgs/main mirrors.

**Hardware / parallelism**
- **6 physical cores / 12 SMT threads** (`nproc` = 12; never pass `detectCores()` to heavy work), 61 GiB RAM,
  1.6 TB free on `/run/media/meng/Seqdata`.
- One level of parallelism, **6 workers**, `OMP/OPENBLAS/MKL_NUM_THREADS=1`.

**Network reality (2026-11-17)**
- Reachable: tuna & aliyun mirrors, `www.ncbi.nlm.nih.gov` **SOFT/eutils endpoints**.
- Blocked: `ftp.ncbi.nlm.nih.gov`, `bioconductor.org`, `cran.r-project.org`, `zenodo.org`, `omnipathdb.org`,
  `ghproxy.net`, `r-universe.dev`.
- Consequence: GEO *sequence/data* downloads are manual, but **GEO metadata can be fetched via SOFT**
  (`acc.cgi?acc=GSM…&targ=self&form=text&view=brief`) — this is how GSE38832 was completed.

**Entry points (relative to `CRC_BAP/`)**
- Pipeline: `./run_all.sh` → `run_pipeline.R`.
- Discovery/QC: `06_Analysis/1_R/1_Rscripts/{00,01,02,03,04,05,07,10,11,12,13,14}_*.R`.
- Layers: `17_proteomics_validation.R`, `19_methylation_analysis.R`, `20_spatial_analysis.R`,
  `21_phospho_analysis.R`, `23_kinase_activity.R`, `24_spatial_niches.R`.
- Outcomes: `22a_build_gse14333_dfs.R` + `22b_dfs_validation.R`, `22c_build_gse38832_dfs.R` +
  `22e_gse38832_endpoints.R`, `22d_dfs_replication.R`, `25_incremental_value.R`.
- Manuscript: `Scripts/build_claims_map.py`, `Scripts/compose_figures.R` (vector montage; needs `pdflatex` with
  `tikz` + `pdftoppm`; options `--mode=uniform|equal`, `--width_mm=180`, `--height_mm=240`, `--floor_pt=7`),
  `Scripts/check_figure_citations.py`, `Scripts/audit_analysis.py`, `Scripts/fix_qc_project_effect.R`,
  `Scripts/fetch_gse38832_metadata.py`, `Scripts/cell2location_check.py`,
  `Scripts/install_missing_packages.R`, `Scripts/export_frozen_hallmark_gmt.R`,
  `Scripts/sort_findings_ledger.R`.
- Environment artefact: `Environment/hallmark_frozen.gmt` (50 Hallmark sets, Ensembl).
- Manifest: `14_collect_plots.R` (run last; idempotent; safe self-copy).

**Data provenance**
- TCGA RNA + 450K: GDC API `https://api.gdc.cancer.gov/` (DR 46.0).
- GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSE…>/suppl/` (manual); metadata via `acc.cgi` SOFT.
- CPTAC-COAD: cBioPortal datahub `coad_cptac_2019` via `media.githubusercontent.com` (LFS content). PDC000116.

**Cross-machine resume**
1. `Rscript Scripts/install_missing_packages.R` (~25 min; needs ≤6 cores).
2. Recreate `sc-analysis` with `--override-channels`; then pin `pyro-ppl==1.8.6` if cell2location is needed.
3. Copy `03_Raw_Data/` (≈12 GB) or re-download manually.
4. Install the figure toolchain: `texlive-latex-recommended` (TikZ) + `poppler-utils` (`pdftoppm`, `pdfinfo`);
   `Scripts/compose_figures.R` stops with an explicit message if either is missing.
4. `./run_all.sh`; then layers 17 → 19 → 20 → 21 → 22a → 22b → 22d → 22e → 23 → 24 → 25 → 14.
5. Expect `plot_manifest.tsv` = 44 entries, `Findings_Ledger.csv` = F01–F15, and
   `12_Manuscript/Figures/` = 16 files.
6. Verify with `python3 Scripts/audit_analysis.py` (22/22) and `python3 Scripts/check_figure_citations.py`.
7. Rebuild the Word manuscript: `pandoc 12_Manuscript/manuscript.md --reference-doc=/run/media/meng/Seqdata/MN/Paper_1/reference/tnr_fangsong_redcode.docx -o 12_Manuscript/manuscript.docx --toc --toc-depth=2`.

---

## 7. Startup Instruction for the Next Session

> **When starting a new session, the very first message must be: "Read handoff.md first."**

Then restate the state before doing anything (stable K=2 axis; protein + phospho reproducible E3; methylation and
spatial negative; **disease-free survival replicates E2 but is not incremental over CMS**; **no mortality
association in any cohort**), confirm which blockers still hold — especially whether `Authors.xlsx` and the target
journal are settled — and resume from the P0 items in §4. If the previous session was interrupted, record what
actually completed before the interruption, then continue.

Standing rules: list any missing dataset for manual download and wait (NCBI FTP is blocked; GEO SOFT may work);
do not install software without explicit approval; never mix disease-free and mortality endpoints; never let a
claim exceed its evidence grade.
