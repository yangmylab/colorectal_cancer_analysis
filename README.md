# A two-ecotype transcriptomic axis in colorectal cancer

Analysis code and derived outputs for the manuscript *"A two-ecotype transcriptomic axis in colorectal cancer
reproduces across assays and is associated with recurrence"* (Qin, Yang, Wang et al.).

The study derives a K=2 partition of Hallmark programme scores in TCGA-COAD/READ, freezes it as a
nearest-centroid classifier, and applies it **without refitting** to CPTAC-COAD proteome and phosphoproteome,
TCGA 450K methylation, 16 Visium sections and four external expression cohorts. Findings were graded before
writing (E3 / E2 / E1 / descriptive / not_supported) and every claim in the manuscript is licensed by one of
them.

## What is here

| Path | Contents |
|---|---|
| `CRC_BAP/06_Analysis/1_R/1_Rscripts/` | The analysis pipeline in execution order (`00_validate` … `25_incremental`) |
| `CRC_BAP/06_Analysis/2_Python/1_PythonScripts/` | Data acquisition helpers |
| `CRC_BAP/Scripts/` | Figure composition, consistency audits, reference-insertion tooling |
| `CRC_BAP/config/pipeline_config.yaml` | All consumed parameters (K range, consensus runs, FDR, seed, workers) |
| `CRC_BAP/04_Processed_Data/` | Intermediate matrices written by the pipeline (regenerated) |
| `CRC_BAP/07_Validation/` | Per-layer validation tables and reports (protein, phospho, methylation, DFS, DSS, incremental value) |
| `CRC_BAP/13_Reports/Findings_Ledger.csv` | The graded findings F01–F15 — the authoritative claim record |
| `CRC_BAP/12_Manuscript/` | Manuscript source, figure plan, claims-to-evidence map, composed figures |
| `CRC_BAP/Logs/DECISION_LOG.tsv` | One row per pipeline decision, with its rationale |

Raw data are **not** in this repository (≈12 GB). They are public and are re-fetched from the accessions listed
in the manuscript: TCGA-COAD/READ (GDC Data Release 46.0), GSE39582, GSE17536, GSE14333, GSE38832, GSE283052,
GSE267401, GSE132465, and CPTAC-COAD (PDC000116, via the cBioPortal datahub).

## Reproducing

1. **R environment**: R 4.6.1 with GSVA 2.6.6, limma 3.68.5, edgeR, survival 3.8.12, SNFtool 2.3.1, nnls 1.6,
   mclust 6.1.3, impute 1.86.0, preprocessCore 1.74.0, pheatmap, KSEAapp 2.0, data.table, ggplot2.
   `Rscript CRC_BAP/Scripts/install_missing_packages.R` installs what is absent.
2. **Python environment**: 3.11.18 with scanpy 1.11.5 and cell2location 0.1.5 (requires `pyro-ppl==1.8.6`;
   1.9.2 breaks its import).
3. **Run order**: `CRC_BAP/run_all.sh` executes the discovery chain; the layer scripts then follow —
   `17_proteomics` → `19_methylation` → `20_spatial` → `21_phospho` → `22a`/`22b`/`22c`/`22d`/`22e` →
   `23_kinase` → `24_spatial_niches` → `25_incremental` → `14_collect_plots`.
4. **Figures**: `Rscript CRC_BAP/Scripts/compose_figures.R` composes the eight composites from
   `12_Manuscript/figures_plan.tsv` at a fixed printed width, and writes a per-panel typography audit.
5. **Checks**: `python3 CRC_BAP/Scripts/audit_analysis.py` (cross-artefact consistency) and
   `python3 CRC_BAP/Scripts/check_figure_citations.py` (every panel cited, in the required order).

Analyses use one level of parallelism with six workers (6 physical cores) and BLAS pinned to one thread per
worker; the seed is `20261117`.

## Licence

See `LICENSE`.
