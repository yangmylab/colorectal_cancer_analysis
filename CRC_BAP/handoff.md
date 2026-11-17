# Handoff

Updated: 2026-11-17 (phosphoproteome + GSE14333 DFS session)

## State
- Real analysis pipeline completed at K=2 (discovery TCGA-COAD/READ, validation GSE39582 + GSE17536,
  single-cell GSE132465, protein CPTAC-COAD, phosphoproteome CPTAC-COAD, methylation TCGA 450K,
  spatial descriptive, DFS GSE14333).
- Final figures: `09_Figures/plots/plot_1` … `plot_34` with `plot_manifest.tsv`.
- Evidence ledger: `13_Reports/Findings_Ledger.csv` F01–F10 (E3: F06, F09; E1: F02, F10;
  not_supported: F01, F03, F04, F07, F08; descriptive: F05).

## Next action
- Manuscript framing (`12_Manuscript/` is empty; `Authors.xlsx` template needs filling).
- Phospho: kinase-activity inference blocked (no offline substrate library / no authorised install).
- DFS is a secondary endpoint: never pool with the OS analysis.
- Spatial layer: bulk classifier does not transfer to spots (ARI 0.042 self-check); descriptive only.

## Key files
- `00_Project_Design/Analysis_Plan.md`
- `02_Metadata/Sample_Annotation.csv`
- `13_Reports/Findings_Ledger.csv`, `13_Reports/Final_Report.md`
- `07_Validation/21_Phospho_Report.md`, `07_Validation/22b_DFS_Report.md`
- `Logs/DECISION_LOG.tsv`

## Caveat
The R/Python environment has drifted (GSVA broken via magick/ImageMagick 7; minfi, SNFtool, GEOquery and
the conda env `sc-analysis` are gone; msigdbr cannot fetch offline). `./run_all.sh` cannot run end-to-end
until those are repaired — see the project-root `handoff.md` §3.8.
