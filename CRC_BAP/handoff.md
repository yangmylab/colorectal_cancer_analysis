# Handoff

Updated: 2026-11-17 08:44:13

## State
- Real analysis pipeline completed at K=2 (discovery TCGA-COAD/READ, validation GSE39582 + GSE17536, single-cell GSE132465, protein CPTAC-COAD, methylation TCGA 450K, spatial descriptive: 16 sections / 43,818 spots).
- Final figures: `09_Figures/plots/plot_1` … `plot_24` with `plot_manifest.tsv`.

## Next action
- Review findings, robustness and validation before manuscript claims.
- Spatial layer: bulk classifier does not transfer to spots (ARI 0.042 self-check); descriptive only — see `07_Validation/Spatial_Report.md`.

## Key files
- `00_Project_Design/Analysis_Plan.md`
- `02_Metadata/Sample_Annotation.csv`
- `13_Reports/Findings_Ledger.csv`
- `Logs/DECISION_LOG.tsv`
