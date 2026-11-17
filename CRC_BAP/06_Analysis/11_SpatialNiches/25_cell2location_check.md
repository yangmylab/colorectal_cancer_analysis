# cell2location vs NNLS concordance (F13 robustness check)

- Sections: GSE283052/GSM8655157_P02, GSE267401/GSM8265214_CTC17M; max_epochs = 100.
- Both methods use the same GSE132465-derived signature matrix (83 markers, 6 lineages).
- Purpose: check that the NNLS fractions underlying F13 are not an artefact of the linear solver.

## Median concordance per cell type

| cell type | Pearson r | Spearman r |
|---|---|---|
| B cells | 0.29 | 0.23 |
| Epithelial cells | 0.12 | 0.10 |
| Mast cells | 0.05 | 0.05 |
| Myeloids | 0.03 | 0.04 |
| Stromal cells | 0.53 | 0.35 |
| T cells | 0.02 | 0.01 |

- Overall median Pearson r = 0.11, Spearman = 0.09.

**The two methods do NOT agree at this training budget.** cell2location returns a nearly uniform fraction across lineages while NNLS returns an epithelial-dominated composition. With only 83 marker genes and 100 epochs the cell2location model is not identifiable, so this check is inconclusive: it neither validates nor refutes the NNLS fractions. Consequence for F13: the programme-level spatial statistics (Moran's I) do not depend on deconvolution, but the E2-stromal co-localisation statistic does, and must be treated as method-dependent until a properly trained cell2location reference is available.

Configuration: 400 spots per section, 100 epochs, signature matrix used directly as the cell_state_df prior (no reference model trained on the single-cell data).
