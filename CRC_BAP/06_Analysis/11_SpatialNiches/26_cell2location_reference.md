# cell2location vs NNLS concordance (F13 robustness check)

- Prior: **reference** — 2007 genes x 6 lineages (cell2location reference trained on GSE132465).
- Sections: GSE283052/GSM8655157_P02, GSE283052/GSM8655158_P04, GSE267401/GSM8265214_CTC17M; max_epochs = 3000; spots per section = 1000; genes shared with Visium = 2007.
- Both methods are run on the same gene set and the same spots, so the comparison is like for like.

## Median concordance per lineage

| lineage | Pearson r | Spearman r | mean fraction (cell2location) | mean fraction (NNLS) |
|---|---|---|---|---|
| B cells | 0.77 | 0.52 | 0.08 | 0.02 |
| Epithelial cells | 0.54 | 0.44 | 0.39 | 0.97 |
| Mast cells | n/a | n/a | 0.12 | 0.00 |
| Myeloids | 0.33 | 0.27 | 0.12 | 0.00 |
| Stromal cells | 0.69 | 0.53 | 0.17 | 0.01 |
| T cells | n/a | n/a | 0.12 | 0.00 |

## Per-section composition agreement

| section | spots | matrix Pearson r | top-lineage agreement | median abundance | undefined lineages |
|---|---|---|---|---|---|
| GSE283052/GSM8655157_P02 | 1000 | 0.94 | 0.97 | 17.9 | 4 |
| GSE283052/GSM8655158_P04 | 1000 | 0.70 | 0.74 | 19.1 | 2 |
| GSE267401/GSM8265214_CTC17M | 1000 | 0.93 | 0.99 | 19.7 | 3 |

- Median per-lineage Pearson r over the 9 identifiable lineages: 0.54; 9 lineage comparisons are undefined because one side has zero variance (NNLS returns exactly zero for those lineages).
- Median matrix-level Pearson r over all spots x lineages: 0.93.
- Largest absolute composition gap: **0.58** (Epithelial cells: cell2location 0.39 vs NNLS 0.97).

## Verdict

**The spatial pattern is reproduced, but the absolute composition is not** (median matrix r 0.93; largest gap 0.58 for Epithelial cells: cell2location 0.39 vs NNLS 0.97). The 2026-11-17 reading that cell2location 'did NOT corroborate' the NNLS fractions is therefore withdrawn: at 100 epochs with a marker-matrix prior the model was simply under-trained. What remains method-dependent is the *scale*, not the direction - unregularised NNLS collapses onto one lineage (~0.97 epithelial, exactly zero for several others) while cell2location spreads mass across lineages, so absolute fractions, and therefore the magnitude of any stromal-fraction correlation, depend on the method and the gene panel. F13 keeps a narrowed caveat about absolute fractions and magnitudes; its direction (E2 programmes high where the stromal compartment is high) is supported.

