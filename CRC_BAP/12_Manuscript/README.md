# Manuscript workspace

Draft v0.1 (2026-11-17), generated from the graded evidence ledger.

## Files

| File | Purpose |
|---|---|
| `manuscript.md` | Main draft: abstract, introduction, results, discussion, methods, figure legends, pending items |
| `claims_to_evidence.tsv` | Every claim mapped to its finding id, evidence grade, datasets, figures, ledger statistic and a citation caveat |
| `manuscript.docx` | Word rendering of `manuscript.md` (generated; regenerate after edits) |

## The rule this manuscript follows

No sentence may assert more than its finding's grade allows. Grades in the ledger
(`13_Reports/Findings_Ledger.csv`, F01–F15):

- **E3** cross-platform convergence across assays — F06 (protein), F09 (phosphoproteome)
- **E2** association replicated in an independent cohort — F11 (disease-free survival, 3 cohorts)
- **E1** exploratory or post hoc — F02 (K=4), F12 (kinase activity, inferred)
- **descriptive** — F05 (single-cell deconvolution), F13 (spatial structure, no outcome)
- **not_supported** — F01, F03, F04 (OS), F07 (methylation), F08 (spatial transfer), F14 (no incremental value), F15 (cancer-specific death)

Two constraints matter most when editing:

1. **F11 must never appear without F14.** The disease-free association replicates (E2) *and* is not incremental
   over CMS. Citing the first without the second would misrepresent the project.
2. **Disease-free endpoints and death endpoints are never pooled.** DFS/RFS (F11) and OS/CSS/DSS (F01, F03, F15)
   come from different endpoint families.

## Regenerating

```bash
# rebuild the claim map after editing the ledger or the section structure
python3 Scripts/build_claims_map.py

# rebuild the Word file
pandoc manuscript.md -o manuscript.docx --toc --toc-depth=2
```

## Pending before submission

1. **Author list** — `Authors.xlsx` (project root) is still an empty template: names, emails, affiliations, ORCIDs.
2. **Journal choice** — drives reference style, word limits and supplementary structure.
3. **Figure assembly** — the 44 plots exist as `09_Figures/plots/plot_N.{pdf,jpg}` with a manifest; a submission
   figure set will need combining/relabelling into main + supplementary.
4. Optional scientific strengthening: a properly trained cell2location reference (GSE132465 labels and TPM matrix
   are on disk) to test whether F13's E2–stromal co-localisation survives a second deconvolution method.
