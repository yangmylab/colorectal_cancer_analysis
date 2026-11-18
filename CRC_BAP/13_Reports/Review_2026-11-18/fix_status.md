# Revision status against the review (2026-11-18)

| ID | Severity | Status | Note |
|---|---|---|---|
| R01 | P0 | **OPEN** | No bibliography yet; needs a verified reference list (~40 entries) |
| R02 | P0 | **DRAFTED** | Statements section added; Funding, Author contributions, ethics approvals and code DOI are `[TODO]` |
| R03 | P0 | **PARTIAL** | Rendered at CRC_PANEL_MAX_W=3.4: Figures 2/4/5/6 now 7.33/7.33/11.25/7.34 pt. Still below floor: Fig 1 3.10 (discovery chain not re-run, label risk), Fig 3 4.01 (20_spatial_analysis aborted on a curl download), S1 5.08, S2 6.35 |
| R04 | P0 | **OPEN** | Group sizes and per-model n still missing |
| R05 | P0 | **OPEN** | Legends still carry no n/test/interval definition |
| R06 | P1 | **FIXED** | Title and Abstract Conclusions now say "associated with" and "no detectable incremental value" |
| R07 | P1 | **FIXED** | Abstract names the three cohorts and GSE14333 explicitly |
| R08 | P1 | **FIXED** | "without overlap" -> "substantial but incomplete separation" |
| R09 | P1 | **REWORDED** | "equally common" -> "similar ... no formal test reported" (no test invented) |
| R10 | P1 | **FIXED** | "Fourteen findings (F01-F15; there is no F10)" |
| R11 | P1 | **OPEN** | Methods expansion (parameters, CV scheme, software versions) still to write |
| R12 | P1 | **CODE FIXED** | Library-size panel now by TCGA project with a Wilcoxon subtitle; takes effect on re-render |
| R13 | P1 | **CODE FIXED** | Heatmap annotation renamed "project"; takes effect on re-render |
| R14 | P1 | **CODE FIXED** | Shared ECOTYPE_COL palette applied to 7 scripts; takes effect on re-render |
| R15 | P1 | **OPEN** | KM at-risk tables / log-rank p not yet added |
| R16 | P1 | **OPEN** | Forest plots not yet annotated |
| R17 | P1 | **OPEN** | Permutation null for label agreement not yet computed |
| R18 | P2 | **OPEN** | Figure 6D y-axis |
| R19 | P2 | **OPEN** | Decision-curve references |
| R20 | P2 | **OPEN** | Per-panel n in legends (folds into R05) |
| R21 | P2 | **CHECKED** | 0.558 vs 0.56 is consistent rounding; p=0.997 is the exact value from the TSV |
| R22 | P2 | **FIXED** | "five of them with unchanged kinase abundance" added |
| R23 | P2 | **OPEN** | Recoverable-feature reporting |
| R24 | P2 | **OPEN** | Figure 6B label offsets |
| R25 | P2 | **OPEN** | Grade tags in Results |
| R26 | P2 | **FIXED** | "947-4,757 spots per section (median 2,670)" |

Next three actions, in order: (1) run the panel re-render and recompose; (2) build the reference list; (3) expand Methods and rewrite the legends with n and tests.
