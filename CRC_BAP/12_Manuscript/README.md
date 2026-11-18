# Manuscript workspace

Draft **v0.2** (2026-11-17): figure set composed, Results rewritten so every panel is cited in order.

## Files

| File | Purpose |
|---|---|
| `manuscript.md` | Main draft: abstract, introduction, Results (panel-by-panel), discussion, methods, figure legends, pending items |
| `manuscript.docx` | Word rendering built with the reference template `MN/Paper_1/reference/tnr_fangsong_redcode.docx` (Times New Roman + FangSong) |
| `figures_plan.tsv` | **Single source of truth** for panel composition: 44 panels → 6 main + 2 supplementary figures, with a `supports_main` column recording which main panel each supplementary panel backs |
| `claims_to_evidence.tsv` | Every claim → finding id, grade, datasets, figures, ledger statistic, citation caveat |
| `Figures/` | The composed figures: `Figure_1` … `Figure_6`, `Supplementary_Figure_1`, `_2` (vector PDF + 300 dpi JPEG), generated here by `Scripts/compose_figures.R` |
| `compose_figures_typography.tsv` / `.md` / `.docx` | Per-panel typography audit written by the composer: authored vs printed size, scale factor, resulting body/title type in pt, and whether each panel meets the legibility floor |

## The rule this manuscript follows

No sentence may assert more than its finding's grade allows (`13_Reports/Findings_Ledger.csv`, F01–F15):

- **E3** cross-platform convergence — F06 (protein), F09 (phosphoproteome)
- **E2** replicated association — F11 (disease-free survival, 3 cohorts)
- **E1** exploratory / post hoc — F02 (K=4), F12 (kinase activity, inferred)
- **descriptive** — F05 (single-cell deconvolution), F13 (spatial structure, no outcome)
- **not_supported** — F01, F03, F04 (OS), F07 (methylation), F08 (spatial transfer), F14 (no incremental value), F15 (cancer-specific death)

Two constraints matter most when editing:

1. **F11 must never appear without F14.** The disease-free association replicates (E2) *and* is not incremental
   over CMS. Citing the first without the second would misrepresent the project.
2. **Disease-free endpoints and mortality endpoints are never pooled.** DFS/RFS (F11) versus OS/CSS/DSS
   (F01, F03, F15).

## Figure composition

| Figure | Panels | Plots | Role |
|---|---|---|---|
| Figure 1 | 6 (A–F) | 6, 7, 8, 9, 1, 13 | the axis in discovery |
| Figure 2 | 6 (A–F) | 16, 17, 19, 28, 29, 31 | reproduction at protein and phosphoproteome |
| Figure 3 | 7 (A–G) | 20, 21, 22, 24, 25, 26, 27 | failure to reproduce (methylation, spatial transfer) |
| Figure 4 | 6 (A–F) | 30, 32, 37, 38, 39, 40 | mechanism: phosphosites, kinase activity, spatial structure |
| Figure 5 | 6 (A–F) | 33, 34, 35, 36, 43, 44 | recurrence replicates, cancer-specific death does not |
| Figure 6 | 5 (A–E) | 10, 14, 15, 41, 42 | OS null, K=4 post hoc, no incremental value over CMS |
| Supplementary Figure 1 | 4 (A–D) | 2, 3, 4, 5 | quality control |
| Supplementary Figure 2 | 4 (A–D) | 11, 12, 18, 23 | diagnostics and negative controls |

Nothing is dropped: all 44 panels appear, and every panel is cited in Results with **first citations in strict
order** (Figure 1A → 1B → … → 6E → Supplementary Figure 1A → … → 2D); later back-references are allowed.
Supplementary panels are cited in a dedicated part of Results and each is explicitly tied to the main panel it
supports (the `supports_main` column of `figures_plan.tsv` records the same mapping). Panels that could overstate
(post-hoc K=4, calibration/decision-curve for a null model, underpowered KMs) sit in the supplement **and are cited
with their limitation stated**.

## Regenerating

```bash
python3 Scripts/build_claims_map.py          # rebuild the claim map from the ledger
Rscript Scripts/compose_figures.R            # recompose the 8 figures from figures_plan.tsv
python3 Scripts/check_figure_citations.py    # assert 44/44 panels cited, strictly in order
python3 Scripts/audit_analysis.py            # project-wide consistency checks
pandoc manuscript.md --reference-doc=/run/media/meng/Seqdata/MN/Paper_1/reference/tnr_fangsong_redcode.docx \
       -o manuscript.docx --toc --toc-depth=2
```

### How the composer sets typography

`Scripts/compose_figures.R` places the panel **PDF pages** (not rasters) with `pdflatex` + TikZ and writes a
true vector PDF plus a 300 dpi JPEG. Its defaults:

| Flag | Default | Meaning |
|---|---|---|
| `--mode` | `uniform` | one scale factor for every panel, so a given font size prints the same in every panel. `equal` keeps equal-width cells (fills the page, but type size then varies with each panel's authored width) |
| `--width_mm` | `180` | printed width of the composite (double-column page) |
| `--height_mm` | `240` | height budget; the layout picks the fewest columns that fit |
| `--floor_pt` | `7` | legibility floor used by the audit |
| `--gutter_mm` / `--letter_strip_mm` / `--letter_pt` | `4` / `4` / `8` | panel spacing, reserved strip above each panel, panel-letter size |

Every run rewrites `compose_figures_typography.{tsv,md,docx}`. **Read it before quoting figure legibility:** the
panels were authored at 6–10 in with ~11 pt type, so at 180 mm they print at 2.8–5.4 pt and *all 44* currently sit
below the 7 pt floor. Reaching it requires re-authoring the panels at print width (~3.4 in, `base_size = 8`), which
means re-running the layer scripts — not a montage change (see P2 in the project `handoff.md`).

## Pending before submission

1. **Author list** — `Authors.xlsx` (project root) is still an empty template: names, emails, affiliations, ORCIDs.
2. **Journal choice** — drives reference style, word limits and supplementary structure.
3. **Panel re-authoring for print legibility** — the montage is now a true vector assembly with one type scale per
   figure and uniform letters, but the panels themselves were drawn at 6–10 in, so their text prints at 2.8–5.4 pt
   (`compose_figures_typography.tsv` lists every panel against the 7 pt floor). Compliant figures need the panels
   re-rendered at ~3.4 in with `base_size = 8`.
4. Optional scientific strengthening: a properly trained cell2location reference (GSE132465 labels and TPM matrix
   are on disk) to test whether F13's E2–stromal co-localisation survives a second deconvolution method.
