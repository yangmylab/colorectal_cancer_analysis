# Strict pre-submission review — CRC two-ecotype manuscript

**Manuscript:** `CRC_BAP/12_Manuscript/manuscript.docx` (draft v0.2, built from `manuscript.md`)
**Figures reviewed:** `CRC_BAP/12_Manuscript/Figures/` — `Figure_1`–`Figure_6`, `Supplementary_Figure_1`, `Supplementary_Figure_2` (44 panels)
**Review date:** 2026-11-18
**Reviewer:** pre-submission review performed on behalf of, and at the request of, the author (no independence claimed)
**Method:** `peer-review` (evidence-bounded review, Location/Observation/Evidence/Why/Fix per finding), `nature-statistics` (P0/P1/P2 severity, n and replication, figure-statistics alignment), `paper-claim-audit` (every number compared with raw result files), dispatched through `family-paper-review`

**Verdict: NOT READY FOR SUBMISSION — major revision required.**
Five blocking (P0) issues, 13 important (P1) issues, 8 polish (P2) issues. The science and the claim discipline are sound; the blockers are infrastructure and reporting, not evidence.

---

## 1. What was checked, and how

| Check | Method | Result |
|---|---|---|
| Numbers vs evidence | All 138 distinct numeric tokens extracted from the manuscript, each high-risk value compared line-by-line with its source TSV or report (`07_Validation/`, `13_Reports/`, `06_Analysis/`) | No fabricated or mis-transcribed number found |
| Claims vs grades | All 14 ledger findings (F01–F15, no F10) mapped to their Results passages | No Results sentence exceeds its assigned grade |
| Ledger integrity | `13_Reports/Findings_Ledger.csv` row count, ids, grades | 14 rows; F10 does not exist (contradicts Methods) |
| Figures vs data | Panel-level comparison with `figures_plan.tsv`, `compose_figures_typography.tsv`; Figure 6D values measured from the rendered panel and compared with `25_incremental_value.tsv` | Panel D faithful to ±0.0001; legends match panels one-to-one |
| Figure legibility | `compose_figures_typography.tsv` against a 7 pt floor | 44/44 panels below the floor (2.76–5.44 pt) |
| Completeness | Reference list, availability, ethics, funding, contributions, conflicts | All absent |

Not assessable from the supplied material: peer-review independence, journal-specific formatting, and whether the raw public data would reproduce under a third party's hands (the scripts and logs suggest they would).

---

## 2. Summary assessment

**Strengths (stated because they are unusual, not as courtesy):**

1. **Claim discipline is genuinely applied.** Every Results sentence about outcomes carries a direction, an effect size and a p value or an explicit null. I looked for a claim exceeding its grade and found none in Results; the overreach that exists is confined to the title and the Abstract Conclusions (R06).
2. **The negative result is the paper, and it is framed as such.** The dissociation — reproducible across assays, not incremental over CMS, no mortality association — is stated plainly in the Discussion and reinforced by "What we do not claim".
3. **Endpoint discipline holds.** Disease-free and mortality endpoints are never pooled; GSE38832's 9 recurrences are explicitly declared untestable rather than forced into a model.
4. **F13 was settled honestly.** A trained-reference cell2location rerun narrowed the caveat (pattern reproduced, magnitude method-dependent) instead of converting a null into a positive.
5. **Reproducibility infrastructure is better than most published work**: scripted pipeline, `DECISION_LOG.tsv`, two automated consistency checks (22/22 and 44/44), and a machine-readable claims-to-evidence map.

**That is exactly why the blockers matter.** A paper this disciplined will be judged by the parts that are currently missing, and right now a desk editor would not send it out.

---

## 3. P0 — blocking issues

### R01. The manuscript has no references at all
- **Location:** whole document (`manuscript.md` 1–309).
- **Observation:** there is no reference list, no in-text citation, no bracketed marker, no "et al." anywhere. Claims that CMS1–4, CRIS and PDS1–3 exist and are the comparators, and that taxonomies "repeatedly failed to reach clinical use", carry no source.
- **Evidence:** `rg "\[[0-9]+\]|\(20[0-9]{2}\)|et al\." manuscript.md` returns nothing; no `## References` heading.
- **Why it matters:** this is an automatic desk-reject at every target journal, and it leaves the paper's central premise (that these taxonomies exist and underperform) unverifiable.
- **Fix:** add a reference list (expect 45–60 entries) covering CMS/CRIS/PDS consensus papers, the source cohorts (TCGA, CPTAC, GSE39582/17536/14333/38832/283052/267401/132465), MSI and CIMP definitions, and the methods cited (SNF, GSVA, nnls, cell2location, KSEA). Then run `citation-audit` to verify each entry exists and supports the sentence it is attached to.

### R02. No availability, ethics, funding, contribution or conflict statements
- **Location:** end of document; the `## Pending before submission` section lists only authors and journal.
- **Observation:** absent: Data Availability, Code Availability, Ethics/consent, Funding, Author Contributions, Competing Interests, Acknowledgements.
- **Why it matters:** all six are mandatory at every candidate journal; two of them (data and code) are also the paper's strongest reproducibility asset and are currently invisible.
- **Fix:** add, with specifics already on disk — TCGA via GDC Data Release 46.0 and its data-use certification; GEO accessions GSE39582, GSE17536, GSE14333, GSE38832, GSE283052, GSE267401, GSE132465; CPTAC via PDC000116/cBioPortal; code as the repository plus a frozen archived release (Zenodo DOI) with the environment files (`Environment/packages_R.tsv`, conda export); ethics statement as secondary analysis of de-identified public data with the relevant data-access approvals.

### R03. Every figure panel prints below the legibility floor
- **Location:** all eight composites, 44/44 panels.
- **Observation:** at the 180 mm submission width the panels render their type at **2.76–5.44 pt**: Figure 3 2.76, Figures 1 and 4 3.10, Figure 2 3.24, Figures 5 and 6 4.76, Supplementary 1 5.08, Supplementary 2 5.44.
- **Evidence:** `12_Manuscript/compose_figures_typography.tsv` (measured, not estimated); the composites are now vector and typographically uniform, so the limit is the *authored* panel width (6–10 in), not the montage.
- **Why it matters:** at 3 pt, axis labels and Hallmark programme names are unreadable in print. A reviewer will ask for the figures to be redone, which delays acceptance by a full cycle.
- **Fix:** re-render each panel at ~3.4 in with a shared theme (`base_size = 8`), then recompose with the existing script; this yields 7–8 pt with no change to the analysis. Tracked as P0-4 in the project handoff.

### R04. Group sizes and model denominators are never reported
- **Location:** Abstract; Results §1–§6; figure legends.
- **Observation:** the number of tumours per ecotype (E1/E2) appears nowhere. Neither does the n behind any model beyond a few cohort totals (GSE14333 226, GSE39582 483, GSE17536 135, discovery 476). Methods promises "complete-case analysis for covariates, with counts reported" — the counts are not in the manuscript.
- **Why it matters:** class balance drives the stability and the chance level of every agreement statistic; without it a reviewer cannot judge whether ARI 0.439 is strong or trivial, and cannot assess overfitting in the nested models.
- **Fix:** add a cohort/analysis table (cohort, platform, n analysed, events, E1/E2 counts, covariates available), and give per-model n in Results or in the table.

### R05. Figure legends carry no statistics
- **Location:** `## Figure legends`, 260–296.
- **Observation:** every legend is a one-line description of panels ("(A) Hallmark programme scores ordered by ecotype."). None states n, the test or model, error-bar or interval definition, correction status, or what a dot/bar represents.
- **Why it matters:** at Nature-family venues, missing n and test definitions in legends are treated as scientific reporting gaps, not formatting. Several panels (Figure 4E/4F, Figure 6D, Supplementary 2A/2B) cannot be interpreted from the legend alone.
- **Fix:** rewrite each legend with: what is plotted, n and its unit, summary convention, test/model, correction, exact p or threshold policy, and whether source data underlie the panel.

---

## 4. P1 — important issues

### R06. Title and Abstract overstate a null result
- **Location:** title (line 1); Abstract Conclusions (31–34); Abstract Results (28–29).
- **Observation:** "adds nothing to CMS" and "added nothing over CMS" describe a test with **HR 1.42, 95% CI 0.93–2.17, p=0.107** and a C-index gain of +0.015. Likewise "predicts recurrence" describes a Cox association that was never tested as a prediction (no held-out prediction, no calibration for DFS).
- **Evidence:** `25_incremental_value.tsv` row "+ MSI + CMS": delta_c 0.01509, hr 1.4182, CI 0.9269–2.1699, p 0.10738.
- **Why it matters:** the paper's thesis is that claims must not exceed evidence; a title that converts a null into "adds nothing" and an association into "predicts" invites exactly the criticism the paper is about. It is also the single easiest attack for a reviewer.
- **Fix:** "…is associated with recurrence but not with death, and shows no detectable incremental value over CMS" in the title; "no detectable incremental value in the one cohort able to test it" in the Abstract. Keep "adds nothing detectable" (Discussion) — that phrasing is already correct.

### R07. The Abstract double-counts the cohorts
- **Location:** Abstract Methods, lines 19–21.
- **Observation:** "three independent cohorts with disease-free endpoints and three with mortality endpoints" reads as six cohorts; there are three (GSE39582, GSE17536, GSE38832), used for both endpoint families, plus GSE14333 for DFS only.
- **Fix:** "three independent cohorts (GSE39582, GSE17536, GSE38832) analysed on the endpoint each can support, plus GSE14333 for disease-free survival".

### R08. Results §1 contradicts its own figure
- **Location:** Results §1, lines 90–92.
- **Observation:** "the fused-view embedding (Figure 1C) separates the two groups **without overlap**". Figure 1C shows two overlapping point clouds with a shared central region.
- **Why it matters:** a factual statement contradicted by the panel it cites is the kind of error reviewers check first, and it undermines trust in the rest.
- **Fix:** "separates the two groups with substantial but incomplete separation (Figure 1C)" — or report a quantitative separation statistic (silhouette, or the classifier's cross-validated accuracy already quoted).

### R09. "Equally common" is asserted without a test
- **Location:** Results §3, lines 118–119.
- **Observation:** "CIMP-high status was equally common in both ecotypes (74.5% E1 versus 84.3% E2)" — a 9.8-percentage-point difference described as equal, with no test for the proportion comparison. Only global methylation has a p value (0.113).
- **Fix:** report the CIMP comparison with a test and CI (or remove "equally common" and say "similar, with no formal test performed"), and align the precision with the source report's "differ only modestly".

### R10. Methods cites a finding that does not exist
- **Location:** Methods `Grading`, line 76.
- **Observation:** "Fifteen findings (F01–F15) were graded". The ledger has **14** rows: F01–F09, F11–F15. There is no F10.
- **Fix:** "Fourteen graded findings (F01–F15, no F10)" or restore F10 with an entry; state which is intended.

### R11. Methods cannot be reproduced from the text
- **Location:** Methods, 56–80.
- **Observation:** six short paragraphs, missing: SNF parameters (neighbours, iterations, views), consensus k-means settings (runs, bootstraps, resampling fraction), the K-selection rule and its threshold, the cross-validation scheme behind "accuracy 0.979", normalisation and any batch handling, FDR procedure per analysis family, software and package versions, random seeds, and how missing clinical covariates were handled beyond "complete-case".
- **Evidence:** the values exist in `DECISION_LOG.tsv` and the scripts but are not in the manuscript.
- **Fix:** expand to a conventional Methods with a parameters table; quote versions (R 4.6.1, GSVA 2.6.6, limma 3.68.5, survival 3.8.12, cell2location 0.1.5, pyro-ppl 1.8.6, scanpy 1.11.5).

### R12. Supplementary Figure 1C cannot support the sentence it is cited for
- **Location:** Results §7, lines 203–204; Supplementary Figure 1C.
- **Observation:** "Library sizes were comparable between groups (Supplementary Figure 1C)". Panel 1C is a **single-group** boxplot (all samples labelled `tumour`); no between-group comparison exists in the panel.
- **Evidence:** `01_bulk_qc.R:122` plots `aes(group, library_size)` and `group` is `"tumour"` for every sample (`01_bulk_qc.R:79`).
- **Fix:** colour panel C by TCGA project (COAD versus READ) and state the test, or drop the panel and the sentence.

### R13. Supplementary Figure 1D labels the project as "batch"
- **Location:** Supplementary Figure 1D legend key; Results §7 lines 197–203.
- **Observation:** the heatmap's column annotation is named **"batch"** with levels TCGA-COAD/TCGA-READ, while the manuscript — correctly, and after a dedicated audit fix — states that no technical batch exists.
- **Why it matters:** it is the exact conflation the QC work was corrected to avoid; a reviewer who reads "batch" will ask why it was not corrected.
- **Fix:** rename the annotation to "project" (or "cohort") in the plotting code and regenerate.

### R14. The same two ecotypes are encoded in four different colour schemes
- **Location:** Figure 1C (salmon/teal), Figure 2A/2D (teal/salmon annotation bars), Figure 3A/3B (grey/red), Figure 5A/5D/5F (black/red), Figure 6A (blue/red), Figure 6E (grey/dark red), Supplementary 2D (black/dark red).
- **Why it matters:** E1 and E2 swap colours between figures, so a reader cannot track a group across the paper — and in Figure 1C versus Figure 3A the *same* colours mean different groups.
- **Fix:** one palette for E1/E2 across all figures (and for cell types within Figure 1F), applied through a shared theme.

### R15. Kaplan–Meier panels lack numbers at risk and the log-rank test
- **Location:** Figures 5A, 5D, 5F, 6A, Supplementary 2C, 2D.
- **Observation:** no at-risk table, no log-rank p (except Supplementary 2C), no n per curve.
- **Fix:** add numbers at risk beneath each curve and annotate the log-rank p; state in the legend whether the quoted HR comes from a Cox model or the log-rank test (Supplementary 2C currently mixes both: figure log-rank p=0.27, text Cox p=0.286).

### R16. Forest plots print no numbers
- **Location:** Figures 5B, 5C, 5E, 6C.
- **Observation:** hazard ratios and CIs are only encoded positionally; no numeric labels per row.
- **Fix:** annotate "HR (95% CI)" per row, with n and events in a side column.

### R17. "Far above chance" is asserted, not quantified
- **Location:** Results §2, lines 106–107.
- **Observation:** "83.3% … imperfect, but far above the chance level implied by two groups" — the chance level is never computed, and with unequal group sizes it is not 50%.
- **Fix:** report the permutation null for label agreement (the project already computes permutation nulls for ARI), or quote the expected agreement under independence from the observed group sizes.

---

## 5. P2 — clarity and polish

| ID | Location | Issue | Fix |
|---|---|---|---|
| R18 | Figure 6D | y-axis truncated at 0.50–0.65, so 0.01–0.03 C-index differences look large | plot ΔC-index per model, or show the full 0.5–1.0 range as an inset |
| R19 | Supplementary 2B | decision curve has no "treat all"/"treat none" reference lines, so net benefit is uninterpretable | add both references, or drop the panel given the null model |
| R20 | Figures 2C vs 2E | denominators differ (96 tumours for mRNA–protein, 97 for the phospho layer) with no panel-specific n | state per-panel n in the legend |
| R21 | Results §2, §4 | ARI quoted as 0.558 in text and 0.56 in the report; p as 0.997 in text and "1" in the report | fix a precision policy and apply it everywhere |
| R22 | Results §4 | F12 records that 5 of 9 kinases are abundance-independent; Results says only "9 kinases at FDR<0.05" | add "5 of which were abundance-independent" |
| R23 | Methods | "recoverable features are reported for each layer" — they are in the layer reports, not the manuscript | either report them in a supplementary table or reword |
| R24 | Figure 6B | p-value labels overlap the K=4 point | offset the labels |
| R25 | Title/Abstract | grades (E3/E2/E1) are used in the Abstract but never shown in Results | tag each Results subsection with its grade, or add a grade column to the claims table |
| R26 | Results §4 | "~2,570 spots per section" — the exact per-section counts are available | quote the range (947–4,757 spots; median 2,670) |

---

## 6. Claim-by-claim audit (F01–F15)

| Finding | Grade | Manuscript location | Reported statistic matches source? | Claim within grade? |
|---|---|---|---|---|
| F01 OS not incremental (discovery) | not_supported | §6 | Yes — C 0.763→0.772, p=0.107 | Yes |
| F02 K=4 OS (post hoc) | E1 | §6 | Yes — C 0.785, p=0.024, labelled post hoc | Yes |
| F03 OS external null | not_supported | §6 | Yes — 0.333 / 0.0785 / meta 1.25, p=0.081 | Yes |
| F04 K=4 not replicated | not_supported | §6 | Yes — qualitative, correct | Yes |
| F05 single-cell / deconvolution | descriptive | §1 | Yes — six lineages, fractions sum to 1.001 | Yes |
| F06 protein reproduction | E3 | §2 | Yes — 83.3%, ARI 0.439, 96 tumours | Yes |
| F07 methylation non-reproduction | not_supported | §3 | Yes — ARI 0.006, 5/19,269, p=0.113 | Yes (see R09) |
| F08 no transfer to spots | not_supported | §3 | Yes — ARI 0.042, ≈5 SD | Yes |
| F09 phospho reproduction | E3 | §2, §4 | Yes — 87.5%, ARI 0.558, p=0.000999 | Yes |
| F11 disease-free replication | E2 | §5 | Yes — 2.37 / 1.50 / 2.51, meta 1.67 | Yes |
| F12 kinase activity | E1 | §4 | Yes — 78 testable, 9 at FDR<0.05, CDK7 d=−1.00 | Yes, but see R22 |
| F13 spatial structure | descriptive | §4 | Yes — I 0.528 vs 0.587, ρ 0.49 | Yes, caveat correctly narrowed |
| F14 no incremental value | not_supported | §6 | Yes — 1.42 (0.93–2.17), p=0.107, V=0.52 | Yes; wording issue only (R06) |
| F15 no cancer-specific death association | not_supported | §5 | Yes — 1.00 (0.46–2.17), p=0.997 | Yes |

---

## 7. Author input needed

These are factual questions only the author can answer; I have not substituted guesses.

1. **Are the E1/E2 counts available for every cohort and model?** (They are derivable from the label files — say the word and I will tabulate them.)
2. **What cross-validation scheme produced accuracy 0.979?** (`DECISION_LOG` records the value, not the scheme: k, folds, stratification.)
3. **Was a test ever run on the CIMP proportions?** If not, I will reword rather than invent one.
4. **Which repository and DOI should the Code Availability statement point to?**
5. **Does the TCGA/CPTAC data-use certification permit the intended statement**, and which approvals should be named in the ethics section?
6. **For Supplementary 1C: shall I re-colour it by project, or remove the panel?**

---

## 8. Reviewer-risk note

If this went out today, the first three reviewer comments would be, in order: (i) "no references"; (ii) "figures are unreadable at print size"; (iii) "the title claims the axis adds nothing to CMS, yet the reported CI (0.93–2.17) is compatible with a clinically meaningful effect — this is absence of evidence, not evidence of absence." Points (i) and (ii) are mechanical to fix; (iii) is a wording change. After those, the sharpest scientific question is the one the paper already answers honestly: whether a single CMS-carrying cohort can bear the incremental-value conclusion at all. Stating that limitation in the Abstract, not only in the Discussion, would pre-empt it.

---

## 9. Reproducibility of this review

- Numbers: extracted programmatically from `manuscript.md` (138 tokens) and compared with `07_Validation/*.{tsv,md}`, `13_Reports/Findings_Ledger.csv`, `06_Analysis/**`.
- Figure values: Figure 6D dots measured from `09_Figures/plots/plot_41.jpg` by self-calibrating on the ten bar tops (fit residual 0.32 px) and compared with `25_incremental_value.tsv`; agreement ±0.0001.
- Legibility: `12_Manuscript/compose_figures_typography.tsv`.
- A machine-readable version of every finding is in `review_findings.tsv` beside this file.
