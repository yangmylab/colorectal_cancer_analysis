# Three-week de-risking pilot plan (Track A)

**Purpose:** before committing to full data acquisition and a 9-12 month analysis, test whether the ecotype framework is stable, reproducible and spatially detectable on a minimal dataset.
**Duration:** 3 weeks (can be compressed to 2 weeks if data are already on disk).
**Cost:** low; uses only open data and existing compute.
**Decision:** go / adjust / pivot, documented in a short pilot report.

---

## 1. Pilot question

> On a minimal but representative dataset, do single-cell-anchored ecotypes form a stable solution that replicates in an independent bulk cohort and shows spatial structure?

If yes, the full protocol is feasible. If no, we adjust the views/K or pivot to Track D before spending months on the full resource.

---

## 2. Minimal dataset

| Role | Dataset | What to use | Why |
|---|---|---|---|
| Discovery (bulk) | TCGA-COAD/READ | RNA-seq + clinical + methylation (optional) | Largest multi-omic discovery set |
| Validation (bulk) | GSE39582 | Expression matrix + survival | Independent, well-annotated |
| Reference (single-cell) | GSE132465 | Processed count matrix + cell annotations | Compact, well-annotated |
| Spatial spot-check | GSE225857 | One primary CRC Visium section | Tests spatial mapping quickly |
| Optional protein anchor | CPTAC colon (subset) | Proteome table for overlapping genes | Tests cross-platform concordance |
| Optional proteomic anchor (new) | PXD079707 (subset) | Processed protein matrix for 50-100 cases | Tests 1,000-CRC proteome integration and WSI linkage |

**Do not download the full resource for the pilot.** Use processed matrices and a subset of TCGA data types (RNA + clinical + optional methylation). For PXD079707, use a processed protein report if available; otherwise process a small DIA subset with DIA-NN/Spectronaut.

---

## 3. Pilot workflow (day-by-day)

### Days 1-2: Setup and QC
1. Create `analysis/` skeleton, `study_config.yaml`, sample sheet.
2. Download the pilot datasets manually (TCGA, GSE39582, GSE132465, GSE225857, CPTAC subset, optional PXD079707 subset); record checksums.
3. Parse TCGA clinical; define CSS/OS endpoints; harmonize stage/MSI/sidedness.
4. QC: sample counts, missingness, PCA, outlier audit.
5. Build the pilot manifest and freeze sample inclusion.

**Output:** `pilot/manifest.tsv`, `pilot/qc_report.html`.

### Days 3-5: Bulk processing and reference signatures
1. TCGA RNA: TPM/log2, filtering, batch check.
2. GSE39582: RMA/quantile, probe re-annotation, ComBat if needed.
3. GSE132465: QC, scVI/Harmony integration, annotation, malignant inference.
4. NMF programs for epithelial, fibroblast, myeloid, T/NK, B/plasma.
5. Export 10-20 robust programs as signatures.

**Output:** `pilot/signatures/`, `pilot/bulk_qc/`.

### Days 6-8: Deconvolution and program scoring
1. Consensus deconvolution (CIBERSORTx/EPIC/quanTIseq/MCP-counter) on TCGA and GSE39582.
2. Score programs (ssGSEA/GSVA/AUCell).
3. Check purity confounding; adjust if needed.
4. Build the program matrix for integration.

**Output:** `pilot/programs/`, `pilot/deconvolution/`.

### Days 9-11: Ecotype discovery and stability
1. MOFA+ on program scores + optional methylation/CNV views.
2. Consensus clustering K=2-6; 200-1,000 bootstraps (pilot can use 200).
3. Compute ARI/PAC/silhouette; select K.
4. Label ecotypes by top programs.
5. Train a simple classifier (nearest centroid) and apply to GSE39582.
6. Compare with CMS if available (CMSclassifier or published labels).

**Output:** `pilot/ecotypes/`, stability plots, replication ARI.

### Days 12-13: Spatial spot-check
1. Process one GSE225857 Visium section.
2. Deconvolve with the pilot signatures (cell2location or RCTD).
3. Score ecotypes per spot; detect niches (CellCharter) or simple neighborhood enrichment.
4. Test whether at least one ecotype score is spatially structured and enriched in a plausible niche.

**Output:** `pilot/spatial/`, one spatial map.

### Day 14: Pilot report and decision
1. Summarize QC, stability, replication, spatial signal, runtime.
2. Apply the go/adjust/pivot criteria (Section 5).
3. Update `study_config.yaml` with the chosen K, views and parameters.
4. Write a 2-page pilot report (`pilot/PILOT_REPORT.md`).

---

## 4. Compute and runtime

| Step | Resource | Expected time |
|---|---|---|
| Data download (pilot) | network + disk | 0.5-1 day |
| Bulk QC | 8 cores, 32 GB | 2-4 h |
| Single-cell integration (GSE132465, ~64k cells) | 16 cores, 64 GB; GPU optional | 4-8 h |
| NMF programs | 8-16 cores | 2-4 h |
| Deconvolution | 8 cores | 2-4 h |
| MOFA+ | 16 cores, 64 GB | 1-3 h |
| Consensus clustering (200 bootstraps) | 16 cores | 3-6 h |
| Spatial processing (one section) | 8-16 cores, 32-64 GB | 2-6 h |
| **Total** | | **~3-5 working days of compute** |

If GPU is unavailable, scVI still runs on CPU; subsample cells or use Harmony for the pilot.

---

## 5. Go / adjust / pivot criteria

### Go (proceed to full protocol)
- TCGA: ≥400 tumors with RNA + clinical; ≥100 CSS events.
- 3-5 ecotypes, each ≥5% of samples.
- Bootstrap ARI >0.70 for the chosen K.
- Replication in GSE39582: ARI >0.50 with discovery.
- At least one ecotype differs in survival after adjustment (nominal p<0.05; exploratory).
- Spatial spot-check: ≥1 ecotype score is spatially structured and enriched in a plausible niche.
- Runtime and memory within budget.

### Adjust (repeat pilot with changes)
- ARI 0.50-0.70: change views (add methylation/proteome), feature selection, or K range; rerun.
- Replication ARI 0.30-0.50: harmonize signatures, use rank-based scores, or restrict to shared genes.
- No spatial signal: try a different spatial dataset/platform or use descriptive spatial analysis only.
- Weak survival signal: focus on stage II/III or MSS subgroups; reconsider endpoint.

### Pivot
- ARI <0.50 after adjustments, or no interpretable ecotypes: switch to Track D (benchmark/resource paper) or a factor-continuum paper.
- No replication in any independent cohort: the ecotype taxonomy is not ready; do not proceed to the full protocol.
- Spatial data unusable: drop spatial claims; target a computational journal with bulk/proteome integration only.

---

## 6. Pilot deliverables

1. `pilot/PILOT_REPORT.md` (2 pages: data, QC, stability, replication, spatial signal, decision).
2. `pilot/manifest.tsv` and checksums.
3. Stability plots and ecotype heatmap.
4. Replication ARI table (TCGA -> GSE39582).
5. One spatial ecotype map.
6. Updated `analysis/config/study_config.yaml` with locked parameters.
7. Decision memo: go/adjust/pivot with rationale.

---

## 7. What the pilot deliberately does NOT test

- Full multi-view integration (methylation + proteome + phospho) at scale.
- All external cohorts.
- Drug repurposing.
- The portal.
- Manuscript-level figures.

These are deferred to the full protocol; the pilot only tests the core feasibility assumptions.

---

## 8. Pilot success statement (for the lab notebook)

> "On TCGA + GSE39582 + GSE132465 + one spatial section, we obtained K ecotypes with bootstrap ARI = X, replicated in GSE39582 with ARI = Y, and observed spatial enrichment for ecotype Z in niche N. Runtime was T. Decision: go/adjust/pivot."

Record this verbatim; it becomes the basis for the full study and the pre-registration.

---

## 9. Immediate checklist

- [ ] Approve pilot datasets and download them manually
- [ ] Approve software installation (see `SOFTWARE.md`)
- [ ] Create `analysis/` skeleton
- [ ] Run days 1-2 QC
- [ ] Freeze pilot manifest
- [ ] Run days 3-13
- [ ] Write pilot report and decision memo
- [ ] If go: register SAP and start full P1
