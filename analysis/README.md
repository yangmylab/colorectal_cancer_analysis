# Analysis skeleton — Track A

**Status:** scaffold only. No data downloaded, no software installed, no analysis run.
**Protocol:** `../PROTOCOL_TrackA.md`
**Pilot:** `../PILOT_PLAN.md`
**Datasets:** `../DATASETS.md`
**Software:** `../SOFTWARE.md`

---

## 1. Purpose

This directory is the working area for the Track A analysis. It is intentionally empty of data and results until the protocol and pilot are approved.

The scaffold provides:
- a single source of truth for parameters (`config/study_config.yaml`);
- a sample-sheet template (`metadata/sample_sheet_template.tsv`);
- a planned module layout (`scripts/README.md`);
- conventions for reproducibility and provenance.

---

## 2. Directory layout

```
analysis/
  README.md
  config/
    study_config.yaml
  metadata/
    sample_sheet_template.tsv
    manifest.tsv              # created after download
    overlap_report.tsv        # created after cohort assembly
  scripts/
    README.md                 # planned module entry points
  data/
    raw/<source>/<accession>/ # user-downloaded, immutable
  processed/
    rna/ methylation/ proteome/ phospho/ cnv/ mutation/
  reference/
    signatures/ programs/ atlas/
  results/
    deconvolution/ programs/ integration/ ecotypes/ biology/
    spatial/ clinical/ drug/
  figures/
  manifests/
  logs/
```

Create subdirectories as needed; never write into `data/raw/`.

---

## 3. How to use the config

- `config/study_config.yaml` is the only place for paths, cohort definitions, QC thresholds, model parameters and seeds.
- Every script must read this file; do not hard-code paths or parameters.
- Freeze the config (git commit + checksum) before outcome modeling.
- Changes after freezing must be recorded in `manifests/config_changelog.tsv` with reason and date.

---

## 4. Sample sheet rules

- One row per sample; one patient may have multiple samples only if explicitly needed (e.g., tumor/normal).
- `sample_id` must be globally unique; `patient_id` must be stable across cohorts.
- Use ISO-like dates and days for time-to-event variables.
- `included` is `TRUE`/`FALSE`; `exclusion_reason` is mandatory when `FALSE`.
- Do not overwrite the template; copy it to `metadata/sample_sheet.tsv`.
- Record platform, assay, modality, batch and source accession for every sample.
- Survival fields: `*_time_days`, `*_event` (1 = event, 0 = censored).
- Purity: record method (ESTIMATE/ABSOLUTE/pathology) and value.
- Stage: record AJCC version; harmonize to a common stage variable in a separate column if needed.

---

## 5. Planned module entry points

See `scripts/README.md` for the complete list. Naming convention:

```
NN_<module>.<ext>        # executable entry point
NN_<module>_utils.<ext>  # helpers
```

Every script must:
1. read `config/study_config.yaml`;
2. write to `results/` or `processed/`, never to `data/raw/`;
3. write a log to `logs/` with start/end time, versions and seed;
4. write a session info/version file;
5. be idempotent (safe to rerun).

---

## 6. Reproducibility rules

- Pin R (`renv.lock`) and Python (`environment.yml`) environments.
- Use containers for the final workflow.
- Fix seeds (config `project.seed`).
- Record checksums for raw and processed files in `manifests/`.
- Keep intermediate objects (`.rds`, `.h5ad`, `.zarr`) with the code version that produced them.
- Do not upgrade packages mid-analysis; use module-specific environments if needed.
- A fresh clone must reproduce the key figures/tables.

---

## 7. What not to do

- Do not download controlled-access data without explicit user action and a DUA.
- Do not merge cohorts without checking patient overlap.
- Do not select features or tune models on the full dataset.
- Do not use a validation cohort that shares patients with discovery.
- Do not make causal or clinical-efficacy claims from these data.
- Do not commit large data files to git; use manifests/checksums and external storage.

---

## 8. Current status

- [ ] Protocol approved
- [ ] Datasets downloaded (Phase 1)
- [ ] Software installed
- [ ] Pilot completed
- [ ] SAP registered
- [ ] Config frozen
