# Metadata Audit Report — CRC Track A

Produced under **G1** of the BAP protocol. This document records the evidence that group assignment is defensible.

- Datasets audited: TCGA-COAD and TCGA-READ (GDC Data Release 46.0)
- Audit date: 2026-11-17
- Annotation table: `Sample_Annotation.csv` (624 rows, 624 donors)
- Machine check: `bap.py check-annotation` and `bap.py gate g1` — **PASS**

---

## 1. Sources consulted

| Source | Retrieved? | What it contributed |
|---|---|---|
| Repository metadata (GDC clinical/case API) | yes | sample type, tissue type, stage, vital status, days_to_death, days_to_last_follow_up, age, sex_at_birth, site of resection |
| Original TCGA publications | partial | TCGA-COAD/READ cohort context; not used to override GDC fields |
| Supplementary sample sheet / series matrix | not needed | GDC case/file metadata used as source of truth |
| Data-internal evidence | yes | sample ID matching, MD5 verification of 701 RNA-seq files, duplicate-patient collapse, expression distribution/QC |

## 2. Final group assignment

| Group | n samples | n donors | Biological replicates | Technical replicates | Definition used |
|---|---|---|---|---|---|
| Primary tumour | 624 | 624 | 624 | 0 (collapsed) | GDC `sample_type = Primary Tumor` |
| Solid tissue normal | 51 files downloaded | — | — | — | not used in discovery |
| Metastatic / recurrent | 3 files | — | — | — | not used in discovery |

The 647 primary-tumour RNA-seq files were collapsed to **624 patients** by choosing, per patient, the sample with the highest total TPM. The 23 additional files are duplicate/alternate samples and are retained in `TCGA/clinical_all.tsv`.

## 3. The ten verification questions

1. **Control / disease / treatment samples:** discovery contains 624 primary tumours; 458 COAD and 166 READ. No controls in the discovery matrix. Treatment status is not consistently available.
2. **Subgroups:** stage I (96), II (196), III (150), IV (64), missing (118). MSI/MMR status is not available in the open GDC clinical table used here.
3. **Paired or unpaired:** unpaired; one primary tumour per patient after collapse.
4. **Longitudinal:** none.
5. **Multiple tissues/cell types per donor:** one primary tumour per donor in the discovery set; normal samples are separate.
6. **Biological replicates per group:** 624 donors. **Technical replicates:** collapsed; no technical replicate retained.
7. **Balance and confounding:** age and sex are available; stage is missing for 118. `batch` is the project (TCGA-COAD vs TCGA-READ), partially confounded with tissue (colon vs rectum). MSI and purity are unavailable.
8. **Ambiguous samples:** 0 critical ambiguities after ID matching; 118 missing stage and 66 missing OS time are recorded as missing, not imputed.
9. **Appropriateness of the reference group:** no healthy control in the discovery matrix; ecotype discovery is unsupervised. Normal samples are not used as a case-control reference.
10. **Do publication and metadata define groups identically?** GDC metadata is the source of truth; no conflicting publication-level group definition was used.

## 4. Conflicts found

| Sample | Field | Conflict | Sources | Resolution |
|---|---|---|---|---|
| — | — | none | GDC case/file metadata | not applicable |

## 5. Confounders and design consequences

| Variable | Association with group | Consequence for the model |
|---|---|---|
| batch/project (COAD/READ) | partial (tissue) | included as batch for QC; survival models use tissue/sidedness where available |
| stage | missing in 118 | complete-case survival model; missingness reported |
| sex | 3 missing | complete-case survival model |
| MSI | unavailable | not modelled |
| purity | unavailable | not modelled; limitation |
| tissue composition | unavailable | limitation; Hallmark scores partly reflect stroma |

## 6. Sex / identity check

- Method: GDC `sex_at_birth` (330 M / 291 F / 3 missing); no genotype-based sex check performed.
- Result: no discordance assessable; missing sex recorded.
- Duplicate detection: 647 primary-tumour files collapsed to 624 patients by total TPM; no near-duplicate expression pairs were analysed.

## 7. Gate statement

- **G1 verdict:** **PASS** (2026-11-17)
- Samples: 624 tumours / 624 donors; 2 project batches; 0 critical ambiguities.
- Unit of inference: biological sample = patient primary tumour.
- Primary contrast: unsupervised ecotype (K=2 primary) with survival association; not a case-control contrast.
- Design statement carried forward to `Analysis_Plan.md` §4.

Signed: Reasonix agent · 2026-11-17
