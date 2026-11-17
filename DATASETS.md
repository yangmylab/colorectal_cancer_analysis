# Dataset download checklist (manual)

**Project:** Colorectal cancer multi-omics strategy
**Prepared:** 2026-11-17
**Status:** NOTHING HAS BEEN DOWNLOADED. Per project standing instruction, all datasets are listed here with source, size estimate and access route for **manual download by the user**. Do not start bulk downloads until the study track and target journal are confirmed.

---

## New resources added in the 2026 scan

The detailed, verified inventory is in **`DATASETS_EXTENDED.md`**. The highest-impact additions are:

| Resource | Content | n | Access |
|---|---|---|---|
| **PXD079707** | DIA proteomics, targeted genomics, WSIs, clinical follow-up | 1,000 treatment-naive CRCs | PRIDE, CC0 |
| **PXD052391** / IPX0002188000 | Proteome + phosphoproteome across progression stages | 148 patients, 435 samples | Proteomics open; WES restricted |
| **IPX0007391001** / PRJCA020890 | Primary/normal/liver-metastasis proteogenomics | 34 patients, 102 samples | iProx/GSA |
| **PXD019103** | Colon cancer proteome, tumor/normal | study-level | PRIDE |
| **GSE241101** | Multiregional primary CRC | 286 | GEO |
| **GSE263030** | cfDNA 5hmC, occult CRC | 602 | GEO |
| **GSE97923** | cfDNA methylation, COAD | 773 | GEO |
| **GSE140973** | Anti-EGFR secondary resistance | 114 | GEO |
| **GSE333907** | Early-onset CRC cohort | 83 | GEO |

**Before downloading:** confirm processed tables vs raw files, metadata completeness, and license. PXD079707 and PXD052391 may require DIA proteomics processing; download a subset first.

---

## 0. Conventions

- Record for every file: accession, release/version date, download date, checksum, license/DUA, and the exact command used.
- Suggested local layout:
  ```
  data/
    raw/<source>/<accession>/
    processed/<modality>/<accession>/
    metadata/
    manifests/
  ```
- Prefer **processed matrices** for single-cell/spatial atlases; download raw FASTQ/BAM only when re-processing is required.
- Check patient overlap between GEO series before using them as independent validation.
- Keep controlled-access data (dbGaP/EGA/Broad Single Cell Portal) separate; the core study is designed to run on open data only.

---

## 1. Discovery cohorts (must have)

### 1.1 TCGA-COAD and TCGA-READ
- **Source:** NCI Genomic Data Commons (GDC), https://portal.gdc.cancer.gov/
- **Content:** RNA-seq, 450K methylation, copy number, WXS/WGS, miRNA-seq, RPPA, H&E diagnostic slides, clinical.
- **Verified counts (GDC API, 2026-11-17):**
  - TCGA-COAD: 461 cases, 29,242 files; TCGA-READ: 172 cases, 10,387 files.
  - Pooled data-type facets: 701 gene-expression quantification, 630 isoform expression, 630 miRNA expression, 2,202 methylation array, 1,890 gene-level copy number, 495 RPPA, 41 ATAC-seq, 1,972 slide-image files (625 diagnostic-slide files).
- **Access:** Open (some germline/controlled files require dbGaP; not needed for this study).
- **Size estimate:** 30-80 GB for selected open data types (RNA/methylation/CNV/clinical/RPPA); full slide images add ~100+ GB.
- **Download:**
  1. Create a cart/project filter in the GDC portal.
  2. Download the `manifest.txt` (and `metadata.json`).
  3. Install and run `gdc-client`:
     ```bash
     gdc-client download -m manifest.txt -d data/raw/tcga
     ```
  4. Verify file counts and MD5 sums against the manifest.
- **Notes:** Use `sample_type` codes: 01 primary tumor, 11 solid tissue normal, 06 metastatic. Exclude duplicates (same patient, multiple aliquots).

### 1.2 CPTAC colon cancer (proteogenomics)
- **Source:** NCI Proteomic Data Commons (PDC), https://pdc.cancer.gov/
- **Study:** CPTAC colon cancer prospective cohort (Vasaikar et al., *Cell* 2019, PMID 31031003). PDC study ID **PDC000116** - confirm the current ID/link on the PDC portal before download.
- **Content:** Tumor and matched normal proteome and phosphoproteome; genomic/transcriptomic layers.
- **n:** 110 prospectively collected colon tumors.
- **Access:** Open (PDC).
- **Size estimate:** 20-60 GB depending on level (raw vs. processed).
- **Download:** PDC portal study page -> "Download" manifest; or use the PDC API/client if available. Prefer processed (`*_proteome_quant`, phosphosite) tables for integration.
- **Notes:** Record the exact PDC data level and normalization; do not mix reprocessed and original tables without harmonization.

---

## 2. Single-cell / single-nucleus reference atlases

| Accession | Content | n | Est. size | Access | Priority |
|---|---|---|---|---|---|
| GSE132465 | scRNA-seq, 23 Korean patients (63,689 cells) | 23 patients | 1-3 GB processed | Open (GEO) | Core |
| GSE144735 | scRNA-seq, 6 Belgian patients (27,414 cells) | 6 patients | 1-3 GB processed | Open | Core |
| GSE178341 | scRNA-seq MMRd/MMRp atlas | 181 samples | 10-30 GB processed | Open | Core |
| GSE166555 | scRNA-seq normal/adjacent tumor | 25 samples | 2-5 GB | Open | Core |
| GSE188711 | scRNA-seq left vs right CRC | 6 samples | 1-3 GB | Open | Optional |
| GSE201348 | scRNA + scATAC, polyp-adenoma-CRC | 72 samples | 20-50 GB | Open | Core |
| GSE231559 | scRNA-seq, MD Anderson CRC Moon Shot | 26 samples | 5-15 GB | Open | Core |
| GSE312983 | snRNA-seq across age continuum | 189 tumors | 10-30 GB | Open | Core (EOCRC) |
| GSE285117 | Neutrophil atlas (4.27M-cell study) | series-level | 5-20 GB | Open | Optional |

**Download (GEO processed files):**
```bash
# Example: list and fetch supplementary files
wget -r -np -nH --cut-dirs=3 \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/
```
- Use the GEO accession page "Supplementary file" links; processed matrices are usually `.txt.gz`, `.h5ad`, `.rds`, or `.mtx.gz`.
- For raw FASTQ (only if needed), use SRA Toolkit (`prefetch`, `fasterq-dump`) or ENA.
- Check the GEO series for **raw data availability** (e.g., GSE132465 states raw data are not provided for privacy; processed matrices are available).

**Notes:**
- Harmonize cell-type labels to a common ontology before integration.
- Check for patient overlap across GSE132465/GSE144735/GSE178341 (different institutions, but verify metadata).
- GSE201348 is the key chromatin resource; confirm whether scATAC fragments are available.

---

## 3. Spatial transcriptomics cohorts

| Accession | Content | n | Est. size | Access | Priority |
|---|---|---|---|---|---|
| GSE225857 | scRNA + spatial, primary CRC and liver metastases | 8 samples | 5-20 GB | Open | Core |
| GSE283052 | Spatial transcriptomics of CRC tissues | 12 samples | 5-20 GB | Open | Core |
| GSE320042 | Spatial ecotypes, high-resolution | 38 samples | 10-40 GB | Open | Core |
| GSE284989 | Spatial genomics, cholesterol/ICI resistance | 22 samples | 5-20 GB | Open | Core |
| GSE280311 | High-definition spatial immune profiling | 12 samples | 5-20 GB | Open | Optional |
| GSE280312 | High-definition spatial immune profiling | 16 samples | 5-20 GB | Open | Optional |
| GSE270745 | Spatial atlas, DMBT1/dysplasia | 6 samples | 2-10 GB | Open | Optional |
| GSE217414 | Liver metastases from CRC | 4 samples | 2-10 GB | Open | Optional |

**Download:** same GEO supplementary-file method; for 10x Visium, spatial folders contain `filtered_feature_bc_matrix.h5`, `positions.csv`, and H&E images. Keep the image files for spatial registration.
**Notes:** Platform mix (Visium, Stereo-seq, Xenium-like) - record platform and spot/bin size; do not pool without platform-aware normalization.

---

## 4. Bulk transcriptomic validation cohorts

| Accession | Content | n | Platform | Access |
|---|---|---|---|---|
| GSE39582 | Colon cancer molecular subtypes/survival | 585 | Affymetrix U133 Plus 2.0 | Open |
| GSE17536 | Moffitt colon cancer metastasis profile | 177 | Affymetrix | Open |
| GSE14333 | 290 primary CRCs | 290 | Affymetrix U133 Plus 2.0 | Open |
| GSE38832 | Colon cancer | 122 | Affymetrix | Open |
| GSE87211 | Rectal cancer +/- controls | 363 | Affymetrix | Open |
| GSE103479 | Stage II/III treated vs untreated | 156 | Affymetrix | Open |

**Download:** GEO supplementary series matrix + platform annotation.
```bash
wget -r -np -nH --cut-dirs=3 \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/
```
**Notes:**
- Re-map probes to current gene annotations; verify sample-level stage/MSI/sidedness from series metadata or associated publications.
- Check for cohort overlap (e.g., GSE14333 and GSE17536 are independent; confirm).
- GSE87211 is rectal-enriched; use for sidedness/rectal-specific validation, not as a pan-colon cohort.

---

## 5. Clinical-genomic validation (cBioPortal)

- **Sources:** cBioPortal (https://www.cbioportal.org/): TCGA PanCancer Atlas, MSK-IMPACT, AACR GENIE.
- **Content:** mutations, copy number, limited clinical outcomes; useful for mutation co-occurrence, real-world metastatic context, and classifier scoring where expression is available.
- **Access:** Open (some studies require registration).
- **Download:** cBioPortal "Download" tab -> data files (`data_clinical_patient.txt`, `data_mutations.txt`, expression if available) or use the cBioPortal API.
- **Notes:** MSK-IMPACT/GENIE have no expression for most samples; use mutation/CNV-based validation only. Do not treat as transcriptomic validation.

---

## 6. Microbiome resources (optional companion)

- **curatedMetagenomicData (Bioconductor):** uniformly processed taxonomic and functional profiles from published cohorts, including CRC case-control studies (e.g., Zeller 2014, Feng 2015, Wirbel 2019, Thomas 2019, Yachida 2019, Yu 2017, Vogtmann 2019).
  ```r
  # After installing curatedMetagenomicData (see SOFTWARE.md)
  library(curatedMetagenomicData)
  # sampleMetadata() lists studies and sample-level covariates
  ```
- **Raw metagenomes:** ENA/SRA via the original publications.
- **Access:** Open; some cohorts have usage terms.
- **Size:** tens of GB depending on selected studies.
- **Notes:** Record profiling pipeline/version; batch-correct with MMUPHin; adjust for BMI, diet, antibiotics, geography.

---

## 7. Metabolomics (optional)

- **Sources:** Metabolomics Workbench (https://www.metabolomicsworkbench.org/), MetaboLights (https://www.ebi.ac.uk/metabolights/).
- **Content:** mostly serum/plasma and mouse models; human CRC tissue metabolomics is sparse.
- **Access:** Open.
- **Use:** targeted validation of predicted metabolites, not discovery.
- **Notes:** Verify species, matrix (tissue vs. plasma), and platform before pooling. There is no large, uniformly processed human CRC tissue metabolome comparable to TCGA/CPTAC.

---

## 8. Imaging and orthogonal resources (optional)

| Resource | Content | Access | Use |
|---|---|---|---|
| TCGA diagnostic slides (GDC) | H&E, 625 diagnostic slides across COAD/READ | Open | Morphology, digital pathology, ecotype correlation |
| NCI Imaging Data Commons (IDC) | Radiology/pathology imaging collections | Open | Optional imaging validation |
| Human Protein Atlas (HPA) | IHC for candidate proteins | Open | Orthogonal protein validation |
| HCA Gut Cell Atlas | Normal colon/rectum single-cell reference | Open | Normal baseline |

---

## 9. Gated/controlled sources (NOT required for the core study)

Listed for completeness; obtain DUAs and do not include in the core plan unless needed:
- dbGaP (TCGA controlled-access germline/raw data)
- EGA (some European CRC cohorts)
- Broad Single Cell Portal (some atlases)
- Synapse/HTAN (check per-dataset terms)

If any of these become necessary, the user must obtain access manually; no controlled data will be downloaded by the assistant.

---

## 10. Download order and pilot plan

**Phase 1 (week 1):**
1. TCGA-COAD/READ selected data types + clinical.
2. GSE39582 (validation).
3. GSE132465 + GSE144735 (reference).
4. GSE225857 + GSE283052 (spatial pilot).

**Phase 2 (week 2-3):**
5. CPTAC colon.
6. GSE178341, GSE201348, GSE312983, GSE231559.
7. Remaining bulk validation cohorts.
8. Remaining spatial datasets.

**Phase 3 (optional):**
9. curatedMetagenomicData cohort subset.
10. Metabolomics/imaging.

**Pilot success criteria:**
- TCGA + GSE39582 + GSE132465 can reproduce a stable 3-5 ecotype solution (bootstrap ARI > 0.7).
- At least one spatial dataset shows ecotype-specific niche enrichment.
- Cross-platform RNA-protein concordance is acceptable after normalization.
If the pilot fails, revisit the ecotype definition before downloading the full set.

---

## 11. Verification checklist (per dataset)

- [ ] Accession and version recorded
- [ ] Download date and command recorded
- [ ] Checksums verified
- [ ] Sample-level metadata extracted (tissue, stage, MSI, sidedness, age, sex, treatment)
- [ ] Patient IDs parsed; duplicates/overlaps flagged
- [ ] Platform and normalization recorded
- [ ] License/DUA recorded
- [ ] Processed vs. raw level documented
- [ ] QC report generated (sample count, missingness, outliers)
- [ ] Manifest added to `metadata/manifest.tsv`

---

## 12. Hard constraints

- **Do not** download controlled-access data without explicit user action and DUA.
- **Do not** assume a dataset's sample count from the accession alone; verify in the series matrix/metadata.
- **Do not** pool cohorts without checking patient overlap and platform.
- **Do not** use a validation cohort that shares patients with discovery.
- Record everything; the paper's reproducibility claim depends on this file.
