# Extended dataset inventory (verified 2026-11-17)

**Purpose:** supplement `DATASETS.md` with additional public resources discovered in the 2026 database scan.
**Counts:** `n` is the GEO `n_samples` field or the cohort size stated in the source/paper, as of 2026-11-17. Counts can change; verify on the accession page before download.
**Status:** NOTHING HAS BEEN DOWNLOADED. Per project standing instruction, all datasets are listed for manual download.

---

## 1. Proteogenomic resources (highest priority)

| Resource | Accession | Content | n | Access | Role |
|---|---|---|---|---|---|
| 1,000 treatment-naive CRC proteogenomic resource | **PXD079707** (PRIDE) | DIA proteomics, targeted genomic profiling, whole-slide histology, clinical follow-up; AI-based histology | 1,000 cases | Open, CC0 | Proteomic discovery/anchor; morphology link |
| CRC progression proteogenomics | **PXD052391** / iProx **IPX0002188000** | Proteome + phosphoproteome across nontumor, intraepithelial neoplasia, infiltration, advanced CRC; WES under restricted HRA007452 | 148 patients, 435 samples | Proteomics open; WES restricted | Adenoma-carcinoma progression axis |
| CRLM proteogenomics | iProx **IPX0007391001** / GSA **PRJCA020890** | Genomics, transcriptomics, proteomics, phosphoproteomics of primary CRC, normal colon and matched liver metastasis | 34 patients, 102 samples | Open (iProx/GSA) | Metastasis axis; C1/C2 proteomic subtypes |
| CPTAC colon cancer | **PDC000116** (confirm) | Proteome + phosphoproteome, matched normal | 110 tumors | Open | Independent proteomic validation |
| Colon cancer proteome (MSKCC) | **PXD019103** | Label-free proteome, tumor/normal | ~30 files | Open | Orthogonal proteomic validation |
| CRC secretome | **PXD078528** | Cancer-cell-intrinsic secretome rewiring | study-level | Open | Secretome/CAF-immune hypotheses |
| BOLD-100 metastatic CRC | **PXD069599** | Proteomics of treated metastatic CRC | study-level | Open | Therapy perturbation context |

**Notes.**
- PXD079707 is the single most important new resource: it provides proteomic depth at population scale and matched H&E. Use it as the proteomic anchor, but differentiate on single-cell/spatial anchoring, reproducibility, morphology, circulating signatures and therapy mapping.
- Confirm the current PDC ID for CPTAC colon before download.
- PXD052391 is ideal for progression; PXD079707 for primary tumors; IPX0007391001 for metastasis.

---

## 2. Bulk transcriptomic cohorts

### 2.1 Classic cohorts (verified counts)

| Accession | n | Platform | Content / role |
|---|---|---|---|
| GSE39582 | 585 | Affymetrix | Colon cancer subtypes/survival; core validation |
| GSE41258 | 390 | Affymetrix | Primary/metastatic CRC; validation |
| GSE68468 | 390 | Affymetrix | Colon cancer molecular dissection; validation |
| GSE87211 | 363 | Affymetrix | Rectal cancer +/- controls; rectal validation |
| GSE14333 | 290 | Affymetrix | 290 primary CRCs; survival |
| GSE44076 | 246 | Affymetrix | Healthy/adjacent/tumor colon; normal baseline |
| GSE17538 | 244 | Affymetrix | Moffitt metastasis profile; survival |
| GSE71187 | 189 | Affymetrix | Embryonic-derived CRC signatures |
| GSE17536 | 177 | Affymetrix | Moffitt colon cancer; survival |
| GSE103479 | 156 | Affymetrix | Stage II/III treated vs untreated; adjuvant chemo |
| GSE21510 | 148 | Affymetrix | CRC tissues; validation |
| GSE37892 | 130 | Affymetrix | Stage II colon cancer; validation |
| GSE72970 | 124 | Affymetrix | Metastatic CRC chemo response |
| GSE38832 | 122 | Affymetrix | Colon cancer; survival |
| GSE44861 | 111 | Affymetrix | Colon cancer tissues; validation |
| GSE33113 | 96 | Affymetrix | Stage II colon cancer |
| GSE37364 | 94 | Affymetrix | Adenoma-carcinoma biopsies |
| GSE26906 | 90 | Affymetrix | Stage II colon cancer |
| GSE41657 | 88 | Affymetrix | Normal/adenoma/adenocarcinoma |
| GSE9348 | 82 | Affymetrix | Early-stage CRC vs healthy |
| GSE39084 | 70 | Affymetrix | Sporadic early-onset CRC |
| GSE29621 | 65 | Affymetrix | Colon cancer mRNA/miRNA |
| GSE8671 | 64 | Affymetrix | Colorectal adenomas |
| GSE24514 | 49 | Affymetrix | MSI CRC vs normal |
| GSE31595 | 37 | Affymetrix | Stage II/III 128-gene signature |

### 2.2 New 2023-2026 cohorts (verified counts)

| Accession | n | Content / role |
|---|---|---|
| GSE241101 | 286 | Multiregional primary CRC; intratumor heterogeneity |
| GSE248381 | 258 | Pathway-level slow-cycling phenotype; subtype validation |
| GSE161158 | 250 | SMAD4-modulated profile; stage II/III DFS |
| GSE156915 | 361 | DNA damage immune response (DDIR); oxaliplatin |
| GSE159216 | 283 | CRC liver metastases; metastasis validation |
| GSE233616 | 188 | TF repertoire; epigenetic/transcriptional subtypes (array + methylation) |
| GSE220150 | 172 | Thailand CRC transcriptomics; ancestry diversity |
| GSE294511 | 169 | Tumor organoids of CRLM; functional validation |
| GSE106584 | 156 | Chemotherapy response stage II/III |
| GSE277272 | 150 | Small extracellular vesicle RNAs |
| GSE236055 | 143 | Colitis-associated neoplastic progression |
| GSE336139 | 128 | CRLM pure desmoplastic vs non-pure |
| GSE164191 | 121 | Blood-based expression CRC vs healthy |
| GSE140973 | 114 | Secondary resistance to anti-EGFR; transcriptional reprogramming |
| GSE178120 | 105 | CRC tissue expression |
| GSE271719 | 98 | CRC and paired normal (Affymetrix X3P) |
| GSE197802 | 88 | Genomic/transcriptional coherent signatures |
| GSE333907 | 83 | Dana-Farber early-onset CRC cohort |
| GSE175433 | 79 | Early-onset vs late-onset immune profiles |
| GSE206613 | 77 | Immune-metabolic rewiring patterns |
| GSE170999 | 76 | Locally advanced rectal cancer by KRAS |
| GSE237249 | 72 | CMS transition during progression |
| GSE313035 | 67 | METIMMOX (oxaliplatin + immunotherapy) |
| GSE233602 | 61 | Colon polyps |
| GSE292601 | 62 | Interferon-stimulated genes and chemotherapy sensitivity |
| GSE339362 | 56 | Total neoadjuvant therapy response |
| GSE315715 | 56 | Circulating mRNA biomarkers |
| GSE290756 | 50 | GATA6 loss and metastasis |
| GSE333300 | 48 | CD4-licensed CD8 T cells in MHC-I-negative MMRd |
| GSE349730 | 48 | LTBP2/collagen I niches in patient-derived models |
| GSE294300 | 36 | scRNA primary CRC and adjacent normal |
| GSE325868 | 30 | Stage-specific immune/molecular signatures |
| GSE330748 | 28 | CD4-licensed CD8 T cells in MHC-I-negative MMRd |
| GSE320038 | 14 | Spatial ecotype study (bulk RNA component) |
| GSE320040 | 7 | Spatial ecotype study (scRNA component) |

**Notes.**
- Avoid using overlapping cohorts (e.g., GSE17536/GSE17538, GSE14333/GSE17536) as independent validation without checking patient overlap.
- GSE87211 is rectal-enriched; use for rectal-specific validation.
- GSE220150 adds ancestry diversity; do not overgeneralize from it.

---

## 3. Single-cell and single-nucleus atlases

| Accession / resource | Content | n | Role |
|---|---|---|---|
| GSE132465 | scRNA-seq, 23 Korean patients (63,689 cells) | 23 patients | Core reference |
| GSE144735 | scRNA-seq, 6 Belgian patients (27,414 cells), core/border/normal | 6 patients | Regional heterogeneity |
| GSE178341 | scRNA-seq MMRd/MMRp atlas | 181 samples | Immune phenotype reference |
| GSE166555 | scRNA-seq normal/adjacent tumor | 25 samples | Normal baseline |
| GSE188711 | scRNA-seq left vs right CRC | 6 samples | Sidedness |
| GSE201348 | scRNA + scATAC, polyp-adenoma-CRC continuum | 72 samples | Progression/chromatin |
| GSE231559 | scRNA-seq, MD Anderson CRC Moon Shot | 26 samples | Independent reference |
| GSE312983 | snRNA-seq across age continuum | 189 tumors | EOCRC |
| GSE285117 | Neutrophil-focused series from the 4.27M-cell study | series-level | Myeloid states |
| GSE294300 | scRNA primary/adjacent | 36 samples | Small independent reference |
| GSE320040 | scRNA component of spatial ecotype study | 7 samples | Ecotype reference |
| HCA Gut Cell Atlas | Normal colon/rectum cell atlas | public | Normal baseline |
| HTAN Data Portal | Human Tumor Atlas Network CRC atlases (check availability) | cohort-dependent | Independent spatial/single-cell validation |
| CROST / SpatialDB | Curated spatial transcriptomics databases | public | Additional spatial sections |

**Notes.**
- Check raw-data availability; some series provide processed matrices only.
- Harmonize cell-type labels before integration.
- Do not treat the 4.27M-cell atlas and GSE285117 as independent.

---

## 4. Spatial transcriptomics and spatial multi-omics

| Accession / resource | Content | n | Platform/role |
|---|---|---|---|
| GSE225857 | scRNA + spatial, primary CRC and liver metastases | 8 samples | Core metastasis niches |
| GSE283052 | Spatial transcriptomics of CRC tissues | 12 samples | Core primary tumors |
| GSE320042 | Spatial ecotypes, high-resolution | 38 samples | Direct ecotype comparison |
| GSE284989 | Spatial genomics, cholesterol/ICI resistance | 22 samples | Therapy resistance |
| GSE280311 | High-definition spatial immune profiling | 12 samples | Immune niches |
| GSE280312 | High-definition spatial immune profiling | 16 samples | Immune niches |
| GSE270745 | Spatial atlas, DMBT1/dysplasia | 6 samples | Early lesions |
| GSE217414 | Liver metastases from CRC | 4 samples | Metastasis |
| GSE320038 | Spatial ecotype study (bulk RNA) | 14 samples | Companion |
| GSE320040 | Spatial ecotype study (scRNA) | 7 samples | Companion |
| HTAN / 10x public / CROST / SpatialDB | Additional CRC spatial sections | cohort-dependent | Replication |

**Platform benchmark:** four subcellular platforms differ in sensitivity, diffusion, segmentation and clustering (PMID 41107232). Record platform, spot/bin size and segmentation method; replicate niches in ≥2 datasets.

---

## 5. Methylation and epigenetics

| Accession / resource | Content | n | Role |
|---|---|---|---|
| TCGA-COAD/READ 450K | Methylation array | 2,202 files | Core epigenome |
| GSE233616 | CRC array + methylation subtypes | 188 | Subtype/epigenome validation |
| GSE220160 | CRC and polyp methylation | 16 | Progression |
| GSE231984 | Methylation, metastatic cancer origin | 791 | Tissue-of-origin deconvolution (optional) |
| GSE231969 | Methylation, metastatic cancer origin | 215 | Optional |
| GSE230193 | Methylation, metastatic cancer origin | 498 | Optional |
| GSE186575 | Plasma MeDIP-seq, GI cancer | 98 | Liquid biopsy |
| GSE124775 | cfDNA methylation (MMBC-seq) | 198 | Liquid biopsy |
| GSE97923 | cfDNA methylation, common cancers (COAD) | 773 | Liquid biopsy |
| GSE263030 | cfDNA 5hmC, occult CRC | 602 | Liquid biopsy |
| GSE86078 | Methylated circulating DNA monitoring | 149 | Liquid biopsy |
| GSE81314 | cfDNA 5hmC | 126 | Liquid biopsy |
| GSE79277 | Methylation haplotype blocks | 134 | Deconvolution |
| GSE208596 | cfDNA methylome monitoring | 95 | Liquid biopsy |
| GSE203071 | cfDNA 5hmC response prediction | 155 | Therapy monitoring |

**Notes.**
- TCGA methylation is the core epigenomic view.
- cfDNA methylation datasets are for the circulating-ecotype extension, not for tissue ecotype discovery.
- GSE231984/GSE231969/GSE230193 are metastatic-origin classifiers; use only for tissue-of-origin deconvolution.

---

## 6. Liquid biopsy and circulating signals

| Accession | n | Content | Role |
|---|---|---|---|
| GSE97923 | 773 | cfDNA methylation, COAD | Circulating methylation |
| GSE263030 | 602 | cfDNA 5hmC, occult CRC | Early detection |
| GSE188627 | 290 | EV-miRNA, metastatic CRC | EV signature |
| GSE174302 | 230 | Plasma cfRNA microbiome | Microbial cfRNA |
| GSE216561 | 113 | cfRNA human/microbial | Circulating multi-omic |
| GSE259422 | 156 | Plasma EV transcriptomics | EV signature |
| GSE290667 | 185 | cfDNA 5hmC, peritoneal metastasis | Metastasis |
| GSE203071 | 155 | cfDNA 5hmC response | Therapy monitoring |
| GSE86078 | 149 | cfDNA methylation monitoring | MRD |
| GSE188433 | 123 | Plasma EV long RNA, CRC/adenoma | Early detection |
| GSE164191 | 121 | Blood expression CRC vs healthy | Blood signature |
| GSE315715 | 56 | Circulating mRNA | Early detection |
| GSE186587 | 127 | Plasma long RNA, GI cancer | Multi-omic plasma |
| GSE186575 | 98 | Plasma MeDIP-seq | Methylation |
| GSE186573 | 86 | Plasma DNA-seq | Multi-omic plasma |
| GSE186595 | 49 | Plasma miRNA-seq | Multi-omic plasma |
| GSE277272 | 150 | Small EV RNAs | EV signature |
| GSE280385 | 53 | Mesenchymal CRC EV cargo | EV signature |
| GSE248663 | 40 | Circulating miRNA, co-clinical trial | Response |
| GSE230027 | 120 | Blood 5hmC, anticoagulant effects | Preanalytical |

**Use:** derive a circulating ecotype signature; test association with tissue ecotypes and outcome. Pre-specify feature sets; do not pool analytes without normalization.

---

## 7. Microbiome and metabolome

| Resource | Content | Role |
|---|---|---|
| curatedMetagenomicData | Uniform taxonomic/functional profiles from published cohorts | Multi-cohort meta-analysis |
| ENA/SRA | Raw metagenomes from landmark CRC studies | Re-profiling if needed |
| GSE174302 | Plasma cfRNA microbiome (230) | Circulating microbial signal |
| GSE216561 | Plasma cfRNA human/microbial (113) | Circulating multi-omic |
| GSE264208 | Microbiome-derived cfRNA (167) | Circulating microbial signal |
| GSE186587 / GSE186575 / GSE186573 / GSE186595 | Plasma multi-omics GI cancer | Integrative plasma |
| Metabolomics Workbench / MetaboLights | Sparse human CRC tissue/serum metabolomics | Targeted validation |
| METASPACE / SpaceM | Spatial metabolomics | Optional |

**Notes:** microbiome work is optional and confounded; keep as a companion or exploratory module. No causal claims without MR/colocalization or experimental validation.

---

## 8. Imaging and pathology

| Resource | Content | Role |
|---|---|---|
| TCGA diagnostic slides (GDC) | 625 diagnostic slides, COAD/READ | H&E morphology |
| PXD079707 WSIs | 1,000 treatment-naive CRC slides | Proteogenomic-matched morphology |
| TCIA / IDC | Radiology/pathology imaging collections | Optional imaging |
| Human Protein Atlas | IHC for candidate proteins | Orthogonal protein validation |
| GTEx colon | Normal colon expression | Normal baseline |

---

## 9. Perturbation, drug and dependency resources

| Resource | Content | Role |
|---|---|---|
| DepMap | CRISPR dependencies, cell-line omics | Ecotype-marker dependencies |
| GDSC | Drug sensitivity | Drug repurposing |
| CTRP | Drug sensitivity | Drug repurposing |
| PRISM | Drug sensitivity | Drug repurposing |
| LINCS / CMap | Connectivity signatures | Signature reversal |
| DREIMT | Immune drug repurposing | Immunotherapy hypotheses |
| Open Targets / DrugBank | Target evidence | Annotation |
| scPerturb / GEARS / CPA | Perturbation prediction | Optional |

---

## 10. Priority tiers

### Tier 1 (core; download first)
- TCGA-COAD/READ (RNA, methylation, CNV, mutation, RPPA, clinical, slides).
- **PXD079707** (1,000-CRC proteogenome).
- CPTAC colon (PDC000116).
- GSE39582, GSE17536, GSE14333, GSE38832, GSE87211, GSE103479.
- GSE132465, GSE144735, GSE178341, GSE201348, GSE312983.
- GSE225857, GSE283052, GSE320042, GSE284989.

### Tier 2 (validation and extension)
- PXD052391 (progression), IPX0007391001 (CRLM).
- GSE241101, GSE159216, GSE248381, GSE161158, GSE233616, GSE106584, GSE140973, GSE333907, GSE313035, GSE339362, GSE175433, GSE237249, GSE233602.
- GSE280311/GSE280312, GSE270745, GSE217414.
- GSE263030, GSE97923, GSE188627, GSE188433, GSE315715, GSE216561.

### Tier 3 (optional/companion)
- curatedMetagenomicData subset.
- Metabolomics Workbench/MetaboLights.
- TCIA/IDC, HPA.
- DepMap/GDSC/CTRP/PRISM/LINCS (no download of raw data; use APIs/tables).

---

## 11. Overlap and confounding checks

- GSE17536 and GSE17538 share the Moffitt cohort; do not use both as independent validation.
- GSE14333 and GSE17536 may overlap; verify patient IDs.
- GSE39582 and GSE14333 are independent (different institutions), but check.
- GSE320038/GSE320040/GSE320042 are components of the same ecotype study; do not treat as independent.
- GSE312983 and GSE285117 are related to the same groups/atlases; check.
- cfDNA/EV cohorts differ by analyte (methylation, 5hmC, miRNA, long RNA); do not pool without normalization.
- PXD079707, PXD052391, CPTAC and IPX0007391001 are independent cohorts, but check for shared patients/institutions.

---

## 12. Download caveats

- Verify `n` on the accession page; GEO counts can include technical replicates or subseries.
- Prefer processed matrices; raw DIA proteomics requires Spectronaut/DIA-NN and significant compute.
- PXD079707 is CC0 but check the file list for processed reports vs raw only.
- PXD052391 WES is restricted (HRA007452); proteomics is open.
- IPX0007391001 is open via iProx/GSA; check access terms.
- Record checksums, versions and licenses in `metadata/manifest.tsv`.

---

## 13. What this changes in the strategy

1. **Proteomic discovery is no longer limited to 110 CPTAC tumors**; PXD079707 provides 1,000 cases.
2. **Progression can be studied with 435 samples across stages** (PXD052391), not only bulk adenoma cohorts.
3. **Metastasis can be studied with matched primary/normal/liver samples** (IPX0007391001, 102 samples).
4. **Morphology can be linked to proteome/ecotypes** using 1,000 WSIs (PXD079707) plus 625 TCGA diagnostic slides.
5. **Circulating ecotypes are feasible** with large cfDNA/EV cohorts.
6. **Therapy-response mapping is feasible** across anti-EGFR, chemo, TNT and ICI contexts.
7. **The main novelty shifts** from "proteomic subtypes" to "cross-platform ecotype framework with single-cell/spatial/morphology/circulating validation".
