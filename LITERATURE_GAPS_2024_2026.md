# Literature and gap matrix (2023-2026)

**Prepared:** 2026-11-17
**Sources:** PubMed / Europe PMC, GEO, GDC, PDC, PRIDE, iProx, journal full texts
**Purpose:** define what is already established, what remains unresolved, and how the proposed Track A study fills each gap.

---

## 1. What has changed in 2023-2026

### 1.1 Atlases are no longer the bottleneck
- A **4.27-million-cell CRC immune atlas** (1,670 samples, plus 266 new patients, spatial and protein profiling) defined immune phenotypes and myeloid/neutrophil programs (PMID 41529665, *Cancer Cell* 2026).
- A **189-tumor single-nucleus RNA-seq atlas across the age continuum** defined early-onset CRC (EOCRC) cell states (GSE312983).
- An **isoform-resolution CRC atlas** integrated short- and long-read single-cell data and identified tumor-specific isoforms and neoepitopes (PMID 39216476, *Cell Genomics* 2024).
- A **polyp-adenoma-carcinoma continuum** with scRNA + scATAC (GSE201348) mapped cell-state and chromatin changes.
- MMRd/MMRp (GSE178341), liver metastasis (GSE225857), and multiple spatial atlases now cover most tissue contexts.

**Implication:** a descriptive "another atlas" paper is no longer novel. The open problems are integration, reproducibility and clinical utility.

### 1.2 Proteogenomics has scaled dramatically
- **CPTAC colon cancer** (110 prospectively collected tumors) established proteomic/phosphoproteomic layers and linked Rb phosphorylation and glycolysis to immune evasion (PMID 31031003, *Cell* 2019).
- A **148-patient, 435-sample progression cohort** profiled nontumor, intraepithelial neoplasia, infiltration and advanced CRC, with proteome and phosphoproteome (PMID 38861363, *Cancer Res* 2024; data PXD052391 / IPX0002188000).
- A **34-patient, 102-sample CRLM proteogenomic cohort** defined C1 (metabolism, poor prognosis) and C2 (RNA function, better prognosis) liver-metastasis proteomic subtypes (PMID 41195591, *Adv Sci* 2026; data IPX0007391001 / PRJCA020890).
- A **1,000 treatment-naive CRC proteogenomic resource** (DIA proteomics, targeted genomics, whole-slide histology, clinical follow-up, AI-based histology) has been deposited in PRIDE as **PXD079707** (submitted 2026; no indexed paper found as of 2026-11-17; CC0). Its own abstract states that current molecular classifications incompletely capture CRC heterogeneity and that a proteomic annotation layer refines existing frameworks.

**Implication:** the proteomic layer is now large enough to anchor a taxonomy, but it is not yet integrated with single-cell, spatial, methylation, liquid-biopsy and therapy-response data.

### 1.3 Spatial biology is maturing but fragmented
- CRC spatial studies now cover primary tumors, liver metastases, TLS, immune exclusion, cholesterol/ICI resistance and neutrophil niches (GSE225857, GSE283052, GSE320042, GSE284989, GSE280311/280312, GSE270745, GSE217414).
- A systematic benchmark of four subcellular spatial platforms (Stereo-seq v1.3, Visium HD FFPE, CosMx 6K, Xenium 5K) on colon adenocarcinoma, HCC and ovarian cancer showed platform-specific differences in sensitivity, diffusion, segmentation and clustering (PMID 41107232, *Nat Commun* 2025).
- Reviews highlight the lack of computational tools for CRC spatial heterogeneity and the need for integration into precision medicine (PMID 39711954; PMID 39075612; PMID 40107672).
- "Spatial ecotype" is now an explicit concept in immune exclusion (PMID 42736558, *Mol Cancer* 2026), but most studies remain single-platform and descriptive.

**Implication:** spatial data can validate ecotypes only if platform effects are modeled and findings replicate across ≥2 platforms.

### 1.4 Molecular classification remains unresolved
- CMS1-4, CRIS A-D and PDS1-3 coexist; a 2024 update documents their rationale and limitations (PMID 38319359).
- A 2026 editorial explicitly asks "why one classification system is no longer enough" (PMID 41607760).
- The 1,000-CRC proteogenomic resource reports that current classifications incompletely capture heterogeneity.
- Multi-omics reviews call for integration and clinical translation (PMID 42192865; PMID 42795017).

**Implication:** a reproducibility benchmark plus a locked, cross-platform classifier is a publishable contribution even without new biology.

### 1.5 Therapy resistance and liquid biopsy are data-rich but disconnected
- Anti-EGFR secondary resistance has a dedicated transcriptional cohort (GSE140973, 114 samples).
- Chemotherapy response cohorts exist (GSE103479, 156; GSE106584, 156; GSE156915, 361 DDIR/oxaliplatin).
- Immunotherapy settings include METIMMOX (GSE313035, 67), total neoadjuvant therapy (GSE339362, 56), and MMRd ICI models (GSE333300, 48; GSE330748, 28).
- ctDNA MRD is clinically maturing (PMID 39675560, *Ann Oncol* 2025), and large cfDNA methylation/5hmC cohorts exist (GSE263030, 602; GSE97923, 773; GSE203071, 155; GSE86078, 149).
- Plasma EV/RNA cohorts exist (GSE188627, 290; GSE188433, 123; GSE259422, 156; GSE315715, 56).

**Implication:** tissue ecotypes have not been systematically linked to circulating signatures or therapy response.

### 1.6 Microbiome and imaging are promising but under-integrated
- CRC microbiome reviews call for spatial multi-omics and AI-driven analytics (PMID 41356485).
- Microbial peptides can be detected directly in CPTAC CRC proteomes (PMID 40263502).
- Plasma microbial cfRNA datasets exist (GSE174302, 230; GSE216561, 113; GSE264208, 167).
- Pathology foundation models (CHIEF, PMID 39232164; MUSK, PMID 39779851) can predict molecular profiles and prognosis from H&E, and the 1,000-CRC resource explicitly enables AI-based histology.
- Reviews describe imaging-to-omics prediction as an emerging field (PMID 41802282).

**Implication:** morphology–molecular and microbiome–host integration are credible extensions, but only with rigorous confounding control.

---

## 2. Gap matrix

| Gap | Evidence | Why it matters | How Track A fills it | Risk / mitigation |
|---|---|---|---|---|
| **G1. Subtype instability and no locked classifier** | CMS/CRIS/PDS coexist and are platform/cohort-dependent (PMID 38319359; PMID 41607760; PMID 42192865); 1,000-CRC resource says current classifications are incomplete | Reviewers and clinicians cannot compare or deploy subtypes | Benchmark CMS/CRIS/PDS across ~20 cohorts and 3-4 platforms; report ARI/NMI, reclassification and survival of discordant cases; release a locked classifier | Platform/batch confounding; use ComBat/MMUPHin, negative controls, sensitivity analyses |
| **G2. Proteome-single-cell-spatial bridge missing** | 1,000-CRC proteome (PXD079707), CPTAC (110), progression (148 patients), CRLM (34 patients), >4M single cells and multiple spatial cohorts exist separately | The field lacks a unified cell-state-to-proteome-to-space mapping | Use single-cell programs to define ecotypes, anchor them in proteome/phosphoproteome, validate in spatial data | Sample overlap and platform mismatch; use program-level (not gene-level) integration and external validation |
| **G3. Spatial validation is single-platform and small** | Platform benchmark shows differences (PMID 41107232); CRC spatial reviews note limited tools (PMID 39711954; PMID 39075612; PMID 40107672) | Niche claims may not replicate | Multi-platform spatial meta-atlas; platform-aware normalization; require replication in ≥2 datasets; report effect sizes/heterogeneity | Small n per dataset; use random-effects meta-analysis and pre-specified niches |
| **G4. Clinical utility is asserted, not demonstrated** | Most subtype papers report associations; MRD/ctDNA reviews emphasize clinical validation (PMID 39675560); multi-omics reviews call for translation (PMID 42192865; PMID 42795017) | Without calibration/DCA, classifiers are not actionable | Primary endpoint CSS; base vs base+ecotype; C-index, time-AUC, calibration, DCA, NRI/IDI; external validation in ≥4 cohorts | Weak incremental value; focus on stage II/III and MSS; report null results |
| **G5. EOCRC is under-integrated** | 189-tumor snRNA atlas (GSE312983); EOCRC bulk cohorts (GSE333907, 83; GSE175433, 79); EOCRC is rising | EOCRC may have distinct biology and therapy implications | Test ecotype distribution by age; define age-associated programs; validate in proteome/spatial data | Confounding by stage/sidedness; stratified and adjusted analyses |
| **G6. Liver metastasis tropism unresolved** | CRLM proteogenomic C1/C2 subtypes (PMID 41195591); single-cell/spatial CRLM (GSE225857; GSE217414; GSE159216; PMID 40107672) | Metastasis is the main cause of death | Map primary ecotypes to matched/semi-matched liver metastases; test tropism and niche remodeling; integrate C1/C2 proteomic subtypes | Few matched primary-metastasis public samples; use multiple cohorts and clearly label matched vs unpaired |
| **G7. Therapy-response ecotypes missing** | Anti-EGFR resistance (GSE140973), METIMMOX (GSE313035), TNT (GSE339362), chemo (GSE103479; GSE106584; GSE156915), MMRd ICI (GSE333300; GSE330748); resistance review (PMID 42901333) | Therapy selection is the main clinical need | Build ecotype-therapy association map; test anti-EGFR, chemo and ICI contexts; label exploratory | Small/public response cohorts; avoid overclaiming; use pre-specified analyses |
| **G8. Microbiome-host integration is correlative** | Microbiome review (PMID 41356485); microbial peptides in CPTAC (PMID 40263502); plasma microbial cfRNA (GSE174302; GSE216561; GSE264208) | Mechanistic and translational claims are weak | Optional companion: batch-aware multi-cohort metagenomics + microbial peptides + confounder control; no causal claims | Confounding (diet, BMI, antibiotics, geography); keep separate or clearly exploratory |
| **G9. Liquid biopsy not linked to tissue ecotypes** | ctDNA MRD review (PMID 39675560); cfDNA methylation/5hmC (GSE263030; GSE97923; GSE203071; GSE86078); EV/RNA (GSE188627; GSE188433; GSE259422; GSE315715) | Non-invasive monitoring is clinically valuable | Derive a circulating ecotype signature; deconvolve cfDNA/EV signals; test association with tissue ecotypes and outcome | Platform/analyte heterogeneity; use pre-specified feature sets and external validation |
| **G10. Morphology-molecular link underexplored** | CHIEF (PMID 39232164), MUSK (PMID 39779851), imaging-to-omics review (PMID 41802282), 1,000-CRC resource enables AI histology (PXD079707) | H&E is universally available; morphology could surrogate molecular ecotypes | Extract pathology foundation-model features from TCGA/PXD079707 WSIs; predict ecotypes; test prognostic value | Slide quality/domain shift; use multiple-instance learning and external slide sets |
| **G11. No CRC-specific methods benchmark** | Spatial benchmark (PMID 41107232); multi-omics reviews (PMID 42192865; PMID 42736558); many single-method papers | Reviewers cannot judge which integration/deconvolution method is appropriate | Benchmark MOFA+/iCluster/SNF/DIABLO, deconvolution methods and spatial niche tools on CRC; report runtime, stability, calibration | Compute-heavy; restrict to pre-specified methods and datasets |
| **G12. Progression/early detection gap** | Progression proteogenomics (PMID 38861363; PXD052391); polyp scRNA/scATAC (GSE201348); adenoma cohorts (GSE8671, 64; GSE37364, 94; GSE41657, 88; GSE233602, 61) | Interception and early detection need stage-resolved targets | Build an adenoma-carcinoma ecotype trajectory; identify stage-specific vulnerabilities; validate in proteome/spatial data | Few matched multi-omic progression samples; use the 148-patient cohort and bulk adenoma cohorts |
| **G13. Data equity/ancestry** | Thailand cohort (GSE220150, 172); most atlases are East Asian/European/North American | Generalizability claims are limited | Multi-ancestry validation; report cohort composition; avoid pan-cancer/general-population claims | Limited ancestry metadata; state limitations explicitly |
| **G14. Causal validation absent** | Most ecotype papers are correlative; MR/colocalization resources are public (GTEx, UKB-PPP, deCODE) | Reviewers ask for causal support | Optional: eQTL/pQTL colocalization/MR for ecotype driver genes; label exploratory | MR assumptions; weak instruments; do not overclaim |

---

## 3. Novelty map (what the paper can claim)

**N1. First cross-platform ecotype taxonomy anchored in the 1,000-CRC proteogenomic resource and single-cell/spatial atlases.**
- The 1,000-CRC resource provides proteomic depth; the single-cell atlases provide cell-state resolution; the spatial cohorts provide tissue context. No study integrates all three into one taxonomy.

**N2. First systematic reproducibility benchmark of CMS/CRIS/PDS across ~20 cohorts and 3-4 platforms, with a locked classifier.**
- This turns a known weakness into a primary result and provides a deployable alternative.

**N3. First multi-platform spatial meta-atlas of ecotype-associated niches.**
- Platform-aware normalization and replication across Visium/Stereo-seq/Xenium/CosMx-like data; niche effect sizes and heterogeneity reported.

**N4. First morphology-molecular ecotype model using pathology foundation models.**
- H&E-derived features predict ecotypes and add prognostic value; tested on TCGA and PXD079707 slides.

**N5. First circulating ecotype signature linking tissue ecotypes to cfDNA/EV RNA.**
- Tests whether ecotypes are detectable non-invasively and whether they track outcome/therapy.

**N6. First ecotype-resolved therapy-response map across anti-EGFR, chemotherapy and ICI contexts.**
- Uses public response cohorts; exploratory but clinically relevant.

**N7. First progression/metastasis ecotype axis using the 148-patient progression cohort and the 34-patient CRLM cohort.**
- Links adenoma-carcinoma progression, primary ecotypes and liver-metastasis niches.

**N8. Open resource: locked classifier, benchmark, portal, code and manifests.**
- Suitable for *Briefings in Bioinformatics* / *CSBJ* framing; clinical translation framing for *J Transl Med*.

**Positioning against PXD079707.** If the 1,000-CRC paper is published, do not compete on "proteomic subtypes". Cite it as the proteomic anchor and differentiate on: single-cell anchoring, cross-platform reproducibility, spatial niche validation, morphology, circulating signatures, therapy mapping and open classifier.

---

## 4. What would be incremental (avoid claiming as novel)

- Another single-cohort CMS/CRIS reclassification.
- Another descriptive single-cell atlas without cross-platform validation.
- Another correlation between a gene and survival in TCGA.
- Another drug-repurposing list without external evidence or mechanistic context.
- Another microbiome association study without confounder control.
- Another spatial study on one section/platform.

---

## 5. Key references (verified)

**Atlases and CRC biology**
- PMID 41529665 - 4.27M-cell CRC immune atlas, *Cancer Cell* 2026
- PMID 39216476 - isoform-resolution CRC atlas, *Cell Genomics* 2024
- PMID 42106149 - MSI/MSS immunopeptidomics, *Mol Cell Proteomics* 2026
- PMID 39137726 - CCL19+ fibroblasts/TLS in CRLM, *Cancer Cell* 2024
- PMID 38335276 - Th17/TWEAK EMT in CRLM, *Cancer Res* 2024
- PMID 40702186 - NNMT CAF immunosuppression, *Nature* 2025
- PMID 40826474 - COL10A1+ CAF, *J Exp Clin Cancer Res* 2025
- PMID 40050872 - CTHRC1+ CAF/WNT5A, *J Transl Med* 2025

**Proteogenomics**
- PMID 31031003 - CPTAC colon proteogenomics, *Cell* 2019
- PMID 38861363 - 148-patient progression proteogenomics, *Cancer Res* 2024 (PXD052391 / IPX0002188000)
- PMID 41195591 - 34-patient CRLM proteogenomics, *Adv Sci* 2026 (IPX0007391001 / PRJCA020890)
- PMID 40263502 - microbial peptides from CPTAC CRC proteomics, *Sci Rep* 2025
- PXD079707 - 1,000 treatment-naive CRC proteogenomic resource, PRIDE 2026 (no indexed paper as of 2026-11-17)

**Classification and reviews**
- PMID 38319359 - molecular pathological classification update, *Virchows Arch* 2024
- PMID 41607760 - "why one classification system is no longer enough", *World J Gastrointest Oncol* 2026
- PMID 42192865 - multi-omics integration in CRC, *Cancers* 2026
- PMID 42795017 - next-generation biomarkers and AI in CRC, *Cancers* 2026
- PMID 39711954 - spatial biology in CRC, *Front Oncol* 2024
- PMID 39075612 - CAF spatial transcriptomics review, *Biomark Res* 2024
- PMID 40107672 - liver metastasis spatial immune review, *J Immunother Cancer* 2025
- PMID 41722207 - single-cell/precision surgery review, *Surg Oncol* 2026
- PMID 39675560 - ctDNA MRD review, *Ann Oncol* 2025
- PMID 41356485 - microbiome review, *Front Microbiol* 2025
- PMID 42901333 - immunotherapy resistance in dMMR CRC, *Crit Rev Oncol Hematol* 2026

**Spatial and AI methods**
- PMID 41107232 - benchmark of four subcellular spatial platforms, *Nat Commun* 2025
- PMID 42736558 - spatial ecotype in immune exclusion, *Mol Cancer* 2026
- PMID 39232164 - CHIEF pathology foundation model, *Nature* 2024
- PMID 39779851 - MUSK vision-language pathology model, *Nature* 2025
- PMID 41802282 - predicting omics from imaging, *Brief Bioinform* 2026
- PMID 41038164 - RAEFISH whole-genome spatial transcriptomics, *Cell* 2025
- PMID 42190664 - SPAC-seq/TARDIS spatial CRISPR screens, *Cell* 2026
- PMID 40488620 - ProCanFDL federated proteomics deep learning, *Cancer Discov* 2025

**New methods (2026)**
- PMID 42020427 - CellNiche contrastive learning for atlas-scale spatial omics
- PMID 42568193 - image-guided spatial omics enhancement
- PMID 42754563 - CoxFormer multimodal generative spatial inference
- PMID 42786276 - DePass paired single-cell/spatial multi-omics integration
- PMID 42986684 - spOT-NMF spatial deconvolution
- PMID 42412825 - Riemannian metric learning for spatial multi-omics alignment
- PMID 42721220 - sosta multicellular anatomical structures from spatial omics
