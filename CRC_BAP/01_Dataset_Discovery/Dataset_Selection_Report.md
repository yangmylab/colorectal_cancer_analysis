# Dataset selection report

## 1. Question-driven selection

The study asks whether CRC ecotypes are reproducible, spatially organised and clinically relevant. Selection is therefore organised by **role**, not by convenience:

- **Discovery:** TCGA-COAD/READ provides multi-omic breadth (RNA, methylation, CNV, mutation, RPPA, clinical) and adequate events.
- **Proteomic anchor:** PXD079707 provides 1,000 treatment-naive CRCs with DIA proteomics and matched H&E; CPTAC (110) provides an independent, well-annotated proteome/phosphoproteome.
- **Progression and metastasis:** PXD052391 (148 patients/435 samples) and IPX0007391001 (34 patients/102 samples) supply stage-resolved and matched primary/liver-metastasis protein data.
- **Reference:** GSE132465, GSE144735, GSE178341, GSE201348, GSE231559, GSE312983 and GSE285117 provide cell states, chromatin and age/EOCRC context.
- **Spatial validation:** GSE225857, GSE283052, GSE320042, GSE284989, GSE280311/GSE280312, GSE270745 and GSE217414 span multiple platforms and tissue contexts.
- **External validation:** GSE39582, GSE17536, GSE14333, GSE38832, GSE87211 and GSE103479 provide independent bulk cohorts with survival/treatment data.
- **Circulating extension:** GSE263030, GSE97923, GSE188627, GSE188433, GSE315715 and GSE216561 provide cfDNA/EV/microbial signals.
- **Orthogonal resources:** DepMap, GDSC, CTRP, PRISM and LINCS provide drug/dependency context; TCGA and PXD079707 WSIs provide morphology.

## 2. Compatibility assessment

| Criterion | Assessment |
|---|---|
| Organism | All human; no cross-species merging in the core |
| Tissue | Primary CRC for discovery; normal/adenoma/metastasis analysed as separate strata |
| Platform | Bulk RNA-seq, microarray, methylation array, MS proteomics, scRNA/snRNA, spatial; never pooled without platform-aware normalisation |
| Gene identifier | Harmonise to GENCODE/Ensembl; re-annotate microarray probes |
| Outcome | CSS/OS/RFS available in TCGA and most GEO cohorts; missing outcomes are excluded from survival analyses, not imputed |
| Treatment | Partial; retained for exploratory therapy analyses only |
| Batch | Cohort/platform are registered; ComBat/MMUPHin/limma used with sensitivity analyses |
| Overlap | GSE17536/GSE17538 and GSE14333/GSE17536 may share patients; check IDs before declaring independence |

## 3. Datasets rejected or deferred

| Dataset / class | Reason |
|---|---|
| Cell lines, organoids, PDX | Not patient tissue; used only for drug/functional context, not ecotype discovery |
| Normal-only atlases | Used as reference, not as disease cohorts |
| Metastatic-only cohorts in discovery | Reserved for the metastasis extension to avoid mixing primary and metastatic biology |
| Controlled-access data (dbGaP/EGA/Broad Single Cell Portal) | Not required for the core; user must obtain DUA if needed |
| Metabolomics | Public human CRC tissue metabolomics is sparse; retained as targeted validation only |
| Microbiome | Confounding is severe; retained as an optional companion, not a core claim |

## 4. Provenance and download rules

- Every dataset is recorded in `Dataset_Inventory.csv` with accession, repository, technology, n, groups, access conditions and limitations.
- Counts are from GEO/GDC/PRIDE metadata as of 2026-11-17; re-verify at download.
- **No dataset has been downloaded by the agent.** Per project standing instruction, the user downloads datasets manually into `03_Raw_Data/` and records checksums in `Logs/`.
- PXD079707 is a new PRIDE submission with no indexed paper; treat its metadata as provisional.
- PXD052391 WES is restricted (HRA007452); only proteomics/phosphoproteomics are used in the core.
- Do not merge cohorts until G1 passes and overlap checks are complete.

## 5. Minimum viable set (if time is limited)

1. TCGA-COAD/READ (RNA + clinical + methylation).
2. GSE39582 (external survival validation).
3. GSE132465 (single-cell reference).
4. One spatial cohort (GSE225857).
5. A PXD079707 subset (protein anchor).

This supports a defensible core paper; the full set supports the extended Track A manuscript.
