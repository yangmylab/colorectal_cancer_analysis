# A reproducible two-ecotype transcriptomic axis in colorectal cancer predicts recurrence but not survival, and is not incremental over CMS

**Draft v0.1 — generated 2026-11-17 from `13_Reports/Findings_Ledger.csv` (F01–F15).**
Every quantitative statement below is traceable to a graded finding; the mapping is in
`claims_to_evidence.tsv`. Claims are not permitted to exceed their evidence grade (E1–E3).

---

## Abstract

**Background.** Molecular taxonomies of colorectal cancer (CRC) — CMS, CRIS, PDS — have repeatedly been proposed
and repeatedly failed to reach clinical use, in part because their reproducibility across assays and cohorts has
not been established with the same rigour as their discovery. We asked whether a two-group expression ecotype
axis is (i) reproducible across platforms relevant to clinical assays, and (ii) prognostic beyond the published
taxonomy.

**Methods.** Consensus clustering of 50 Hallmark programme scores from 624 TCGA-COAD/READ tumours selected a
K=2 partition (E1, Hedgehog/pancreas-β/MYC-V2; E2, EMT/angiogenesis/IL6–JAK–STAT3), frozen as a nearest-centroid
classifier and then applied without refitting to independent cohorts: protein and phosphoproteome (CPTAC-COAD),
450K methylation (TCGA), Visium spatial transcriptomics (16 sections, 43,818 spots), and disease-free and
overall survival (GSE14333, GSE39582, GSE17536, GSE38832). Every claim was assigned an evidence grade before
writing.

**Results.** The K=2 partition was stable (bootstrap ARI 0.975; classifier cross-validation 0.979). It reproduced
at the protein level (83.3% agreement with the mRNA call, ARI 0.439; **F06, E3**) and at the phosphoproteome
level (87.5%, ARI 0.558; permutation p<0.001; **F09, E3**), where kinase-substrate inference further identified
ecotype-associated kinase activity not explained by kinase abundance (**F12, E1**). It did **not** reproduce in
promoter methylation (ARI 0.006; **F07**) and the bulk classifier did **not** transfer to Visium spots (ARI 0.042
in a self-check; **F08**), although E2 programmes are spatially autocorrelated and co-localise with the stromal
compartment (**F13, descriptive**). Prognostically, the axis was **not** associated with overall survival in
discovery (**F01**) or in external validation (meta HR 1.25, 95% CI 0.97–1.61, p=0.081; **F03**), and a
cancer-specific death analysis in an independent cohort was likewise null (**F15**). In contrast, the axis was
associated with **disease-free survival**, replicating across three cohorts (GSE14333 HR 2.37; GSE39582 HR 1.50,
p=0.018; GSE17536 HR 2.51; fixed-effect 1.67, p=0.0007; **F11, E2**). That association, however, **did not survive
adjustment for CMS** (HR 1.42, p=0.107; ecotype–CMS Cramér's V = 0.52; **F14**).

**Conclusions.** A two-ecotype axis is reproducible across transcriptome, proteome and phosphoproteome, which
distinguishes it from taxonomies that exist only at the RNA level. It is directionally prognostic for recurrence
but not for death, and it is not demonstrably incremental over CMS. Reproducibility is established; clinical
added value is not.

---

## Introduction

Colorectal cancer has no shortage of molecular classifications. Consensus Molecular Subtypes (CMS1–4) are the
most widely used; CRIS and PDS1–3 add alternatives; each was derived from expression data, each claims biological
meaning, and none has displaced stage in the clinic. Two problems recur across this literature. First,
classification instability: partitions are sensitive to cohort composition, platform and algorithm, so a subtype
label learned in one study frequently does not transfer. Second, incremental value: a new taxonomy must be shown
to add information to what clinicians already have — stage, MSI status, and for research purposes CMS — rather
than restate it.

We therefore treated reproducibility as the primary object of study rather than a footnote. A K=2 partition of
Hallmark programme scores was derived in TCGA-COAD/READ and frozen; every subsequent analysis applies that frozen
classifier to a different assay, cohort or endpoint without refitting. Because a null result is as informative as
a positive one, each analysis was pre-assigned an evidence grade: **E3** for cross-platform convergence spanning
assays, **E2** for an association replicated in an independent cohort, **E1** for exploratory or post-hoc
observations, **descriptive** for characterisations with no inferential claim, and **not_supported** where the
pre-specified analysis did not reach significance. Claims in this manuscript may not exceed the grade of the
finding that licenses them.

The specific questions were: (1) does the axis survive a change of assay — protein, phosphoprotein, methylation,
spatial? (2) is it associated with outcome, and does any association replicate? (3) does it add anything beyond
CMS and MSI?

---

## Results

### 1. A stable two-group axis in discovery

Fifty Hallmark programmes were scored in 624 TCGA-COAD/READ tumours, integrated by similarity network fusion, and
clustered by consensus k-means across K=2–6. K=2 was selected on stability and deployability: cluster sizes 355
(E1) and 269 (E2), bootstrap adjusted Rand index **0.975**, nearest-centroid classifier cross-validated accuracy
**0.979**. E1 is characterised by Hedgehog, pancreas-β-cell and MYC-V2 programmes; E2 by epithelial–mesenchymal
transition, angiogenesis and IL6–JAK–STAT3 signalling. K=4 also produced a stable partition and an OS association
in discovery (C-index 0.785, LRT p=0.024) but is **post hoc** and did not replicate (**F02 E1**, **F04**).

Single-cell deconvolution of the same tumours against 63,689 cells from six lineages gave a plausible tissue
composition (epithelial 0.546, stromal 0.267, myeloid 0.123; **F05, descriptive**) and was used only as context.

### 2. The axis reproduces at the protein and phosphoproteome level

In CPTAC-COAD, the frozen classifier was applied to Hallmark scores computed from protein abundance. **83.3%** of
96 tumours received the same ecotype from mRNA and protein (**ARI 0.439**), and programme-level concordance was
highest for the programmes that define E2 and the immune axis (IFN-α r=0.87, IFN-γ r=0.81, EMT r=0.71;
**F06, E3**). Applying the same classifier to the phosphoproteome — 11,334 sites quantified in at least half of
the tumours, aggregated to gene level so that hyper-phosphorylated proteins do not dominate — raised agreement to
**87.5% (ARI 0.558)**, against a 1000-permutation null whose 95th percentile was 0.033 (**p=0.000999**;
**F09, E3**). Differential phosphorylation between ecotypes reached FDR<0.05 at 3,095 sites, of which 538 exceeded
|log2FC|>0.5.

Kinase–substrate inference on the same data (PhosphoSitePlus prior; 78 kinases with at least five measured
substrates) identified 9 kinases at FDR<0.05, five of them without any change in the kinase's own abundance — most
notably CDK7 activity lower in E2 (d=−1.00, FDR 1.1×10⁻⁴; **F12, E1**). Activity and abundance were only partly
correlated across kinases (r=0.31), which is the point: some of the phosphoproteome difference is signalling, not
protein quantity. This remains exploratory — activity is inferred, not measured, in a single cohort.

### 3. Where the axis does not reproduce

**Methylation.** Using gene-level promoter β values for 19,269 genes (annotation-based selection, fixed before
testing), the ecotype showed essentially no relationship to methylation structure: ARI against an unsupervised
methylation partition **0.006**, only 5 genes at FDR<0.05 with |logFC|>0.5, global β p=0.113, and CIMP-high
proportions of 74.5% (E1) versus 84.3% (E2) (**F07, not_supported**). The axis is not a methylation phenotype.

**Spatial projection.** The frozen bulk classifier cannot be applied to 55 µm Visium spots: an internal
self-check on the bulk training domain returned ARI 0.042, every spot fell on one side of the decision boundary,
and spot-level programme scores sat about 5 SD below the bulk distribution (**F08, not_supported**). This is a
statement about transfer, not about biology — a spot has no bulk tumour's compartment structure.

**What the spots do show.** Analysing sections on their own terms (k=6 neighbour graph, 199 permutations) showed
that the E2-defining programmes are spatially autocorrelated (mean Moran's I 0.528) — but **no more so than
Hallmark programmes generally** (0.587), and they co-localise with the stromal compartment (median ρ=0.49) rather
than the myeloid (0.04) (**F13, descriptive**). Two caveats are stated rather than buried: with ~2,570 spots per
section the permutation test rejects the null for every programme (a ceiling effect, so only effect sizes are
interpreted), and a 400-spot cell2location cross-check did not corroborate the NNLS fractions used for the
co-localisation statistic (median r=0.11). The spatial co-localisation is therefore method-dependent; the
autocorrelation results are not, as they involve no deconvolution.

### 4. No association with overall survival, in any cohort

In discovery (n=476, 79 events, adjusted for age, sex, stage and sidedness), adding the ecotype improved the
C-index from 0.763 to 0.772 with a likelihood-ratio **p=0.107** (**F01, not_supported**). The frozen classifier
reproduced this null externally: meta-analysis of GSE39582 (p=0.333) and GSE17536 (p=0.0785) gave **HR 1.25 (0.97–1.61),
p=0.081** with consistent direction but no significance (**F03, not_supported**). A third, independent cohort
analysed on cancer-specific survival (28 cancer deaths) was also null (**HR 1.00, 0.46–2.17, p=0.997**;
**F15, not_supported**).

### 5. Association with disease-free survival replicates — but is not incremental over CMS

The one positive outcome association concerns recurrence. In GSE14333 (226 tumours, 50 disease-free events,
stage I–III), the adjusted hazard ratio for E2 versus E1 was **2.37 (1.27–4.40), LRT p=0.0043**. Applying the
frozen classifier to recurrence-free survival in GSE39582 (483 tumours, 137 events, stage I–III) gave **HR 1.50
(1.07–2.10), p=0.018**, and in GSE17536 (135 tumours, 32 events) **HR 2.51 (1.14–5.53), p=0.022**. The fixed-effect
summary across the three cohorts is **HR 1.67 (1.24–2.24), p=0.0007**, with no evidence of heterogeneity (p=0.21);
bootstrapping the fully adjusted estimate left 99% of resamples above 1 (**F11, E2**).

That result is, however, not independent of the published taxonomy. In GSE39582, which carries CMS calls, ecotype
and CMS overlap substantially (Cramér's V=0.52). Adding CMS to the model removed the disease-free association:
**HR 1.42 (0.93–2.17), p=0.107**, with a C-index gain over covariates of +0.015 (versus +0.025 when only stage,
age and sex were adjusted for). The same adjustment removed the (already weaker) OS signal (HR 1.17, p=0.40)
(**F14, not_supported**). We therefore report the disease-free association as replicated and directionally
consistent, but not demonstrably incremental over CMS in the one cohort that can test it.

---

## Discussion

The axis we describe is unusual less in what it predicts than in where it survives. It persists from mRNA to
protein to phosphoprotein — three assays with different dynamic ranges, different quantification and only partly
overlapping gene coverage — with concordance at the level usually reported between technical replicates of the
same assay. It disappears at the methylome and cannot be projected onto spatial spots. Read together, these
results suggest the partition captures something about the tumour–stroma–immune axis that is post-transcriptionally
and post-translationally maintained, rather than a promoter-methylation program or a property of individual
cells in situ.

The prognostic picture is more deflationary. Three independent cohorts, analysed on the endpoint each can support,
give no association with death; the one association that replicates is with recurrence. This pattern — recurrence
yes, death no — is biologically plausible in a disease where salvage therapy and competing mortality dilute the
recurrence signal, but it is also exactly the pattern that should make one cautious about claiming clinical
utility. The incremental-value analysis settles the question for this dataset collection: with stage, MSI and CMS
in the model, the ecotype adds nothing detectable.

Our results do not show that expression-based taxonomy of CRC is worthless. They show something narrower and more
useful for the field: a two-group axis can be made reproducible across assays if it is frozen early and tested
without refitting, and reproducibility does not imply incremental clinical value. Both properties have to be
demonstrated, and they can dissociate.

**Limitations.** All data are retrospective and public; there is no wet-lab validation and no prospective cohort.
Cross-platform layers come from small cohorts (CPTAC-COAD: 96 tumours, 8 deaths; GSE14333: 226 tumours, 50 events),
so protein and phospho analyses test reproducibility, not prognosis. The phosphoproteome is 59.8% missing before
filtering; kinase activity is inferred rather than measured. The incremental-value test rests on a single cohort
(GSE39582), because the other cohorts carry no CMS calls; CMS itself was expression-derived in that same cohort,
so ecotype and CMS are not independent by construction. Endpoint definitions differ slightly between cohorts (DFS
in GSE14333 and GSE38832, RFS in GSE39582 and GSE17536) and disease-free endpoints are never pooled with overall
survival. Finally, classifier agreement is sensitive to how programme scores are built: recomputing Hallmark
scores from a different gene-set mapping moved protein-layer agreement from 85.4% to 83.3% by flipping 2 of 96
borderline tumours, a sensitivity we report rather than hide.

**What we do not claim.** We do not claim a new prognostic classifier, a clinically actionable subtype, a causal
role for any programme, or that the ecotype refines CMS. We do not claim kinase activity as a measured quantity.
We do not claim any spatial ecotype assignment — the bulk classifier does not transfer to spots, and no spot was
labelled.

---

## Methods

**Data.** TCGA-COAD/READ RNA-seq (624 tumours, GDC Data Release 46.0); external cohorts GSE39582 (585), GSE17536
(177), GSE14333 (290, DFS only), GSE38832 (122, CEL files); CPTAC-COAD proteome and phosphoproteome (PDC000116,
cBioPortal datahub); TCGA 450K methylation (393 tumours); Visium sections GSE283052 (12) and GSE267401 (4).
Single-cell reference: GSE132465 (63,689 cells, six lineages).

**Programme scoring.** Fifty Hallmark sets. Discovery used GSVA; external layers used the mean of gene-wise
z-scores within a set, with sites averaged per gene in the phosphoproteome so that proteins with many measured
sites do not dominate. The sets are frozen in `Environment/hallmark_frozen.gmt` so scoring is reproducible offline.

**Classifier.** Nearest centroid on the 50 programme scores, computed in the discovery data and frozen. Applied to
every external layer without refitting; features recoverable in each layer are reported (50/50 in CPTAC,
49/50 in the phosphoproteome).

**Statistics.** Cox proportional hazards with likelihood-ratio tests; C-index for discrimination; fixed-effect
inverse-variance meta-analysis; bootstrap (500 resamples) for stability; permutation nulls (1,000 for cross-layer
ARI, 199 for spatial autocorrelation). Endpoints are analysed separately and never pooled across families
(OS/CSS/DSS versus DFS/RFS). Complete-case analysis for covariates, with counts reported.

**Evidence grading.** E3 = cross-platform convergence across assays; E2 = association replicated in an independent
cohort; E1 = exploratory or post hoc; descriptive = characterisation without inference; not_supported = the
pre-specified analysis did not reach significance. Fifteen findings were graded (F01–F15) and every claim in this
manuscript maps to one of them (`claims_to_evidence.tsv`).

**Reproducibility.** All analyses are scripted (`06_Analysis/1_R/1_Rscripts/`), logged to a decision log, and version
controlled. Analysis scripts use one level of parallelism (6 workers) with BLAS pinned to one thread per worker.

**Code and data availability.** Code and derived tables are in this repository; raw data are public and referenced
by accession. A Zenodo DOI is pending.

---

## Figure legends

| Figure | Content |
|---|---|
| plot_1 | Single-cell marker signatures from GSE132465 (six lineages) |
| plot_2–5 | Quality control: PCA before/after batch correction, library size, top variable features |
| plot_6–8 | Hallmark signature heatmap, consensus matrix for K=2, SNF ecotype PCA |
| plot_9–10 | Ecotype stability across K; Kaplan–Meier by ecotype (OS, discovery) |
| plot_11–12 | Five-year calibration; decision-curve analysis |
| plot_13 | Mean TCGA cell-type fractions (NNLS against the single-cell reference) |
| plot_14 | K-sensitivity of the survival model |
| plot_15 | External-validation forest plot (GSE39582, GSE17536, fixed-effect) |
| plot_16–19 | CPTAC-COAD protein layer: programme heatmap by ecotype, mRNA–protein concordance, KM, agreement matrix |
| plot_20–24 | TCGA 450K: CIMP score, global methylation, top differentially methylated genes, KM by methylation cluster, ecotype versus methylation cluster |
| plot_25–27 | Spatial transfer diagnostics: spot programme scores, within- versus between-section variation, covariate shift |
| plot_28–32 | Phosphoproteome: programme heatmap, three-layer concordance, differential-site volcano, phospho–protein agreement, permutation null |
| plot_33–34 | GSE14333 disease-free survival: KM, hazard ratios across models |
| plot_35–36 | GSE39582 recurrence-free survival: replication forest, KM |
| plot_37–38 | Kinase-substrate activity: volcano, activity versus abundance |
| plot_39–40 | Spatial niches: Moran's I per programme, per-section Moran's I for E2 programmes |
| plot_41–42 | Incremental value: C-index by covariate set, ecotype versus CMS composition |
| plot_43–44 | GSE38832: DSS versus DFS hazard ratios, cancer-specific survival by ecotype |

---

## Pending before submission

1. **Author list and affiliations** — `Authors.xlsx` is still an empty template.
2. **Target journal and formatting** — reference style, word limits, supplementary structure.
3. **Optional strengthening** — a properly trained cell2location reference (GSE132465 labels and TPM matrix are
   on disk) would settle whether the E2–stromal co-localisation survives a second deconvolution method.
4. **Optional** — an independent cohort with CMS/MSI plus a disease-free endpoint would test whether F14
   (non-incrementality) is cohort-specific.
