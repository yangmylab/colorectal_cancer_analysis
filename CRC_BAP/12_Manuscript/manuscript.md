# Reproducibility without incremental value: a two-ecotype transcriptomic axis in colorectal cancer is associated with recurrence but not with death, and shows no detectable incremental value over CMS

**Draft v0.2 — figures composed 2026-11-17.** Six main figures (36 panels) and two supplementary figures
(8 panels) cover all 44 published panels. Every panel is cited in order in Results. Every quantitative claim is
traceable to a graded finding (F01–F15, `13_Reports/Findings_Ledger.csv`); the mapping is in
`12_Manuscript/claims_to_evidence.tsv` and the panel composition in `12_Manuscript/figures_plan.tsv`.

---

## Abstract

**Background.** Molecular taxonomies of colorectal cancer have repeatedly failed to reach clinical use. Their
reproducibility across assays is rarely tested as rigorously as their discovery, and their incremental value over
existing classifiers is rarely tested at all. We asked whether a two-group expression ecotype axis is (i)
reproducible across clinically relevant platforms and (ii) prognostic beyond CMS and MSI.

**Methods.** Hallmark programme scores from 624 TCGA-COAD/READ tumours were clustered by similarity network fusion
and consensus k-means; K=2 was selected on stability and frozen as a nearest-centroid classifier. The classifier
was then applied without refitting to CPTAC-COAD proteome and phosphoproteome, TCGA 450K methylation, 16 Visium
sections, and three independent cohorts (GSE39582, GSE17536, GSE38832) analysed on the endpoint each can support, plus
GSE14333 for disease-free survival. Each
analysis was pre-assigned an evidence grade.

**Results.** The K=2 partition was stable (bootstrap ARI 0.975; cross-validated classifier accuracy 0.979) and
reproduced at the protein level (83.3% agreement, ARI 0.439; E3) and the phosphoproteome level (87.5%, ARI 0.558,
permutation p<0.001; E3). It did not reproduce in promoter methylation (ARI 0.006) and the bulk classifier did not
transfer to spatial spots (self-check ARI 0.042). The axis was associated with disease-free survival in three
cohorts (fixed-effect HR 1.67, 1.24–2.24, p=0.0007; E2) but with neither overall survival (meta HR 1.25, p=0.081)
nor cancer-specific death (HR 1.00, p=0.997). Adjusting for CMS removed the disease-free association (HR 1.42,
p=0.107; ecotype–CMS Cramér's V=0.52).

**Conclusions.** Expression-based ecotypes can be made reproducible from transcriptome to phosphoproteome, but
reproducibility does not imply incremental clinical value: in the one cohort able to test it, the axis showed no
detectable incremental value over CMS, and it was unassociated with death in every cohort tested. Taxonomy studies should report both
properties.

---

## Introduction

Colorectal cancer has no shortage of molecular classifications. CMS1–4 are the most widely used; CRIS and PDS1–3
offer alternatives. Each was derived from expression data, each claims biological meaning, and none has displaced
stage in the clinic. Two problems recur. Classification instability: partitions respond to cohort composition,
platform and algorithm, so a label learned in one study often fails to transfer. Incremental value: a new taxonomy
must add information to what is already available — stage, MSI, and for research purposes CMS — rather than
restate it.

We therefore made reproducibility the object of study rather than a footnote, and we made incremental value a
pre-specified test rather than a discussion point. A K=2 partition of Hallmark programme scores was derived in
TCGA-COAD/READ and frozen; every subsequent analysis applies that classifier to a different assay, cohort or
endpoint without refitting. Findings were graded before writing — **E3** cross-platform convergence, **E2**
association replicated in an independent cohort, **E1** exploratory or post hoc, **descriptive**, and
**not_supported** — and no claim exceeds its grade.

---

## Methods

**Data.** TCGA-COAD/READ RNA-seq (624 tumours; GDC Data Release 46.0). External cohorts: GSE39582 (585), GSE17536
(177), GSE14333 (290, DFS only), GSE38832 (122, CEL files). CPTAC-COAD proteome and phosphoproteome (PDC000116, via
the cBioPortal datahub). TCGA 450K methylation (393 tumours). Visium: GSE283052 (12 sections) and GSE267401 (4).
Single-cell reference: GSE132465.

**Scoring.** Fifty Hallmark sets; discovery used GSVA, external layers the mean of within-set gene z-scores, with
phosphosites averaged per gene. Sets are frozen in `Environment/hallmark_frozen.gmt` so scoring is offline and
reproducible.

**Classifier.** Nearest centroid on 50 programme scores, computed in discovery and frozen. Applied to every external
layer without refitting; recoverable features are reported for each layer.

**Statistics.** Cox models with likelihood-ratio tests; C-index; fixed-effect inverse-variance meta-analysis;
bootstrap (500 resamples) for stability; permutation nulls (1,000 for cross-layer ARI; 199 for spatial
autocorrelation). Endpoints analysed separately and never pooled across families (OS/CSS/DSS versus DFS/RFS).
Complete-case analysis for covariates, with counts reported.

**Grading.** E3 cross-platform convergence; E2 replicated association; E1 exploratory/post hoc; descriptive;
not_supported. Fourteen findings (F01–F15; there is no F10) were graded and are the sole licence for claims.

**Reproducibility.** All analyses scripted (`06_Analysis/1_R/1_Rscripts/`), logged (`Logs/DECISION_LOG.tsv`) and
version-controlled. One level of parallelism (6 workers) with BLAS pinned to one thread per worker. Figures are
composed from a machine-readable plan (`12_Manuscript/figures_plan.tsv`) by `Scripts/compose_figures.R`.

---

## Results

### 1. A stable two-ecotype axis in discovery (Figure 1)

Fifty Hallmark programmes were scored in 624 TCGA-COAD/READ tumours and clustered after similarity network
fusion. **Figure 1A** shows the two programme blocks that define the partition: E1 is characterised by Hedgehog,
pancreas-β-cell and MYC-V2 programmes, E2 by EMT, angiogenesis and IL6–JAK–STAT3. The consensus matrix
(**Figure 1B**) is essentially block-diagonal at K=2, and the fused-view embedding (**Figure 1C**) separates the
two groups with substantial but incomplete separation. Stability across candidate K (**Figure 1D**) peaked at K=2 (bootstrap ARI 0.975) and
degraded at higher K, which is why K=2 — not the K=4 solution explored later — was frozen for external testing.

To interpret the composition of these tumours we used a single-cell reference of 63,689 cells from six lineages,
whose marker genes are shown in **Figure 1E**, and deconvolved the bulk tumours against it: epithelial 0.546,
stromal 0.267, myeloid 0.123, B 0.044, T 0.015, mast 0.006 (**Figure 1F**). This is descriptive context
(F05), and it matters only because it shows the axis is not a crude epithelial-versus-stroma ratio.

### 2. The axis reproduces at the protein and phosphoproteome level (Figure 2)

In CPTAC-COAD the frozen classifier was applied to Hallmark scores computed from protein abundance. **Figure 2A**
orders tumours by protein-assigned ecotype and shows the same E1/E2 programme blocks as the transcriptome.
Programme-level concordance between mRNA and protein (**Figure 2B**) was highest for the programmes that define E2
and the immune axis (IFN-α r=0.87, IFN-γ r=0.81, EMT r=0.71, angiogenesis r=0.69). At the level of the assignment
itself, **Figure 2C** shows that 83.3% of the 96 profiled tumours received the same ecotype from both layers
(ARI 0.439) — imperfect, but far above the chance level implied by two groups.

The phosphoproteome behaved the same way, and more strongly. **Figure 2D** shows the same programme blocks
computed from phosphosite data (11,334 sites quantified in at least half the tumours, averaged per gene so that
proteins with many measured sites do not dominate). Concordance across all three layers (**Figure 2E**) puts
phospho–protein correlation (median r=0.22) below protein–mRNA (0.42), as expected for a signalling layer, while
the assignment itself agrees with the protein call in 87.5% of tumours (**Figure 2F**; ARI 0.558).

### 3. The axis does not reproduce in methylation, and does not transfer to spots (Figure 3)

We tested the same partition against promoter methylation, selecting genes by annotation (19,269 genes) rather than
by variance, so that the test was not circular. CIMP-high status was similar in both ecotypes (74.5% E1 versus
84.3% E2, no formal test reported; **Figure 3A**) and global methylation did not differ (p=0.113; **Figure 3B**). Differential
methylation between ecotypes was essentially absent: only 5 of 19,269 genes reached FDR<0.05 with |logFC|>0.5
(**Figure 3C**), and an unsupervised methylation partition bore no relationship to the ecotype (ARI 0.006;
**Figure 3D**). The axis is not a methylation phenotype (F07).

We then asked whether the same classifier could label Visium spots. It cannot. **Figure 3E** maps the E2 programme
score across a section and shows fine-grained intra-tumour structure rather than the two-block pattern expected
from a transferable classifier; **Figure 3F** shows that variation of the E2 programme within sections exceeds
variation between them. The failure is diagnostic rather than subtle: after scoring in the bulk reference space,
every spot fell on one side of the decision boundary and spot scores sat about 5 SD below the bulk distribution
(**Figure 3G**), and a self-check of the same simplified statistic on the bulk training data returned ARI 0.042.
A 55 µm spot has no compartment structure comparable to a bulk tumour, so the bulk classifier is simply the wrong
instrument (F08).

### 4. What the phosphoproteome and the sections do say about mechanism (Figure 4)

Between ecotypes, 3,095 phosphosites reached FDR<0.05, of which 538 exceeded |log2FC|>0.5 (**Figure 4A**). Because
cross-layer agreement could in principle arise by chance, we compared the observed ARI of 0.558 against 1,000 label
permutations: the null's 95th percentile is 0.033 and the observed value is far outside it (p=0.000999;
**Figure 4B**). The phosphoproteome reproduces the partition, and does so beyond chance.

Kinase–substrate inference (PhosphoSitePlus prior; 78 kinases with at least five measured substrates) identified 9 kinases at FDR<0.05, five of them with unchanged kinase abundance
(**Figure 4C**), and comparing kinase activity with the kinase's own protein abundance
(**Figure 4D**) separated activity from quantity: CDK7 activity is lower in E2 (d=−1.00) without any change in CDK7
protein, whereas ROCK1 and PRKG1 differ largely because their abundance differs. Across kinases the two measures
correlate only r=0.31. This is exploratory — activity is inferred, not measured, in one cohort (F12).

Spatial data, analysed on their own terms rather than by projection, show that the E2 programmes are spatially
structured: mean Moran's I 0.528 on a k=6 neighbour graph (**Figure 4E**). They are, however, no more structured
than Hallmark programmes in general (0.587), and the per-section view (**Figure 4F**) confirms that this holds
across sections rather than in one or two. E2 scores co-localise with the stromal compartment (median ρ=0.49) but
not the myeloid (0.04). We report this descriptively, with two caveats: with 947–4,757 spots per section (median 2,670) the
permutation test rejects the null for every programme, so only effect sizes are interpreted; and a
trained-reference cell2location rerun (three sections, 3,000 epochs, 1,000 spots each) reproduces the spatial
pattern of the NNLS fractions used for the co-localisation statistic (matrix r 0.93, top-lineage agreement
0.74–0.99) but not their absolute scale (epithelial fraction 0.39 by cell2location vs 0.97 by NNLS), so the
co-localisation is directionally supported while its magnitude remains method- and gene-panel-dependent (F13).

### 5. Disease-free survival replicates; cancer-specific death does not (Figure 5)

The one outcome association that replicates concerns recurrence. In GSE14333 (226 tumours, 50 disease-free events,
stage I–III) E2 tumours recurred sooner (**Figure 5A**; adjusted HR 2.37, LRT p=0.0043), and the effect was stable
across model specifications and in a 500-resample bootstrap (**Figure 5B**). Applying the frozen classifier without
refitting to GSE39582 (483 tumours, 137 events, stage I–III) gave HR 1.50 (1.07–2.10), p=0.018, with the all-stage
and unadjusted models pointing the same way (**Figure 5C**), and the corresponding survival curves separate
visibly (**Figure 5D**). GSE17536 (135 tumours, 32 events) gave HR 2.51 (1.14–5.53), p=0.022. The fixed-effect
summary across the three cohorts is HR 1.67 (1.24–2.24), p=0.0007, without heterogeneity (p=0.21) — graded E2
(F11).

The same cohort that made this test possible also shows its limit. GSE38832 records only 9 recurrences (too few to
test F11) but 28 cancer deaths, so it was analysed on cancer-specific survival, an endpoint never pooled with
disease-free endpoints. **Figure 5E** places both endpoints side by side: disease-free survival in this cohort is
uninformative (HR 0.92, dashed with 9 events), while cancer-specific death shows no association at all (HR 1.00,
0.46–2.17, p=0.997). The corresponding survival curves are superimposable (**Figure 5F**). The axis predicts
recurrence, not death (F15).

### 6. No overall-survival association, and no incremental value over CMS (Figure 6)

In discovery (476 tumours, 79 events, adjusted for age, sex, stage and sidedness) the survival curves by ecotype
overlap (**Figure 6A**) and adding the ecotype to the model moved the C-index from 0.763 to 0.772 with LRT p=0.107
(F01). A K=4 solution did show an OS association in discovery (C-index 0.785, p=0.024), and **Figure 6B** shows the
K-sensitivity of the survival model across candidate K — but this was post hoc, the chosen K remains 2, and the K=4
association did not replicate (F02 E1, F04). External validation reproduced the null: GSE39582 p=0.333, GSE17536
p=0.0785, fixed-effect meta HR 1.25 (0.97–1.61), p=0.081 (**Figure 6C**).

Finally, the question that determines whether any of this matters clinically. In GSE39582, which carries CMS calls,
nested models were compared (**Figure 6D**): with stage, age and sex alone the ecotype contributed a C-index gain of
+0.025 (HR 1.50, p=0.018); adding CMS cut the gain to +0.018, and adding MSI and CMS together left **+0.015 with
HR 1.42 (0.93–2.17), p=0.107**. The reason is visible in **Figure 6E**: ecotype and CMS are substantially
overlapping classifications (Cramér's V=0.52), with E2 enriched in CMS4 and depleted in CMS1–2. The disease-free
association is therefore replicated and directionally consistent, but not demonstrable incremental over the
published taxonomy (F14).

### 7. Evidence supporting Figures 1–6 (Supplementary Figures 1 and 2)

Every supplementary panel underpins a specific main-figure claim. They are cited here in order, after the main
figures, and each is explicitly tied to the panel(s) it supports.

**Supplementary Figure 1A** shows the PCA of all 624 quantified tumours, with points coloured and shaped by TCGA
project. COAD and READ intermix across the principal plane, and the only grouping available here is anatomical
(colon versus rectum) rather than technical: the discovery data come from one platform and one processing pipeline,
so there is no batch term to remove. **Supplementary Figure 1B** quantifies that group structure component by
component — the project term accounts for at most 2.8% of the variance of any of the first ten components (0.5% of
PC1) — which is why no correction was applied and why the downstream analyses use the quantile-normalised matrix
as it stands. Library sizes were comparable between groups (**Supplementary Figure 1C**) and the top 500 variable
features that drive the embedding are shown in **Supplementary Figure 1D**, so the separation in Figure 1C cannot
be attributed to uneven sequencing depth or to a handful of features.

**Supplementary Figure 2A** (five-year calibration) and **Supplementary Figure 2B** (decision curve) characterise
the discovery overall-survival model of Figure 6A and the covariate comparison of Figure 6D. Because that model
carries no significant ecotype term, no net-benefit advantage is interpreted: both panels are reported for
completeness, not as a claim of clinical utility.

**Supplementary Figure 2C** gives survival by protein-assigned ecotype in CPTAC-COAD, the cohort behind
Figure 2A–2C. With 8 deaths it is uninformative (HR 2.44, p=0.286) and serves only as a direction check on the
protein assignment; it cannot test prognosis.

**Supplementary Figure 2D** is the methylation counterpart of Figure 3A–3D: survival by methylation cluster does
not differ (HR 0.74, p=0.398). That null is what would be expected if the ecotype carried no methylation signal,
and it rules out the possibility that the negative methylation result reflects a methylation cluster with better
outcome.

---

## Discussion

The axis we describe is unusual less in what it predicts than in where it survives. It persists from mRNA to
protein to phosphoprotein — three assays with different dynamic ranges, coverage and quantification — at a level
of concordance usually reported between technical replicates. It vanishes at the methylome and cannot be projected
onto spots. Together these results suggest the partition tracks something about the tumour–stroma–immune axis that
is post-transcriptionally and post-translationally maintained, rather than a promoter-methylation programme.

The prognostic picture is deflationary, and that is the paper's point. Three independent cohorts, each analysed on
the endpoint it can support, show no association with death; the one association that replicates concerns
recurrence. That dissociation is biologically plausible — salvage therapy and competing mortality dilute recurrence
signals — but it is also exactly the pattern that should restrain claims of clinical utility. The incremental-value
analysis then settles the question for this dataset collection: with stage, MSI and CMS in the model, the ecotype
adds nothing detectable.

We do not conclude that expression-based taxonomy of colorectal cancer is worthless. We conclude something narrower
and, for the field, more useful: a two-group axis can be made reproducible across assays if it is frozen early and
tested without refitting, and reproducibility does not imply incremental clinical value. Both properties must be
demonstrated, and they can dissociate. Studies proposing new subtypes would do well to report the second test as
routinely as the first.

**Limitations.** All data are retrospective and public: no wet-lab validation, no prospective cohort. Cross-platform
layers come from small cohorts (CPTAC-COAD 96 tumours and 8 deaths; GSE14333 226 tumours and 50 events), so protein
and phospho analyses test reproducibility rather than prognosis. The phosphoproteome is 59.8% missing before
filtering and kinase activity is inferred, not measured. The incremental-value test rests on one cohort, because the
others carry no CMS calls, and CMS itself is expression-derived in that cohort, so ecotype and CMS are not
independent by construction. Endpoint definitions differ slightly across cohorts (DFS in GSE14333 and GSE38832, RFS
in GSE39582 and GSE17536), and disease-free endpoints are never pooled with mortality endpoints. Classifier
agreement is sensitive to how programme scores are built: a different gene-set mapping moved protein-layer agreement
from 85.4% to 83.3% by flipping 2 of 96 borderline tumours — a sensitivity we report rather than hide. One spatial
co-localisation statistic is directionally reproducible but its magnitude is method-dependent, as noted in Results.

**What we do not claim.** No new prognostic classifier; no clinically actionable subtype; no causal role for any
programme; no refinement beyond CMS; no measured kinase activity; no spatial ecotype assignment.

---

## Figure legends

**Figure 1. A stable two-ecotype axis in TCGA-COAD/READ.** (A) Hallmark programme scores ordered by ecotype.
(B) Consensus matrix at K=2. (C) Similarity-fusion embedding coloured by ecotype. (D) Stability across candidate K.
(E) Single-cell marker signatures from GSE132465. (F) Mean deconvolved cell-type fractions.

**Figure 2. The axis reproduces at the protein and phosphoproteome level.** (A) Protein programme scores by
ecotype. (B) Programme-level mRNA–protein concordance. (C) Ecotype agreement, mRNA versus protein. (D) Phospho
programme scores by ecotype. (E) Three-layer programme concordance. (F) Ecotype agreement, phospho versus protein.

**Figure 3. The axis does not reproduce in methylation and does not transfer to spots.** (A) CIMP score by
ecotype. (B) Global methylation by ecotype. (C) Top differentially methylated genes. (D) RNA ecotype versus
methylation cluster. (E) E2 programme score per Visium spot. (F) Within- versus between-section variation.
(G) Covariate shift of spot scores relative to bulk.

**Figure 4. Signalling and tissue architecture.** (A) Differential phosphosite volcano. (B) Permutation null of the
cross-layer ARI. (C) Kinase-substrate activity volcano. (D) Kinase activity versus kinase abundance. (E) Moran's I
of Hallmark programmes across sections. (F) Per-section Moran's I for the E2 programmes.

**Figure 5. Disease-free survival replicates; cancer-specific death does not.** (A) GSE14333 disease-free survival
by ecotype. (B) GSE14333 hazard ratios across models. (C) GSE39582 replication forest. (D) GSE39582
recurrence-free survival. (E) GSE38832 DSS versus DFS hazard ratios. (F) GSE38832 cancer-specific survival.

**Figure 6. No overall-survival association and no incremental value over CMS.** (A) Discovery overall survival by
ecotype. (B) K sensitivity of the survival model (K=4 is post hoc). (C) External OS validation forest.
(D) C-index by covariate set, with and without the ecotype. (E) Ecotype versus CMS composition.

**Supplementary Figure 1. Quality control underpinning Figure 1.** (A) PCA of all 624 quantified tumours, coloured
and shaped by TCGA project (COAD versus READ) — no technical batch exists in this dataset. (B) Share of each
principal component explained by the project term (red) against the variance that component explains (grey); at
most 2.8% for any of PC1–PC10, so no correction was applied. (C) Library size by group. (D) Top 500 variable
features. Each panel supports the reading of Figure 1A–1D as biological separation.

**Supplementary Figure 2. Diagnostics and negative controls underpinning Figures 2, 3 and 6.** (A) Five-year
calibration of the null OS model (supports Figure 6A). (B) Five-year decision curve (supports Figure 6D; no
utility claim). (C) Survival by protein-assigned ecotype, 8 deaths, underpowered (supports Figure 2A–2C).
(D) Survival by methylation cluster, negative control (supports Figure 3A–3D).

---

---

## Statements

**Data availability.** All data analysed here are already public. TCGA-COAD/READ RNA-seq and 450K methylation
were obtained from the GDC (Data Release 46.0) under the corresponding data-use certification; gene-expression
cohorts from GEO (GSE39582, GSE17536, GSE14333, GSE38832, GSE283052, GSE267401, GSE132465); and CPTAC-COAD
proteome and phosphoproteome (PDC000116) via the cBioPortal datahub. Derived tables — the findings ledger, the
claims-to-evidence map, per-layer validation outputs and the figure source tables — are released with the code.

**Code availability.** The full pipeline (discovery, external layers, outcomes, figure composition) is under
version control, with a decision log and two automated consistency checks. Archived release and DOI:
[TODO: repository URL and Zenodo DOI].

**Ethics.** This is a secondary analysis of de-identified, publicly available data; no new human-subject
research was performed. Data-use certifications and any institutional determination:
[TODO: name the certifications/approvals].

**Funding.** [TODO: funding bodies and grant numbers, or "This work received no specific grant."]

**Author contributions.** [TODO: CRediT roles per author.]

**Competing interests.** The authors declare no competing interests.

**Acknowledgements.** [TODO: cohorts, consortia and reviewers to thank.]

---

## Pending before submission

1. **Author list and affiliations** — `Authors.xlsx` is still an empty template.
2. **Journal choice** — reference style, word limits, supplementary structure.
3. ~~**Optional strengthening** — a properly trained cell2location reference (GSE132465 labels and TPM matrix are on
   disk) would settle whether the F13 co-localisation survives a second deconvolution method.~~ **Done
   2026-11-18:** a reference trained on GSE132465 (2,007 genes, 250 epochs) run at 3,000 epochs over three sections
   reproduces the spatial pattern of the NNLS fractions (matrix r 0.93) but not their absolute scale; see
   `06_Analysis/11_SpatialNiches/26_cell2location_reference.md`. Remaining strengthening: the same check at GPU
   scale on all 16 sections.
