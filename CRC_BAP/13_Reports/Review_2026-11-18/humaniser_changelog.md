# AI-accent removal — change log (2026-11-18)

Scope: `12_Manuscript/manuscript.md`, rebuilt to `manuscript.docx`. Method: `academic-humanizer`
(Layer 1 general tells, Layer 2 academic tells, Layer 3 preservation rules, Layer 4 claim-evidence).
No figure, supplementary-figure or citation was reordered, added or removed.

## Invariants verified after the edit

| Check | Result |
|---|---|
| Numeric tokens (441) before vs after | identical multiset |
| Figure/supplementary citations (53) | identical sequence |
| First-citation order | unchanged: 1A-1F, S1A-S1D, 2A-2F, S2C, 3A-3G, S2D, 4A-4F, 5A-5F, 6A-6E, S2A, S2B |
| `check_figure_citations.py` | PASS |
| `audit_analysis.py` | 22/22 |

## Patterns removed

**Em-dashes: 24 → 0.** Recast with commas, colons, parentheses or full stops. Examples:
- "what is already available — stage, MSI, and for research purposes CMS — rather than restate it"
  → "what is already available (stage, MSI, and, for research purposes, CMS) rather than restate it"
- "Findings were graded before writing — E3 … not_supported — and no claim exceeds its grade"
  → "Findings were graded before writing: … . No claim exceeds its grade."
- "K=2 — not the K=4 solution explored later — was frozen"
  → "K=2, and not the K=4 solution explored later, was frozen"
- "Derived tables — the findings ledger, … — are released with the code"
  → "Derived tables (the findings ledger, …) are released with the code"

**Formulaic openers.** "Quality control supports that reading." removed (the sentence was redundant with the
paragraph's closing clause). "A negative control closes the alternative reading." → "A methylation-cluster
control addresses the opposite explanation." "Two diagnostics support reading this null as a null rather than a
mis-specified model" → "Two diagnostics check that the null reflects the data rather than a mis-specified model."
"Finally, the question that determines whether any of this matters clinically." → "The clinical question is
whether any of this adds value."

**Claims matched to evidence.**
- "The axis predicts recurrence, not death (F15)" → "The axis is associated with recurrence but not with death
  (F15)" (also closes review R06's distinction between association and prediction).
- "far above the chance level implied by two groups" → "well above the agreement expected under a two-group
  chance model".
- "simply the wrong instrument" → "the wrong instrument"; "exactly the pattern" → "the pattern".

**Self-congratulatory framing.** "a sensitivity we report rather than hide" → "we report that sensitivity".

**Negative antithesis.** "We do not conclude that expression-based taxonomy … is worthless. We conclude
something narrower and, for the field, more useful:" → "Our conclusion is narrower:".

**Clause-stacked sentences split.** The grading list, the dissociation sentence, and the qualification of the
protein-assignment result.

## Left deliberately

33 sentences remain above 30 words and 12 uses of "rather than", because where they survive they carry a real
contrast (e.g. "anatomical rather than technical") rather than padding. Trimming those further trades precision
for rhythm, so it is left for the author to judge.
