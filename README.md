# Proximity — Spatial Immune Profiling for Melanoma

**HIMPact Hacks '26 submission** · Christine Li, solo

> Live demo: `website/index.html` (open directly, or deploy to GitHub Pages)
> Analysis pipeline: `analysis/`

## The problem

Standard melanoma pathology reports **how many** CD8+ cytotoxic T-cells are present in a tumor sample. But a T-cell walled off from the tumor by cancer-associated fibroblast stroma cannot kill the cancer cells it can't reach — **proximity, not just count, determines whether the immune system can functionally engage the tumor.**

**Research question:** How does the average nearest-neighbour (NN) distance between CD8+ T-cells and malignant melanocytes correlate with tumor thickness (Breslow depth) and patient survival outcomes?

### ITN framing
- **Importance** — melanoma drives a disproportionate share of skin-cancer mortality; immunotherapy response varies widely in ways standard cell-count pathology doesn't explain.
- **Tractability** — mIF imaging already exists in research pipelines; nearest-neighbour spatial statistics (the point-pattern "G-function") are established math. The gap is a fast, interpretable layer between coordinates and a clinical recommendation.
- **Neglectedness** — tools like Squidpy compute spatial statistics for computational biologists, but few translate that into an auditable, plain-language readout a practicing oncologist can use directly.

## What's in this repo

```
himpact/
├── data/
│   └── generate_synthetic_tme.py     # synthetic cell-coordinate generator
├── analysis/
│   ├── spatial_analysis.py           # NN distance, stats, correlation, recommendation logic
│   ├── make_visuals.py               # distribution + density heatmaps, correlation plots
│   └── make_case_visuals.py          # per-case demo images
├── website/
│   ├── index.html                    # the tool's clinician-facing front end
│   └── assets/                       # generated figures
└── docs/
    └── (literature review, project proposal — see appendices)
```

## The method

1. **Nearest-neighbour distance** — for every CD8+ T-cell, the Euclidean distance to the nearest melanocyte, computed with a KD-tree (`scipy.spatial.cKDTree`) for O(n log n) queries.
2. **Aggregate statistics per patient** — mean and standard deviation of NN distance. Mean captures typical proximity; standard deviation captures whether exclusion is uniform or patchy across the tumor.
3. **Density vs. distribution** — distribution shows *where* cells sit (raw scatter); density shows *how concentrated* they are (2D heatmap). The two can diverge in ways a scatter plot alone hides.
4. **Correlation** — Pearson r (plus Spearman's ρ as a non-parametric check) between mean NN distance and Breslow depth / survival outcome.
5. **Recommendation logic** — a small, explicit rule set (not a black-box model) flags cases as Immune-Infiltrated / Intermediate / Immune-Excluded, using a 50 μm functional-engagement threshold drawn from the literature review. Auditability was a deliberate design choice: a physician can see exactly why a case was flagged, which matters both for trust and for any future regulatory review.

## Why synthetic data

Real HTAN / TCGA-SKCM multiplexed IF data requires registered access, multi-GB downloads, and QuPath-based cell segmentation — none of which fit a hackathon weekend. Instead, `generate_synthetic_tme.py` builds 60 virtual patients whose CD8+ T-cell exclusion distance scales with a simulated Breslow depth, **using the exact stromal-exclusion mechanism described in the literature review** (CAFs physically barrier T-cells from tumor cells, with the effect intensifying in thicker tumors), plus realistic per-patient noise and heterogeneity so the resulting correlation is representative rather than trivially perfect.

**Results on synthetic data (n = 60):**

| Relationship | Statistic | p-value |
|---|---|---|
| NN distance vs. Breslow depth | Pearson r = 0.53 | p = 1.1 × 10⁻⁵ |
| NN distance vs. Breslow depth (robustness) | Spearman ρ = 0.49 | p = 8.2 × 10⁻⁵ |
| NN distance vs. simulated survival | Pearson r = −0.10 | p = 0.44 (not significant) |

The survival result is reported as-is, not tuned away — a real validation cohort with more patients and longer follow-up would be needed to properly power that comparison.

## Limitations

- Synthetic, not real, patient data — correlation strength on real biopsies will differ.
- Survival correlation is weak and not statistically significant at this sample size.
- The 50 μm engagement threshold is drawn from the literature review, not independently validated against outcomes.
- No cell-segmentation step — the pipeline assumes pre-extracted coordinates, not raw microscopy images.
- Only two cell types modeled (CD8+ T-cells, melanocytes); real melanoma TME involves CAFs, TAMs, and Tregs (see literature review), deliberately excluded here for tractability.

## Where this goes next

- Extend to a multi-cell-type spatial graph (add Tregs, TAMs, CAFs) for a fuller exclusion picture.
- Empirically validate the 50 μm threshold rather than adopting it directly from prior literature.
- Integrate QuPath cell-segmentation output so the pipeline runs on raw mIF slides.
- Apply the same spatial framework to other solid tumors where T-cell exclusion is prognostically relevant.

## Path to clinical use (directional, not regulatory advice)

| Phase | Milestone | Est. timeline |
|---|---|---|
| 0 | Retrospective validation on a real HTAN/TCGA-SKCM cohort | 3–6 mo |
| 1 | Multi-site concordance study vs. pathologist manual scoring | 6–12 mo |
| 2 | SaMD risk classification & pre-submission meeting (FDA) / IVDR classification (EU) | 3–6 mo |
| 3 | Clinical validation study sized for regulatory submission | 12–18 mo |
| 4 | 510(k) / De Novo (US) or IVDR technical file (EU) | 6–12 mo review |
| 5 | Post-market surveillance | ongoing |

## Illustrative financial model

| Year | Focus | Est. revenue |
|---|---|---|
| Y1 | Research-use-only licensing to academic pathology labs | $0–15K |
| Y2 | Pilot deployments, grant-funded validation studies | $15K–75K |
| Y3 | Post-clearance per-case SaaS licensing to pathology labs | $150K–500K |

Directional planning estimates for a pre-revenue research tool, not a funded projection.

## Appendices

- **Website code:** `website/` (this repo)
- **Tool / analysis code:** `data/`, `analysis/` (this repo)
- **Data:** fully synthetic, generated by `data/generate_synthetic_tme.py` — no real patient data used, so no de-identification step is applicable (would be required under HIPAA's 18-identifier Safe Harbor standard for any real cohort data used in future work).
- **Literature review & project proposal:** `docs/`

## AI use disclosure

Claude (Anthropic) was used to: draft the synthetic data generation methodology and analysis code, build the website front end, and structure this writeup, based on the author's original research question, literature review, and project direction. All modeling assumptions, thresholds, and limitations were reviewed and are disclosed above.
