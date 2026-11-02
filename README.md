# DLBCL-PRISM

**P**rognostic **R**esource for **I**ntegrated **S**urvival **M**eta-analysis in diffuse large B-cell lymphoma.

An interactive, serverless web app (Shiny + [Shinylive](https://posit-dev.github.io/r-shinylive/) / webR) that lets anyone explore gene-level prognostic meta-analysis across 24 non-redundant public DLBCL cohorts (6,058 patients) directly in the browser.

## Analyses

| Tab | What it does |
|---|---|
| **Overview** | Summary statistics, non-redundant cohort catalog, clinical-data availability heatmap |
| **Meta forest** | Per-gene forest plot (REML random effects with Hartung-Knapp-Sidik-Jonkman CI), leave-one-out, exclude cohorts on the fly, PNG/PDF/CSV export |
| **Discovery** | Volcano filtered by minimum number of cohorts *k*, FDR, effect size, I², direction consistency; hover tooltips; top-hits summary forest; results table. Unannotated and multi-gene probe identifiers are excluded |
| **Clinical associations** | Gene z-score vs COO, IPI, stage, ECOG, LDH, age, sex, histology, treatment, early failure on R-CHOP-based therapy; Wilcoxon/Kruskal/Spearman + cohort-adjusted linear model; per-cohort Hedges' g meta-analysis; all-variable summary table with FDR |
| **Survival** | Kaplan-Meier by within-cohort median/tertiles/quartiles, cohort-stratified Cox (optionally adjusted for IPI, COO, age, sex, stage), subgroup forest with interaction tests |
| **Treatment & response** | Early treatment failure on R-CHOP-based therapy (event < 24 months, POD24-like): per-gene boxplot with Wilcoxon test, per-cohort logistic forest pooled by REML-HKSJ, cohort- and IPI-adjusted OR; Kaplan-Meier within R-CHOP-based therapy (regimen, endpoint, COO and IPI filters; cohort-stratified log-rank and Cox, 24/60-month estimates) |
| **Pathways** | Pre-ranked GSEA (fgsea; MSigDB Hallmark and Reactome) of the pooled meta-analysis statistic for any survival analysis or early failure; top-pathway bar chart, enrichment plot with clickable leading-edge genes, table |
| **Microenvironment** | Per-patient cell-type / program scores from canonical markers (T, CD8, Treg, NK, macrophages, M2, DC, neutrophils, IFN-γ, checkpoints, MHC-II, matrix/CAF, endothelium, proliferation, B-cell); pooled prognostic value, correlation with the selected gene, association with COO, early failure on R-CHOP-based therapy (boxplot and per-cohort OR forest), Kaplan-Meier by score groups, per-cohort forest |
| **Signature** | Multi-gene score with **leave-one-cohort-out** weights, validation forest, KM tertiles, C-index vs IPI with LRT |

Every plot can be downloaded as PNG (300 dpi) or PDF.

## Repository layout

```
app/app.R            Shiny app (runs in desktop R or in the browser)
prep/prep_data.R     Builds docs/data/ from the raw cohort files
prep/prep_extra.R    HKSJ pooled results and genome-wide early failure on R-CHOP-based therapy
prep/prep_tme_gsea.R GSEA and microenvironment scores (needs fgsea, msigdbr)
docs/                Static site served by GitHub Pages (Shinylive export + data/)
docs/data/           Compact data: clinical, pooled results, gene chunks (fetched on demand)
build.R              Re-exports the app to docs/ with shinylive
```

## Build & deploy

```r
# 1. (only when the data change) rebuild the web data, from the DLBCL-PRISM folder
Rscript prep/prep_data.R
Rscript prep/prep_extra.R
Rscript prep/prep_tme_gsea.R
# 2. export the app to docs/ (keeps docs/data/)
Rscript build.R
# 3. preview locally
Rscript -e 'httpuv::runStaticServer("docs", port = 8008)'
```

Push to GitHub, then **Settings → Pages → Deploy from a branch → `main` / `/docs`**.
The app will be at `https://<user>.github.io/<repo>/`.

Run in desktop R (no export needed): `shiny::runApp("app")`; data are read from `docs/data/`.

## Data notes (important for the manuscript)

* Expression is z-scored **within each cohort**; HRs are per 1 SD. Values are stored at 0.1-SD resolution.
* Survival times were harmonized to months (12 cohorts were reported in days).
* Treatment response categories were not available; early failure (< 24 months) is used as a proxy.
* LDH reported in U/L without ULN (GSE117556) is set to missing.
