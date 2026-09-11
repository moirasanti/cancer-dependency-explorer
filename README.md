# Cancer Dependency Explorer

**Which molecular features are associated with dependency on a selected cancer target, and which exploratory molecular predictors are reproducible when the data are considered jointly?**

Cancer Dependency Explorer is a reproducible R workflow that integrates DepMap CRISPR dependency with cancer lineage, expression, damaging mutations and copy number. The integrated modelling uses elastic net so conditional or suppressor relationships can emerge without needing to pass a marginal association threshold first.

The workflow produces data quality summaries, target self-associations, analyses on pre-specified hypotheses, within-lineage follow-up, whole-genome analyses, results on an integrated model with performance evaluated with cross-validation, and a self-contained HTML report.

## Run a target assessment

Create the environment:

```bash
conda env create -f environment.yml
conda activate cancer-dependency-explorer
```

Clone the repository in your local (see structure in **Repository structure**), download the files listed in [data/README.md](data/README.md) and place them in `data/raw/`, then run:

```bash
Rscript run_analysis.R \
  --target SOX10 \
  --release "DepMap Public 26Q1" \
  --raw-dir data/raw \
  --biomarkers config/SOX10_biomarkers.yml
```

`DEPMAP_RAW_DIR` and `DEPMAP_RELEASE` can replace `--raw-dir` and `--release`.

The YAML file is optional (see below **Pre-specified hypotheses**). When `--biomarkers` is omitted, the workflow automatically uses `config/<TARGET>_biomarkers.yml` if it exists. A supplied file must follow that filename convention, and its `target` value must match `--target`.

## Pre-specified hypotheses

A pre-specified hypothesis is a molecular feature nominated before modelling from biological knowledge in the YAML configuration. The designation does not assert an association or assign the feature preferential model weight. The files look like this:

```yaml
target: SOX10
expression:
  - MITF
mutation:
  - NF1
copy_number:
  - MITF
```

Every YAML feature is first evaluated pan-cancer in the same modality wide univariate screen as every other gene. In the integrated model, an eligible YAML feature receives the same penalty and stability requirements as every other molecular feature. Eligibility still applies:

- Expression and copy number require at least 80% coverage among models with any data from the corresponding assay, plus nonzero variance. Models without that assay are retained for integrated modelling and receive training fold median imputation; the report shows overall assay coverage separately.
- Damaging mutations require at least 3 altered and 3 wild type models.

Eligible damaging-mutation associations are tested with R's built-in unequal-variance Welch t-test. The reported mean difference, 95% confidence interval and P value come from the same test, followed by modality-wide FDR correction.

##### What is the purpose of this YAML file?

To have a highlighted section when some biological pre-specified hypotheses are known, and quickly follow up on their results.

The target's eligible expression, damaging mutation and copy number follow the same modelling rule. They are reported as target characterisation instead of pre-specified hypotheses. Only cancer lineage indicators are unpenalized in the integrated model.

After the pan-cancer results, eligible target/pre-specified-hypothesis associations are estimated separately within each lineage containing at least 10 complete models. Mutation associations additionally require 3 altered and 3 wild type models in that lineage. FDR is corrected jointly across all estimable feature by lineage tests in this family.

### SOX10 worked example

The retained worked example assesses SOX10, a lineage dependency in melanoma. Its pre-specified expression hypotheses represent the melanocytic SOX10/MITF transcriptional state (MITF, PAX3 and TFAP2A), an upstream SOX10 regulator (TYRO3), and receptor tyrosine kinase features associated with SOX10/MITF state or adaptation (ERBB3 and AXL). NF1, CDKN2A and PTEN damaging mutations represent loss of function melanoma contexts. MITF and EP300 copy number test reported genomic relationships with the SOX10 programme.

These hypotheses are motivated by experimental studies showing that SOX10 depletion restricts melanoma cell proliferation and alters MITF and cell-cycle regulation, that PAX3/SOX10/MITF cooperate in melanoma transcriptional control, and that SOX10/MITF state is linked to receptor tyrosine kinase programmes. YAML designation makes these explicit hypotheses but does not guarantee eligibility, association or selection.

## Integrated exploratory model

Every eligible expression, damaging mutation and copy number feature enters one integrated molecular search with penalty factor 1. Target characterisation, pre-specified hypothesis and genome wide candidate are biological source labels only, they do not affect model treatment. Cancer lineage indicators enter with penalty factor 0 to adjust for histological context. This means that lineage will always be present in the model, whilst other features may receive penalisation factors to drive their coefficients to zero. A feature with little individual association can become informative after lineage and other molecular data are considered jointly. Elastic net handles correlated predictors during that joint fitting.

The model uses `glmnet` with:

- Gaussian family because Chronos dependency is continuous.
- `alpha = 0.5`, combining lasso sparsity with ridge-like handling of correlated biological signals.
- The one-standard-error lambda, favouring the most regularized compact model whose inner-CV error remains within one standard error of the minimum.

Preprocessing is learned inside each training fold. Continuous values use training fold median imputation and standardization. Mutation values use training fold modal imputation, binary encoding, and then centering and scaling from the training fold mutation prevalence and SD. This gives all molecular predictors a comparable one-SD scale. Training derived transformations are applied unchanged to held-out models.

The pipeline runs 3 repeats of lineage stratified 5-fold outer CV, with lineage-stratified 5-fold inner CV for lambda selection. A molecular predictor is stable only if it:

- Is selected in at least 60% of the 15 outer models.
- Has the same coefficient direction in at least 80% of models in which it is selected.
- Is nonzero in the full data elastic net fit.

Zero stable molecular predictors is a valid result. In that case the pipeline reports up to ten highest frequency nonzero genome wide near misses as explicitly unstable exploratory candidates.

## Result terminology

- **Univariate association:** descriptive one feature at a time pan-cancer result.
- **Within-lineage association:** eligible target/pre-specified-hypothesis association estimated separately within one lineage; this is descriptive and is not a formal test that effects differ between lineages.
- **Target characterisation:** an eligible target expression, damaging mutation or copy number feature.
- **Pre-specified hypothesis:** an eligible YAML-listed molecular feature nominated before modelling.
- **Genome wide candidate:** any other eligible molecular feature.
- **Lineage adjustment:** the only predictor class included with zero penalty; it allows different lineage baselines but does not include feature-by-lineage interactions.
- **Stable exploratory molecular predictor:** any molecular predictor meeting all frequency, direction and full fit criteria.
- **Near miss:** nonzero in at least one outer model but failing the stability criteria; it is not retained.
- **Predictive performance:** nested-CV held-out Pearson r, R² and root mean squared error (RMSE) for the whole procedure, an observed-versus-predicted diagnostic for each repeat, and mean absolute error (MAE), RMSE and within lineage R² for lineages with at least 10 models. Pearson r measures whether predictions and observations move together. Prediction error R² compares squared errors with assigning every model the observed mean, RMSE reports error size in Chronos units.

## Workflow – What is happening?

1. Validate model identifiers, feature availability, coverage, variance and mutation group sizes.
2. Run descriptive pan-cancer expression, damaging mutation and copy number screens.
3. Report target self-associations and YAML hypotheses, then follow eligible ones within lineages.
4. Build one predictor matrix containing all eligible molecular features plus lineage adjustment.
5. Run repeated nested CV with fold specific preprocessing and elastic net selection.
6. Report designation, compact coefficients, stability, held out performance, figures and the HTML assessment.

## Outputs

Each run writes a target-specific result tree and a self-contained report. The principal integrated-model interfaces are:

```text
results/<TARGET>/tables/12_feature_designation.csv
results/<TARGET>/tables/13_integrated_model_coefficients.csv
results/<TARGET>/tables/14_selection_stability.csv
results/<TARGET>/tables/15_nested_cv_performance.csv
results/<TARGET>/tables/16_nested_cv_predictions.csv  # generated locally, excluded from Git
results/<TARGET>/tables/17_lineage_cv_performance.csv
results/<TARGET>/intermediate/integrated_model.rds
reports/<TARGET>_target_assessment.html
```

## Repository structure

```text
.
├── R/util.R
├── config/<TARGET>_biomarkers.yml
├── data/README.md  # data dictionary to work with
├── reports/target_assessment.Rmd  # report script
├── results/<TARGET>/
├── scripts/
|    ├── 01.prepare.depmap.data.R
|    ├── 02.dependency.landscape.R
|    ├── 03.molecular.associations.R
|    ├── 04.multivariable.model.R
|    └── 05.generate.figures.R
├── environment.yml
└── run_analysis.R
```

## Limitations

- DepMap associations can reflect lineage composition or unmeasured confounding.
- Nested CV reduces selection/evaluation leakage but is not an independent external validation.
- Correlated features may share signal, so selection frequencies can be distributed across substitutes.
- Cell line knockout phenotypes may not reproduce in tumours, normal tissues or other experiments.

Chronos scores become more negative as knockout has a stronger effect on cell fitness. See Dempster *et al.*, [Chronos: a cell population dynamics model of CRISPR experiments](https://doi.org/10.1186/s13059-021-02540-7).

The SOX10 hypothesis rationale draws on experimental work on [SOX10 loss and melanoma proliferation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3803156/), [PAX3/SOX10/MITF regulation in melanoma](https://pmc.ncbi.nlm.nih.gov/articles/PMC2979310/), [SOX10 addiction and copy-number gain](https://pmc.ncbi.nlm.nih.gov/articles/PMC5540806/), and the [MITF-low/AXL-high adaptive state](https://pmc.ncbi.nlm.nih.gov/articles/PMC4428333/).

The example data use [DepMap Public 26Q1](https://depmap.org/portal/data_page/?tab=currentRelease). Raw matrices, model-level intermediates and cell-line-level prediction tables remain excluded from version control.
