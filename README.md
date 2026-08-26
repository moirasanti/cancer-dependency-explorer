# Cancer Dependency Explorer

**Which molecular features predict dependency on a selected cancer target?**

Cancer Dependency Explorer is a reproducible R workflow for investigating why cancer cell lines differ in their dependency on a selected target gene. The selected gene’s CRISPR dependency score is treated as the outcome. The pipeline tests whether differences in dependency are associated with cancer lineage, expression and copy number of the target, damaging mutations elsewhere in the genome, and genome-wide copy-number features.
The workflow produces data-quality summaries, univariate association screens, an interpretable multivariable linear regression with HC3-adjusted uncertainty, standardised figures and a self-contained HTML report.

#### Worked example: MDM2

The included example treats the MDM2 CRISPR dependency score as the outcome and asks:

> Why are some cancer cell lines more dependent on MDM2 than others?

To address this question, the pipeline tests whether stronger or weaker MDM2 dependency is associated with cancer lineage, MDM2 expression, MDM2 copy number, or damaging mutations in other genes.
Using DepMap Public 26Q1, the worked example identifies TP53 mutation status as a major molecular feature associated with MDM2 dependency. It also demonstrates how prior biological hypotheses and genome-wide discoveries can be evaluated alongside target expression, copy number and lineage, while reducing redundancy among highly correlated predictors.
These associations are intended to generate hypotheses. They do not demonstrate causality or establish that the selected target is therapeutically actionable.

## Run a target assessment

Create the environment:

```bash
conda env create -f environment.yml
conda activate cancer-dependency-explorer
```

Download the files listed in [data/README.md](data/README.md). These are the raw files from the portal, and they need processing. Save these files in data/raw/ before continuing with the pipeline.
This workflow has chained the numbered scripts into one workflow given the genes of interest (``--target MDM2`` or any target of interest), run:

```bash
Rscript run_analysis.R \
  --target MDM2 \
  --raw-dir /path/to/depmap/26Q1
```

`DEPMAP_RAW_DIR` can replace `--raw-dir`. 

An optional YAML file can pre-specify biologically motivated biomarkers for inclusion in the model before data-driven candidates are selected:

```bash
Rscript run_analysis.R \
  --target MDM2 \
  --raw-dir /path/to/depmap/26Q1 \
  --biomarkers config/MDM2_biomarkers.yml
```

The configuration format is:

```yaml
target: MDM2
mutation:
  - TP53
copy_number: []
```

or

```yaml
target: MDM2
mutation:
  - TP53
copy_number:
  - gene1
  - gene2
```

The ``target: value`` inside the file must match ``--target``; this prevents accidentally applying hypotheses intended for another target. No configuration is required. Any target-only run works when the gene is present in the dependency, expression and copy-number matrices.

The pipeline validates the pre-specified biomarkers, performs the genome-wide screens, and excludes automatic candidates that are highly correlated with target features, configured biomarkers or previously selected candidates. (Without a YAML file the tool remains fully functional can select up to two significant, nonredundant mutation and two copy-number biomarkers automatically).

## What the workflow does

1. Validates gene availability, default omics entries, model identifiers, coverage, missingness and variance.
2. Describes the target dependency distribution and lineage specificity.
3. Tests target expression and target copy number, then screens genome-wide damaging mutations and copy-number features.
4. Fits an ordinary least-squares linear regression containing lineage, target expression, target copy number and a small nonredundant biomarker set.
5. Keeps those fitted coefficients and predictions unchanged, but uses HC3-adjusted standard errors for confidence intervals and P values (and therefore FDR). Fixed-feature 10-fold cross-validation estimates held-out performance.
6. Generates a self-contained HTML report.

### Analysis rules

- Eligible biomarkers supplied in the optional YAML file are considered before automatic discoveries.
- A mutation association is tested only when at least 10 models carry a damaging mutation in that candidate gene and at least 10 have no damaging-mutation call for it.
- Copy-number features require at least 80% coverage among copy-number-assayed models and nonzero variance.
- Automatic features require a BH FDR below 0.05 within their mutation or copy-number screen.
- At most two automatic features per genomic modality enter the model.
- If an automatic candidate has absolute correlation of at least 0.8 with a retained predictor, it stays in the univariate tables but is excluded from the multivariable linear regression to avoid unstable coefficients.
- Numeric predictors are standardised so their model coefficients are easier to compare.
- The ordinary linear-regression coefficients are reported with HC3-adjusted uncertainty.
- Fixed-feature 10-fold cross-validation reports held-out prediction error (RMSE) and explained variation (R²).

#### Why use HC3-adjusted uncertainty?

Ordinary linear-regression confidence intervals are most reliable when prediction errors have similar variability across cell lines. That assumption may not be realistic in DepMap because:

- Some cancer lineages are more variable than others.
- Mutation groups can be small and unequal.
- Some cell lines have unusual combinations of expression, copy number and mutations.
- The regression may predict some groups more accurately than others.

HC3 does not refit the regression or change its coefficients and predictions. It recalculates the standard errors, confidence intervals and P values so that uncertainty is more cautious when unusual cell lines or uneven prediction errors are present. Model FDR values are calculated from these HC3-adjusted P values.

> **In plain language:** HC3 asks the model to be more cautious when a small number of unusual cell lines could be driving an apparently precise result.

Chronos scores become more negative as knockout has a stronger effect on cell fitness. See Dempster *et al.*, [Chronos: a cell population dynamics model of CRISPR experiments](https://doi.org/10.1186/s13059-021-02540-7).

## Repository structure

```text
.
├── R/util.R
├── config/MDM2_biomarkers.yml
├── data/README.md
├── reports/
│   ├── target_assessment.Rmd
│   └── MDM2_target_assessment.html
├── results/<TARGET>/
│   ├── figures/
│   └── tables/
├── scripts/
│   ├── 01.prepare.depmap.data.R
│   ├── 02.dependency.landscape.R
│   ├── 03.molecular.associations.R
│   ├── 04.multivariable.model.R
│   └── 05.generate.figures.R
├── tests/run_tests.R
├── environment.yml
└── run_analysis.R
```

Each target creates its own output files, repeated runs don't overwrite analyses of other genes. Full association objects and model-level data are reproducible intermediates and remain Git-ignored; only compact top-ranked tables, figures and reports are intended for version control.

## Worked example and interpretation

The checked-in MDM2 report is generated from **DepMap Public 26Q1** using the optional TP53 biomarker configuration. Target labels, columns, paths, feature screens and reporting statements are generated automatically at runtime.

Please note: Associations should be treated as hypothesis-generating. A feature can remain associated after adjustment without being causal, and a genetic dependency is not equivalent to a druggable or safe therapeutic target.

## Important limitations

- DepMap associations are observational and can reflect lineage composition or unmeasured confounding.
- Copy-number features from the same chromosomal event can be highly correlated; one representative may be selected for modelling.
- Discovery and model evaluation use the same DepMap release. HC3-adjusted intervals and fixed-feature cross-validation do not provide independent validation.
- The multivariable linear regression uses complete cases, which may change cohort composition.
- Cell-line knockout phenotypes may not reproduce in tumours, normal tissues or pharmacological experiments.
- Damaging-mutation calls use the supplied binary/ordinal DepMap matrix and do not model allele-specific function.

## Data source

The example uses [DepMap Public 26Q1](https://depmap.org/portal/data_page/?tab=currentRelease). Downloaded datasets retain their own data-use terms and requested citations. Raw and model-level data are intentionally excluded from this MIT-licensed software repository.
