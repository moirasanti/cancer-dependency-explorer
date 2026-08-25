# Cancer Dependency Explorer

**Which molecular features predict dependency on a selected cancer target?**

Cancer Dependency Explorer is a compact, reproducible R workflow that integrates CRISPR gene-effect scores with cancer lineage, expression, copy number and damaging somatic mutations. Give it a target gene and it produces the same quality checks, association screens, interpretable multivariable model, figures and HTML assessment report.

The included worked example asks what predicts **MDM2 dependency**. In DepMap Public 26Q1, TP53 mutation status, MDM2 expression and genomic context provide an interpretable translational example while also demonstrating why lineage adjustment and correlated-feature filtering matter.

> This repository asks *what predicts dependency on target X?* It does not perform the inverse genome-wide screen *which dependencies are selective in models with mutation of tumour suppressor X?*

## Run a target assessment

Create the environment:

```bash
conda env create -f environment.yml
conda activate cancer-dependency-explorer
```

Download the files listed in [data/README.md](data/README.md), then run:

```bash
Rscript run_analysis.R \
  --target MDM2 \
  --raw-dir /path/to/depmap/26Q1
```

`DEPMAP_RAW_DIR` can replace `--raw-dir`. An optional YAML file can force biological hypotheses into the model before automatic discoveries are considered:

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

No configuration is required. Any target-only run works when the gene is present in the dependency, expression and copy-number matrices.

## What the workflow does

1. Validates gene availability, default omics entries, model identifiers, coverage, missingness and variance.
2. Describes the target dependency distribution and lineage specificity.
3. Tests target expression and target copy number, then screens genome-wide damaging mutations and copy-number features.
4. Fits a lineage-adjusted linear model containing target expression, target copy number and a small nonredundant biomarker set.
5. Reports HC3 robust confidence intervals and fixed-feature 10-fold cross-validation.
6. Generates a self-contained HTML report with a short interpretation after every major result.

Automatic mutation features require at least 10 altered and 10 reference models. Copy-number features require at least 80% coverage among copy-number-assayed models and nonzero variance. Tests are corrected within modality using Benjamini–Hochberg FDR. At most two automatic features per genomic modality enter the model. If an automatic candidate has absolute correlation of at least 0.8 with a retained predictor, it stays in the univariate tables but is excluded from the multivariable model to avoid unstable coefficients.

Chronos scores become more negative as knockout has a stronger effect on cell fitness. See Dempster *et al.*, [Chronos: a cell population dynamics model of CRISPR experiments](https://doi.org/10.1186/s13059-021-02540-7).

## Repository structure

```text
.
├── R/utils.R
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

Each target writes to its own output directory, so repeated runs do not overwrite analyses of other genes. Full association objects and model-level data are reproducible intermediates and remain Git-ignored; only compact top-ranked tables, figures and reports are intended for version control.

## Worked example and interpretation

The checked-in MDM2 report is generated from **DepMap Public 26Q1** using the optional TP53 biomarker configuration. It is an example of the tool, not a special code path: target labels, columns, paths, feature screens and narrative sentences are generated at runtime.

Associations should be treated as hypothesis-generating. A feature can remain associated after adjustment without being causal, and a genetic dependency is not equivalent to a druggable or safe therapeutic target.

## Important limitations

- DepMap associations are observational and can reflect lineage composition or unmeasured confounding.
- Copy-number features from the same chromosomal event can be highly correlated; one representative may be selected for modelling.
- Discovery and model evaluation use the same DepMap release. Robust intervals and fixed-feature cross-validation do not provide independent validation.
- The multivariable model uses complete cases, which may change cohort composition.
- Cell-line knockout phenotypes may not reproduce in tumours, normal tissues or pharmacological experiments.
- Damaging-mutation calls use the supplied binary/ordinal DepMap matrix and do not model allele-specific function.

## Data source

The example uses [DepMap Public 26Q1](https://depmap.org/portal/data_page/?tab=currentRelease). Downloaded datasets retain their own data-use terms and requested citations. Raw and model-level data are intentionally excluded from this MIT-licensed software repository.

## Tests

```bash
Rscript tests/run_tests.R
```

The tests cover gene parsing, input validation, minimum mutation groups, copy-number coverage, redundant-feature detection and deterministic cross-validation.
