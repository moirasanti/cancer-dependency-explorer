# DepMap input data

Download these files from the [DepMap data portal](https://depmap.org/portal/data_page/?tab=allData) for the same release and place them in one directory:

```text
CRISPRGeneEffect.csv
Model.csv
OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
OmicsCNGeneWGS.csv
OmicsSomaticMutationsMatrixDamaging.csv
```

The worked example is pinned to **DepMap Public 26Q1**. DepMap updates its pipelines and model collection between releases, so results from another release may differ even when filenames are compatible.

Point the workflow to the download directory without copying the matrices into the repository:

```bash
DEPMAP_RAW_DIR=/path/to/depmap/26Q1 \
  Rscript run_analysis.R --target MDM2
```

Raw matrices and model-level intermediates are excluded from Git because they are large and retain the terms attached to their DepMap downloads. Before redistributing any data-derived artifact, review the licence and requested citation displayed for every contributing file.
