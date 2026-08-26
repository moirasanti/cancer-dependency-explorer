## DESCRIPTION ################################################################
## PREPARE DEPMAP DATA
# This stage builds the model-level dataset used throughout the target
# assessment. It reads the selected target's CRISPR dependency, expression and
# copy-number measurements, then joins them to broad OncoTree lineage
# annotations using the stable DepMap ModelID.
#
# Expected objects from the interactive setup or run_analysis.R:
#   - target: requested HGNC gene symbol.
#   - input.files: named paths to the dependency, expression, copy-number and
#     model-annotation files. These are either provided by the user or derived
#     from the DEPMAP_RAW_DIR environment variable.
#   - table.dir and intermediate.dir: target-specific output directories.
#
# The script runs some sanity checks before analysis.
# Saves compact coverage and quality-control tables, records session
# information, and saves analysis.data as an (.git)ignored intermediate
# for the subsequent numbered scripts.

## LOAD DATA ###################################################################
# Read the target dependency, expression and copy-number measurements.
dependency.data <- read_depmap_gene_data(input.files[["dependency"]], target)
setnames(dependency.data, target, "dependency")

expression.data <- read_depmap_gene_data(
  input.files[["expression"]], target, default.only = TRUE
)
setnames(expression.data, target, "target_expression")

copy.number.data <- read_depmap_gene_data(
  input.files[["copy_number"]], target, default.only = TRUE
)
setnames(copy.number.data, target, "target_copy_number")

# Read the model annotations needed for interpretation and lineage adjustment.
model.columns <- c(
  "ModelID",
  "StrippedCellLineName",
  "OncotreeLineage"
)
available.model.columns <- intersect(
  model.columns,
  names(fread(input.files[["model"]], nrows = 0L, check.names = FALSE))
)
annotation.data <- fread(input.files[["model"]], select = available.model.columns)
annotation.data <- unique(annotation.data, by = "ModelID")
if (!"OncotreeLineage" %in% names(annotation.data)) stop("Model.csv lacks OncotreeLineage.")
setnames(
  annotation.data,
  c("StrippedCellLineName", "OncotreeLineage"),
  c("cell_line", "lineage"),
  skip_absent = TRUE
)

# Join modalities to the dependency cohort using stable DepMap model IDs.
analysis.data <- merge(dependency.data, annotation.data, by = "ModelID", all.x = TRUE)
analysis.data <- merge(analysis.data, expression.data, by = "ModelID", all.x = TRUE)
analysis.data <- merge(analysis.data, copy.number.data, by = "ModelID", all.x = TRUE)

## Sanity check:
# Stop before analysis when identifiers or core modality coverage are inadequate.
if (anyDuplicated(analysis.data$ModelID)) stop("ModelID is not unique after modality joins.")
if (sum(!is.na(analysis.data$dependency)) < 50L) {
  stop("Fewer than 50 dependency observations are available for ", target, ".")
}
if (sum(!is.na(analysis.data$lineage)) < 50L) {
  stop("Insufficient lineage annotation for ", target, ".")
}

## MAIN ANALYSIS ###############################################################
## Summarise continuous associations between target dependency and its own expression or copy number.
message("Screening continuous associations...")
coverage.summary <- data.frame(
  modality = c("dependency", "lineage", "expression", "copy_number"),
  observed = c(
    sum(!is.na(analysis.data$dependency)),
    sum(!is.na(analysis.data$lineage)),
    sum(!is.na(analysis.data$target_expression)),
    sum(!is.na(analysis.data$target_copy_number))
  ),
  total_dependency_models = nrow(analysis.data)
)
coverage.summary$proportion <- coverage.summary$observed /
  coverage.summary$total_dependency_models

quality.checks <- data.frame(
  check = c(
    "Target in CRISPR matrix", "Unique ModelID after multiomic joins",
    "At least 50 dependency observations", "At least 50 models with lineage annotation",
    "Expression has nonzero variance", "Copy number has nonzero variance"
  ),
  passed = c(
    TRUE,
    !anyDuplicated(analysis.data$ModelID),
    coverage.summary$observed[coverage.summary$modality == "dependency"] >= 50L,
    coverage.summary$observed[coverage.summary$modality == "lineage"] >= 50L,
    sd(analysis.data$target_expression, na.rm = TRUE) > 0,
    sd(analysis.data$target_copy_number, na.rm = TRUE) > 0
  ),
  detail = c(
    target,
    as.character(nrow(analysis.data)),
    as.character(coverage.summary$observed[coverage.summary$modality == "dependency"]),
    as.character(coverage.summary$observed[coverage.summary$modality == "lineage"]),
    sprintf("SD %.3f", sd(analysis.data$target_expression, na.rm = TRUE)),
    sprintf("SD %.3f", sd(analysis.data$target_copy_number, na.rm = TRUE))
  )
)

# Write compact summaries and retain the model-level table only as an ignored intermediate.
fwrite(coverage.summary, file.path(table.dir, "01_modality_coverage.csv"))
fwrite(quality.checks, file.path(table.dir, "02_data_quality_checks.csv"))
saveRDS(analysis.data, file.path(intermediate.dir, "analysis_data.rds"))

session.information <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session.information, file.path(table.dir, "03_session_info.txt"))

message("Prepared ", nrow(analysis.data), " dependency models.")
