# Build a compact model-level table for the requested target. Wide genomic
# matrices are screened later and are never written as model-level CSV files.

dependency <- read_gene_data(input_files[["dependency"]], target)
setnames(dependency, target, "dependency")

expression <- read_gene_data(input_files[["expression"]], target, default_only = TRUE)
setnames(expression, target, "target_expression")

copy_number <- read_gene_data(input_files[["copy_number"]], target, default_only = TRUE)
setnames(copy_number, target, "target_copy_number")

model_columns <- c(
  "ModelID", "StrippedCellLineName", "OncotreeLineage",
  "OncotreePrimaryDisease", "OncotreeSubtype"
)
available_model_columns <- intersect(
  model_columns,
  names(fread(input_files[["model"]], nrows = 0L, check.names = FALSE))
)
annotation <- fread(input_files[["model"]], select = available_model_columns)
annotation <- unique(annotation, by = "ModelID")
if (!"OncotreeLineage" %in% names(annotation)) stop("Model.csv lacks OncotreeLineage.")
setnames(annotation, c("StrippedCellLineName", "OncotreeLineage"), c("cell_line", "lineage"), skip_absent = TRUE)

analysis_data <- merge(dependency, annotation, by = "ModelID", all.x = TRUE)
analysis_data <- merge(analysis_data, expression, by = "ModelID", all.x = TRUE)
analysis_data <- merge(analysis_data, copy_number, by = "ModelID", all.x = TRUE)

if (anyDuplicated(analysis_data$ModelID)) stop("ModelID is not unique after modality joins.")
if (sum(!is.na(analysis_data$dependency)) < 50L) stop("Fewer than 50 dependency observations are available for ", target, ".")
if (sum(!is.na(analysis_data$lineage)) < 50L) stop("Insufficient lineage annotation for ", target, ".")

coverage <- data.frame(
  modality = c("dependency", "lineage", "expression", "copy_number"),
  observed = c(
    sum(!is.na(analysis_data$dependency)), sum(!is.na(analysis_data$lineage)),
    sum(!is.na(analysis_data$target_expression)), sum(!is.na(analysis_data$target_copy_number))
  ),
  total_dependency_models = nrow(analysis_data)
)
coverage$proportion <- coverage$observed / coverage$total_dependency_models

quality_checks <- data.frame(
  check = c(
    "Target present in CRISPR matrix", "Unique ModelID after joins",
    "At least 50 dependency observations", "At least 50 lineage annotations",
    "Expression has nonzero variance", "Copy number has nonzero variance"
  ),
  passed = c(
    TRUE, !anyDuplicated(analysis_data$ModelID),
    coverage$observed[coverage$modality == "dependency"] >= 50L,
    coverage$observed[coverage$modality == "lineage"] >= 50L,
    sd(analysis_data$target_expression, na.rm = TRUE) > 0,
    sd(analysis_data$target_copy_number, na.rm = TRUE) > 0
  ),
  detail = c(
    target, as.character(nrow(analysis_data)),
    as.character(coverage$observed[coverage$modality == "dependency"]),
    as.character(coverage$observed[coverage$modality == "lineage"]),
    sprintf("SD %.3f", sd(analysis_data$target_expression, na.rm = TRUE)),
    sprintf("SD %.3f", sd(analysis_data$target_copy_number, na.rm = TRUE))
  )
)

fwrite(coverage, file.path(table_dir, "01_modality_coverage.csv"))
fwrite(quality_checks, file.path(table_dir, "02_data_quality_checks.csv"))
saveRDS(analysis_data, file.path(intermediate_dir, "analysis_data.rds"))
session_information <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_information, file.path(table_dir, "03_session_info.txt"))

message("Prepared ", nrow(analysis_data), " dependency models.")
