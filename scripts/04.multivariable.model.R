# Select a small, nonredundant genomic feature set and estimate an interpretable
# lineage-adjusted model. Optional configured biomarkers are considered first.

analysis_data <- readRDS(file.path(intermediate_dir, "analysis_data.rds"))
mutation_results <- readRDS(file.path(intermediate_dir, "mutation_associations.rds"))
copy_number_results <- readRDS(file.path(intermediate_dir, "copy_number_associations.rds"))
configured <- load_biomarker_config(biomarker_file, target)

eligible_configured <- function(genes, results, modality) {
  if (!length(genes)) return(character())
  absent <- setdiff(genes, results$gene)
  if (length(absent)) stop("Configured ", modality, " genes absent from the matrix: ", paste(absent, collapse = ", "))
  ineligible <- genes[!results$eligible[match(genes, results$gene)]]
  if (length(ineligible)) {
    warning("Configured ", modality, " biomarkers failed coverage/group filters and were omitted: ",
            paste(ineligible, collapse = ", "))
  }
  setdiff(genes, ineligible)
}

curated_mutation <- eligible_configured(configured$mutation, mutation_results, "mutation")
curated_copy_number <- eligible_configured(configured$copy_number, copy_number_results, "copy-number")

automatic_mutation_pool <- mutation_results$gene[
  mutation_results$eligible & !is.na(mutation_results$fdr) & mutation_results$fdr < 0.05
]
automatic_mutation_pool <- head(setdiff(automatic_mutation_pool, curated_mutation), 100L)
automatic_copy_number_pool <- copy_number_results$gene[
  copy_number_results$eligible & !is.na(copy_number_results$fdr) & copy_number_results$fdr < 0.05
]
automatic_copy_number_pool <- head(
  setdiff(automatic_copy_number_pool, c(target, curated_copy_number)), 100L
)

mutation_to_read <- unique(c(curated_mutation, automatic_mutation_pool))
copy_number_to_read <- unique(c(curated_copy_number, automatic_copy_number_pool))

mutation_values <- if (length(mutation_to_read)) {
  read_feature_matrix(input_files[["mutation"]], mutation_to_read, default_only = TRUE)
} else data.frame(ModelID = analysis_data$ModelID)
copy_number_values <- if (length(copy_number_to_read)) {
  read_feature_matrix(input_files[["copy_number"]], copy_number_to_read, default_only = TRUE)
} else data.frame(ModelID = analysis_data$ModelID)

mutation_values <- mutation_values[match(analysis_data$ModelID, mutation_values$ModelID), , drop = FALSE]
copy_number_values <- copy_number_values[match(analysis_data$ModelID, copy_number_values$ModelID), , drop = FALSE]

feature_vector <- function(modality, gene) {
  if (modality == "mutation") as.integer(mutation_values[[gene]] > 0) else copy_number_values[[gene]]
}

retained_vectors <- list(
  target_expression = analysis_data$target_expression,
  target_copy_number = analysis_data$target_copy_number
)
selection <- data.frame(
  modality = character(), gene = character(), source = character(),
  selected = logical(), maximum_absolute_correlation = numeric(), reason = character(),
  stringsAsFactors = FALSE
)

record_candidate <- function(modality, gene, source, enforce_redundancy = TRUE) {
  vector <- feature_vector(modality, gene)
  maximum <- maximum_absolute_correlation(vector, retained_vectors)
  redundant <- enforce_redundancy && !is.na(maximum) && maximum >= 0.80
  selected <- !redundant
  reason <- if (redundant) "Excluded: |correlation| >= 0.80 with a retained predictor" else "Retained"
  selection <<- rbind(selection, data.frame(
    modality = modality, gene = gene, source = source, selected = selected,
    maximum_absolute_correlation = maximum, reason = reason
  ))
  if (selected) retained_vectors[[paste(modality, gene, sep = "__")]] <<- vector
  selected
}

for (gene in curated_mutation) record_candidate("mutation", gene, "configured", FALSE)
for (gene in curated_copy_number) record_candidate("copy_number", gene, "configured", FALSE)

selected_automatic_mutation <- 0L
for (gene in automatic_mutation_pool) {
  if (selected_automatic_mutation >= 2L) break
  if (record_candidate("mutation", gene, "automatic", TRUE)) {
    selected_automatic_mutation <- selected_automatic_mutation + 1L
  }
}
selected_automatic_copy_number <- 0L
for (gene in automatic_copy_number_pool) {
  if (selected_automatic_copy_number >= 2L) break
  if (record_candidate("copy_number", gene, "automatic", TRUE)) {
    selected_automatic_copy_number <- selected_automatic_copy_number + 1L
  }
}

selected <- selection[selection$selected, , drop = FALSE]
model_data <- data.frame(
  dependency = analysis_data$dependency,
  z_target_expression = as.numeric(scale(analysis_data$target_expression)),
  z_target_copy_number = as.numeric(scale(analysis_data$target_copy_number)),
  lineage = analysis_data$lineage,
  stringsAsFactors = FALSE
)

predictor_labels <- c(
  z_target_expression = paste0(target, " expression (per SD)"),
  z_target_copy_number = paste0(target, " copy number (per SD)")
)
selected_terms <- character()
if (nrow(selected)) {
  for (row in seq_len(nrow(selected))) {
    modality <- selected$modality[[row]]
    gene <- selected$gene[[row]]
    raw_vector <- retained_vectors[[paste(modality, gene, sep = "__")]]
    if (modality == "mutation") {
      term <- make.names(paste0("mutation__", gene))
      model_data[[term]] <- raw_vector
      predictor_labels[[term]] <- paste0(gene, " damaging mutation")
    } else {
      term <- make.names(paste0("z_copy_number__", gene))
      model_data[[term]] <- as.numeric(scale(raw_vector))
      predictor_labels[[term]] <- paste0(gene, " copy number (per SD)")
    }
    selected_terms <- c(selected_terms, term)
  }
}

lineage_counts <- sort(table(model_data$lineage), decreasing = TRUE)
retained_lineages <- names(lineage_counts[lineage_counts >= 20L])
model_data$lineage_model <- ifelse(model_data$lineage %in% retained_lineages, model_data$lineage, "Other")
reference_lineage <- names(sort(table(model_data$lineage_model), decreasing = TRUE))[[1L]]
model_data$lineage_model <- relevel(factor(model_data$lineage_model), ref = reference_lineage)

model_terms <- c("z_target_expression", "z_target_copy_number", selected_terms, "lineage_model")
model_formula <- reformulate(model_terms, response = "dependency")
complete <- complete.cases(model_data[, c("dependency", model_terms), drop = FALSE])
complete_model_data <- model_data[complete, , drop = FALSE]
if (nrow(complete_model_data) < 100L) stop("Fewer than 100 complete observations are available for the multivariable model.")

model <- lm(model_formula, data = complete_model_data)
robust_coefficients <- hc3_tidy(model)
robust_coefficients$label <- robust_coefficients$term
for (term in names(predictor_labels)) {
  robust_coefficients$label[robust_coefficients$term == term] <- predictor_labels[[term]]
}
robust_coefficients$label[grepl("^lineage_model", robust_coefficients$term)] <- paste0(
  "Lineage: ", sub("^lineage_model", "", robust_coefficients$term[grepl("^lineage_model", robust_coefficients$term)]),
  " vs ", reference_lineage
)
robust_coefficients$fdr <- p.adjust(robust_coefficients$p_value, method = "BH")

cv_performance <- fixed_feature_cv(complete_model_data, model_formula)
model_performance <- data.frame(
  n = nobs(model), predictors = model$rank - 1L,
  adjusted_r_squared = summary(model)$adj.r.squared,
  residual_sd = summary(model)$sigma,
  cv_folds = cv_performance$folds,
  cv_rmse = cv_performance$rmse,
  cv_r_squared = cv_performance$r_squared,
  reference_lineage = reference_lineage
)
model_diagnostics <- data.frame(
  metric = c("Maximum leverage", "Condition number", "Maximum absolute studentized residual"),
  value = c(max(hatvalues(model)), kappa(model.matrix(model)), max(abs(rstudent(model)), na.rm = TRUE))
)

fwrite(selection, file.path(table_dir, "10_biomarker_selection.csv"))
fwrite(robust_coefficients, file.path(table_dir, "11_multivariable_model_coefficients.csv"))
fwrite(model_performance, file.path(table_dir, "12_model_performance.csv"))
fwrite(model_diagnostics, file.path(table_dir, "13_model_diagnostics.csv"))
saveRDS(list(model = model, coefficients = robust_coefficients, data = complete_model_data,
             labels = predictor_labels, selection = selection),
        file.path(intermediate_dir, "multivariable_model.rds"))

message("Fitted model using ", nobs(model), " complete observations and ", model$rank - 1L, " coefficients.")
