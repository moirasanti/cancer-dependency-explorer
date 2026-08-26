## DESCRIPTION ################################################################
## SELECT FEATURES AND FIT THE MULTIVARIABLE MODEL
# This stage determines which associations remain informative when considered
# together. It combines the target's expression and copy number with broad
# lineage and a small set of genomic biomarkers selected from the univariate
# screens in script 03.
#
# Expected objects from the interactive setup or run_analysis.R:
#   - target and biomarker.file: the requested gene and optional YAML priors.
#   - input.files: named paths used to retrieve candidate genomic features.
#   - table.dir and intermediate.dir: target-specific output directories.
#   - analysis and association RDS files produced by scripts 01 and 03.
#
# See README.md under "Analysis rules" for selection and validation details.


# Load the prepared cohort and complete genome-wide association results.
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))
mutation.results <- readRDS(file.path(intermediate.dir, "mutation_associations.rds"))
copy.number.results <- readRDS(file.path(intermediate.dir, "copy_number_associations.rds"))

# Read optional prior hypotheses from YAML, or use empty lists when none are supplied.
configured.biomarkers <- list(mutation = character(), copy_number = character())
if (nzchar(biomarker.file)) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for biomarker configuration.")
  }
  biomarker.config <- yaml::read_yaml(biomarker.file)
  if (!is.null(biomarker.config$target) && toupper(biomarker.config$target) != target) {
    stop("Biomarker configuration target is ", biomarker.config$target, ", not ", target, ".")
  }
  if (!is.null(biomarker.config$mutation)) {
    configured.biomarkers$mutation <- unique(toupper(unlist(biomarker.config$mutation)))
  }
  if (!is.null(biomarker.config$copy_number)) {
    configured.biomarkers$copy_number <- unique(toupper(unlist(biomarker.config$copy_number)))
  }
}

# Validate configured mutation biomarkers against matrix presence and group sizes.
configured.mutation <- configured.biomarkers$mutation
if (length(configured.mutation)) {
  absent.mutation <- setdiff(configured.mutation, mutation.results$gene)
  if (length(absent.mutation)) {
    stop("Configured mutation genes absent from the matrix: ", paste(absent.mutation, collapse = ", "))
  }
  ineligible.mutation <- configured.mutation[
    !mutation.results$eligible[match(configured.mutation, mutation.results$gene)]
  ]
  if (length(ineligible.mutation)) {
    warning(
      "Configured mutation biomarkers failed coverage/group filters and were omitted: ",
      paste(ineligible.mutation, collapse = ", ")
    )
  }
  configured.mutation <- setdiff(configured.mutation, ineligible.mutation)
}

# Validate configured copy-number biomarkers against matrix presence and coverage.
configured.copy.number <- configured.biomarkers$copy_number
if (length(configured.copy.number)) {
  absent.copy.number <- setdiff(configured.copy.number, copy.number.results$gene)
  if (length(absent.copy.number)) {
    stop(
      "Configured copy-number genes absent from the matrix: ",
      paste(absent.copy.number, collapse = ", ")
    )
  }
  ineligible.copy.number <- configured.copy.number[
    !copy.number.results$eligible[match(configured.copy.number, copy.number.results$gene)]
  ]
  if (length(ineligible.copy.number)) {
    warning(
      "Configured copy-number biomarkers failed coverage/group filters and were omitted: ",
      paste(ineligible.copy.number, collapse = ", ")
    )
  }
  configured.copy.number <- setdiff(configured.copy.number, ineligible.copy.number)
}

# Define ranked automatic pools after FDR filtering.
automatic.mutation.pool <- mutation.results$gene[
  mutation.results$eligible & !is.na(mutation.results$fdr) & mutation.results$fdr < 0.05
]
automatic.mutation.pool <- head(
  setdiff(automatic.mutation.pool, configured.mutation),
  100L
)

automatic.copy.number.pool <- copy.number.results$gene[
  copy.number.results$eligible &
    !is.na(copy.number.results$fdr) &
    copy.number.results$fdr < 0.05
]
automatic.copy.number.pool <- head(
  setdiff(automatic.copy.number.pool, c(target, configured.copy.number)),
  100L
)

# Read only the configured and candidate genomic features needed for selection.
mutation.genes.to.read <- unique(c(configured.mutation, automatic.mutation.pool))
copy.number.genes.to.read <- unique(c(configured.copy.number, automatic.copy.number.pool))

if (length(mutation.genes.to.read)) {
  mutation.values <- as.data.frame(
    read_depmap_gene_data(
      input.files[["mutation"]],
      mutation.genes.to.read,
      default.only = TRUE,
      modality = "mutation"
    ),
    check.names = FALSE
  )
} else {
  mutation.values <- data.frame(ModelID = analysis.data$ModelID)
}

if (length(copy.number.genes.to.read)) {
  copy.number.values <- as.data.frame(
    read_depmap_gene_data(
      input.files[["copy_number"]],
      copy.number.genes.to.read,
      default.only = TRUE,
      modality = "continuous"
    ),
    check.names = FALSE
  )
} else {
  copy.number.values <- data.frame(ModelID = analysis.data$ModelID)
}

mutation.values <- mutation.values[
  match(analysis.data$ModelID, mutation.values$ModelID),
  , drop = FALSE
]
copy.number.values <- copy.number.values[
  match(analysis.data$ModelID, copy.number.values$ModelID),
  , drop = FALSE
]

# Start with target expression and target copy number as retained predictors.
retained.vectors <- list(
  target.expression = analysis.data$target_expression,
  target.copy.number = analysis.data$target_copy_number
)
selection.results <- data.frame(
  modality = character(),
  gene = character(),
  source = character(),
  selected = logical(),
  maximum_absolute_correlation = numeric(),
  reason = character(),
  stringsAsFactors = FALSE
)

# Retain configured biomarkers before evaluating automatic discoveries.
configured.candidates <- rbind(
  data.frame(
    modality = rep("mutation", length(configured.mutation)),
    gene = configured.mutation
  ),
  data.frame(
    modality = rep("copy_number", length(configured.copy.number)),
    gene = configured.copy.number
  )
)
if (nrow(configured.candidates)) {
  for (candidate.index in seq_len(nrow(configured.candidates))) {
    candidate.modality <- configured.candidates$modality[[candidate.index]]
    candidate.gene <- configured.candidates$gene[[candidate.index]]
    if (candidate.modality == "mutation") {
      candidate.vector <- as.integer(mutation.values[[candidate.gene]] > 0)
    } else {
      candidate.vector <- copy.number.values[[candidate.gene]]
    }
    maximum.correlation <- calculate_maximum_absolute_correlation(
      candidate.vector,
      retained.vectors
    )
    selection.results <- rbind(selection.results, data.frame(
      modality = candidate.modality,
      gene = candidate.gene,
      source = "configured",
      selected = TRUE,
      maximum_absolute_correlation = maximum.correlation,
      reason = "Retained"
    ))
    retained.vectors[[paste(candidate.modality, candidate.gene, sep = "__")]] <-
      candidate.vector
  }
}

# Select up to two nonredundant automatic mutation predictors.
selected.automatic.mutation <- 0L
for (candidate.gene in automatic.mutation.pool) {
  if (selected.automatic.mutation >= 2L) break

  candidate.vector <- as.integer(mutation.values[[candidate.gene]] > 0)
  maximum.correlation <- calculate_maximum_absolute_correlation(
    candidate.vector,
    retained.vectors
  )
  candidate.redundant <-
    !is.na(maximum.correlation) && maximum.correlation >= 0.80
  candidate.selected <- !candidate.redundant
  candidate.reason <- if (candidate.redundant) {
    "Excluded: |correlation| >= 0.80 with a retained predictor"
  } else {
    "Retained"
  }

  selection.results <- rbind(selection.results, data.frame(
    modality = "mutation",
    gene = candidate.gene,
    source = "automatic",
    selected = candidate.selected,
    maximum_absolute_correlation = maximum.correlation,
    reason = candidate.reason
  ))
  if (candidate.selected) {
    retained.vectors[[paste("mutation", candidate.gene, sep = "__")]] <- candidate.vector
    selected.automatic.mutation <- selected.automatic.mutation + 1L
  }
}

# Select up to two nonredundant automatic copy-number predictors.
selected.automatic.copy.number <- 0L
for (candidate.gene in automatic.copy.number.pool) {
  if (selected.automatic.copy.number >= 2L) break

  candidate.vector <- copy.number.values[[candidate.gene]]
  maximum.correlation <- calculate_maximum_absolute_correlation(
    candidate.vector,
    retained.vectors
  )
  candidate.redundant <-
    !is.na(maximum.correlation) && maximum.correlation >= 0.80
  candidate.selected <- !candidate.redundant
  candidate.reason <- if (candidate.redundant) {
    "Excluded: |correlation| >= 0.80 with a retained predictor"
  } else {
    "Retained"
  }

  selection.results <- rbind(selection.results, data.frame(
    modality = "copy_number",
    gene = candidate.gene,
    source = "automatic",
    selected = candidate.selected,
    maximum_absolute_correlation = maximum.correlation,
    reason = candidate.reason
  ))
  if (candidate.selected) {
    retained.vectors[[paste("copy_number", candidate.gene, sep = "__")]] <- candidate.vector
    selected.automatic.copy.number <- selected.automatic.copy.number + 1L
  }
}

# Construct readable model terms from the selected genomic biomarkers.
selected.biomarkers <- selection.results[selection.results$selected, , drop = FALSE]
model.data <- data.frame(
  dependency = analysis.data$dependency,
  z_target_expression = as.numeric(scale(analysis.data$target_expression)),
  z_target_copy_number = as.numeric(scale(analysis.data$target_copy_number)),
  lineage = analysis.data$lineage,
  stringsAsFactors = FALSE
)

predictor.labels <- c(
  z_target_expression = paste0(target, " expression (for a 1 SD increase)"),
  z_target_copy_number = paste0(target, " copy number (for a 1 SD increase)")
)
selected.terms <- character()
if (nrow(selected.biomarkers)) {
  for (biomarker.index in seq_len(nrow(selected.biomarkers))) {
    biomarker.modality <- selected.biomarkers$modality[[biomarker.index]]
    biomarker.gene <- selected.biomarkers$gene[[biomarker.index]]
    raw.vector <- retained.vectors[[paste(biomarker.modality, biomarker.gene, sep = "__")]]

    if (biomarker.modality == "mutation") {
      model.term <- make.names(paste0("mutation__", biomarker.gene))
      model.data[[model.term]] <- raw.vector
      predictor.labels[[model.term]] <- paste0(biomarker.gene, " damaging mutation")
    } else {
      model.term <- make.names(paste0("z_copy_number__", biomarker.gene))
      model.data[[model.term]] <- as.numeric(scale(raw.vector))
      predictor.labels[[model.term]] <- paste0(
        biomarker.gene,
        " copy number (for a 1 SD increase)"
      )
    }
    selected.terms <- c(selected.terms, model.term)
  }
}

# Collapse sparsely represented lineages and choose the largest group as reference.
lineage.counts <- sort(table(model.data$lineage), decreasing = TRUE)
retained.lineages <- names(lineage.counts[lineage.counts >= 20L])
model.data$lineage_model <- ifelse(
  model.data$lineage %in% retained.lineages,
  model.data$lineage,
  "Other"
)
reference.lineage <- names(sort(table(model.data$lineage_model), decreasing = TRUE))[[1L]]
model.data$lineage_model <- relevel(
  factor(model.data$lineage_model),
  ref = reference.lineage
)

# Fit the complete-case lineage-adjusted linear model.
model.terms <- c(
  "z_target_expression",
  "z_target_copy_number",
  selected.terms,
  "lineage_model"
)
model.formula <- reformulate(model.terms, response = "dependency")
complete.rows <- complete.cases(model.data[, c("dependency", model.terms), drop = FALSE])
complete.model.data <- model.data[complete.rows, , drop = FALSE]
if (nrow(complete.model.data) < 100L) {
  stop("Fewer than 100 complete observations are available for the multivariable model.")
}
dependency.model <- lm(model.formula, data = complete.model.data)

# Calculate HC3 robust standard errors and confidence intervals explicitly.
model.coefficients <- coef(dependency.model)
retained.coefficients <- !is.na(model.coefficients)
model.matrix.data <- model.matrix(dependency.model)[, retained.coefficients, drop = FALSE]
cross.product <- crossprod(model.matrix.data)
bread.matrix <- try(solve(cross.product), silent = TRUE)
if (inherits(bread.matrix, "try-error")) bread.matrix <- qr.solve(cross.product)

adjusted.residuals <- residuals(dependency.model) /
  pmax(1 - hatvalues(dependency.model), 1e-8)
meat.matrix <- crossprod(model.matrix.data * adjusted.residuals)
robust.covariance <- bread.matrix %*% meat.matrix %*% bread.matrix
robust.standard.error <- sqrt(diag(robust.covariance))
robust.estimate <- model.coefficients[retained.coefficients]
robust.statistic <- robust.estimate / robust.standard.error
residual.degrees.freedom <- df.residual(dependency.model)
robust.p.value <- 2 * pt(
  abs(robust.statistic),
  df = residual.degrees.freedom,
  lower.tail = FALSE
)
confidence.critical <- qt(0.975, df = residual.degrees.freedom)

robust.coefficients <- data.frame(
  term = names(robust.estimate),
  estimate = unname(robust.estimate),
  std_error = robust.standard.error,
  statistic = robust.statistic,
  p_value = robust.p.value,
  conf_low = robust.estimate - confidence.critical * robust.standard.error,
  conf_high = robust.estimate + confidence.critical * robust.standard.error,
  stringsAsFactors = FALSE
)

# Replace internal term names with labels suitable for tables and figures.
robust.coefficients$label <- robust.coefficients$term
for (model.term in names(predictor.labels)) {
  robust.coefficients$label[robust.coefficients$term == model.term] <-
    predictor.labels[[model.term]]
}
lineage.coefficients <- grepl("^lineage_model", robust.coefficients$term)
robust.coefficients$label[lineage.coefficients] <- paste0(
  "Lineage: ",
  sub("^lineage_model", "", robust.coefficients$term[lineage.coefficients]),
  " vs ",
  reference.lineage
)
robust.coefficients$fdr <- p.adjust(robust.coefficients$p_value, method = "BH")

# Estimate fixed-feature performance with deterministic 10-fold cross-validation.
set.seed(20260825L)
cross.validation.folds <- 10L
fold.id <- sample(rep(seq_len(cross.validation.folds), length.out = nrow(complete.model.data)))
observed.values <- predicted.values <- rep(NA_real_, nrow(complete.model.data))
for (fold in seq_len(cross.validation.folds)) {
  training.data <- complete.model.data[fold.id != fold, , drop = FALSE]
  test.data <- complete.model.data[fold.id == fold, , drop = FALSE]
  fold.model <- lm(model.formula, data = training.data)
  predicted.values[fold.id == fold] <- suppressWarnings(
    predict(fold.model, newdata = test.data)
  )
  observed.values[fold.id == fold] <- test.data$dependency
}
prediction.complete <- complete.cases(observed.values, predicted.values)
cross.validation.rmse <- sqrt(mean(
  (observed.values[prediction.complete] - predicted.values[prediction.complete])^2
))
cross.validation.r.squared <- 1 -
  sum((observed.values[prediction.complete] - predicted.values[prediction.complete])^2) /
  sum(
    (observed.values[prediction.complete] - mean(observed.values[prediction.complete]))^2
  )

# Summarise model fit, fixed-feature performance and basic diagnostics.
model.performance <- data.frame(
  n = nobs(dependency.model),
  predictors = dependency.model$rank - 1L,
  adjusted_r_squared = summary(dependency.model)$adj.r.squared,
  residual_sd = summary(dependency.model)$sigma,
  cv_folds = cross.validation.folds,
  cv_rmse = cross.validation.rmse,
  cv_r_squared = cross.validation.r.squared,
  reference_lineage = reference.lineage
)
model.diagnostics <- data.frame(
  metric = c(
    "Maximum leverage",
    "Condition number",
    "Maximum absolute studentized residual"
  ),
  value = c(
    max(hatvalues(dependency.model)),
    kappa(model.matrix(dependency.model)),
    max(abs(rstudent(dependency.model)), na.rm = TRUE)
  )
)

# Write compact results and retain the fitted model as an ignored intermediate.
fwrite(selection.results, file.path(table.dir, "10_biomarker_selection.csv"))
fwrite(
  robust.coefficients,
  file.path(table.dir, "11_multivariable_model_coefficients.csv")
)
fwrite(model.performance, file.path(table.dir, "12_model_performance.csv"))
fwrite(model.diagnostics, file.path(table.dir, "13_model_diagnostics.csv"))
saveRDS(
  list(
    model = dependency.model,
    coefficients = robust.coefficients,
    data = complete.model.data,
    labels = predictor.labels,
    selection = selection.results
  ),
  file.path(intermediate.dir, "multivariable_model.rds")
)

message(
  "Fitted model using ", nobs(dependency.model),
  " complete observations and ", dependency.model$rank - 1L,
  " coefficients."
)
