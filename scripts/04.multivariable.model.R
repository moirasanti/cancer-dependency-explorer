## DESCRIPTION ################################################################
## FIT AN EXPLORATORY INTEGRATED ELASTIC-NET MODEL
# This stage combines eligible target features, pre-specified hypotheses,
# lineage and all other eligible genome-wide molecular features. Lineage is the
# only unpenalized adjustment; every molecular feature competes under the same
# elastic-net penalty.
# Repeated nested cross-validation estimates held-out performance and selection
# stability without using held-out outcomes during preprocessing or selection.

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for integrated exploratory modelling.")
}

## LOAD AND VALIDATE ###########################################################
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))
expression.results <- readRDS(file.path(intermediate.dir, "expression_associations.rds"))
mutation.results <- readRDS(file.path(intermediate.dir, "mutation_associations.rds"))
copy.number.results <- readRDS(file.path(intermediate.dir, "copy_number_associations.rds"))

configured.hypotheses <- read_biomarker_configuration(biomarker.file, target)
configured.hypotheses <- deduplicate_target_hypotheses(
  configured.hypotheses, target
)
requested.configured.hypotheses <- configured.hypotheses

configured.expression <- validate_configured_features(
  configured.hypotheses$expression,
  expression.results,
  "expression",
  "coverage/variance"
)
configured.mutation <- validate_configured_features(
  configured.hypotheses$mutation,
  mutation.results,
  "mutation",
  "altered/reference group size"
)
configured.copy.number <- validate_configured_features(
  configured.hypotheses$copy_number,
  copy.number.results,
  "copy number",
  "coverage/variance"
)

target.expression.eligible <- isTRUE(expression.results$eligible[
  match(target, expression.results$gene)
])
target.mutation.eligible <- is_eligible_mutation_feature(mutation.results, target)
target.copy.number.eligible <- isTRUE(copy.number.results$eligible[
  match(target, copy.number.results$gene)
])
if (!target.expression.eligible) {
  warning("Target expression was omitted because basic eligibility failed.")
}
if (!target.mutation.eligible) {
  warning(
    "Target damaging mutation was omitted because ", target,
    " did not have at least 10 altered and 10 wild-type models."
  )
}
if (!target.copy.number.eligible) {
  warning("Target copy number was omitted because basic eligibility failed.")
}

eligible.expression <- expression.results$gene[expression.results$eligible]
eligible.mutation <- mutation.results$gene[mutation.results$eligible]
eligible.copy.number <- copy.number.results$gene[copy.number.results$eligible]

## READ ALL ELIGIBLE FEATURES #################################################
message("Reading eligible genome-wide features for integrated modelling...")
expression.data <- read_aligned_genomic_matrix(
  input.files[["expression"]], analysis.data$ModelID, modality = "expression"
)
expression.matrix <- expression.data$matrix[
  , match(eligible.expression, expression.data$symbols), drop = FALSE
]
colnames(expression.matrix) <- eligible.expression
rm(expression.data)
gc()

mutation.data <- read_aligned_genomic_matrix(
  input.files[["mutation"]], analysis.data$ModelID, modality = "mutation"
)
mutation.matrix <- mutation.data$matrix[
  , match(eligible.mutation, mutation.data$symbols), drop = FALSE
]
mutation.matrix <- ifelse(is.na(mutation.matrix), NA_real_, mutation.matrix > 0)
colnames(mutation.matrix) <- eligible.mutation
rm(mutation.data)
gc()

copy.number.data <- read_aligned_genomic_matrix(
  input.files[["copy_number"]], analysis.data$ModelID, modality = "copy_number"
)
copy.number.matrix <- copy.number.data$matrix[
  , match(eligible.copy.number, copy.number.data$symbols), drop = FALSE
]
colnames(copy.number.matrix) <- eligible.copy.number
rm(copy.number.data)
gc()

## LINEAGE-SPECIFIC FOLLOW-UP #################################################
# This descriptive follow-up is separate from elastic-net selection.
lineage.candidates <- rbind(
  make_feature_family_rows(
    "target_characterization",
    c("expression", "mutation", "copy_number"),
    rep(target, 3L)
  ),
  make_feature_family_rows(
    "prespecified_hypothesis", "expression",
    requested.configured.hypotheses$expression
  ),
  make_feature_family_rows(
    "prespecified_hypothesis", "mutation",
    requested.configured.hypotheses$mutation
  ),
  make_feature_family_rows(
    "prespecified_hypothesis", "copy_number",
    requested.configured.hypotheses$copy_number
  )
)
lineage.specific.results <- vector("list", nrow(lineage.candidates))
for (candidate.index in seq_len(nrow(lineage.candidates))) {
  candidate <- lineage.candidates[candidate.index, ]
  result.table <- switch(
    candidate$modality,
    expression = expression.results,
    mutation = mutation.results,
    copy_number = copy.number.results
  )
  result.index <- match(candidate$gene, result.table$gene)
  candidate.eligible <-
    !is.na(result.index) && isTRUE(result.table$eligible[[result.index]])
  if (!candidate.eligible) next

  candidate.vector <- switch(
    candidate$modality,
    expression = expression.matrix[, candidate$gene],
    mutation = mutation.matrix[, candidate$gene],
    copy_number = copy.number.matrix[, candidate$gene]
  )
  lineage.result <- fit_lineage_specific_associations(
    analysis.data$dependency,
    candidate.vector,
    analysis.data$lineage,
    continuous = candidate$modality != "mutation",
    minimum.lineage.size = 20L,
    minimum.mutation.group.size = 10L
  )
  if (!nrow(lineage.result)) next
  candidate.label <- switch(
    candidate$modality,
    expression = paste(candidate$gene, "expression"),
    mutation = paste(candidate$gene, "damaging mutation"),
    copy_number = paste(candidate$gene, "copy number")
  )
  lineage.specific.results[[candidate.index]] <- cbind(
    data.frame(
      family = candidate$family,
      modality = candidate$modality,
      gene = candidate$gene,
      feature = candidate.label,
      estimate_type = if (candidate$modality == "mutation") {
        "Altered minus reference"
      } else {
        "Chronos change per 1 overall SD increase"
      },
      stringsAsFactors = FALSE
    ),
    lineage.result
  )
}
lineage.specific.associations <- rbindlist(lineage.specific.results, fill = TRUE)
lineage.specific.associations[, fdr := NA_real_]
estimable.lineage.rows <-
  lineage.specific.associations$estimable &
  !is.na(lineage.specific.associations$p_value)
lineage.specific.associations$fdr[estimable.lineage.rows] <- p.adjust(
  lineage.specific.associations$p_value[estimable.lineage.rows], method = "BH"
)
setorder(lineage.specific.associations, fdr, p_value, feature, lineage)

## BUILD THE INTEGRATED PREDICTOR MATRIX #####################################
make_feature_metadata <- function(modality, genes, target.gene, configured) {
  if (!length(genes)) return(data.frame())
  role <- ifelse(
    genes == target.gene,
    "target_characterization",
    ifelse(
      genes %in% configured,
      "prespecified_hypothesis",
      "genome_wide_candidate"
    )
  )
  label <- switch(
    modality,
    expression = paste(genes, "expression"),
    mutation = paste(genes, "damaging mutation"),
    copy_number = paste(genes, "copy number")
  )
  data.frame(
    term = paste(modality, genes, sep = "__"),
    modality = modality,
    gene = genes,
    role = role,
    label = label,
    stringsAsFactors = FALSE
  )
}

expression.metadata <- make_feature_metadata(
  "expression", eligible.expression,
  if (target.expression.eligible) target else "", configured.expression
)
mutation.metadata <- make_feature_metadata(
  "mutation", eligible.mutation,
  if (target.mutation.eligible) target else "", configured.mutation
)
copy.number.metadata <- make_feature_metadata(
  "copy_number", eligible.copy.number,
  if (target.copy.number.eligible) target else "", configured.copy.number
)
molecular.metadata <- rbind(
  expression.metadata, mutation.metadata, copy.number.metadata
)

colnames(expression.matrix) <- expression.metadata$term
colnames(mutation.matrix) <- mutation.metadata$term
colnames(copy.number.matrix) <- copy.number.metadata$term
molecular.matrix <- cbind(expression.matrix, mutation.matrix, copy.number.matrix)
rm(expression.matrix, mutation.matrix, copy.number.matrix)
gc()

continuous.columns <- which(
  molecular.metadata$modality %in% c("expression", "copy_number")
)
mutation.columns <- which(molecular.metadata$modality == "mutation")
molecular.columns <- seq_len(nrow(molecular.metadata))

# Collapse rare lineages using outcome-independent cohort counts.
lineage.counts <- sort(table(analysis.data$lineage), decreasing = TRUE)
retained.lineages <- names(lineage.counts[lineage.counts >= 20L])
lineage.model <- ifelse(
  analysis.data$lineage %in% retained.lineages,
  analysis.data$lineage,
  "Other"
)
reference.lineage <- names(sort(table(lineage.model), decreasing = TRUE))[[1L]]
lineage.model <- relevel(factor(lineage.model), ref = reference.lineage)
lineage.matrix <- model.matrix(~ lineage.model)[, -1L, drop = FALSE]
colnames(lineage.matrix) <- paste0(
  "lineage__", sub("^lineage.model", "", colnames(lineage.matrix))
)
lineage.metadata <- data.frame(
  term = colnames(lineage.matrix),
  modality = "lineage",
  gene = sub("^lineage__", "", colnames(lineage.matrix)),
  role = "lineage_adjustment",
  label = paste0(
    "Lineage: ", sub("^lineage__", "", colnames(lineage.matrix)),
    " vs ", reference.lineage
  ),
  stringsAsFactors = FALSE
)
all.metadata <- rbind(molecular.metadata, lineage.metadata)
penalty.factor <- c(
  rep(1, nrow(molecular.metadata)),
  rep(0, ncol(lineage.matrix))
)

analysis.rows <- which(
  !is.na(analysis.data$dependency) & !is.na(analysis.data$lineage)
)
if (length(analysis.rows) < 100L) {
  stop("Fewer than 100 complete outcome/lineage observations are available.")
}
outcome <- analysis.data$dependency

## REPEATED NESTED CROSS-VALIDATION ###########################################
elastic.alpha <- 0.5
outer.folds <- 5L
inner.folds <- 5L
outer.repeats <- 3L
selection.tolerance <- 1e-10
outer.model.count <- outer.folds * outer.repeats
outer.coefficients <- matrix(
  0,
  nrow = nrow(all.metadata),
  ncol = outer.model.count,
  dimnames = list(all.metadata$term, NULL)
)
nested.performance <- vector("list", outer.repeats)
nested.predictions <- vector("list", outer.repeats)
outer.model.index <- 0L

message(
  "Running ", outer.repeats, " repeats of ", outer.folds,
  "-fold nested elastic-net validation across ",
  length(molecular.columns), " equally penalized molecular features..."
)
for (repeat.index in seq_len(outer.repeats)) {
  repeat.fold.id <- make_stratified_folds(
    analysis.data$lineage[analysis.rows],
    folds = outer.folds,
    seed = 20260825L + repeat.index
  )
  repeat.prediction <- rep(NA_real_, length(analysis.rows))

  for (outer.fold in seq_len(outer.folds)) {
    outer.model.index <- outer.model.index + 1L
    test.position <- which(repeat.fold.id == outer.fold)
    train.position <- which(repeat.fold.id != outer.fold)
    training.rows <- analysis.rows[train.position]
    test.rows <- analysis.rows[test.position]

    preprocessor <- fit_integrated_preprocessor(
      molecular.matrix, training.rows, continuous.columns, mutation.columns
    )
    training.molecular <- apply_integrated_preprocessor(
      molecular.matrix, training.rows, preprocessor
    )
    test.molecular <- apply_integrated_preprocessor(
      molecular.matrix, test.rows, preprocessor
    )
    training.matrix <- cbind(training.molecular, lineage.matrix[training.rows, ])
    test.matrix <- cbind(test.molecular, lineage.matrix[test.rows, ])
    inner.fold.id <- make_stratified_folds(
      analysis.data$lineage[training.rows],
      folds = inner.folds,
      seed = 20260825L + repeat.index * 100L + outer.fold
    )

    fitted.cv <- glmnet::cv.glmnet(
      x = training.matrix,
      y = outcome[training.rows],
      family = "gaussian",
      alpha = elastic.alpha,
      foldid = inner.fold.id,
      penalty.factor = penalty.factor,
      standardize = FALSE,
      type.measure = "mse"
    )
    repeat.prediction[test.position] <- as.numeric(predict(
      fitted.cv, newx = test.matrix, s = "lambda.1se"
    ))
    fold.coefficients <- as.matrix(coef(fitted.cv, s = "lambda.1se"))
    outer.coefficients[, outer.model.index] <-
      fold.coefficients[all.metadata$term, 1L]

    rm(
      training.molecular, test.molecular, training.matrix, test.matrix,
      fitted.cv, fold.coefficients
    )
    gc()
  }

  observed <- outcome[analysis.rows]
  repeat.rmse <- sqrt(mean((observed - repeat.prediction)^2))
  repeat.pearson.r <- cor(
    observed, repeat.prediction,
    use = "complete.obs", method = "pearson"
  )
  repeat.r.squared <- 1 -
    sum((observed - repeat.prediction)^2) /
    sum((observed - mean(observed))^2)
  nested.performance[[repeat.index]] <- data.frame(
    repeat_id = repeat.index,
    n = length(observed),
    outer_folds = outer.folds,
    inner_folds = inner.folds,
    alpha = elastic.alpha,
    lambda_rule = "lambda.1se",
    rmse = repeat.rmse,
    pearson_r = repeat.pearson.r,
    r_squared = repeat.r.squared
  )
  nested.predictions[[repeat.index]] <- data.frame(
    ModelID = analysis.data$ModelID[analysis.rows],
    lineage = as.character(analysis.data$lineage[analysis.rows]),
    repeat_id = repeat.index,
    outer_fold = repeat.fold.id,
    observed = observed,
    predicted = repeat.prediction,
    residual = observed - repeat.prediction,
    stringsAsFactors = FALSE
  )
}
nested.performance <- rbindlist(nested.performance)
nested.predictions <- rbindlist(nested.predictions)

# Summarize genuinely held-out prediction accuracy within sufficiently large
# cancer lineages. R-squared is calculated relative to the observed mean within
# each lineage, so zero means no improvement over that lineage-specific mean.
lineage.performance.by.repeat <- nested.predictions[
  , .(
    n = .N,
    mae = mean(abs(residual)),
    rmse = sqrt(mean(residual^2)),
    r_squared = if (sum((observed - mean(observed))^2) > 0) {
      1 - sum(residual^2) / sum((observed - mean(observed))^2)
    } else {
      NA_real_
    }
  ),
  by = .(repeat_id, lineage)
][n >= 20L]
lineage.performance <- lineage.performance.by.repeat[
  , .(
    n = unique(n)[1L],
    mae_mean = mean(mae),
    mae_sd = sd(mae),
    rmse_mean = mean(rmse),
    rmse_sd = sd(rmse),
    r_squared_mean = mean(r_squared),
    r_squared_sd = sd(r_squared)
  ),
  by = lineage
]
setorder(lineage.performance, -r_squared_mean, lineage)

## MOLECULAR STABILITY AND FULL-DATA FIT ######################################
molecular.outer.coefficients <- outer.coefficients[
  molecular.columns, , drop = FALSE
]
stability.summary <- summarize_elastic_net_stability(
  molecular.outer.coefficients,
  tolerance = selection.tolerance
)
selected.count <- stability.summary$selected_count
selection.frequency <- stability.summary$selection_frequency
direction.concordance <- stability.summary$direction_concordance
median.selected.coefficient <- stability.summary$median_selected_coefficient

# Full-data discovery fit supplies the third stability requirement.
full.preprocessor <- fit_integrated_preprocessor(
  molecular.matrix, analysis.rows, continuous.columns, mutation.columns
)
full.molecular <- apply_integrated_preprocessor(
  molecular.matrix, analysis.rows, full.preprocessor
)
full.matrix <- cbind(full.molecular, lineage.matrix[analysis.rows, ])
full.fold.id <- make_stratified_folds(
  analysis.data$lineage[analysis.rows],
  folds = inner.folds,
  seed = 20260825L
)
full.discovery.cv <- glmnet::cv.glmnet(
  x = full.matrix,
  y = outcome[analysis.rows],
  family = "gaussian",
  alpha = elastic.alpha,
  foldid = full.fold.id,
  penalty.factor = penalty.factor,
  standardize = FALSE,
  type.measure = "mse"
)
full.discovery.coefficients <- as.matrix(coef(
  full.discovery.cv, s = "lambda.1se"
))
full.discovery.coefficient <-
  full.discovery.coefficients[molecular.metadata$term, 1L]

stable.molecular <-
  selection.frequency >= 0.60 &
  !is.na(direction.concordance) &
  direction.concordance >= 0.80 &
  abs(full.discovery.coefficient) > selection.tolerance

selection.stability <- molecular.metadata
selection.stability <- cbind(
  selection.stability,
  stability.summary[, c(
    "selected_count", "models_evaluated", "selection_frequency",
    "positive_count", "negative_count", "direction_concordance",
    "median_selected_coefficient"
  )]
)
selection.stability$full_discovery_coefficient <- full.discovery.coefficient
selection.stability$stable <- stable.molecular
selection.stability$near_miss <- FALSE
genome.wide.rows <- selection.stability$role == "genome_wide_candidate"
if (!any(stable.molecular) && any(selected.count[genome.wide.rows] > 0L)) {
  near.miss.order <- order(
    -selection.frequency,
    -direction.concordance,
    -abs(median.selected.coefficient),
    selection.stability$term
  )
  near.miss.order <- near.miss.order[
    genome.wide.rows[near.miss.order] &
      selected.count[near.miss.order] > 0L
  ]
  selection.stability$near_miss[head(near.miss.order, 10L)] <- TRUE
}
selection.stability <- selection.stability[
  order(
    -selection.stability$stable,
    -selection.stability$selection_frequency,
    -selection.stability$direction_concordance,
    -abs(selection.stability$median_selected_coefficient)
  ),
]

# The one full-data fit retains the complete eligible predictor matrix. The
# coefficient interface is a reporting view containing only molecular features
# that also passed the repeated-fold stability criteria; it is not a refit.
stable.molecular.columns <- molecular.columns[stable.molecular]
final.cv <- full.discovery.cv
final.model <- full.discovery.cv$glmnet.fit
final.lambda <- full.discovery.cv$lambda.1se
final.coefficients <- selection.stability[
  selection.stability$stable,
  c(
    "term", "modality", "gene", "role", "label",
    "full_discovery_coefficient", "selection_frequency",
    "direction_concordance", "stable"
  ),
  drop = FALSE
]
names(final.coefficients)[
  names(final.coefficients) == "full_discovery_coefficient"
] <- "coefficient"

## FEATURE DESIGNATION ########################################################
target.designation <- data.frame(
  modality = c("expression", "mutation", "copy_number"),
  feature = target,
  source = "target_characterization",
  model_entry = ifelse(
    c(
      target.expression.eligible,
      target.mutation.eligible,
      target.copy.number.eligible
    ),
    "Entered", "Omitted"
  ),
  predictors_entered = as.integer(c(
    target.expression.eligible,
    target.mutation.eligible,
    target.copy.number.eligible
  )),
  penalty_factor = ifelse(c(
    target.expression.eligible,
    target.mutation.eligible,
    target.copy.number.eligible
  ), 1, NA_real_),
  reason = c(
    if (target.expression.eligible) "Entered with the common molecular penalty" else "Omitted: expression eligibility failed",
    if (target.mutation.eligible) "Entered with the common molecular penalty" else "Omitted: fewer than 10 altered or 10 wild-type models",
    if (target.copy.number.eligible) "Entered with the common molecular penalty" else "Omitted: copy-number eligibility failed"
  )
)
make_hypothesis_designation <- function(modality, requested, retained) {
  if (!length(requested)) return(data.frame())
  selected <- requested %in% retained
  data.frame(
    modality = modality,
    feature = requested,
    source = "prespecified_hypothesis",
    model_entry = ifelse(selected, "Entered", "Omitted"),
    predictors_entered = as.integer(selected),
    penalty_factor = ifelse(selected, 1, NA_real_),
    reason = ifelse(
      selected,
      "Entered with the common molecular penalty",
      "Omitted: basic modality eligibility failed"
    )
  )
}
hypothesis.designation <- rbind(
  make_hypothesis_designation(
    "expression", requested.configured.hypotheses$expression,
    configured.expression
  ),
  make_hypothesis_designation(
    "mutation", requested.configured.hypotheses$mutation,
    configured.mutation
  ),
  make_hypothesis_designation(
    "copy_number", requested.configured.hypotheses$copy_number,
    configured.copy.number
  )
)
lineage.designation <- data.frame(
  modality = "annotation",
  feature = "Cancer lineage",
  source = "lineage_adjustment",
  model_entry = "Entered",
  predictors_entered = ncol(lineage.matrix),
  penalty_factor = 0,
  reason = "Entered without penalization to adjust for histological context"
)
genome.wide.designation <- do.call(rbind, lapply(
  c("expression", "mutation", "copy_number"),
  function(modality.name) {
    modality.rows <-
      selection.stability$modality == modality.name &
      selection.stability$role == "genome_wide_candidate"
    data.frame(
      modality = modality.name,
      feature = "Remaining eligible genome-wide features",
      source = "genome_wide_candidate",
      model_entry = "Entered",
      predictors_entered = sum(modality.rows),
      penalty_factor = 1,
      reason = "Entered with the same penalty as target and pre-specified-hypothesis features"
    )
  }
))
feature.designation <- rbind(
  target.designation,
  hypothesis.designation,
  genome.wide.designation,
  lineage.designation
)

## WRITE OUTPUTS ##############################################################
fwrite(
  lineage.specific.associations,
  file.path(table.dir, "11_lineage_specific_associations.csv")
)
fwrite(feature.designation, file.path(table.dir, "12_feature_designation.csv"))
fwrite(
  final.coefficients,
  file.path(table.dir, "13_integrated_model_coefficients.csv")
)
fwrite(selection.stability, file.path(table.dir, "14_selection_stability.csv"))
fwrite(nested.performance, file.path(table.dir, "15_nested_cv_performance.csv"))
fwrite(nested.predictions, file.path(table.dir, "16_nested_cv_predictions.csv"))
fwrite(
  lineage.performance,
  file.path(table.dir, "17_lineage_cv_performance.csv")
)

unlink(file.path(table.dir, c(
  "11_lineage_adjusted_associations.csv",
  "12_lineage_specific_associations.csv",
  "12_biomarker_selection.csv",
  "13_biomarker_selection.csv",
  "13_multivariable_model_coefficients.csv",
  "14_multivariable_model_coefficients.csv",
  "14_model_performance.csv",
  "15_model_performance.csv",
  "15_model_diagnostics.csv",
  "16_model_diagnostics.csv"
)))
unlink(file.path(intermediate.dir, "multivariable_model.rds"))
saveRDS(
  list(
    model = final.model,
    cv_model = final.cv,
    lambda = final.lambda,
    alpha = elastic.alpha,
    coefficients = final.coefficients,
    stability = selection.stability,
    performance = nested.performance,
    predictions = nested.predictions,
    lineage_performance = lineage.performance,
    preprocessor = full.preprocessor,
    molecular_columns = colnames(molecular.matrix),
    stable_molecular_columns = molecular.metadata$term[stable.molecular],
    predictor_metadata = all.metadata,
    penalty_factor = penalty.factor,
    lineage_levels = levels(lineage.model),
    reference_lineage = reference.lineage
  ),
  file.path(intermediate.dir, "integrated_model.rds")
)

message(
  "Integrated model searched ", length(molecular.columns),
  " equally penalized molecular features and identified ",
  sum(stable.molecular), " stable molecular predictors."
)
