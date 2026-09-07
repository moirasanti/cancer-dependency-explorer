#!/usr/bin/env Rscript

source("R/util.R")

# Gene symbols are recovered from both DepMap and duplicate-column formats.
stopifnot(identical(
  extract_gene_symbol(c("TP53 (7157)", "A1BG..1")),
  c("TP53", "A1BG")
))

# Small synthetic DepMap data exercise selective and aligned matrix readers.
test.file <- tempfile(fileext = ".csv")
test.data <- data.frame(
  ModelID = c("ACH-1", "ACH-2"),
  IsDefaultEntryForModel = c("Yes", "No"),
  check.names = FALSE
)
test.data[["TP53 (7157)"]] <- c(1, 0)
test.data[["MDM2 (4193)"]] <- c(-1.2, -0.3)
fwrite(test.data, test.file)

selected.data <- read_depmap_gene_data(
  test.file,
  c("TP53", "MDM2"),
  default.only = TRUE
)
stopifnot(
  nrow(selected.data) == 1L,
  identical(names(selected.data), c("ModelID", "TP53", "MDM2"))
)

aligned.data <- read_aligned_genomic_matrix(test.file, c("ACH-2", "ACH-1"))
stopifnot(
  identical(aligned.data$symbols, c("TP53", "MDM2")),
  is.na(aligned.data$matrix[1, 1]),
  aligned.data$matrix[2, 1] == 1,
  aligned.data$matrix[2, 2] == -1.2
)

# Duplicate mutation symbols retain an alteration from any source column;
# continuous duplicates use all available values and report disagreements.
duplicate.file <- tempfile(fileext = ".csv")
duplicate.data <- data.frame(
  ModelID = paste0("ACH-", 1:4),
  check.names = FALSE
)
duplicate.data[["DUP (1)"]] <- c(0, 0, NA, NA)
duplicate.data[["DUP..1"]] <- c(1, 0, 1, NA)
fwrite(duplicate.data, duplicate.file)

mutation.messages <- capture.output(
  mutation.duplicate <- read_aligned_genomic_matrix(
    duplicate.file,
    duplicate.data$ModelID,
    modality = "mutation"
  ),
  type = "message"
)
stopifnot(
  identical(mutation.duplicate$symbols, "DUP"),
  isTRUE(all.equal(
    as.numeric(mutation.duplicate$matrix[, "DUP"]),
    c(1, 0, 1, NA)
  )),
  any(grepl("classified as altered", mutation.messages, fixed = TRUE))
)

copy.number.messages <- capture.output(
  copy.number.duplicate <- read_aligned_genomic_matrix(
    duplicate.file,
    duplicate.data$ModelID,
    modality = "copy_number"
  ),
  type = "message"
)
stopifnot(
  isTRUE(all.equal(
    as.numeric(copy.number.duplicate$matrix[, "DUP"]),
    c(0.5, 0, 1, NA)
  )),
  any(grepl("used the mean", copy.number.messages, fixed = TRUE))
)

selected.messages <- capture.output(
  selected.duplicate <- read_depmap_gene_data(
    duplicate.file,
    "DUP",
    modality = "mutation"
  ),
  type = "message"
)
stopifnot(
  isTRUE(all.equal(selected.duplicate$DUP, c(1, 0, 1, NA))),
  any(grepl("classified as altered", selected.messages, fixed = TRUE))
)

# Statistical helpers return stable values for valid and degenerate inputs.
spearman.result <- calculate_safe_spearman(1:10, 10:1)
stopifnot(
  spearman.result[["n"]] == 10,
  abs(spearman.result[["estimate"]] + 1) < 1e-12
)

# The mutation helper uses one built-in Welch test for both its P value and CI.
altered.scores <- c(-1.4, -1.1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2)
wild.type.scores <- seq(-0.3, 0.3, length.out = 25L)
expected.welch <- t.test(altered.scores, wild.type.scores, var.equal = FALSE)
observed.welch <- calculate_safe_welch(altered.scores, wild.type.scores)
ineligible.welch <- calculate_safe_welch(altered.scores[1:9], wild.type.scores)
stopifnot(
  isTRUE(all.equal(observed.welch[["standard_error"]], expected.welch$stderr)),
  isTRUE(all.equal(
    unname(observed.welch[c("conf_low", "conf_high")]),
    as.numeric(expected.welch$conf.int)
  )),
  isTRUE(all.equal(observed.welch[["p_value"]], expected.welch$p.value)),
  observed.welch[["eligible"]] == 1,
  ineligible.welch[["eligible"]] == 0
)

# Lineage follow-up models recover feature associations within each lineage.
set.seed(20260828L)
adjusted.lineage <- rep(c("A", "B", "C"), each = 40L)
adjusted.feature <- rnorm(120L)
adjusted.outcome <-
  0.7 * as.numeric(scale(adjusted.feature)) +
  rep(c(-1, 0, 1), each = 40L) +
  rnorm(120L, sd = 0.2)
lineage.specific.association <- fit_lineage_specific_associations(
  adjusted.outcome,
  adjusted.feature,
  adjusted.lineage,
  continuous = TRUE,
  minimum.lineage.size = 20L
)
stopifnot(
  nrow(lineage.specific.association) == 3L,
  all(lineage.specific.association$estimable),
  all(lineage.specific.association$p_value < 0.05)
)

# Logical damaging-mutation inputs are encoded as numeric 0/1 so the fitted
# coefficient retains the common feature_model term name.
logical.mutation.association <- fit_lineage_specific_associations(
  outcome = c(rep(0, 20), rep(-1, 20)),
  feature = c(rep(FALSE, 20), rep(TRUE, 20)),
  lineage = rep("A", 40),
  continuous = FALSE,
  minimum.lineage.size = 20L,
  minimum.mutation.group.size = 10L
)
stopifnot(
  nrow(logical.mutation.association) == 1L,
  logical.mutation.association$estimable,
  isTRUE(all.equal(logical.mutation.association$estimate, -1))
)

# Block-wise Spearman screening matches direct rank correlations with ties and
# feature-specific missingness, and safely handles a constant feature.
screen.matrix <- cbind(
  increasing = c(1, 2, 2, 4, 5, NA),
  decreasing = c(6, 5, NA, 3, 2, 1),
  constant = rep(1, 6)
)
screen.outcome <- c(1, 3, 2, 5, 4, 6)
block.screen <- calculate_block_spearman(
  screen.matrix,
  screen.outcome,
  block.size = 2L
)
direct.rho <- vapply(seq_len(ncol(screen.matrix)), function(column.index) {
  complete.rows <- complete.cases(screen.matrix[, column.index], screen.outcome)
  suppressWarnings(cor(
    screen.matrix[complete.rows, column.index],
    screen.outcome[complete.rows],
    method = "spearman"
  ))
}, numeric(1))
stopifnot(
  isTRUE(all.equal(block.screen$spearman_rho[1:2], direct.rho[1:2])),
  identical(block.screen$n, c(5L, 5L, 6L)),
  all(is.finite(block.screen$p_value[1:2])),
  is.na(block.screen$spearman_rho[[3L]])
)

stopifnot(identical(
  classify_correlation_strength(c(NA, 0.29, -0.30, 0.49, -0.50)),
  c(NA_character_, "small", "moderate", "moderate", "large")
))

# Empty YAML modalities create valid zero-row feature-family tables.
empty.hypothesis.rows <- make_feature_family_rows(
  "prespecified_hypothesis", "expression", character()
)
target.family.rows <- make_feature_family_rows(
  "target_characterization",
  c("expression", "mutation", "copy_number"),
  rep("EGFR", 3L)
)
stopifnot(
  identical(
    names(empty.hypothesis.rows),
    c("family", "modality", "gene")
  ),
  nrow(empty.hypothesis.rows) == 0L,
  identical(
    target.family.rows$modality,
    c("expression", "mutation", "copy_number")
  ),
  all(target.family.rows$gene == "EGFR")
)

# Fold assignment is deterministic, lineage-stratified and assigns singleton
# strata rather than accidentally sampling an unrelated row index.
fold.strata <- c(rep("A", 20), rep("B", 12), "singleton")
fold.one <- make_stratified_folds(fold.strata, folds = 5L, seed = 42L)
fold.two <- make_stratified_folds(fold.strata, folds = 5L, seed = 42L)
stopifnot(
  identical(fold.one, fold.two),
  all(fold.one %in% 1:5),
  all(table(fold.one[fold.strata == "A"]) >= 4L)
)

# Imputation and scaling are learned only from training rows and then applied
# unchanged to held-out observations.
preprocess.matrix <- cbind(
  continuous = c(1, 2, 3, 100),
  mutation = c(0, 0, 1, 1)
)
preprocess.fit <- fit_integrated_preprocessor(
  preprocess.matrix, 1:3, continuous.columns = 1L, mutation.columns = 2L
)
preprocess.test <- apply_integrated_preprocessor(
  preprocess.matrix, 4L, preprocess.fit
)
stopifnot(
  preprocess.fit$imputation[[1L]] == 2,
  preprocess.fit$imputation[[2L]] == 0,
  preprocess.fit$center[[1L]] == 2,
  isTRUE(all.equal(preprocess.fit$center[[2L]], 1 / 3)),
  preprocess.test[[1L]] > 90,
  isTRUE(all.equal(
    preprocess.test[[2L]],
    (1 - preprocess.fit$center[[2L]]) / preprocess.fit$scale[[2L]]
  ))
)

# A mutation that is constant in one training fold remains transformable, and
# its held-out value cannot influence the training-derived center or scale.
constant.fold.matrix <- cbind(
  continuous = 1:5,
  mutation = c(0, 0, 0, 0, 1)
)
constant.fold.fit <- fit_integrated_preprocessor(
  constant.fold.matrix, 1:4, continuous.columns = 1L, mutation.columns = 2L
)
constant.fold.test <- apply_integrated_preprocessor(
  constant.fold.matrix, 5L, constant.fold.fit
)
stopifnot(
  constant.fold.fit$center[[2L]] == 0,
  constant.fold.fit$scale[[2L]] == 1,
  constant.fold.test[[2L]] == 1
)

# Stability summaries are deterministic and distinguish frequency from
# direction concordance for correlated or substitutable signals.
stability.matrix <- rbind(
  consistent = c(0.4, 0.2, 0, 0.1, 0.3),
  inconsistent = c(-0.2, 0.1, 0, -0.3, 0.2),
  absent = rep(0, 5)
)
stability.one <- summarize_elastic_net_stability(stability.matrix)
stability.two <- summarize_elastic_net_stability(stability.matrix)
stopifnot(
  identical(stability.one, stability.two),
  stability.one["consistent", "selection_frequency"] == 0.8,
  stability.one["consistent", "direction_concordance"] == 1,
  stability.one["inconsistent", "direction_concordance"] == 0.5,
  is.na(stability.one["absent", "direction_concordance"])
)

if (requireNamespace("glmnet", quietly = TRUE)) {
  # A suppressor pair can be jointly selected despite weak marginal effects.
  set.seed(20260901L)
  suppressor.n <- 500L
  shared.signal <- rnorm(suppressor.n, sd = 5)
  suppressor.one <- shared.signal + rnorm(suppressor.n)
  suppressor.two <- shared.signal + rnorm(suppressor.n)
  suppressor.y <- 2 * (suppressor.one - suppressor.two) + rnorm(suppressor.n, sd = 0.2)
  suppressor.x <- cbind(suppressor_one = suppressor.one, suppressor_two = suppressor.two)
  suppressor.folds <- make_stratified_folds(rep(c("A", "B"), each = 250), 5L, 7L)
  suppressor.fit <- glmnet::cv.glmnet(
    suppressor.x, suppressor.y, family = "gaussian", alpha = 0.5,
    foldid = suppressor.folds, penalty.factor = c(1, 1),
    standardize = TRUE, type.measure = "mse"
  )
  suppressor.coef <- as.matrix(coef(suppressor.fit, s = "lambda.1se"))
  stopifnot(
    all(abs(cor(suppressor.x, suppressor.y)) < 0.30),
    all(abs(suppressor.coef[c("suppressor_one", "suppressor_two"), 1]) > 0)
  )

  # Correlated predictors and a null outcome produce deterministic selections;
  # a zero-feature molecular result remains permissible under lambda.1se.
  set.seed(20260902L)
  null.x <- matrix(rnorm(300 * 20), 300, 20)
  null.y <- rnorm(300)
  null.folds <- make_stratified_folds(rep(c("A", "B", "C"), each = 100), 5L, 9L)
  null.fit.one <- glmnet::cv.glmnet(
    null.x, null.y, family = "gaussian", alpha = 0.5,
    foldid = null.folds, standardize = TRUE, type.measure = "mse"
  )
  null.fit.two <- glmnet::cv.glmnet(
    null.x, null.y, family = "gaussian", alpha = 0.5,
    foldid = null.folds, standardize = TRUE, type.measure = "mse"
  )
  null.coef.one <- as.matrix(coef(null.fit.one, s = "lambda.1se"))[-1, 1]
  null.coef.two <- as.matrix(coef(null.fit.two, s = "lambda.1se"))[-1, 1]
  stopifnot(
    identical(null.coef.one, null.coef.two),
    sum(abs(null.coef.one) > 1e-10) == 0L
  )

  # Target, pre-specified-hypothesis and genome-wide molecular features share
  # one penalty and can all be reduced to zero; only lineage is unpenalized.
  set.seed(20260903L)
  equal.penalty.n <- 300L
  lineage.adjustment <- rep(c(0, 1), each = equal.penalty.n / 2L)
  equal.penalty.x <- cbind(
    target_characterization = rnorm(equal.penalty.n),
    prespecified_hypothesis = rnorm(equal.penalty.n),
    genome_wide_candidate = rnorm(equal.penalty.n),
    lineage_adjustment = lineage.adjustment
  )
  equal.penalty.y <- 2 * lineage.adjustment + rnorm(equal.penalty.n, sd = 0.1)
  equal.penalty.vector <- c(1, 1, 1, 0)
  equal.penalty.fit <- glmnet::cv.glmnet(
    equal.penalty.x, equal.penalty.y,
    family = "gaussian", alpha = 0.5,
    foldid = null.folds,
    penalty.factor = equal.penalty.vector,
    standardize = FALSE, type.measure = "mse"
  )
  equal.penalty.coefficients <- as.matrix(coef(
    equal.penalty.fit, s = "lambda.1se"
  ))[-1L, 1L]
  stopifnot(
    identical(equal.penalty.vector, c(1, 1, 1, 0)),
    all(equal.penalty.coefficients[1:3] == 0),
    equal.penalty.coefficients[[4L]] != 0
  )
}

# YAML expression hypotheses are normalized while legacy files remain valid.
legacy.config.file <- tempfile(fileext = ".yml")
expression.config.file <- tempfile(fileext = ".yml")
writeLines(c(
  "target: MDM2",
  "mutation:",
  "  - TP53",
  "copy_number: []"
), legacy.config.file)
writeLines(c(
  "target: MDM2",
  "expression:",
  "  - egfr",
  "  - EGFR",
  "mutation: []",
  "copy_number: []"
), expression.config.file)
legacy.config <- read_biomarker_configuration(legacy.config.file, "MDM2")
expression.config <- read_biomarker_configuration(expression.config.file, "MDM2")
stopifnot(
  !length(legacy.config$expression),
  identical(legacy.config$mutation, "TP53"),
  identical(expression.config$expression, "EGFR")
)

# Configured expression hypotheses bypass marginal FDR/magnitude thresholds but
# still fail fast for unknown genes and omit basic eligibility failures.
configured.expression.results <- data.frame(
  gene = c("FORCED", "LOW_COVERAGE", "ZERO_VARIANCE"),
  eligible = c(TRUE, FALSE, FALSE),
  fdr = c(0.90, 0.001, 0.001),
  spearman_rho = c(0.01, 0.80, 0.80)
)
stopifnot(identical(
  validate_configured_features(
    "FORCED",
    configured.expression.results,
    "expression",
    "coverage/variance"
  ),
  "FORCED"
))
ineligible.warning <- character()
retained.expression <- withCallingHandlers(
  validate_configured_features(
    c("LOW_COVERAGE", "ZERO_VARIANCE"),
    configured.expression.results,
    "expression",
    "coverage/variance"
  ),
  warning = function(condition) {
    ineligible.warning <<- conditionMessage(condition)
    invokeRestart("muffleWarning")
  }
)
unknown.expression <- try(
  validate_configured_features(
    "UNKNOWN",
    configured.expression.results,
    "expression",
    "coverage/variance"
  ),
  silent = TRUE
)
stopifnot(
  !length(retained.expression),
  any(grepl("coverage/variance", ineligible.warning, fixed = TRUE)),
  inherits(unknown.expression, "try-error")
)

duplicate.messages <- capture.output(
  deduplicated.config <- deduplicate_target_hypotheses(
    list(
      expression = c("MDM2", "EGFR"),
      mutation = "MDM2",
      copy_number = "MDM2"
    ),
    "MDM2"
  ),
  type = "message"
)
stopifnot(
  identical(deduplicated.config$expression, "EGFR"),
  !length(deduplicated.config$mutation),
  !length(deduplicated.config$copy_number),
  length(duplicate.messages) == 3L
)

# Core target mutation follows the same established eligibility flag.
mutation.eligibility <- data.frame(
  gene = c("ELIGIBLE", "RARE"),
  eligible = c(TRUE, FALSE)
)
stopifnot(
  is_eligible_mutation_feature(mutation.eligibility, "ELIGIBLE"),
  !is_eligible_mutation_feature(mutation.eligibility, "RARE"),
  !is_eligible_mutation_feature(mutation.eligibility, "ABSENT")
)

# The shared plotting helper returns a ggplot2 theme object.
stopifnot(inherits(dependency_plot_theme(), "theme"))

# Molecular features are screened once and the report reuses those results for
# target characterization, pre-specified hypotheses and discovery displays.
molecular.script.text <- paste(
  readLines("scripts/03.molecular.associations.R"),
  collapse = "\n"
)
figure.script.text <- paste(
  readLines("scripts/05.generate.figures.R"),
  collapse = "\n"
)
model.script.text <- paste(
  readLines("scripts/04.multivariable.model.R"),
  collapse = "\n"
)
run.script.text <- paste(readLines("run_analysis.R"), collapse = "\n")
readme.text <- paste(readLines("README.md"), collapse = "\n")
gitignore.text <- paste(readLines(".gitignore"), collapse = "\n")
report.lines <- readLines("reports/target_assessment.Rmd")
report.text <- paste(report.lines, collapse = "\n")
raw.div.opens <- sum(grepl("^<div(?: |>)", report.lines, perl = TRUE))
raw.div.closes <- sum(grepl("^</div>$", report.lines))
raw.details.opens <- sum(grepl("^<details(?: |>)", report.lines, perl = TRUE))
raw.details.closes <- sum(grepl("^</details>$", report.lines))
performance.section <- sub(
  ".*# Performance on held-out cell lines",
  "# Performance on held-out cell lines",
  report.text
)
performance.section <- sub(
  "## Observed versus predicted dependency.*",
  "",
  performance.section
)
performance.metric.positions <- vapply(
  c("**Pearson r**", "**R-squared**", "**RMSE**"),
  function(metric) regexpr(metric, performance.section, fixed = TRUE)[[1L]],
  integer(1)
)
report.headings <- c(
  "# Target characterisation",
  "# Pre-specified hypotheses",
  "# Lineage-specific molecular associations",
  "# Genomic associations",
  "# Feature designation for integrated exploratory multiomic modelling",
  "# Integrated exploratory model",
  "# Performance on held-out cell lines"
)
report.heading.positions <- vapply(
  report.headings,
  function(heading) regexpr(heading, report.text, fixed = TRUE)[[1L]],
  integer(1)
)
stopifnot(
  !grepl("calculate_safe_spearman(", molecular.script.text, fixed = TRUE),
  grepl("calculate_safe_welch(", molecular.script.text, fixed = TRUE),
  grepl("t.test(altered_scores, wild_type_scores", paste(readLines("R/util.R"), collapse = "\n"), fixed = TRUE),
  !grepl("1.96 * mutation.standard.error", molecular.script.text, fixed = TRUE),
  grepl("07_target_associations.csv", molecular.script.text, fixed = TRUE),
  grepl("special_features.rds", molecular.script.text, fixed = TRUE),
  !grepl("!expression.associations$gene %in% special.features", molecular.script.text, fixed = TRUE),
  !grepl("!mutation.associations$gene %in% special.features", molecular.script.text, fixed = TRUE),
  !grepl("!copy.number.associations$gene %in% special.features", molecular.script.text, fixed = TRUE),
  grepl('method = "lm"', figure.script.text, fixed = TRUE),
  !grepl('method = "loess"', figure.script.text, fixed = TRUE),
  !grepl("%in% special.features", figure.script.text, fixed = TRUE),
  grepl("abs(copy.number.plot.data$correlation) >= 0.30", figure.script.text, fixed = TRUE),
  grepl("|Pearson r| ≥ 0.30", figure.script.text, fixed = TRUE),
  grepl('paste("Spearman correlation with", target, "Chronos dependency")', figure.script.text, fixed = TRUE),
  grepl('"Mean", target, "Chronos difference: altered minus wild type"', figure.script.text, fixed = TRUE),
  grepl('title = "Top copy number associations"', figure.script.text, fixed = TRUE),
  grepl('paste("Pearson correlation with", target, "Chronos dependency")', figure.script.text, fixed = TRUE),
  !grepl('title = "Top copy-number associations"', figure.script.text, fixed = TRUE),
  !grepl("altered minus reference", figure.script.text, fixed = TRUE),
  !grepl("lineage_adjusted_associations", report.text, fixed = TRUE),
  grepl('--release LABEL', run.script.text, fixed = TRUE),
  grepl('release = release', run.script.text, fixed = TRUE),
  grepl('message("DepMap release: ", release)', run.script.text, fixed = TRUE),
  grepl('paste0(target, "_biomarkers.yml")', run.script.text, fixed = TRUE),
  grepl('Biomarker file must be named', run.script.text, fixed = TRUE),
  grepl('biomarker_file = if (nzchar(biomarker.file))', run.script.text, fixed = TRUE),
  regexpr('rmarkdown::render(', run.script.text, fixed = TRUE)[[1L]] <
    regexpr('capture.output(sessionInfo())', run.script.text, fixed = TRUE)[[1L]],
  !grepl('capture.output(sessionInfo())', paste(readLines("scripts/01.prepare.depmap.data.R"), collapse = "\n"), fixed = TRUE),
  grepl('results/*/tables/16_nested_cv_predictions.csv', gitignore.text, fixed = TRUE),
  grepl('--release "DepMap Public 26Q1"', readme.text, fixed = TRUE),
  grepl('params$biomarker_file', report.text, fixed = TRUE),
  !grepl('paste0("config/", target, "_biomarkers.yml")', report.text, fixed = TRUE),
  grepl('rather than across all dependency models', report.text, fixed = TRUE),
  grepl('among models with any data from the corresponding assay', readme.text, fixed = TRUE),
  grepl('Nature Genetics 57, 522–529 (2025)', report.text, fixed = TRUE),
  identical(raw.div.opens, raw.div.closes),
  identical(raw.details.opens, raw.details.closes),
  grepl("</details>\n\n---\n\n# Data quality", report.text, fixed = TRUE),
  grepl("predictions, residuals and RMSE remain directly interpretable in Chronos units", report.text, fixed = TRUE),
  !grepl("fit_lineage_adjusted_association", model.script.text, fixed = TRUE),
  grepl("11_lineage_specific_associations.csv", model.script.text, fixed = TRUE),
  grepl("12_feature_designation.csv", model.script.text, fixed = TRUE),
  grepl("13_integrated_model_coefficients.csv", model.script.text, fixed = TRUE),
  grepl("14_selection_stability.csv", model.script.text, fixed = TRUE),
  grepl("15_nested_cv_performance.csv", model.script.text, fixed = TRUE),
  grepl("pearson_r = repeat.pearson.r", model.script.text, fixed = TRUE),
  grepl("16_nested_cv_predictions.csv", model.script.text, fixed = TRUE),
  grepl("17_lineage_cv_performance.csv", model.script.text, fixed = TRUE),
  grepl('family = "gaussian"', model.script.text, fixed = TRUE),
  grepl("elastic.alpha <- 0.5", model.script.text, fixed = TRUE),
  grepl('s = "lambda.1se"', model.script.text, fixed = TRUE),
  grepl("penalty.factor = penalty.factor", model.script.text, fixed = TRUE),
  grepl("rep(1, nrow(molecular.metadata))", model.script.text, fixed = TRUE),
  grepl("rep(0, ncol(lineage.matrix))", model.script.text, fixed = TRUE),
  grepl('"prespecified_hypothesis"', model.script.text, fixed = TRUE),
  grepl('"genome_wide_candidate"', model.script.text, fixed = TRUE),
  !grepl("get_qualifying_continuous_features", model.script.text, fixed = TRUE),
  !grepl("select_ranked_automatic_features", model.script.text, fixed = TRUE),
  !grepl("calculate_maximum_absolute_correlation", model.script.text, fixed = TRUE),
  !grepl("biologically_forced", model.script.text, fixed = TRUE),
  !grepl('role == "automatic"', model.script.text, fixed = TRUE),
  grepl("integrated_model.rds", figure.script.text, fixed = TRUE),
  grepl("08_observed_vs_predicted.png", figure.script.text, fixed = TRUE),
  grepl("Held-out Pearson r", figure.script.text, fixed = TRUE),
  grepl("names(repeat.labels) <- names(repeat.pearson.r)", figure.script.text, fixed = TRUE),
  grepl("Each point is one cell line. The dashed diagonal", figure.script.text, fixed = TRUE),
  grepl("09_lineage_cv_performance.png", figure.script.text, fixed = TRUE),
  grepl("## Observed versus predicted dependency", report.text, fixed = TRUE),
  grepl("The three panels show the three complete cross-validation repeats", report.text, fixed = TRUE),
  grepl("larger distances increase RMSE and reduce R-squared", report.text, fixed = TRUE),
  !grepl("Points far from the line are prediction errors", report.text, fixed = TRUE),
  !grepl("systematic compression toward the centre", report.text, fixed = TRUE),
  !grepl("predictions cover a narrower range", report.text, fixed = TRUE),
  !grepl("underestimates the strength of the most extreme", report.text, fixed = TRUE),
  grepl("## Performance by lineage", report.text, fixed = TRUE),
  all(performance.metric.positions > 0L),
  identical(performance.metric.positions, sort(performance.metric.positions)),
  grepl("The mean held-out Pearson r was", performance.section, fixed = TRUE),
  grepl("diff(pearson.range) <= 0.05", performance.section, fixed = TRUE),
  grepl("diff(r2.range) <= 0.05", performance.section, fixed = TRUE),
  grepl("diff(rmse.range) <= 0.05 * dependency.sd", performance.section, fixed = TRUE),
  grepl("similar across repeats", performance.section, fixed = TRUE),
  grepl("varied across repeats", performance.section, fixed = TRUE),
  !grepl("How to read these values", report.text, fixed = TRUE),
  !grepl("Overall, this is", performance.section, fixed = TRUE),
  !grepl("Overall, the model had", performance.section, fixed = TRUE),
  !grepl("conf_low", tail(readLines("scripts/05.generate.figures.R"), 100L), fixed = TRUE),
  !grepl("HC3", report.text, fixed = TRUE),
  !grepl("ordinary linear regression", report.text, fixed = TRUE),
  !grepl("multivariable_model_coefficients", report.text, fixed = TRUE),
  !grepl("<summary><strong>View feature-designation decisions", report.text, fixed = TRUE),
  !grepl("\\b(prior|priors|predefined hypothesis|predefined hypotheses)\\b", report.text, ignore.case = TRUE),
  !grepl("\\b(prior|priors|predefined hypothesis|predefined hypotheses)\\b", paste(readLines("README.md"), collapse = "\n"), ignore.case = TRUE),
  all(report.heading.positions > 0L),
  identical(report.heading.positions, sort(report.heading.positions))
)

# All utility function names use underscore nomenclature.
utility.lines <- readLines("R/util.R")
function.lines <- grep("^[a-z][a-z0-9_]* <- function", utility.lines, value = TRUE)
function.names <- sub(" <- function.*$", "", function.lines)
stopifnot(all(grepl("^[a-z][a-z0-9]*(_[a-z0-9]+)+$", function.names)))

unlink(c(test.file, duplicate.file, legacy.config.file, expression.config.file))
cat("All unit tests passed.\n")
