## DESCRIPTION ################################################################
## SCREEN MOLECULAR ASSOCIATIONS
# This stage tests each molecular feature once for association with dependency
# on the selected target. The genome-wide result objects include the target's
# own features and any YAML-specified hypotheses; those rows are reused later
# for dedicated report sections and excluded only from discovery displays.
#
# Expected objects from the interactive setup or run_analysis.R:
#   - target: requested HGNC gene symbol.
#   - input.files: named paths to the mutation and copy-number matrices.
#   - table.dir and intermediate.dir: target-specific output directories.
#   - analysis_data.rds: the model-level dataset created by script 01.
#
# Expression features are screened with block-wise Spearman correlations and
# require at least 80% coverage and nonzero variance. Mutation features are
# tested by comparing mean dependency between altered and
# reference models using Welch stats. Candidate genes require at least 10
# models in each group (defective vs proficient).
# Copy-number features are screened with correlations and require at least
# 80% coverage and nonzero variance.
# FDR correction is applied separately within each genome-wide modality.
# Complete association objects are saved for feature selection. Compact
# leading-result tables are written for the report.

## LOAD DATA ###################################################################
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))

# Record features that need dedicated reporting. They remain part of their
# complete modality-wide screens, FDR families, plots and leading-result tables.
configured.for.reporting <- read_biomarker_configuration(biomarker.file, target)
configured.for.reporting <- deduplicate_target_hypotheses(
  configured.for.reporting,
  target
)
special.features <- list(
  target = target,
  configured = configured.for.reporting,
  expression = unique(c(target, configured.for.reporting$expression)),
  mutation = unique(c(target, configured.for.reporting$mutation)),
  copy_number = unique(c(target, configured.for.reporting$copy_number))
)
saveRDS(
  special.features,
  file.path(intermediate.dir, "special_features.rds")
)

## MAIN ANALYSIS ###############################################################
# Screen genome-wide expression in bounded blocks to limit ranking temporaries.
message("Screening genome-wide expression...")
expression.data <- read_aligned_genomic_matrix(
  input.files[["expression"]],
  analysis.data$ModelID,
  modality = "expression"
)
expression.matrix <- expression.data$matrix
expression.outcome <- analysis.data$dependency
expression.complete <- !is.na(expression.outcome)
expression.matrix <- expression.matrix[expression.complete, , drop = FALSE]
expression.outcome <- expression.outcome[expression.complete]

expression.assayed <- rowSums(!is.na(expression.matrix)) > 0L
expression.matrix <- expression.matrix[expression.assayed, , drop = FALSE]
expression.outcome <- expression.outcome[expression.assayed]
expression.associations <- calculate_block_spearman(
  expression.matrix,
  expression.outcome,
  block.size = 500L
)
expression.associations$coverage <-
  expression.associations$n / length(expression.outcome)
expression.associations$correlation_strength <- classify_correlation_strength(
  expression.associations$spearman_rho
)
expression.associations$eligible <-
  expression.associations$coverage >= 0.80 &
  is.finite(expression.associations$sd) &
  expression.associations$sd > 0 &
  is.finite(expression.associations$spearman_rho) &
  abs(expression.associations$spearman_rho) < 1
expression.associations$fdr <- NA_real_
expression.associations$fdr[expression.associations$eligible] <- p.adjust(
  expression.associations$p_value[expression.associations$eligible],
  method = "BH"
)
expression.associations <- expression.associations[
  order(
    expression.associations$fdr,
    -abs(expression.associations$spearman_rho)
  ),
]
saveRDS(
  expression.associations,
  file.path(intermediate.dir, "expression_associations.rds")
)
fwrite(
  head(expression.associations[
    expression.associations$eligible,
  ], 50L),
  file.path(table.dir, "08_top_expression_associations.csv")
)
rm(expression.data, expression.matrix)
gc()

# Load and align the damaging-mutation matrix to models with target dependency data.
message("Screening damaging mutations...")
mutation.data <- read_aligned_genomic_matrix(
  input.files[["mutation"]],
  analysis.data$ModelID,
  modality = "mutation"
)
mutation.matrix <- mutation.data$matrix
mutation.outcome <- analysis.data$dependency
mutation.complete <- !is.na(mutation.outcome)
mutation.matrix <- mutation.matrix[mutation.complete, , drop = FALSE]
mutation.outcome <- mutation.outcome[mutation.complete]

# Define altered and reference groups for every mutation feature.
altered.matrix <- mutation.matrix > 0
altered.matrix[is.na(mutation.matrix)] <- NA
altered.n <- colSums(altered.matrix, na.rm = TRUE)
mutation.observed.n <- colSums(!is.na(altered.matrix))
reference.n <- mutation.observed.n - altered.n

# Calculate group means for every mutation gene.
altered.sum <- colSums(altered.matrix * mutation.outcome, na.rm = TRUE)
reference.sum <- colSums((!altered.matrix) * mutation.outcome, na.rm = TRUE)
altered.mean <- altered.sum / altered.n
reference.mean <- reference.sum / reference.n

# Let R's built-in t.test calculate mutually consistent Welch P values and
# confidence intervals for genes meeting the established 10/10 group rule.
mutation.candidate <- altered.n >= 10L & reference.n >= 10L
mutation.test.results <- matrix(
  NA_real_,
  nrow = ncol(altered.matrix),
  ncol = 5L,
  dimnames = list(NULL, c(
    "standard_error", "conf_low", "conf_high", "p_value", "eligible"
  ))
)
for (gene.index in which(mutation.candidate)) {
  mutation.test.results[gene.index, ] <- calculate_safe_welch(
    mutation.outcome[which(altered.matrix[, gene.index] %in% TRUE)],
    mutation.outcome[which(altered.matrix[, gene.index] %in% FALSE)]
  )
}
mutation.standard.error <- mutation.test.results[, "standard_error"]
mutation.conf.low <- mutation.test.results[, "conf_low"]
mutation.conf.high <- mutation.test.results[, "conf_high"]
mutation.p.value <- mutation.test.results[, "p_value"]
mutation.eligible <- mutation.test.results[, "eligible"] == 1
mutation.eligible[is.na(mutation.eligible)] <- FALSE

# Assemble, correct and rank the genome-wide mutation results.
mutation.associations <- data.frame(
  gene = mutation.data$symbols,
  altered_n = altered.n,
  reference_n = reference.n,
  altered_mean_dependency = altered.mean,
  reference_mean_dependency = reference.mean,
  effect = altered.mean - reference.mean,
  std_error = mutation.standard.error,
  conf_low = mutation.conf.low,
  conf_high = mutation.conf.high,
  p_value = mutation.p.value,
  eligible = mutation.eligible,
  stringsAsFactors = FALSE
)
mutation.associations$fdr <- NA_real_
mutation.associations$fdr[mutation.eligible] <- p.adjust(
  mutation.associations$p_value[mutation.eligible],
  method = "BH"
)
mutation.associations <- mutation.associations[
  order(mutation.associations$fdr, mutation.associations$effect),
]

saveRDS(
  mutation.associations,
  file.path(intermediate.dir, "mutation_associations.rds")
)
fwrite(
  head(mutation.associations[
    mutation.associations$eligible,
  ], 50L),
  file.path(table.dir, "09_top_mutation_associations.csv")
)
rm(mutation.data, mutation.matrix, altered.matrix)
gc()

## Load and align copy-number measurements to the target dependency cohort.
message("Screening copy-number features...")
copy.number.data <- read_aligned_genomic_matrix(
  input.files[["copy_number"]],
  analysis.data$ModelID,
  modality = "copy_number"
)
copy.number.matrix <- copy.number.data$matrix
copy.number.outcome <- analysis.data$dependency
copy.number.complete <- !is.na(copy.number.outcome)
copy.number.matrix <- copy.number.matrix[copy.number.complete, , drop = FALSE]
copy.number.outcome <- copy.number.outcome[copy.number.complete]

# Restrict coverage calculations to models with a copy-number assay.
copy.number.assayed <- rowSums(!is.na(copy.number.matrix)) > 0L
copy.number.matrix <- copy.number.matrix[copy.number.assayed, , drop = FALSE]
copy.number.outcome <- copy.number.outcome[copy.number.assayed]
copy.number.n <- colSums(!is.na(copy.number.matrix))
copy.number.coverage <- copy.number.n / length(copy.number.outcome)

# Calculate genome-wide copy-number correlations and their P values.
copy.number.means <- colMeans(copy.number.matrix, na.rm = TRUE)
copy.number.centered <- sweep(copy.number.matrix, 2L, copy.number.means, "-")
copy.number.variance <- colSums(copy.number.centered^2, na.rm = TRUE) /
  pmax(copy.number.n - 1, 1)
copy.number.correlation <- as.numeric(
  cor(copy.number.outcome, copy.number.matrix, use = "pairwise.complete.obs")
)
copy.number.eligible <-
  copy.number.coverage >= 0.80 &
  copy.number.variance > 0 &
  is.finite(copy.number.correlation) &
  abs(copy.number.correlation) < 1

copy.number.statistic <- copy.number.correlation * sqrt(
  pmax(copy.number.n - 2, 1) /
    pmax(1 - copy.number.correlation^2, 1e-12)
)
copy.number.p.value <- rep(NA_real_, length(copy.number.correlation))
copy.number.p.value[copy.number.eligible] <- 2 * pt(
  abs(copy.number.statistic[copy.number.eligible]),
  df = copy.number.n[copy.number.eligible] - 2,
  lower.tail = FALSE
)

# Assemble, correct and rank the genome-wide copy-number results.
copy.number.associations <- data.frame(
  gene = copy.number.data$symbols,
  n = copy.number.n,
  coverage = copy.number.coverage,
  sd = sqrt(copy.number.variance),
  correlation = copy.number.correlation,
  p_value = copy.number.p.value,
  eligible = copy.number.eligible,
  stringsAsFactors = FALSE
)
copy.number.associations$correlation_strength <- classify_correlation_strength(
  copy.number.associations$correlation
)
copy.number.associations$fdr <- NA_real_
copy.number.associations$fdr[copy.number.eligible] <- p.adjust(
  copy.number.associations$p_value[copy.number.eligible],
  method = "BH"
)
copy.number.associations <- copy.number.associations[
  order(copy.number.associations$fdr, -abs(copy.number.associations$correlation)),
]

saveRDS(
  copy.number.associations,
  file.path(intermediate.dir, "copy_number_associations.rds")
)
fwrite(
  head(copy.number.associations[
    copy.number.associations$eligible,
  ], 50L),
  file.path(table.dir, "10_top_copy_number_associations.csv")
)

# Reuse the modality-wide results to summarise all three target features. This
# table performs no additional statistical tests.
target.expression <- expression.associations[
  expression.associations$gene == target,
]
target.mutation <- mutation.associations[
  mutation.associations$gene == target,
]
if (!nrow(target.mutation)) {
  target.mutation <- data.frame(
    altered_n = 0L,
    reference_n = 0L,
    effect = NA_real_,
    p_value = NA_real_,
    fdr = NA_real_,
    eligible = FALSE
  )
}
target.copy.number <- copy.number.associations[
  copy.number.associations$gene == target,
]
target.associations <- rbind(
  data.frame(
    feature = paste(target, "expression"),
    modality = "expression",
    test = "Spearman correlation",
    n = as.character(target.expression$n),
    estimate = target.expression$spearman_rho,
    estimate_type = "Spearman rho",
    p_value = target.expression$p_value,
    fdr = target.expression$fdr,
    eligible = target.expression$eligible
  ),
  data.frame(
    feature = paste(target, "damaging mutation"),
    modality = "mutation",
    test = "Welch mean comparison",
    n = sprintf(
      "%d altered / %d reference",
      target.mutation$altered_n,
      target.mutation$reference_n
    ),
    estimate = target.mutation$effect,
    estimate_type = "Altered minus reference",
    p_value = target.mutation$p_value,
    fdr = target.mutation$fdr,
    eligible = target.mutation$eligible
  ),
  data.frame(
    feature = paste(target, "copy number"),
    modality = "copy_number",
    test = "Pearson correlation",
    n = as.character(target.copy.number$n),
    estimate = target.copy.number$correlation,
    estimate_type = "Pearson r",
    p_value = target.copy.number$p_value,
    fdr = target.copy.number$fdr,
    eligible = target.copy.number$eligible
  )
)
fwrite(
  target.associations,
  file.path(table.dir, "07_target_associations.csv")
)

rm(copy.number.data, copy.number.matrix, copy.number.centered)
gc()

message("Completed genome-wide expression, mutation and copy-number screens.")
