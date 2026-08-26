## DESCRIPTION ################################################################
## SCREEN MOLECULAR ASSOCIATIONS
# This stage tests which molecular features are associated with dependency on
# the selected target. It first evaluates the target's own expression and copy
# number, then performs genome-wide screens of damaging mutations and gene-level
# copy number.
#
# Expected objects from the interactive setup or run_analysis.R:
#   - target: requested HGNC gene symbol.
#   - input.files: named paths to the mutation and copy-number matrices.
#   - table.dir and intermediate.dir: target-specific output directories.
#   - analysis_data.rds: the model-level dataset created by script 01.
#
# Mutation features are tested by comparing mean dependency between altered and
# reference models using Welch stats. Candidate genes require at least 10
# models in each group (defective vs proficient).
# Copy-number features are screened with correlations and require at least
# 80% coverage and nonzero variance.
# FDR correction is applied within each modality.
# Complete association objects are saved for feature selection. Compact
# leading-result tables are written for the report.

## LOAD DATA ###################################################################
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))

## MAIN ANALYSIS ###############################################################
# Test the target's own expression and copy number against its dependency score.
expression.test <- calculate_safe_spearman(
  analysis.data$target_expression,
  analysis.data$dependency
)
copy.number.test <- calculate_safe_spearman(
  analysis.data$target_copy_number,
  analysis.data$dependency
)

continuous.associations <- data.frame(
  feature = c(paste0(target, " expression"), paste0(target, " copy number")),
  modality = c("expression", "copy_number"),
  n = c(expression.test[["n"]], copy.number.test[["n"]]),
  spearman_rho = c(expression.test[["estimate"]], copy.number.test[["estimate"]]),
  p_value = c(expression.test[["p_value"]], copy.number.test[["p_value"]])
)
fwrite(
  continuous.associations,
  file.path(table.dir, "07_target_continuous_associations.csv")
)

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

# Calculate vectorised Welch-test statistics for every eligible mutation gene.
altered.sum <- colSums(altered.matrix * mutation.outcome, na.rm = TRUE)
reference.sum <- colSums((!altered.matrix) * mutation.outcome, na.rm = TRUE)
altered.mean <- altered.sum / altered.n
reference.mean <- reference.sum / reference.n

altered.squared.sum <- colSums(altered.matrix * mutation.outcome^2, na.rm = TRUE)
reference.squared.sum <- colSums((!altered.matrix) * mutation.outcome^2, na.rm = TRUE)
altered.variance <- (altered.squared.sum - altered.n * altered.mean^2) /
  pmax(altered.n - 1, 1)
reference.variance <- (reference.squared.sum - reference.n * reference.mean^2) /
  pmax(reference.n - 1, 1)

mutation.standard.error <- sqrt(
  altered.variance / altered.n + reference.variance / reference.n
)
mutation.statistic <- (altered.mean - reference.mean) / mutation.standard.error
mutation.degrees.freedom <-
  (altered.variance / altered.n + reference.variance / reference.n)^2 /
  (
    (altered.variance / altered.n)^2 / pmax(altered.n - 1, 1) +
      (reference.variance / reference.n)^2 / pmax(reference.n - 1, 1)
  )

mutation.eligible <-
  altered.n >= 10L & reference.n >= 10L & is.finite(mutation.statistic)
mutation.p.value <- rep(NA_real_, length(mutation.statistic))
mutation.p.value[mutation.eligible] <- 2 * pt(
  abs(mutation.statistic[mutation.eligible]),
  df = mutation.degrees.freedom[mutation.eligible],
  lower.tail = FALSE
)

# Assemble, correct and rank the genome-wide mutation results.
mutation.associations <- data.frame(
  gene = mutation.data$symbols,
  altered_n = altered.n,
  reference_n = reference.n,
  altered_mean_dependency = altered.mean,
  reference_mean_dependency = reference.mean,
  effect = altered.mean - reference.mean,
  std_error = mutation.standard.error,
  conf_low = altered.mean - reference.mean - 1.96 * mutation.standard.error,
  conf_high = altered.mean - reference.mean + 1.96 * mutation.standard.error,
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
  head(mutation.associations[mutation.associations$eligible, ], 50L),
  file.path(table.dir, "08_top_mutation_associations.csv")
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
  head(copy.number.associations[copy.number.associations$eligible, ], 50L),
  file.path(table.dir, "09_top_copy_number_associations.csv")
)

rm(copy.number.data, copy.number.matrix, copy.number.centered)
gc()

message("Completed genome-wide molecular association screens.")
