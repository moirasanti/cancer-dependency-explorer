## DESCRIPTION ################################################################
## CHARACTERISE THE DEPENDENCY LANDSCAPE
# This stage describes how strongly cell lines depend on the selected target
# before molecular predictors are considered. It summarises the pan-cancer
# Chronos distribution and asks whether dependency differs across broad lineages.
#
# Expected objects from the interactive setup or run_analysis.R:
#   - target: requested HGNC gene symbol.
#   - table.dir and intermediate.dir: target-specific output directories.
#   - analysis_data.rds: the model-level dataset created by script 01.
#
# Only lineages represented by at least 10 models enter the descriptive lineage
# analysis. A Kruskal-Wallis test assesses the global lineage difference.
# Lineage vs rest Wilcoxon tests provide effect sizes and FDR-adjusted
# comparisons. The script writes the distribution and lineage tables and saves
# lineage.data for plotting.

## LOAD DATA ###################################################################
# Summarise the pan-cancer dependency distribution.
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))
observed.dependency <- analysis.data$dependency[!is.na(analysis.data$dependency)]

dependency.summary <- data.frame(
  target = target,
  n = length(observed.dependency),
  mean = mean(observed.dependency),
  sd = sd(observed.dependency),
  minimum = min(observed.dependency),
  q25 = unname(quantile(observed.dependency, 0.25)),
  median = median(observed.dependency),
  q75 = unname(quantile(observed.dependency, 0.75)),
  maximum = max(observed.dependency),
  proportion_below_minus_0_5 = mean(observed.dependency < -0.5),
  proportion_below_minus_1 = mean(observed.dependency < -1)
)

# Retain lineages with enough models for stable descriptive comparisons.
lineage.data <- analysis.data[
  complete.cases(analysis.data[, c("dependency", "lineage")]),
]
lineage.counts <- table(lineage.data$lineage)
eligible.lineages <- names(lineage.counts[lineage.counts >= 10L])
lineage.data <- lineage.data[lineage %in% eligible.lineages]
if (length(eligible.lineages) < 2L) {
  stop("Fewer than two lineages have at least 10 models.")
}

## MAIN ANALYSIS ###############################################################
# Test whether the dependency distribution differs globally across lineages.
# Use Kruskal-Wallis test non-parametric alternative to one-way ANOVA.
global.lineage.test <- kruskal.test(dependency ~ lineage, data = lineage.data)

# Compare each eligible lineage with all remaining eligible models.
lineage.results <- vector("list", length(eligible.lineages))
for (lineage.index in seq_along(eligible.lineages)) {
  current.lineage <- eligible.lineages[[lineage.index]]
  inside.scores <- lineage.data[lineage == current.lineage, dependency]
  outside.scores <- lineage.data[lineage != current.lineage, dependency]
  lineage.test <- suppressWarnings(wilcox.test(inside.scores, outside.scores, exact = FALSE))

  lineage.results[[lineage.index]] <- data.frame(
    lineage = current.lineage,
    n = length(inside.scores),
    mean = mean(inside.scores),
    median = median(inside.scores),
    rest_n = length(outside.scores),
    rest_median = median(outside.scores),
    median_difference = median(inside.scores) - median(outside.scores),
    p_value = lineage.test$p.value
  )
}
lineage.summary <- rbindlist(lineage.results)
lineage.summary[, fdr := p.adjust(p_value, method = "BH")]
setorder(lineage.summary, median_difference)

# Write distribution and lineage results for plotting and reporting.
fwrite(dependency.summary, file.path(table.dir, "04_dependency_summary.csv"))
fwrite(
  data.frame(
    test = "Kruskal-Wallis",
    n = nrow(lineage.data),
    groups = length(eligible.lineages),
    statistic = unname(global.lineage.test$statistic),
    p_value = global.lineage.test$p.value
  ),
  file.path(table.dir, "05_lineage_global_test.csv")
)
fwrite(lineage.summary, file.path(table.dir, "06_lineage_associations.csv"))
saveRDS(lineage.data, file.path(intermediate.dir, "lineage_data.rds"))

message("Tested ", length(eligible.lineages), " lineages with at least 10 models.")
