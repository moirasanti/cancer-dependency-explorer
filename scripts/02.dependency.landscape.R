# Describe the dependency outcome before testing molecular predictors.

analysis_data <- readRDS(file.path(intermediate_dir, "analysis_data.rds"))
observed_dependency <- analysis_data$dependency[!is.na(analysis_data$dependency)]

dependency_summary <- data.frame(
  target = target,
  n = length(observed_dependency),
  mean = mean(observed_dependency),
  sd = sd(observed_dependency),
  minimum = min(observed_dependency),
  q25 = unname(quantile(observed_dependency, 0.25)),
  median = median(observed_dependency),
  q75 = unname(quantile(observed_dependency, 0.75)),
  maximum = max(observed_dependency),
  proportion_below_minus_0_5 = mean(observed_dependency < -0.5),
  proportion_below_minus_1 = mean(observed_dependency < -1)
)

lineage_data <- analysis_data[complete.cases(analysis_data[, c("dependency", "lineage")]), ]
lineage_counts <- table(lineage_data$lineage)
eligible_lineages <- names(lineage_counts[lineage_counts >= 10L])
lineage_data <- lineage_data[lineage %in% eligible_lineages]
if (length(eligible_lineages) < 2L) stop("Fewer than two lineages have at least 10 models.")

global_lineage <- kruskal.test(dependency ~ lineage, data = lineage_data)
lineage_summary <- rbindlist(lapply(eligible_lineages, function(group) {
  inside <- lineage_data[lineage == group, dependency]
  outside <- lineage_data[lineage != group, dependency]
  test <- suppressWarnings(wilcox.test(inside, outside, exact = FALSE))
  data.frame(
    lineage = group, n = length(inside), mean = mean(inside), median = median(inside),
    rest_n = length(outside), rest_median = median(outside),
    median_difference = median(inside) - median(outside), p_value = test$p.value
  )
}))
lineage_summary[, fdr := p.adjust(p_value, method = "BH")]
setorder(lineage_summary, median_difference)

fwrite(dependency_summary, file.path(table_dir, "04_dependency_summary.csv"))
fwrite(
  data.frame(test = "Kruskal-Wallis", n = nrow(lineage_data), groups = length(eligible_lineages),
             statistic = unname(global_lineage$statistic), p_value = global_lineage$p.value),
  file.path(table_dir, "05_lineage_global_test.csv")
)
fwrite(lineage_summary, file.path(table_dir, "06_lineage_associations.csv"))
saveRDS(lineage_data, file.path(intermediate_dir, "lineage_data.rds"))

message("Tested ", length(eligible_lineages), " lineages with at least 10 models.")
