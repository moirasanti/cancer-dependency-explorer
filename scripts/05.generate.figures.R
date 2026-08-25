# Produce the same compact figure set for every requested target.

analysis_data <- readRDS(file.path(intermediate_dir, "analysis_data.rds"))
lineage_data <- readRDS(file.path(intermediate_dir, "lineage_data.rds"))
mutation_results <- readRDS(file.path(intermediate_dir, "mutation_associations.rds"))
copy_number_results <- readRDS(file.path(intermediate_dir, "copy_number_associations.rds"))
model_result <- readRDS(file.path(intermediate_dir, "multivariable_model.rds"))

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(file.path(figure_dir, filename), plot, width = width, height = height, dpi = 300, bg = "white")
}

distribution_plot <- ggplot(analysis_data[!is.na(dependency)], aes(dependency)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "#28666e", color = "white") +
  geom_density(color = "#7c2e41", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 0), linetype = c("dotted", "dashed"), color = "grey35") +
  labs(
    title = paste(target, "dependency across DepMap models"),
    subtitle = "More-negative Chronos scores indicate stronger dependency",
    x = paste(target, "Chronos gene-effect score"), y = "Density"
  ) + theme_dependency()
save_plot(distribution_plot, "01_dependency_distribution.png")

lineage_order <- lineage_data[, .(median_dependency = median(dependency)), by = lineage][order(median_dependency), lineage]
lineage_data$lineage <- factor(lineage_data$lineage, levels = lineage_order)
lineage_plot <- ggplot(lineage_data, aes(lineage, dependency)) +
  geom_boxplot(fill = "#8ab0ab", outlier.alpha = 0.35, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  coord_flip() +
  labs(title = paste(target, "dependency by cancer lineage"), x = NULL, y = "Chronos score") +
  theme_dependency()
save_plot(lineage_plot, "02_lineage_dependency.png", height = max(5, length(lineage_order) * 0.24))

expression_plot <- ggplot(
  analysis_data[complete.cases(analysis_data[, c("target_expression", "dependency")])],
  aes(target_expression, dependency)
) +
  geom_point(alpha = 0.35, size = 1.4, color = "#28666e") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#7c2e41") +
  labs(
    title = paste(target, "expression and dependency"),
    x = paste0(target, " expression (log2 TPM + 1)"), y = "Chronos score"
  ) + theme_dependency()
save_plot(expression_plot, "03_expression_dependency.png")

mutation_plot_data <- head(mutation_results[mutation_results$eligible, ], 12L)
mutation_plot_data$gene <- factor(mutation_plot_data$gene, levels = rev(mutation_plot_data$gene))
mutation_plot <- ggplot(mutation_plot_data, aes(effect, gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2, color = "#28666e") +
  geom_point(aes(color = fdr < 0.05), size = 2.3) +
  scale_color_manual(values = c(`TRUE` = "#7c2e41", `FALSE` = "#8a8a8a"), guide = "none") +
  labs(
    title = "Top damaging-mutation associations",
    subtitle = "Negative effects indicate stronger dependency in altered models",
    x = "Mean Chronos difference: altered minus reference", y = NULL
  ) + theme_dependency()
save_plot(mutation_plot, "04_mutation_associations.png")

copy_plot_data <- head(copy_number_results[copy_number_results$eligible, ], 12L)
copy_plot_data$gene <- factor(copy_plot_data$gene, levels = rev(copy_plot_data$gene))
copy_number_plot <- ggplot(copy_plot_data, aes(correlation, gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = correlation, yend = gene), color = "#8ab0ab") +
  geom_point(aes(color = fdr < 0.05), size = 2.4) +
  scale_color_manual(values = c(`TRUE` = "#7c2e41", `FALSE` = "#8a8a8a"), guide = "none") +
  labs(
    title = "Top copy-number associations",
    subtitle = "Pearson correlations are used for the scalable genome-wide screen",
    x = "Correlation with Chronos dependency", y = NULL
  ) + theme_dependency()
save_plot(copy_number_plot, "05_copy_number_associations.png")

coefficient_data <- model_result$coefficients[
  model_result$coefficients$term != "(Intercept)" &
    !grepl("^lineage_model", model_result$coefficients$term), , drop = FALSE
]
coefficient_data$label <- factor(coefficient_data$label, levels = rev(coefficient_data$label))
coefficient_plot <- ggplot(coefficient_data, aes(estimate, label)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2, color = "#28666e") +
  geom_point(size = 2.5, color = "#7c2e41") +
  labs(
    title = paste(target, "multivariable molecular associations"),
    subtitle = "Coefficients with HC3 robust 95% confidence intervals",
    x = "Change in Chronos dependency", y = NULL
  ) + theme_dependency()
save_plot(coefficient_plot, "06_multivariable_coefficients.png")

message("Wrote six standardized figures.")
