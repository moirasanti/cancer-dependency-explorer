## DESCRIPTION ################################################################
## GENERATE STANDARDIZED FIGURES
# Produce the same compact figure set for every requested target.

# Define the complete figure palette in one place.
plot.colours <- c(
  primary = "#28666e",
  primary.light = "#a9ced1",
  primary.dark = "#173d42",
  accent = "#e98389",
  secondary = "#8ab0ab",
  neutral = "#8a8a8a",
  threshold = "grey35",
  baseline = "grey45",
  reference = "grey50",
  background = "white"
)

## LOAD DATA ###################################################################
# Load the prepared data and results used by the six figures.
analysis.data <- readRDS(file.path(intermediate.dir, "analysis_data.rds"))
lineage.data <- readRDS(file.path(intermediate.dir, "lineage_data.rds"))
mutation.results <- readRDS(file.path(intermediate.dir, "mutation_associations.rds"))
copy.number.results <- readRDS(file.path(intermediate.dir, "copy_number_associations.rds"))
model.result <- readRDS(file.path(intermediate.dir, "multivariable_model.rds"))

## MAIN ########################################################################
# Plot the pan-cancer dependency distribution.
distribution.plot <- ggplot(analysis.data[!is.na(dependency)], aes(dependency)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 35,
    fill = plot.colours[["primary"]],
    color = plot.colours[["background"]],
    alpha = 0.75
  ) +
  geom_density(color = plot.colours[["accent"]], linewidth = 0.8) +
  geom_vline(
    xintercept = c(-0.5, 0),
    linetype = c("dotted", "dashed"),
    color = plot.colours[["threshold"]]
  ) +
  labs(
    title = paste(target, "dependency across DepMap models"),
    subtitle = "More negative Chronos scores indicate stronger dependency",
    x = paste(target, "Chronos gene-effect score"),
    y = "Density"
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "01_dependency_distribution.png"),
  distribution.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Order lineages by median dependency and estimate one density curve per lineage.
lineage.order <- lineage.data[
  , .(median.dependency = median(dependency)), by = lineage
][order(median.dependency), lineage]
lineage.data$lineage <- factor(lineage.data$lineage, levels = lineage.order)
lineage.palette <- setNames(
  grDevices::colorRampPalette(c(
    plot.colours[["primary.dark"]],
    plot.colours[["primary"]],
    plot.colours[["primary.light"]]
  ))(length(lineage.order)),
  lineage.order
)

dependency.range <- range(lineage.data$dependency, na.rm = TRUE)
ridge.data.list <- vector("list", length(lineage.order))

for (lineage.index in seq_along(lineage.order)) {
  lineage.name <- lineage.order[[lineage.index]]
  lineage.scores <- lineage.data$dependency[lineage.data$lineage == lineage.name]
  lineage.density <- density(
    lineage.scores,
    from = dependency.range[[1L]],
    to = dependency.range[[2L]],
    n = 256,
    na.rm = TRUE
  )
  lineage.density$y[c(1L, length(lineage.density$y))] <- 0
  ridge.baseline <- lineage.index

  ridge.data.list[[lineage.index]] <- data.frame(
    dependency = lineage.density$x,
    baseline = ridge.baseline,
    ridge.height = ridge.baseline +
      1.15 * lineage.density$y / max(lineage.density$y),
    lineage = lineage.name
  )
}
ridge.data <- rbindlist(ridge.data.list)

# Overlap translucent density ridges to compare lineage-specific distributions.
lineage.plot <- ggplot(
  ridge.data,
  aes(
    x = dependency,
    ymin = baseline,
    ymax = ridge.height,
    group = lineage,
    fill = lineage
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = plot.colours[["baseline"]]
  ) +
  geom_ribbon(
    alpha = 0.62,
    color = plot.colours[["background"]],
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = lineage.palette,
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq_along(lineage.order),
    labels = lineage.order,
    expand = expansion(mult = c(0.01, 0.06))
  ) +
  labs(
    title = paste(target, "dependency by cancer lineage"),
    subtitle = "Stronger median dependencies are shown at the bottom",
    x = "Chronos score",
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "02_lineage_dependency.png"),
  lineage.plot,
  width = 7,
  height = max(6, length(lineage.order) * 0.28),
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Show the target expression-dependency relationship and fitted linear trend.
expression.plot.data <- analysis.data[
  complete.cases(analysis.data[, c("target_expression", "dependency")]),
]
expression.plot <- ggplot(
  expression.plot.data,
  aes(target_expression, dependency)
) +
  geom_point(
    alpha = 0.35,
    size = 1.4,
    color = plot.colours[["primary"]]
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = plot.colours[["accent"]]
  ) +
  labs(
    title = paste(target, "expression and dependency"),
    x = paste0(target, " expression (log2 TPM + 1)"),
    y = "Chronos score"
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "03_expression_dependency.png"),
  expression.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Display the leading damaging-mutation effects and confidence intervals.
mutation.plot.data <- head(
  mutation.results[mutation.results$eligible, ],
  12L
)
mutation.plot.data$gene <- factor(
  mutation.plot.data$gene,
  levels = rev(mutation.plot.data$gene)
)
mutation.plot <- ggplot(mutation.plot.data, aes(effect, gene)) +
  geom_vline(
    xintercept = 0,
    color = plot.colours[["reference"]],
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    height = 0.2,
    color = plot.colours[["primary"]]
  ) +
  geom_point(aes(color = fdr < 0.05), size = 2.3) +
  scale_color_manual(
    values = c(
      `TRUE` = plot.colours[["accent"]],
      `FALSE` = plot.colours[["neutral"]]
    ),
    guide = "none"
  ) +
  labs(
    title = "Top damaging-mutation associations",
    subtitle = "Negative effects indicate stronger dependency in altered models",
    x = "Mean Chronos difference: altered minus reference",
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "04_mutation_associations.png"),
  mutation.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Display the leading genome-wide copy-number correlations.
copy.number.plot.data <- head(
  copy.number.results[copy.number.results$eligible, ],
  12L
)
copy.number.plot.data$gene <- factor(
  copy.number.plot.data$gene,
  levels = rev(copy.number.plot.data$gene)
)
copy.number.significant <- copy.number.plot.data[
  !is.na(copy.number.plot.data$fdr) & copy.number.plot.data$fdr < 0.05,
]
copy.number.not.significant <- copy.number.plot.data[
  is.na(copy.number.plot.data$fdr) | copy.number.plot.data$fdr >= 0.05,
]
copy.number.plot <- ggplot(copy.number.plot.data, aes(correlation, gene)) +
  geom_vline(
    xintercept = 0,
    color = plot.colours[["reference"]],
    linetype = "dashed"
  ) +
  geom_segment(
    data = copy.number.significant,
    aes(x = 0, xend = correlation, yend = gene),
    color = plot.colours[["secondary"]]
  ) +
  geom_segment(
    data = copy.number.not.significant,
    aes(x = 0, xend = correlation, yend = gene),
    color = plot.colours[["neutral"]]
  ) +
  geom_point(
    data = copy.number.significant,
    color = plot.colours[["accent"]],
    size = 2.4
  ) +
  geom_point(
    data = copy.number.not.significant,
    color = plot.colours[["neutral"]],
    size = 2.4
  ) +
  labs(
    title = "Top copy-number associations",
    subtitle = "Pearson correlations used for the genome-wide screen",
    x = "Correlation with Chronos dependency",
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "05_copy_number_associations.png"),
  copy.number.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Show adjusted molecular coefficients without the intercept or lineage contrasts.
coefficient.plot.data <- model.result$coefficients[
  model.result$coefficients$term != "(Intercept)" &
    !grepl("^lineage_model", model.result$coefficients$term),
  , drop = FALSE
]
coefficient.plot.data <- coefficient.plot.data[
  order(
    is.na(coefficient.plot.data$fdr),
    coefficient.plot.data$fdr,
    coefficient.plot.data$p_value
  ),
]
coefficient.plot.data$plot.label <- sub(
  " \\(for a 1 SD increase\\)$",
  "\n(1 SD higher)",
  coefficient.plot.data$label
)
coefficient.plot.data$plot.label <- factor(
  coefficient.plot.data$plot.label,
  levels = rev(coefficient.plot.data$plot.label)
)
coefficient.plot.significant <- coefficient.plot.data[
  !is.na(coefficient.plot.data$fdr) & coefficient.plot.data$fdr < 0.05,
]
coefficient.plot.not.significant <- coefficient.plot.data[
  is.na(coefficient.plot.data$fdr) | coefficient.plot.data$fdr >= 0.05,
]
coefficient.plot <- ggplot(coefficient.plot.data, aes(estimate, plot.label)) +
  geom_vline(
    xintercept = 0,
    color = plot.colours[["reference"]],
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    height = 0.2,
    color = plot.colours[["primary"]]
  ) +
  geom_point(
    data = coefficient.plot.significant,
    size = 2.5,
    color = plot.colours[["accent"]]
  ) +
  geom_point(
    data = coefficient.plot.not.significant,
    size = 2.5,
    color = plot.colours[["neutral"]]
  ) +
  labs(
    title = paste("What predicts", target, "dependency after adjustment?"),
    subtitle = "Ordered by FDR; pink points pass the 5% FDR threshold",
    x = paste0(
      "Adjusted change in ", target, " Chronos score\n",
      "← Stronger dependency                 Weaker dependency →"
    ),
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "06_multivariable_coefficients.png"),
  coefficient.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

message("Wrote six standardized figures.")
