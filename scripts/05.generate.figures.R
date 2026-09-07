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
expression.results <- readRDS(file.path(intermediate.dir, "expression_associations.rds"))
mutation.results <- readRDS(file.path(intermediate.dir, "mutation_associations.rds"))
copy.number.results <- readRDS(file.path(intermediate.dir, "copy_number_associations.rds"))
model.result <- readRDS(file.path(intermediate.dir, "integrated_model.rds"))

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

# Show the target expression-dependency relationship and its linear trend.
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
    subtitle = "Line and band show the unadjusted linear fit and 95% confidence interval",
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

# Display the leading genome-wide expression correlations.
genome.expression.plot.data <- head(
  expression.results[
    expression.results$eligible,
  ],
  12L
)
genome.expression.plot.data$gene <- factor(
  genome.expression.plot.data$gene,
  levels = rev(genome.expression.plot.data$gene)
)
genome.expression.plot.data$notable <-
  !is.na(genome.expression.plot.data$fdr) &
  genome.expression.plot.data$fdr < 0.05 &
  abs(genome.expression.plot.data$spearman_rho) >= 0.30
genome.expression.plot <- ggplot(
  genome.expression.plot.data,
  aes(spearman_rho, gene)
) +
  geom_vline(
    xintercept = 0,
    color = plot.colours[["reference"]],
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = 0, xend = spearman_rho, yend = gene, color = notable)
  ) +
  geom_point(aes(color = notable), size = 2.4) +
  scale_color_manual(
    values = c(
      `TRUE` = plot.colours[["accent"]],
      `FALSE` = plot.colours[["neutral"]]
    ),
    guide = "none"
  ) +
  labs(
    title = "Top genome-wide expression associations",
    subtitle = "Pink points pass 5% FDR and |Spearman rho| ≥ 0.30",
    x = paste("Spearman correlation with", target, "Chronos dependency"),
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "04_expression_associations.png"),
  genome.expression.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Display the leading damaging-mutation effects and confidence intervals.
mutation.plot.data <- head(
  mutation.results[
    mutation.results$eligible,
  ],
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
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    orientation = "y",
    width = 0.2,
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
    subtitle = "Points show mean differences and bars show Welch 95% confidence intervals",
    x = paste(
      "Mean", target, "Chronos difference: altered minus wild type"
    ),
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "05_mutation_associations.png"),
  mutation.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Display the leading genome-wide copy-number correlations.
copy.number.plot.data <- head(
  copy.number.results[
    copy.number.results$eligible,
  ],
  12L
)
copy.number.plot.data$gene <- factor(
  copy.number.plot.data$gene,
  levels = rev(copy.number.plot.data$gene)
)
copy.number.plot.data$notable <-
  !is.na(copy.number.plot.data$fdr) &
  copy.number.plot.data$fdr < 0.05 &
  abs(copy.number.plot.data$correlation) >= 0.30
copy.number.plot <- ggplot(copy.number.plot.data, aes(correlation, gene)) +
  geom_vline(
    xintercept = 0,
    color = plot.colours[["reference"]],
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = 0, xend = correlation, yend = gene, color = notable)
  ) +
  geom_point(aes(color = notable), size = 2.4) +
  scale_color_manual(
    values = c(
      `TRUE` = plot.colours[["accent"]],
      `FALSE` = plot.colours[["neutral"]]
    ),
    guide = "none"
  ) +
  labs(
    title = "Top copy number associations",
    subtitle = "Pink points pass 5% FDR and |Pearson r| ≥ 0.30",
    x = paste("Pearson correlation with", target, "Chronos dependency"),
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "06_copy_number_associations.png"),
  copy.number.plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Show only stable molecular predictors from the single full-data fit. Their
# biological source changes the visual annotation, not their model penalty.
coefficient.plot.data <- model.result$coefficients
if (nrow(coefficient.plot.data)) {
  coefficient.plot.data <- coefficient.plot.data[
    order(abs(coefficient.plot.data$coefficient), decreasing = TRUE),
  ]
  coefficient.plot.data$plot.label <- factor(
    coefficient.plot.data$label,
    levels = rev(coefficient.plot.data$label)
  )
  coefficient.plot.data$predictor.source <- unname(c(
    target_characterization = "Target characterisation",
    prespecified_hypothesis = "Pre-specified hypothesis",
    genome_wide_candidate = "Genome-wide candidate"
  )[coefficient.plot.data$role])
  coefficient.plot <- ggplot(
    coefficient.plot.data,
    aes(
      coefficient, plot.label,
      colour = predictor.source, shape = predictor.source
    )
  ) +
    geom_vline(
      xintercept = 0,
      color = plot.colours[["reference"]],
      linetype = "dashed"
    ) +
    geom_segment(
      aes(x = 0, xend = coefficient, yend = plot.label),
      linewidth = 0.45
    ) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = c(
      "Target characterisation" = plot.colours[["primary"]],
      "Pre-specified hypothesis" = plot.colours[["accent"]],
      "Genome-wide candidate" = plot.colours[["secondary"]]
    )) +
    labs(
      title = paste("Stable integrated predictors of", target, "dependency"),
      subtitle = "Standardized elastic-net coefficients; biological source does not alter penalization",
      x = paste0(
        "Elastic-net coefficient for ", target, " Chronos score\n",
        "← Stronger dependency                 Weaker dependency →"
      ),
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    dependency_plot_theme()
  ggsave(
    file.path(figure.dir, "07_integrated_model_coefficients.png"),
    coefficient.plot,
    width = 7,
    height = max(5, nrow(coefficient.plot.data) * 0.38 + 1.5),
    dpi = 300,
    bg = plot.colours[["background"]]
  )
} else {
  unlink(file.path(figure.dir, "07_integrated_model_coefficients.png"))
}
unlink(file.path(figure.dir, "07_multivariable_coefficients.png"))

# Display every held-out prediction without averaging across repeats. Faceting
# makes the repeat-to-repeat consistency visible while the diagonal shows
# perfect calibration.
prediction.plot.data <- model.result$predictions
repeat.pearson.r <- vapply(
  split(prediction.plot.data, prediction.plot.data$repeat_id),
  function(repeat.data) cor(
    repeat.data$observed, repeat.data$predicted,
    use = "complete.obs", method = "pearson"
  ),
  numeric(1)
)
repeat.labels <- sprintf(
  "Repeat %s\nHeld-out Pearson r = %.3f",
  names(repeat.pearson.r), repeat.pearson.r
)
names(repeat.labels) <- names(repeat.pearson.r)
prediction.plot.data$repeat_label <- factor(
  repeat.labels[as.character(prediction.plot.data$repeat_id)],
  levels = repeat.labels
)
prediction.range <- range(
  c(prediction.plot.data$observed, prediction.plot.data$predicted),
  finite = TRUE
)
prediction.plot <- ggplot(
  prediction.plot.data,
  aes(observed, predicted)
) +
  geom_abline(
    slope = 1, intercept = 0,
    colour = plot.colours[["reference"]], linetype = "dashed"
  ) +
  geom_point(
    colour = plot.colours[["primary"]], alpha = 0.35, size = 0.9
  ) +
  facet_wrap(~ repeat_label, nrow = 1L) +
  coord_equal(xlim = prediction.range, ylim = prediction.range) +
  labs(
    title = paste("Held-out predictions of", target, "dependency"),
    subtitle = "Each point is one cell line. The dashed diagonal represents a perfect prediction",
    x = "Observed Chronos score",
    y = "Predicted Chronos score"
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "08_observed_vs_predicted.png"),
  prediction.plot,
  width = 9,
  height = 3.8,
  dpi = 300,
  bg = plot.colours[["background"]]
)

# Compare within-lineage predictive performance. Error bars describe variation
# across the three outer-CV repeats, not coefficient uncertainty.
lineage.performance.plot.data <- model.result$lineage_performance
lineage.performance.plot.data$lineage <- factor(
  lineage.performance.plot.data$lineage,
  levels = rev(lineage.performance.plot.data$lineage)
)
lineage.performance.plot <- ggplot(
  lineage.performance.plot.data,
  aes(r_squared_mean, lineage)
) +
  geom_vline(
    xintercept = 0,
    colour = plot.colours[["reference"]], linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = r_squared_mean - r_squared_sd,
      xmax = r_squared_mean + r_squared_sd
    ),
    orientation = "y",
    width = 0.2,
    colour = plot.colours[["secondary"]]
  ) +
  geom_point(
    aes(size = n, colour = r_squared_mean > 0),
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = c(
      "TRUE" = plot.colours[["primary"]],
      "FALSE" = plot.colours[["neutral"]]
    ),
    guide = "none"
  ) +
  scale_size_continuous(name = "Cell lines", range = c(2, 5)) +
  labs(
    title = paste("Held-out predictive performance by cancer lineage"),
    subtitle = "Within-lineage R² averaged across three repeats; bars show ±1 SD",
    x = "Held-out R² (0 = no improvement over the lineage mean)",
    y = NULL
  ) +
  dependency_plot_theme()
ggsave(
  file.path(figure.dir, "09_lineage_cv_performance.png"),
  lineage.performance.plot,
  width = 8,
  height = max(5, nrow(lineage.performance.plot.data) * 0.28 + 1.8),
  dpi = 300,
  bg = plot.colours[["background"]]
)

message("Wrote standardized figures.")
