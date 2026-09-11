suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Remove Entrez identifiers and duplicate-column suffixes from DepMap headers.
extract_gene_symbol <- function(x) {
  x <- sub("\\s*\\([^)]*\\)$", "", x)
  sub("\\.\\.[0-9]+\\.?$", "", x)
}

# Read selected genes without loading an entire genome-wide matrix into memory.
# Duplicate mutation columns are combined as any-altered; duplicate continuous
# columns are averaged across the available measurements for each model.
read_depmap_gene_data <- function(
  file,
  genes,
  default.only = FALSE,
  modality = c("continuous", "mutation")
) {
  modality <- match.arg(modality)
  file.header <- names(fread(file, nrows = 0L, check.names = FALSE))
  header.symbols <- extract_gene_symbol(file.header)
  gene.column.index <- which(header.symbols %in% genes)
  found.symbols <- unique(header.symbols[gene.column.index])
  missing.genes <- setdiff(genes, found.symbols)

  if (length(missing.genes)) {
    stop("Genes absent from ", basename(file), ": ", paste(missing.genes, collapse = ", "))
  }

  metadata.column.index <- which(
    file.header %in% c("ModelID", "IsDefaultEntryForModel")
  )
  if (!"ModelID" %in% file.header) {
    metadata.column.index <- unique(c(1L, metadata.column.index))
  }

  gene.data <- fread(
    file,
    select = unique(c(metadata.column.index, gene.column.index)),
    check.names = FALSE
  )
  if (!"ModelID" %in% names(gene.data)) setnames(gene.data, 1L, "ModelID")

  if (default.only && "IsDefaultEntryForModel" %in% names(gene.data)) {
    default.rows <- as.character(gene.data$IsDefaultEntryForModel) %in% c("Yes", "TRUE", "1")
    gene.data <- gene.data[default.rows]
  }

  # Build one output column per requested symbol. When several source columns
  # map to the same symbol, retain all observed information rather than silently
  # keeping the first occurrence.
  selected.data <- data.table(ModelID = gene.data$ModelID)
  selected.symbols <- extract_gene_symbol(names(gene.data))
  for (gene in genes) {
    matching.columns <- which(selected.symbols == gene)
    gene.values <- as.matrix(gene.data[, matching.columns, with = FALSE])
    storage.mode(gene.values) <- "double"

    if (ncol(gene.values) == 1L) {
      collapsed.values <- gene.values[, 1L]
    } else if (modality == "mutation") {
      observed.values <- rowSums(!is.na(gene.values))
      altered.values <- rowSums(gene.values > 0, na.rm = TRUE)
      reference.values <- rowSums(gene.values <= 0, na.rm = TRUE)
      collapsed.values <- as.numeric(altered.values > 0L)
      collapsed.values[observed.values == 0L] <- NA_real_
      inconsistent.rows <-
        observed.values > 1L & altered.values > 0L & reference.values > 0L
      if (any(inconsistent.rows)) {
        message(
          "Duplicate mutation values disagreed for ", gene, " in ",
          sum(inconsistent.rows),
          " model(s); classified as altered when any value was > 0."
        )
      }
    } else {
      observed.values <- rowSums(!is.na(gene.values))
      collapsed.values <- rowMeans(gene.values, na.rm = TRUE)
      collapsed.values[observed.values == 0L] <- NA_real_
      maximum.input <- gene.values
      minimum.input <- gene.values
      maximum.input[is.na(maximum.input)] <- -Inf
      minimum.input[is.na(minimum.input)] <- Inf
      observed.range <-
        apply(maximum.input, 1L, max) - apply(minimum.input, 1L, min)
      inconsistent.rows <-
        observed.values > 1L & observed.range > 1e-8
      if (any(inconsistent.rows)) {
        message(
          "Duplicate continuous values disagreed for ", gene, " in ",
          sum(inconsistent.rows),
          " model(s); used the mean of available values."
        )
      }
    }
    selected.data[[gene]] <- collapsed.values
  }

  unique(selected.data, by = "ModelID")
}

# Read and align a genome-wide matrix while resolving duplicate gene symbols.
# Mutation duplicates use any-altered logic. Continuous copy-number duplicates
# use the row mean, which does not preferentially retain gains over losses.
read_aligned_genomic_matrix <- function(
  file,
  analysis.ids,
  modality = c("copy_number", "expression", "mutation"),
  comparison.tolerance = 1e-8
) {
  modality <- match.arg(modality)
  genomic.data <- fread(file, check.names = FALSE, showProgress = interactive())
  if (!"ModelID" %in% names(genomic.data)) setnames(genomic.data, 1L, "ModelID")

  if ("IsDefaultEntryForModel" %in% names(genomic.data)) {
    default.rows <- as.character(genomic.data$IsDefaultEntryForModel) %in% c("Yes", "TRUE", "1")
    genomic.data <- genomic.data[default.rows]
  }
  genomic.data <- unique(genomic.data, by = "ModelID")

  row.index <- match(analysis.ids, genomic.data$ModelID)
  metadata.columns <- c(
    "V1", "ModelID", "SequencingID", "ModelConditionID",
    "IsDefaultEntryForModel", "IsDefaultEntryForMC"
  )
  gene.column.index <- which(!names(genomic.data) %in% metadata.columns)
  gene.columns <- names(genomic.data)[gene.column.index]
  gene.symbols <- extract_gene_symbol(gene.columns)

  genomic.matrix <- as.matrix(
    genomic.data[row.index, gene.column.index, with = FALSE]
  )
  storage.mode(genomic.matrix) <- "double"

  # Resolve every duplicated symbol at model level. Stay silent when duplicate
  # values agree; report only conflicts and the rule used to resolve them.
  duplicated.genes <- unique(gene.symbols[duplicated(gene.symbols)])
  if (length(duplicated.genes)) {
    for (duplicate.index in seq_along(duplicated.genes)) {
      duplicate.gene <- duplicated.genes[[duplicate.index]]
      source.index <- which(gene.symbols == duplicate.gene)
      duplicate.values <- genomic.matrix[, source.index, drop = FALSE]
      observed.values <- rowSums(!is.na(duplicate.values))

      if (modality == "mutation") {
        altered.values <- rowSums(duplicate.values > 0, na.rm = TRUE)
        reference.values <- rowSums(duplicate.values <= 0, na.rm = TRUE)
        collapsed.values <- as.numeric(altered.values > 0L)
        collapsed.values[observed.values == 0L] <- NA_real_
        inconsistent.rows <-
          observed.values > 1L & altered.values > 0L & reference.values > 0L
        if (any(inconsistent.rows)) {
          message(
            "Duplicate mutation values disagreed for ", duplicate.gene, " in ",
            sum(inconsistent.rows),
            " model(s); classified as altered when any value was > 0."
          )
        }
      } else {
        collapsed.values <- rowMeans(duplicate.values, na.rm = TRUE)
        collapsed.values[observed.values == 0L] <- NA_real_
        maximum.input <- duplicate.values
        minimum.input <- duplicate.values
        maximum.input[is.na(maximum.input)] <- -Inf
        minimum.input[is.na(minimum.input)] <- Inf
        maximum.observed <- apply(maximum.input, 1L, max)
        minimum.observed <- apply(minimum.input, 1L, min)
        observed.range <- maximum.observed - minimum.observed
        inconsistent.rows <-
          observed.values > 1L & observed.range > comparison.tolerance
        if (any(inconsistent.rows)) {
          message(
            "Duplicate ", gsub("_", "-", modality),
            " values disagreed for ", duplicate.gene,
            " in ", sum(inconsistent.rows),
            " model(s); used the mean of available values."
          )
        }
      }

      # Replace the first occurrence with the resolved value; later duplicate
      # columns are removed after every duplicate symbol has been checked.
      genomic.matrix[, source.index[[1L]]] <- collapsed.values
    }
  }

  unique.columns <- !duplicated(gene.symbols)
  resolved.matrix <- genomic.matrix[, unique.columns, drop = FALSE]
  resolved.symbols <- gene.symbols[unique.columns]
  colnames(resolved.matrix) <- resolved.symbols

  list(
    matrix = resolved.matrix,
    symbols = resolved.symbols
  )
}

# Return a consistent result when a Spearman test cannot be estimated safely.
calculate_safe_spearman <- function(x, y) {
  complete.rows <- complete.cases(x, y)
  if (
    sum(complete.rows) < 3L ||
    length(unique(x[complete.rows])) < 2L ||
    length(unique(y[complete.rows])) < 2L
  ) {
    return(c(n = sum(complete.rows), estimate = NA_real_, p_value = NA_real_))
  }

  correlation.test <- suppressWarnings(
    cor.test(x[complete.rows], y[complete.rows], method = "spearman", exact = FALSE)
  )
  c(
    n = sum(complete.rows),
    estimate = unname(correlation.test$estimate),
    p_value = correlation.test$p.value
  )
}

# Compare altered and wild-type dependency scores with a built-in Welch test.
# Return an explicitly ineligible result when either group is too small or the
# test cannot be estimated, for example because both groups are constant.
calculate_safe_welch <- function(
  altered_scores,
  wild_type_scores,
  minimum_group_size = 3L
) {
  altered_scores <- altered_scores[is.finite(altered_scores)]
  wild_type_scores <- wild_type_scores[is.finite(wild_type_scores)]
  if (
    length(altered_scores) < minimum_group_size ||
      length(wild_type_scores) < minimum_group_size
  ) {
    return(c(
      standard_error = NA_real_, conf_low = NA_real_,
      conf_high = NA_real_, p_value = NA_real_, eligible = 0
    ))
  }

  welch.test <- tryCatch(
    t.test(altered_scores, wild_type_scores, var.equal = FALSE),
    error = function(condition) NULL
  )
  if (is.null(welch.test)) {
    return(c(
      standard_error = NA_real_, conf_low = NA_real_,
      conf_high = NA_real_, p_value = NA_real_, eligible = 0
    ))
  }

  c(
    standard_error = unname(welch.test$stderr),
    conf_low = unname(welch.test$conf.int[[1L]]),
    conf_high = unname(welch.test$conf.int[[2L]]),
    p_value = welch.test$p.value,
    eligible = 1
  )
}

# Screen a numeric feature matrix against one outcome using bounded blocks.
# Spearman coefficients are Pearson correlations of within-feature ranks. The
# P values use the usual large-sample t approximation, avoiding thousands of
# separate cor.test objects while preserving pairwise missing-value handling.
calculate_block_spearman <- function(feature.matrix, outcome, block.size = 500L) {
  feature.matrix <- as.matrix(feature.matrix)
  storage.mode(feature.matrix) <- "double"
  if (nrow(feature.matrix) != length(outcome)) {
    stop("Feature matrix rows must match the outcome length.")
  }
  if (is.null(colnames(feature.matrix))) {
    colnames(feature.matrix) <- paste0("feature_", seq_len(ncol(feature.matrix)))
  }
  if (!is.numeric(block.size) || length(block.size) != 1L || block.size < 1L) {
    stop("block.size must be a positive integer.")
  }
  block.size <- as.integer(block.size)

  feature.count <- ncol(feature.matrix)
  observed.n <- integer(feature.count)
  feature.sd <- correlation <- p.value <- rep(NA_real_, feature.count)

  for (block.start in seq.int(1L, feature.count, by = block.size)) {
    block.end <- min(feature.count, block.start + block.size - 1L)
    for (feature.index in block.start:block.end) {
      feature <- feature.matrix[, feature.index]
      complete.rows <- complete.cases(feature, outcome)
      observed.n[[feature.index]] <- sum(complete.rows)
      if (observed.n[[feature.index]] < 3L) next

      observed.feature <- feature[complete.rows]
      observed.outcome <- outcome[complete.rows]
      feature.sd[[feature.index]] <- sd(observed.feature)
      if (
        !is.finite(feature.sd[[feature.index]]) ||
        feature.sd[[feature.index]] == 0 ||
        length(unique(observed.outcome)) < 2L
      ) {
        next
      }

      rho <- cor(
        rank(observed.feature, ties.method = "average"),
        rank(observed.outcome, ties.method = "average"),
        method = "pearson"
      )
      correlation[[feature.index]] <- rho
      if (!is.finite(rho) || abs(rho) >= 1) next

      statistic <- rho * sqrt(
        (observed.n[[feature.index]] - 2) / pmax(1 - rho^2, 1e-12)
      )
      p.value[[feature.index]] <- 2 * pt(
        abs(statistic),
        df = observed.n[[feature.index]] - 2,
        lower.tail = FALSE
      )
    }
  }

  data.frame(
    gene = colnames(feature.matrix),
    n = observed.n,
    sd = feature.sd,
    spearman_rho = correlation,
    p_value = p.value,
    stringsAsFactors = FALSE
  )
}

# Read and normalize optional modality-specific pre-specified hypotheses.
read_biomarker_configuration <- function(file, target) {
  configured <- list(
    expression = character(),
    mutation = character(),
    copy_number = character()
  )
  if (!nzchar(file)) return(configured)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for biomarker configuration.")
  }

  biomarker.config <- yaml::read_yaml(file)
  configured.target <- biomarker.config$target
  if (
    !is.null(configured.target) &&
    toupper(trimws(as.character(configured.target))) != target
  ) {
    stop("Biomarker configuration target is ", configured.target, ", not ", target, ".")
  }

  for (modality in names(configured)) {
    values <- biomarker.config[[modality]]
    if (is.null(values)) next
    values <- toupper(trimws(as.character(unlist(values))))
    configured[[modality]] <- unique(values[nzchar(values)])
  }
  configured
}

# Remove configured entries already represented by target-characterisation features.
deduplicate_target_hypotheses <- function(configured, target) {
  for (modality in names(configured)) {
    if (!target %in% configured[[modality]]) next
    message(
      "Configured ", modality, " hypothesis ", target,
      " duplicates a core target feature and was de-duplicated."
    )
    configured[[modality]] <- setdiff(configured[[modality]], target)
  }
  configured
}

# Validate configured genes without applying marginal association thresholds.
validate_configured_features <- function(
  genes,
  results,
  modality,
  eligibility.description
) {
  if (!length(genes)) return(genes)
  absent <- setdiff(genes, results$gene)
  if (length(absent)) {
    stop(
      "Configured ", modality, " genes absent from the matrix: ",
      paste(absent, collapse = ", ")
    )
  }
  ineligible <- genes[!results$eligible[match(genes, results$gene)]]
  if (length(ineligible)) {
    warning(
      "Pre-specified ", modality, " hypotheses failed ",
      eligibility.description, " filters and were omitted: ",
      paste(ineligible, collapse = ", ")
    )
  }
  setdiff(genes, ineligible)
}

# Build consistently typed feature-family rows, including an empty feature list.
make_feature_family_rows <- function(family, modality, genes) {
  genes <- as.character(genes)
  data.frame(
    family = rep(as.character(family), length.out = length(genes)),
    modality = rep(as.character(modality), length.out = length(genes)),
    gene = genes,
    stringsAsFactors = FALSE
  )
}

# Apply the report's descriptive absolute-correlation categories.
classify_correlation_strength <- function(correlation) {
  strength <- rep(NA_character_, length(correlation))
  strength[is.finite(correlation) & abs(correlation) < 0.30] <- "small"
  strength[is.finite(correlation) & abs(correlation) >= 0.30] <- "moderate"
  strength[is.finite(correlation) & abs(correlation) >= 0.50] <- "large"
  strength
}

# Determine whether one mutation feature meets the established 10/10 rule.
is_eligible_mutation_feature <- function(results, gene) {
  feature.index <- match(gene, results$gene)
  !is.na(feature.index) && isTRUE(results$eligible[[feature.index]])
}

# Calculate HC3-adjusted coefficient uncertainty for an already fitted model.
calculate_hc3_coefficients <- function(model) {
  model.coefficients <- coef(model)
  retained.coefficients <- !is.na(model.coefficients)
  model.matrix.data <- model.matrix(model)[, retained.coefficients, drop = FALSE]
  cross.product <- crossprod(model.matrix.data)
  bread.matrix <- try(solve(cross.product), silent = TRUE)
  if (inherits(bread.matrix, "try-error")) bread.matrix <- qr.solve(cross.product)

  adjusted.residuals <- residuals(model) /
    pmax(1 - hatvalues(model), 1e-8)
  meat.matrix <- crossprod(model.matrix.data * adjusted.residuals)
  robust.covariance <- bread.matrix %*% meat.matrix %*% bread.matrix
  robust.standard.error <- sqrt(diag(robust.covariance))
  robust.estimate <- model.coefficients[retained.coefficients]
  robust.statistic <- robust.estimate / robust.standard.error
  residual.degrees.freedom <- df.residual(model)
  robust.p.value <- 2 * pt(
    abs(robust.statistic),
    df = residual.degrees.freedom,
    lower.tail = FALSE
  )
  confidence.critical <- qt(0.975, df = residual.degrees.freedom)

  data.frame(
    term = names(robust.estimate),
    estimate = unname(robust.estimate),
    std_error = robust.standard.error,
    statistic = robust.statistic,
    p_value = robust.p.value,
    conf_low = robust.estimate - confidence.critical * robust.standard.error,
    conf_high = robust.estimate + confidence.critical * robust.standard.error,
    stringsAsFactors = FALSE
  )
}

# Estimate a molecular association separately within each named lineage.
fit_lineage_specific_associations <- function(
  outcome,
  feature,
  lineage,
  continuous = TRUE,
  minimum.lineage.size = 10L,
  minimum.mutation.group.size = 3L
) {
  complete.rows <- complete.cases(outcome, feature, lineage)
  feature.values <- feature[complete.rows]
  if (!continuous) feature.values <- as.numeric(feature.values > 0)
  analysis.frame <- data.frame(
    dependency = outcome[complete.rows],
    feature_model = feature.values,
    lineage = as.character(lineage[complete.rows]),
    stringsAsFactors = FALSE
  )
  if (continuous && nrow(analysis.frame)) {
    feature.sd <- sd(analysis.frame$feature_model)
    if (is.finite(feature.sd) && feature.sd > 0) {
      analysis.frame$feature_model <- as.numeric(scale(
        analysis.frame$feature_model
      ))
    }
  }
  lineage.names <- sort(unique(analysis.frame$lineage))
  lineage.results <- lapply(lineage.names, function(lineage.name) {
    lineage.frame <- analysis.frame[
      analysis.frame$lineage == lineage.name,
      ,
      drop = FALSE
    ]
    altered.n <- if (continuous) NA_integer_ else {
      sum(lineage.frame$feature_model > 0)
    }
    reference.n <- if (continuous) NA_integer_ else {
      sum(lineage.frame$feature_model == 0)
    }
    empty.result <- function(reason) {
      data.frame(
        lineage = lineage.name,
        n = nrow(lineage.frame),
        altered_n = altered.n,
        reference_n = reference.n,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        estimable = FALSE,
        reason = reason,
        stringsAsFactors = FALSE
      )
    }
    if (nrow(lineage.frame) < minimum.lineage.size) {
      return(empty.result(sprintf(
        "Not estimable: fewer than %d complete models",
        minimum.lineage.size
      )))
    }
    if (length(unique(lineage.frame$feature_model)) < 2L) {
      return(empty.result("Not estimable: feature has no within-lineage variation"))
    }
    if (
      !continuous &&
        (altered.n < minimum.mutation.group.size ||
          reference.n < minimum.mutation.group.size)
    ) {
      return(empty.result(sprintf(
        "Not estimable: fewer than %d altered or %d reference models",
        minimum.mutation.group.size,
        minimum.mutation.group.size
      )))
    }
    lineage.model <- lm(dependency ~ feature_model, data = lineage.frame)
    lineage.coefficients <- calculate_hc3_coefficients(lineage.model)
    feature.coefficient <- lineage.coefficients[
      lineage.coefficients$term == "feature_model",
      ,
      drop = FALSE
    ]
    if (nrow(feature.coefficient) != 1L) {
      return(empty.result(
        "Not estimable: the within-lineage feature coefficient could not be fitted"
      ))
    }
    data.frame(
      lineage = lineage.name,
      n = nrow(lineage.frame),
      altered_n = altered.n,
      reference_n = reference.n,
      estimate = feature.coefficient$estimate,
      std_error = feature.coefficient$std_error,
      statistic = feature.coefficient$statistic,
      p_value = feature.coefficient$p_value,
      conf_low = feature.coefficient$conf_low,
      conf_high = feature.coefficient$conf_high,
      estimable = TRUE,
      reason = "Estimated",
      stringsAsFactors = FALSE
    )
  })
  if (!length(lineage.results)) return(data.frame())
  do.call(rbind, lineage.results)
}

# Assign deterministic cross-validation folds while balancing a grouping label.
make_stratified_folds <- function(strata, folds = 5L, seed = 1L) {
  folds <- as.integer(folds)
  if (folds < 2L) stop("folds must be at least 2.")
  strata <- as.character(strata)
  strata[is.na(strata) | !nzchar(strata)] <- "Missing"
  fold.id <- integer(length(strata))
  set.seed(seed)
  for (stratum in sort(unique(strata))) {
    stratum.rows <- which(strata == stratum)
    # Index through sample.int so one-member strata are never interpreted by
    # sample() as a request to draw from 1:stratum.rows.
    stratum.rows <- stratum.rows[sample.int(length(stratum.rows))]
    fold.id[stratum.rows] <- rep(seq_len(folds), length.out = length(stratum.rows))
  }
  if (any(fold.id == 0L)) stop("Stratified fold assignment left rows unassigned.")
  fold.id
}

# Learn fold-specific imputation and scaling without inspecting held-out rows.
fit_integrated_preprocessor <- function(
  feature.matrix,
  training.rows,
  continuous.columns,
  mutation.columns
) {
  feature.matrix <- as.matrix(feature.matrix)
  training.matrix <- feature.matrix[training.rows, , drop = FALSE]
  feature.count <- ncol(training.matrix)
  continuous.columns <- as.integer(continuous.columns)
  mutation.columns <- as.integer(mutation.columns)
  if (!setequal(c(continuous.columns, mutation.columns), seq_len(feature.count))) {
    stop("Every molecular feature must be classified as continuous or mutation.")
  }

  imputation <- numeric(feature.count)
  if (length(continuous.columns)) {
    imputation[continuous.columns] <- vapply(
      continuous.columns,
      function(column.index) {
        value <- median(training.matrix[, column.index], na.rm = TRUE)
        if (is.finite(value)) value else 0
      },
      numeric(1)
    )
  }
  if (length(mutation.columns)) {
    imputation[mutation.columns] <- vapply(
      mutation.columns,
      function(column.index) {
        observed <- training.matrix[, column.index]
        observed <- observed[!is.na(observed)]
        if (!length(observed)) 0 else as.numeric(mean(observed > 0) >= 0.5)
      },
      numeric(1)
    )
  }

  completed.training <- training.matrix
  missing.index <- which(is.na(completed.training), arr.ind = TRUE)
  if (nrow(missing.index)) {
    completed.training[missing.index] <-
      imputation[missing.index[, "col"]]
  }
  if (length(mutation.columns)) {
    completed.training[, mutation.columns] <-
      as.numeric(completed.training[, mutation.columns, drop = FALSE] > 0)
  }

  # All molecular predictors are placed on the same training-fold SD scale so
  # their common elastic-net penalty factor represents comparable shrinkage.
  center <- colMeans(completed.training)
  centered <- sweep(completed.training, 2L, center, "-")
  scale.value <- sqrt(
    colSums(centered^2) / pmax(nrow(centered) - 1L, 1L)
  )
  invalid.scale <- !is.finite(scale.value) | scale.value <= 0
  scale.value[invalid.scale] <- 1

  list(
    imputation = imputation,
    center = center,
    scale = scale.value,
    continuous_columns = continuous.columns,
    mutation_columns = mutation.columns,
    feature_names = colnames(feature.matrix)
  )
}

# Apply training-derived preprocessing to training or held-out molecular data.
apply_integrated_preprocessor <- function(feature.matrix, rows, preprocessor) {
  transformed <- as.matrix(feature.matrix[rows, , drop = FALSE])
  missing.index <- which(is.na(transformed), arr.ind = TRUE)
  if (nrow(missing.index)) {
    transformed[missing.index] <-
      preprocessor$imputation[missing.index[, "col"]]
  }
  if (length(preprocessor$mutation_columns)) {
    transformed[, preprocessor$mutation_columns] <-
      as.numeric(transformed[, preprocessor$mutation_columns] > 0)
  }
  transformed <- sweep(transformed, 2L, preprocessor$center, "-")
  transformed <- sweep(transformed, 2L, preprocessor$scale, "/")
  colnames(transformed) <- preprocessor$feature_names
  transformed
}

# Summarize repeated-model nonzero selection and coefficient directions.
summarize_elastic_net_stability <- function(
  coefficient.matrix,
  tolerance = 1e-10
) {
  coefficient.matrix <- as.matrix(coefficient.matrix)
  selected <- abs(coefficient.matrix) > tolerance
  selected.count <- rowSums(selected)
  positive.count <- rowSums(coefficient.matrix > tolerance)
  negative.count <- rowSums(coefficient.matrix < -tolerance)
  direction.concordance <- ifelse(
    selected.count > 0,
    pmax(positive.count, negative.count) / selected.count,
    NA_real_
  )
  median.selected.coefficient <- vapply(
    seq_len(nrow(coefficient.matrix)),
    function(feature.index) {
      selected.values <- coefficient.matrix[
        feature.index, selected[feature.index, ]
      ]
      if (!length(selected.values)) 0 else median(selected.values)
    },
    numeric(1)
  )
  data.frame(
    selected_count = selected.count,
    models_evaluated = ncol(coefficient.matrix),
    selection_frequency = selected.count / ncol(coefficient.matrix),
    positive_count = positive.count,
    negative_count = negative.count,
    direction_concordance = direction.concordance,
    median_selected_coefficient = median.selected.coefficient,
    row.names = rownames(coefficient.matrix),
    check.names = FALSE
  )
}

# Apply one visual style across all generated figures.
dependency_plot_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}
