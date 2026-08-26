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
  modality = c("copy_number", "mutation"),
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
            "Duplicate copy-number values disagreed for ", duplicate.gene,
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

# Compare a candidate with predictors already retained for the model.
calculate_maximum_absolute_correlation <- function(candidate, retained) {
  if (!length(retained)) return(NA_real_)

  correlations <- rep(NA_real_, length(retained))
  for (retained.index in seq_along(retained)) {
    existing <- retained[[retained.index]]
    complete.rows <- complete.cases(candidate, existing)
    if (
      sum(complete.rows) < 3L ||
      sd(candidate[complete.rows]) == 0 ||
      sd(existing[complete.rows]) == 0
    ) {
      next
    }
    correlations[[retained.index]] <- abs(
      cor(candidate[complete.rows], existing[complete.rows], method = "pearson")
    )
  }

  if (all(is.na(correlations))) NA_real_ else max(correlations, na.rm = TRUE)
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
