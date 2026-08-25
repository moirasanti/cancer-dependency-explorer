suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

required_files <- c(
  dependency = "CRISPRGeneEffect.csv",
  expression = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  copy_number = "OmicsCNGeneWGS.csv",
  mutation = "OmicsSomaticMutationsMatrixDamaging.csv",
  model = "Model.csv"
)

gene_symbol <- function(x) {
  x <- sub("\\s*\\([^)]*\\)$", "", x)
  sub("\\.\\.[0-9]+\\.?$", "", x)
}

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

validate_inputs <- function(raw_dir) {
  paths <- file.path(raw_dir, required_files)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing required DepMap files:\n", paste0("- ", missing, collapse = "\n"))
  }
  setNames(paths, names(required_files))
}

header_info <- function(file) {
  header <- names(fread(file, nrows = 0L, check.names = FALSE))
  data.frame(column = header, symbol = gene_symbol(header), stringsAsFactors = FALSE)
}

find_gene_column <- function(file, gene) {
  info <- header_info(file)
  hits <- info$column[info$symbol == gene]
  if (!length(hits)) stop("Gene ", gene, " is absent from ", basename(file), ".")
  if (length(hits) > 1L) warning("Multiple columns found for ", gene, "; using the first.")
  hits[[1L]]
}

read_gene_data <- function(file, genes, default_only = FALSE) {
  info <- header_info(file)
  gene_columns <- info$column[match(genes, info$symbol, nomatch = 0L)]
  found_symbols <- gene_symbol(gene_columns)
  missing <- setdiff(genes, found_symbols)
  if (length(missing)) {
    stop("Genes absent from ", basename(file), ": ", paste(missing, collapse = ", "))
  }

  metadata <- intersect(c("ModelID", "IsDefaultEntryForModel"), info$column)
  if (!"ModelID" %in% metadata) metadata <- unique(c(info$column[[1L]], metadata))
  out <- fread(file, select = unique(c(metadata, gene_columns)), check.names = FALSE)
  if (!"ModelID" %in% names(out)) setnames(out, 1L, "ModelID")
  if (default_only && "IsDefaultEntryForModel" %in% names(out)) {
    keep <- as.character(out$IsDefaultEntryForModel) %in% c("Yes", "TRUE", "1")
    out <- out[keep]
    out[, IsDefaultEntryForModel := NULL]
  }
  gene_names <- names(out)[gene_symbol(names(out)) %in% genes]
  setnames(out, gene_names, gene_symbol(gene_names))
  unique(out, by = "ModelID")
}

safe_spearman <- function(x, y) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(c(n = sum(keep), estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  c(n = sum(keep), estimate = unname(test$estimate), p_value = test$p.value)
}

mutation_group_eligible <- function(values, minimum_group_n = 10L) {
  values <- values[!is.na(values)] > 0
  sum(values) >= minimum_group_n && sum(!values) >= minimum_group_n
}

copy_number_eligible <- function(values, total_n = length(values), minimum_coverage = 0.80) {
  observed <- values[!is.na(values)]
  length(observed) / total_n >= minimum_coverage && length(observed) > 1L && sd(observed) > 0
}

maximum_absolute_correlation <- function(candidate, retained) {
  if (!length(retained)) return(NA_real_)
  correlations <- vapply(retained, function(existing) {
    keep <- complete.cases(candidate, existing)
    if (sum(keep) < 3L || sd(candidate[keep]) == 0 || sd(existing[keep]) == 0) return(NA_real_)
    abs(cor(candidate[keep], existing[keep], method = "pearson"))
  }, numeric(1))
  if (all(is.na(correlations))) NA_real_ else max(correlations, na.rm = TRUE)
}

format_p <- function(x) {
  ifelse(is.na(x), "not estimable", ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x)))
}

load_biomarker_config <- function(path, target) {
  empty <- list(mutation = character(), copy_number = character())
  if (is.null(path) || !nzchar(path)) return(empty)
  if (!file.exists(path)) stop("Biomarker configuration does not exist: ", path)
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required for biomarker configuration.")
  config <- yaml::read_yaml(path)
  if (!is.null(config$target) && toupper(config$target) != target) {
    stop("Biomarker configuration target is ", config$target, ", not ", target, ".")
  }
  list(
    mutation = unique(toupper(unlist(config$mutation %||% character()))),
    copy_number = unique(toupper(unlist(config$copy_number %||% character())))
  )
}

read_feature_matrix <- function(file, genes, default_only = TRUE) {
  if (!length(genes)) return(data.frame(ModelID = character()))
  as.data.frame(read_gene_data(file, genes, default_only = default_only), check.names = FALSE)
}

hc3_tidy <- function(model) {
  coefficients <- coef(model)
  retained <- !is.na(coefficients)
  x <- model.matrix(model)[, retained, drop = FALSE]
  bread <- tryCatch(solve(crossprod(x)), error = function(e) qr.solve(crossprod(x)))
  adjusted_residuals <- residuals(model) / pmax(1 - hatvalues(model), 1e-8)
  meat <- crossprod(x * adjusted_residuals)
  covariance <- bread %*% meat %*% bread
  standard_error <- sqrt(diag(covariance))
  estimate <- coefficients[retained]
  statistic <- estimate / standard_error
  degrees_freedom <- df.residual(model)
  p_value <- 2 * pt(abs(statistic), df = degrees_freedom, lower.tail = FALSE)
  critical <- qt(0.975, df = degrees_freedom)
  data.frame(
    term = names(estimate), estimate = unname(estimate), std_error = standard_error,
    statistic = statistic, p_value = p_value,
    conf_low = estimate - critical * standard_error,
    conf_high = estimate + critical * standard_error,
    stringsAsFactors = FALSE
  )
}

fixed_feature_cv <- function(data, formula, folds = 10L, seed = 20260825L) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(folds), length.out = nrow(data)))
  observed <- predicted <- rep(NA_real_, nrow(data))
  outcome <- all.vars(formula)[[1L]]
  for (fold in seq_len(folds)) {
    train <- data[fold_id != fold, , drop = FALSE]
    test <- data[fold_id == fold, , drop = FALSE]
    fit <- lm(formula, data = train)
    predicted[fold_id == fold] <- suppressWarnings(predict(fit, newdata = test))
    observed[fold_id == fold] <- test[[outcome]]
  }
  keep <- complete.cases(observed, predicted)
  rmse <- sqrt(mean((observed[keep] - predicted[keep])^2))
  r_squared <- 1 - sum((observed[keep] - predicted[keep])^2) /
    sum((observed[keep] - mean(observed[keep]))^2)
  data.frame(folds = folds, n = sum(keep), rmse = rmse, r_squared = r_squared, seed = seed)
}

parse_cli <- function(args) {
  values <- list(target = NULL, raw_dir = Sys.getenv("DEPMAP_RAW_DIR", unset = ""), biomarkers = "")
  help_requested <- any(args %in% c("-h", "--help"))
  if (help_requested) return(c(values, list(help = TRUE)))
  index <- 1L
  while (index <= length(args)) {
    flag <- args[[index]]
    if (!flag %in% c("--target", "--raw-dir", "--biomarkers")) stop("Unknown argument: ", flag)
    if (index == length(args)) stop("Missing value after ", flag)
    key <- sub("^--", "", flag)
    key <- gsub("-", "_", key)
    values[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  values$help <- FALSE
  values
}

cli_help <- function() {
  paste(
    "Cancer Dependency Explorer",
    "",
    "Usage:",
    "  Rscript run_analysis.R --target GENE [--raw-dir PATH] [--biomarkers FILE]",
    "",
    "--target       HGNC gene symbol whose CRISPR dependency will be assessed",
    "--raw-dir      directory containing the DepMap CSV files (or set DEPMAP_RAW_DIR)",
    "--biomarkers   optional YAML file containing mutation and copy_number gene lists",
    sep = "\n"
  )
}

assert_valid_cli <- function(options) {
  if (is.null(options$target) || !nzchar(options$target)) stop("--target is required.\n\n", cli_help())
  if (!nzchar(options$raw_dir)) stop("Provide --raw-dir or set DEPMAP_RAW_DIR.")
  options$target <- toupper(trimws(options$target))
  if (!grepl("^[A-Z0-9][A-Z0-9._-]*$", options$target)) stop("Invalid gene symbol: ", options$target)
  options
}

theme_dependency <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(), plot.title.position = "plot",
      plot.title = element_text(face = "bold"), legend.position = "bottom"
    )
}
