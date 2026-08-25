# Screen target expression, target copy number and genome-wide genomic features.
# Wide matrices are loaded one at a time to keep peak memory bounded.

analysis_data <- readRDS(file.path(intermediate_dir, "analysis_data.rds"))

expression_test <- safe_spearman(analysis_data$target_expression, analysis_data$dependency)
copy_number_test <- safe_spearman(analysis_data$target_copy_number, analysis_data$dependency)
continuous_associations <- data.frame(
  feature = c(paste0(target, " expression"), paste0(target, " copy number")),
  modality = c("expression", "copy_number"),
  n = c(expression_test[["n"]], copy_number_test[["n"]]),
  spearman_rho = c(expression_test[["estimate"]], copy_number_test[["estimate"]]),
  p_value = c(expression_test[["p_value"]], copy_number_test[["p_value"]])
)
fwrite(continuous_associations, file.path(table_dir, "07_target_continuous_associations.csv"))

read_wide_aligned <- function(file, analysis_ids) {
  wide <- fread(file, check.names = FALSE, showProgress = interactive())
  if (!"ModelID" %in% names(wide)) setnames(wide, 1L, "ModelID")
  if ("IsDefaultEntryForModel" %in% names(wide)) {
    keep <- as.character(wide$IsDefaultEntryForModel) %in% c("Yes", "TRUE", "1")
    wide <- wide[keep]
  }
  wide <- unique(wide, by = "ModelID")
  row_index <- match(analysis_ids, wide$ModelID)
  metadata <- c(
    "V1", "ModelID", "SequencingID", "ModelConditionID",
    "IsDefaultEntryForModel", "IsDefaultEntryForMC"
  )
  gene_columns <- setdiff(names(wide), metadata)
  symbols <- gene_symbol(gene_columns)
  keep_unique <- !duplicated(symbols)
  gene_columns <- gene_columns[keep_unique]
  symbols <- symbols[keep_unique]
  matrix <- as.matrix(wide[row_index, ..gene_columns])
  storage.mode(matrix) <- "double"
  rm(wide)
  list(matrix = matrix, symbols = symbols)
}

screen_mutations <- function(file, data) {
  loaded <- read_wide_aligned(file, data$ModelID)
  matrix <- loaded$matrix
  y <- data$dependency
  valid_y <- !is.na(y)
  matrix <- matrix[valid_y, , drop = FALSE]
  y <- y[valid_y]
  altered <- matrix > 0
  altered[is.na(matrix)] <- NA
  n1 <- colSums(altered, na.rm = TRUE)
  observed <- colSums(!is.na(altered))
  n0 <- observed - n1
  sum1 <- colSums(altered * y, na.rm = TRUE)
  sum0 <- colSums((!altered) * y, na.rm = TRUE)
  mean1 <- sum1 / n1
  mean0 <- sum0 / n0
  squared1 <- colSums(altered * y^2, na.rm = TRUE)
  squared0 <- colSums((!altered) * y^2, na.rm = TRUE)
  var1 <- (squared1 - n1 * mean1^2) / pmax(n1 - 1, 1)
  var0 <- (squared0 - n0 * mean0^2) / pmax(n0 - 1, 1)
  standard_error <- sqrt(var1 / n1 + var0 / n0)
  statistic <- (mean1 - mean0) / standard_error
  degrees_freedom <- (var1 / n1 + var0 / n0)^2 /
    ((var1 / n1)^2 / pmax(n1 - 1, 1) + (var0 / n0)^2 / pmax(n0 - 1, 1))
  eligible <- n1 >= 10L & n0 >= 10L & is.finite(statistic)
  p_value <- rep(NA_real_, length(statistic))
  p_value[eligible] <- 2 * pt(abs(statistic[eligible]), df = degrees_freedom[eligible], lower.tail = FALSE)
  out <- data.frame(
    gene = loaded$symbols, altered_n = n1, reference_n = n0,
    altered_mean_dependency = mean1, reference_mean_dependency = mean0,
    effect = mean1 - mean0, std_error = standard_error,
    conf_low = mean1 - mean0 - 1.96 * standard_error,
    conf_high = mean1 - mean0 + 1.96 * standard_error,
    p_value = p_value, eligible = eligible, stringsAsFactors = FALSE
  )
  out$fdr <- NA_real_
  out$fdr[eligible] <- p.adjust(out$p_value[eligible], method = "BH")
  out[order(out$fdr, out$effect), ]
}

screen_copy_number <- function(file, data) {
  loaded <- read_wide_aligned(file, data$ModelID)
  matrix <- loaded$matrix
  y <- data$dependency
  valid_y <- !is.na(y)
  matrix <- matrix[valid_y, , drop = FALSE]
  y <- y[valid_y]
  assayed <- rowSums(!is.na(matrix)) > 0L
  matrix <- matrix[assayed, , drop = FALSE]
  y <- y[assayed]
  n <- colSums(!is.na(matrix))
  coverage <- n / length(y)
  means <- colMeans(matrix, na.rm = TRUE)
  centered <- sweep(matrix, 2L, means, "-")
  variance <- colSums(centered^2, na.rm = TRUE) / pmax(n - 1, 1)
  correlation <- as.numeric(cor(y, matrix, use = "pairwise.complete.obs"))
  eligible <- coverage >= 0.80 & variance > 0 & is.finite(correlation) & abs(correlation) < 1
  statistic <- correlation * sqrt(pmax(n - 2, 1) / pmax(1 - correlation^2, 1e-12))
  p_value <- rep(NA_real_, length(correlation))
  p_value[eligible] <- 2 * pt(abs(statistic[eligible]), df = n[eligible] - 2, lower.tail = FALSE)
  out <- data.frame(
    gene = loaded$symbols, n = n, coverage = coverage, sd = sqrt(variance),
    correlation = correlation, p_value = p_value, eligible = eligible,
    stringsAsFactors = FALSE
  )
  out$fdr <- NA_real_
  out$fdr[eligible] <- p.adjust(out$p_value[eligible], method = "BH")
  out[order(out$fdr, -abs(out$correlation)), ]
}

message("Screening damaging mutations...")
mutation_associations <- screen_mutations(input_files[["mutation"]], analysis_data)
saveRDS(mutation_associations, file.path(intermediate_dir, "mutation_associations.rds"))
fwrite(head(mutation_associations[mutation_associations$eligible, ], 50L),
       file.path(table_dir, "08_top_mutation_associations.csv"))
gc()

message("Screening copy-number features...")
copy_number_associations <- screen_copy_number(input_files[["copy_number"]], analysis_data)
saveRDS(copy_number_associations, file.path(intermediate_dir, "copy_number_associations.rds"))
fwrite(head(copy_number_associations[copy_number_associations$eligible, ], 50L),
       file.path(table_dir, "09_top_copy_number_associations.csv"))
gc()

message("Completed genome-wide molecular association screens.")
