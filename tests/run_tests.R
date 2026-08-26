#!/usr/bin/env Rscript

source("R/util.R")

# Gene symbols are recovered from both DepMap and duplicate-column formats.
stopifnot(identical(
  extract_gene_symbol(c("TP53 (7157)", "A1BG..1")),
  c("TP53", "A1BG")
))

# Small synthetic DepMap data exercise selective and aligned matrix readers.
test.file <- tempfile(fileext = ".csv")
test.data <- data.frame(
  ModelID = c("ACH-1", "ACH-2"),
  IsDefaultEntryForModel = c("Yes", "No"),
  check.names = FALSE
)
test.data[["TP53 (7157)"]] <- c(1, 0)
test.data[["MDM2 (4193)"]] <- c(-1.2, -0.3)
fwrite(test.data, test.file)

selected.data <- read_depmap_gene_data(
  test.file,
  c("TP53", "MDM2"),
  default.only = TRUE
)
stopifnot(
  nrow(selected.data) == 1L,
  identical(names(selected.data), c("ModelID", "TP53", "MDM2"))
)

aligned.data <- read_aligned_genomic_matrix(test.file, c("ACH-2", "ACH-1"))
stopifnot(
  identical(aligned.data$symbols, c("TP53", "MDM2")),
  is.na(aligned.data$matrix[1, 1]),
  aligned.data$matrix[2, 1] == 1,
  aligned.data$matrix[2, 2] == -1.2
)

# Duplicate mutation symbols retain an alteration from any source column;
# continuous duplicates use all available values and report disagreements.
duplicate.file <- tempfile(fileext = ".csv")
duplicate.data <- data.frame(
  ModelID = paste0("ACH-", 1:4),
  check.names = FALSE
)
duplicate.data[["DUP (1)"]] <- c(0, 0, NA, NA)
duplicate.data[["DUP..1"]] <- c(1, 0, 1, NA)
fwrite(duplicate.data, duplicate.file)

mutation.messages <- capture.output(
  mutation.duplicate <- read_aligned_genomic_matrix(
    duplicate.file,
    duplicate.data$ModelID,
    modality = "mutation"
  ),
  type = "message"
)
stopifnot(
  identical(mutation.duplicate$symbols, "DUP"),
  isTRUE(all.equal(
    as.numeric(mutation.duplicate$matrix[, "DUP"]),
    c(1, 0, 1, NA)
  )),
  any(grepl("classified as altered", mutation.messages, fixed = TRUE))
)

copy.number.messages <- capture.output(
  copy.number.duplicate <- read_aligned_genomic_matrix(
    duplicate.file,
    duplicate.data$ModelID,
    modality = "copy_number"
  ),
  type = "message"
)
stopifnot(
  isTRUE(all.equal(
    as.numeric(copy.number.duplicate$matrix[, "DUP"]),
    c(0.5, 0, 1, NA)
  )),
  any(grepl("used the mean", copy.number.messages, fixed = TRUE))
)

selected.messages <- capture.output(
  selected.duplicate <- read_depmap_gene_data(
    duplicate.file,
    "DUP",
    modality = "mutation"
  ),
  type = "message"
)
stopifnot(
  isTRUE(all.equal(selected.duplicate$DUP, c(1, 0, 1, NA))),
  any(grepl("classified as altered", selected.messages, fixed = TRUE))
)

# Statistical helpers return stable values for valid and degenerate inputs.
spearman.result <- calculate_safe_spearman(1:10, 10:1)
stopifnot(
  spearman.result[["n"]] == 10,
  abs(spearman.result[["estimate"]] + 1) < 1e-12
)

candidate <- 1:20
retained.predictors <- list(
  nearly.same = candidate + rep(c(-0.01, 0.01), 10),
  unrelated = rep(c(0, 1), 10)
)
stopifnot(
  calculate_maximum_absolute_correlation(candidate, retained.predictors) >= 0.99,
  is.na(calculate_maximum_absolute_correlation(candidate, list(constant = rep(1, 20))))
)

# The shared plotting helper returns a ggplot2 theme object.
stopifnot(inherits(dependency_plot_theme(), "theme"))

# All utility function names use underscore nomenclature.
utility.lines <- readLines("R/util.R")
function.lines <- grep("^[a-z][a-z0-9_]* <- function", utility.lines, value = TRUE)
function.names <- sub(" <- function.*$", "", function.lines)
stopifnot(
  length(function.names) == 6L,
  all(grepl("^[a-z][a-z0-9]*(_[a-z0-9]+)+$", function.names))
)

unlink(c(test.file, duplicate.file))
cat("All unit tests passed.\n")
