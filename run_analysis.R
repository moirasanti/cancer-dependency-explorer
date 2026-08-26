#!/usr/bin/env Rscript

source("R/util.R")

# Read the supported command-line arguments without adding a CLI dependency.
command.arguments <- commandArgs(trailingOnly = TRUE)
command.options <- list(
  target = NULL,
  raw.dir = Sys.getenv("DEPMAP_RAW_DIR", unset = ""),
  biomarker.file = ""
)

help.text <- paste(
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

if (any(command.arguments %in% c("-h", "--help"))) {
  cat(help.text, "\n")
  quit(status = 0L)
}

argument.index <- 1L
while (argument.index <= length(command.arguments)) {
  argument.flag <- command.arguments[[argument.index]]
  if (!argument.flag %in% c("--target", "--raw-dir", "--biomarkers")) {
    stop("Unknown argument: ", argument.flag)
  }
  if (argument.index == length(command.arguments)) stop("Missing value after ", argument.flag)

  argument.value <- command.arguments[[argument.index + 1L]]
  if (argument.flag == "--target") command.options$target <- argument.value
  if (argument.flag == "--raw-dir") command.options$raw.dir <- argument.value
  if (argument.flag == "--biomarkers") command.options$biomarker.file <- argument.value
  argument.index <- argument.index + 2L
}

# Validate the requested target and input locations before starting the analysis.
if (is.null(command.options$target) || !nzchar(command.options$target)) {
  stop("--target is required.\n\n", help.text)
}
if (!nzchar(command.options$raw.dir)) stop("Provide --raw-dir or set DEPMAP_RAW_DIR.")

target <- toupper(trimws(command.options$target))
if (!grepl("^[A-Z0-9][A-Z0-9._-]*$", target)) stop("Invalid gene symbol: ", target)

raw.dir <- normalizePath(command.options$raw.dir, mustWork = TRUE)
biomarker.file <- if (nzchar(command.options$biomarker.file)) {
  normalizePath(command.options$biomarker.file, mustWork = TRUE)
} else {
  ""
}

required.files <- c(
  dependency = "CRISPRGeneEffect.csv",
  expression = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  copy_number = "OmicsCNGeneWGS.csv",
  mutation = "OmicsSomaticMutationsMatrixDamaging.csv",
  model = "Model.csv"
)
input.files <- file.path(raw.dir, required.files)
names(input.files) <- names(required.files)
missing.files <- input.files[!file.exists(input.files)]
if (length(missing.files)) {
  stop("Missing required DepMap files:\n", paste0("- ", missing.files, collapse = "\n"))
}

# Create one independent output tree for the requested target.
target.dir <- file.path("results", target)
table.dir <- file.path(target.dir, "tables")
figure.dir <- file.path(target.dir, "figures")
intermediate.dir <- file.path(target.dir, "intermediate")
report.file <- file.path("reports", paste0(target, "_target_assessment.html"))

output.directories <- c(table.dir, figure.dir, intermediate.dir, "reports")
for (output.directory in output.directories) {
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
}

# Run the numbered stages in order within this analysis environment.
message("Cancer Dependency Explorer: ", target)
message("DepMap inputs: ", raw.dir)
analysis.scripts <- file.path("scripts", sprintf("%02d.%s.R", 1:5, c(
  "prepare.depmap.data", "dependency.landscape", "molecular.associations",
  "multivariable.model", "generate.figures"
)))
for (analysis.script in analysis.scripts) {
  message("\nRunning ", analysis.script)
  source(analysis.script, local = environment())
}

# Render a self-contained report after all tables and figures are available.
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the report.")
}
report.target.dir <- normalizePath(target.dir, mustWork = TRUE)
report.parameters <- list(
  target = target,
  target_dir = report.target.dir,
  release = "DepMap Public 26Q1"
)
rmarkdown::render(
  "reports/target_assessment.Rmd",
  output_file = basename(report.file),
  output_dir = dirname(report.file),
  params = report.parameters,
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)
message("\nCompleted: ", report.file)
