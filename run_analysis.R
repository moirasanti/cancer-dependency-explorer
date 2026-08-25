#!/usr/bin/env Rscript

source("R/utils.R")

options <- parse_cli(commandArgs(trailingOnly = TRUE))
if (isTRUE(options$help)) {
  cat(cli_help(), "\n")
  quit(status = 0L)
}
options <- assert_valid_cli(options)

target <- options$target
raw_dir <- normalizePath(options$raw_dir, mustWork = TRUE)
biomarker_file <- if (nzchar(options$biomarkers)) normalizePath(options$biomarkers, mustWork = TRUE) else ""
input_files <- validate_inputs(raw_dir)
target_dir <- file.path("results", target)
table_dir <- file.path(target_dir, "tables")
figure_dir <- file.path(target_dir, "figures")
intermediate_dir <- file.path(target_dir, "intermediate")
report_file <- file.path("reports", paste0(target, "_target_assessment.html"))
invisible(lapply(c(table_dir, figure_dir, intermediate_dir, "reports"), ensure_directory))

message("Cancer Dependency Explorer: ", target)
message("DepMap inputs: ", raw_dir)
scripts <- file.path("scripts", sprintf("%02d.%s.R", 1:5, c(
  "prepare.depmap.data", "dependency.landscape", "molecular.associations",
  "multivariable.model", "generate.figures"
)))
for (script in scripts) {
  message("\nRunning ", script)
  source(script, local = environment())
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Package 'rmarkdown' is required to render the report.")
report_target_dir <- normalizePath(target_dir, mustWork = TRUE)
report_parameters <- list(target = target, target_dir = report_target_dir, release = "DepMap Public 26Q1")
rmarkdown::render(
  "reports/target_assessment.Rmd",
  output_file = basename(report_file), output_dir = dirname(report_file),
  params = report_parameters,
  envir = new.env(parent = globalenv()), quiet = TRUE
)
message("\nCompleted: ", report_file)
