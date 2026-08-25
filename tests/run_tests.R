#!/usr/bin/env Rscript

source("R/utils.R")

stopifnot(identical(gene_symbol(c("TP53 (7157)", "A1BG..1")), c("TP53", "A1BG")))
stopifnot(mutation_group_eligible(c(rep(0, 10), rep(1, 10))))
stopifnot(!mutation_group_eligible(c(rep(0, 19), 1)))
stopifnot(copy_number_eligible(c(rep(1, 8), 2, NA), total_n = 10, minimum_coverage = 0.8))
stopifnot(!copy_number_eligible(c(rep(1, 7), rep(NA, 3)), total_n = 10, minimum_coverage = 0.8))

x <- 1:20
retained <- list(nearly_same = x + rep(c(-0.01, 0.01), 10), unrelated = rep(c(0, 1), 10))
stopifnot(maximum_absolute_correlation(x, retained) >= 0.99)
stopifnot(is.na(maximum_absolute_correlation(x, list(constant = rep(1, 20)))))

parsed <- parse_cli(c("--target", "mdm2", "--raw-dir", "data/raw"))
stopifnot(parsed$target == "mdm2", parsed$raw_dir == "data/raw")
invalid_gene <- tryCatch({
  assert_valid_cli(list(target = "bad gene", raw_dir = ".", biomarkers = "")); FALSE
}, error = function(e) TRUE)
stopifnot(invalid_gene)

set.seed(1)
cv_data <- data.frame(y = rnorm(100), x = rnorm(100))
first <- fixed_feature_cv(cv_data, y ~ x)
second <- fixed_feature_cv(cv_data, y ~ x)
stopifnot(identical(first, second))

cat("All unit tests passed.\n")
