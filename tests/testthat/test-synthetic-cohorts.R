# These fixtures do not require devtools::load_all() or access to real cohorts.
source(testthat::test_path("..", "fixtures", "synthetic-cohorts.R"), local = TRUE)
source(testthat::test_path("..", "fixtures", "check-synthetic-cohorts.R"), local = TRUE)

synthetic_root <- testthat::test_path("..", "..", "inst", "extdata", "test_data", "synthetic_cohorts")

test_that("synthetic NGS and Myeloseq inputs are coherent and reproducible", {
  for (package in c("GenomicRanges", "Rsamtools", "jsonlite")) skip_if_not_installed(package)
  expect_true(run_synthetic_fixture_checks(
    synthetic_root, testthat::test_path("..", "fixtures", "generate-synthetic-cohorts.R")
  ))
})

test_that("graph recipes construct real gGnome objects when available", {
  skip_if_not_installed("gGnome")
  skip_if_not_installed("data.table")
  output <- tempfile("synthetic-graphs-")
  on.exit(unlink(output, recursive = TRUE))
  paths <- materialize_synthetic_graphs(synthetic_root, output)
  expect_length(paths, 3L)
  for (path in paths) {
    graph <- readRDS(path)
    expect_s3_class(graph, "gGraph")
    expect_equal(length(graph$nodes$gr), 3L)
    expect_equal(graph$nodes$gr$cn, c(2, 6, 2))
  }
})

test_that("synthetic inputs can enter the legacy Cohort wrapper", {
  skip_if_not_installed("Skilift")
  skip_if_not_installed("data.table")
  for (profile in c("ngs", "myeloseq")) {
    inputs <- data.table::as.data.table(synthetic_case_inputs(synthetic_root, profile))
    cases <- fixture_read_table(file.path(synthetic_root, profile, "cases.tsv"))
    for (type in unique(cases$cohort_type)) {
      selected <- inputs[which(cases$cohort_type == type)]
      cohort <- suppressWarnings(Skilift::Cohort$new(
        selected, col_mapping = list(pair = "pair"),
        cohort_type = type, merge_tumor_type_db = FALSE
      ))
      expect_identical(cohort$inputs$pair, selected$pair)
      expect_true(all(file.exists(cohort$inputs$oncokb_snv)))
      expect_true(all(file.exists(cohort$inputs$multiplicity)))
    }
  }
})

test_that("existing OncoKB collection accepts both synthetic profiles", {
  skip_if_not_installed("Skilift")
  for (profile in c("ngs", "myeloseq")) {
    inputs <- synthetic_case_inputs(synthetic_root, profile)
    for (i in 1:2) {
      result <- Skilift:::collect_oncokb(inputs$oncokb_snv[i], inputs$multiplicity[i], verbose = FALSE)
      expect_equal(nrow(result), c(3L, 2L)[i])
      expect_identical(result$gene, c("MOCK1", "MOCK2", "MOCK3")[seq_len(nrow(result))])
      expect_equal(result$VAF, c(0.25, 0.2, 0.05)[seq_len(nrow(result))])
      if (profile == "myeloseq") expect_true(all(is.na(result$segment_cn)))
    }
  }
})
