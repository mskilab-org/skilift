# Standalone checks: no installed Skilift/devtools/testthat or private data needed.
# Rscript tests/fixtures/check-synthetic-cohorts.R [fixture-root]

run_synthetic_fixture_checks <- function(root, generator) {
  validate_synthetic_cohorts(root)
  source(generator, local = TRUE)
  temporary <- tempfile("skilift-synthetic-")
  repeated <- tempfile("skilift-synthetic-repeat-")
  on.exit(unlink(c(temporary, repeated), recursive = TRUE))
  generate_synthetic_cohorts(temporary)
  generate_synthetic_cohorts(repeated)
  validate_synthetic_cohorts(temporary)
  paths <- function(dir) sort(list.files(dir, recursive = TRUE, all.files = TRUE))
  fixture_assert(identical(paths(root), paths(temporary)), "Regeneration changed the file inventory")
  hashes <- function(dir) unname(tools::md5sum(file.path(dir, paths(dir))))
  fixture_assert(identical(hashes(temporary), hashes(repeated)), "Regeneration is not byte-for-byte reproducible in the same toolchain")
  # Committed RDS/BGZF bytes can vary with R/Bioconductor/htslib versions.
  # Compare their meaning, while the two fresh builds above must match exactly.
  contents <- function(path) {
    if (grepl("\\.rds$", path)) {
      gr <- readRDS(path)
      return(list(ranges = as.data.frame(gr), seqinfo = as.data.frame(GenomicRanges::seqinfo(gr))))
    }
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    readLines(con, warn = FALSE)
  }
  for (path in paths(root)) {
    if (grepl("\\.tbi$", path)) next # Both corpora's tabix indexes were exercised above.
    fixture_assert(isTRUE(all.equal(contents(file.path(root, path)), contents(file.path(temporary, path)))),
                   paste("Regenerated fixture changed:", path))
  }
  refusal <- tryCatch(generate_synthetic_cohorts(temporary), error = identity)
  fixture_assert(inherits(refusal, "error"), "Generator must refuse a non-empty destination")
  check_synthetic_parallel_reads(root)
  # Prove that validation detects a broken relationship, not just readable files.
  input <- synthetic_case_inputs(temporary, "myeloseq")$multiplicity[1]
  gr <- readRDS(input)
  gr$VAF[1] <- 0.99
  saveRDS(gr, input, version = 2)
  failure <- tryCatch(validate_synthetic_cohorts(temporary), error = identity)
  fixture_assert(inherits(failure, "error") && grepl("VAF", conditionMessage(failure)), "Validation failed to detect inconsistent VAF")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), mustWork = TRUE)
  directory <- dirname(script)
  source(file.path(directory, "synthetic-cohorts.R"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript check-synthetic-cohorts.R [fixture-root]")
  root <- if (length(args)) args[1] else file.path(directory, "../../inst/extdata/test_data/synthetic_cohorts")
  run_synthetic_fixture_checks(root, file.path(directory, "generate-synthetic-cohorts.R"))
  files <- list.files(root, full.names = TRUE, recursive = TRUE, all.files = TRUE)
  cat(sprintf("PASS: six synthetic Cases; %d files; %d bytes.\n", length(files), sum(file.info(files)$size)))
  cat("Schema, links, provenance, regeneration, two-worker reads, and invalid-data detection passed.\n")
}
