test_that("parquet preserves quoted multiline annotations and rejects stale content", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  old <- options(skilift_provenance_format = NULL)
  on.exit(options(old), add = TRUE)
  table <- data.table::data.table(
    variant_summary = c('Contains ""hot spots""\nand a second line', "plain", NA_character_),
    position = c(1L, 2L, NA_integer_), af = c(0.25, NA_real_, 0),
    reported = c(TRUE, FALSE, NA)
  )
  path <- Skilift::write_provenance_table(table, root, "case", "collected")
  expect_match(path, "collected\\.parquet$")
  restored <- Skilift::read_provenance_table(root, "case", "collected")
  expect_identical(as.data.frame(restored), as.data.frame(table))
  manifest <- jsonlite::read_json(file.path(root, "case", "manifest.json"))
  expect_equal(manifest$schemaVersion, 1L)
  expect_equal(manifest$caseId, "case")
  expect_equal(manifest$stages$collected$rows, nrow(table))
  expect_identical(unlist(manifest$stages$collected$columns, use.names = FALSE), names(table))
  expect_equal(manifest$stages$collected$format, "parquet")

  recorded <- manifest$stages$collected$sha256
  writeBin(charToRaw("corrupt cache"), path)
  actual <- digest::digest(file = path, algo = "sha256")
  error <- tryCatch(Skilift::read_provenance_table(root, "case", "collected"), error = identity)
  expect_s3_class(error, "provenance_table_hash_mismatch")
  expect_match(conditionMessage(error), "pair case, stage collected", fixed = TRUE)
  expect_match(conditionMessage(error), recorded, fixed = TRUE)
  expect_match(conditionMessage(error), actual, fixed = TRUE)
})

test_that("curated publication preserves the collected cache and its manifest entry", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  table <- data.table::data.table(value = c(1L, 2L))
  collected <- Skilift::write_provenance_table(table, root, "case", "collected")
  manifest_path <- Skilift::provenance_table_paths(root, "case", "curated")$manifest
  original <- jsonlite::read_json(manifest_path)$stages$collected
  curated <- Skilift::write_provenance_table(table[1L], root, "case", "curated",
    inputs = list(collected = collected), parameters = list(reported_only = TRUE))
  manifest <- jsonlite::read_json(manifest_path)
  expect_identical(manifest$stages$collected, original)
  expect_false(identical(collected, curated))
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case", "collected")),
    as.data.frame(table))
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case")),
    as.data.frame(table[1L]))
  expect_equal(manifest$stages$curated$inputs$collected$sha256, original$sha256)
  expect_true(manifest$stages$curated$parameters$reported_only)
  expect_error(Skilift::provenance_table_paths(root, "case", "unknown"), "arg")
})

test_that("list columns and explicit RDS mode round trip through manifest dispatch", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  old <- options(skilift_provenance_format = NULL)
  on.exit(options(old), add = TRUE)
  table <- data.table::data.table(gene = c("A", "B"),
    transcript_coverage = list(c(TX1 = 0.5, TX2 = 1), numeric()))
  path <- Skilift::write_provenance_table(table, root, "case", "collected")
  expect_match(path, "collected\\.rds$")
  manifest <- jsonlite::read_json(file.path(root, "case", "manifest.json"))
  expect_equal(manifest$stages$collected$format, "rds")
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case", "collected")),
    as.data.frame(table))

  options(skilift_provenance_format = "rds")
  atomic <- data.table::data.table(value = 1L)
  path <- Skilift::write_provenance_table(atomic, root, "case", "curated")
  expect_match(path, "curated\\.rds$")
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case")), as.data.frame(atomic))
})

test_that("conventional caches can be read without a manifest", {
  root <- tempfile()
  dir.create(file.path(root, "case"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE))
  table <- data.table::data.table(value = 1L)
  saveRDS(table, file.path(root, "case", "collected.rds"))
  arrow::write_parquet(table, file.path(root, "case", "curated.parquet"))
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case", "collected")),
    as.data.frame(table))
  expect_identical(as.data.frame(Skilift::read_provenance_table(root, "case")), as.data.frame(table))
})

test_that("production lineage chains raw input through both table stages into metadata", {
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  raw <- file.path(root, "raw.rds")
  table <- data.table::data.table(value = 1L)
  saveRDS(table, raw)
  collected <- Skilift::write_provenance_table(table, root, "case", "collected",
    inputs = list(annotation = raw), parameters = list(filter = "PASS"))
  curate <- function() {
    .provenance_context <- Skilift:::provenance_context(list(pair = "case", collected = collected))
    value <- Skilift::read_provenance_table(root, "case", "collected")
    Skilift::write_provenance_table(value, root, "case", "curated", inputs = list(collected = collected))
  }
  curated <- curate()
  case_dir <- file.path(root, "data", "case")
  dir.create(case_dir, recursive = TRUE)
  downstream <- function() {
    .provenance_context <- Skilift:::provenance_context(list(pair = "case", provenance_table = curated))
    value <- Skilift::read_provenance_table(root, "case")
    Skilift:::provenance_write_json(value, file.path(case_dir, "filtered.events.json"))
  }
  downstream()
  jsonlite::write_json(list(list(pair = "case")), file.path(case_dir, "metadata.json"))
  Skilift:::provenance_finalize(file.path(root, "data"))
  metadata <- jsonlite::read_json(file.path(case_dir, "metadata.json"))
  event <- metadata[[1L]]$provenance$artifacts[[1L]]
  expect_equal(event$sources[[1L]]$path, curated)
  curated_receipt <- event$sources[[1L]]$production
  expect_equal(curated_receipt$sources[[1L]]$path, collected)
  collected_receipt <- curated_receipt$sources[[1L]]$production
  expect_equal(collected_receipt$sources[[1L]]$path, raw)
  expect_equal(collected_receipt$parameters$filter, "PASS")
  expect_equal(collected_receipt$parameters$provenance_stage, "collected")
})

test_that("the provenance consumer preserves legacy event projection and JSON content", {
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  chrom_sizes <- system.file("extdata", "hg19.regularChr.chrom.sizes", package = "gGnome")
  expect_true(file.exists(chrom_sizes))
  withr::local_envvar(c(DEFAULT_GENOME = chrom_sizes))
  fixture <- system.file("extdata", "test_data", "oncotable_test_data", "new_oncotable",
    "unit_oncotable.rds", package = "Skilift")
  graph <- system.file("extdata", "test_data", "oncotable_test_data", "complex.rds", package = "Skilift")
  table <- base::readRDS(fixture)
  # This historical fixture serialized copy counts as strings; current collectors
  # produce numeric values. Keep every event class and repair only that fixture type.
  table$value <- as.numeric(table$value)
  # Its complex rows are per-type counts and predate the footprint column.
  # Reconstruct their full footprints from the paired graph's original events.
  graph_events <- base::readRDS(graph)$meta$events
  footprints <- graph_events[, .(footprint = paste(footprint, collapse = ","), count = .N), by = type]
  complex_rows <- which(table$source == "complex")
  matched <- match(table$type[complex_rows], footprints$type)
  expect_false(anyNA(matched))
  expect_equal(table$value[complex_rows], as.numeric(footprints$count[matched]))
  table$footprint <- NA_character_
  table$footprint[complex_rows] <- footprints$footprint[matched]
  typed_fixture <- file.path(root, "typed_legacy.rds")
  saveRDS(table, typed_fixture)
  Skilift::write_provenance_table(table, root, "case", "curated",
    inputs = list(legacy_fixture = fixture, paired_graph = graph))
  legacy_json <- file.path(root, "legacy.json")
  table_json <- file.path(root, "table.json")
  legacy <- Skilift::create_filtered_events("case", typed_fixture, graph, legacy_json, return_table = TRUE)
  current <- Skilift::create_filtered_events_from_table("case", root, jabba_gg = graph,
    out_file = table_json, return_table = TRUE)
  expect_equal(as.data.frame(current), as.data.frame(legacy))
  expect_true(all(c("SNV", "AMP", "HOMDEL") %in% current$vartype))
  expect_setequal(current$vartype[current$type == "Complex SV"],
    c("chromothripsis", "rigma", "pyrgo", "bfb"))
  expect_true(any(grepl("fusion", current$type, ignore.case = TRUE)))
  expect_identical(jsonlite::read_json(table_json), jsonlite::read_json(legacy_json))
})

test_that("the additive collector preserves table shape and applies field transforms", {
  fixture <- system.file("extdata", "test_data", "oncotable_test_data", "oncokb_snvcn.maf", package = "Skilift")
  gencode <- Skilift:::process_gencode(system.file(
    "extdata", "test_data", "test_gencode_v29lift37.rds", package = "Skilift"))
  cytoband <- Skilift:::process_cytoband(system.file("extdata", "data", "cytoband.rds", package = "Skilift"))
  args <- list(pair = "case", oncokb_snv = fixture, gencode = gencode, cytoband = cytoband, verbose = FALSE)
  legacy <- do.call(Skilift::oncotable, args)
  current <- do.call(Skilift::create_provenance_table, args)
  expect_identical(as.data.frame(current), as.data.frame(legacy))
  transformed <- do.call(Skilift::create_provenance_table, c(args, list(fields = list(
    Skilift::annotation_field("tier", transform = function(x) x + 10L)
  ))))
  expect_identical(transformed$tier, legacy$tier + 10L)
  expect_identical(names(transformed), names(legacy))
})

test_that("cohort event lifting does not swallow cache hash mismatches", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  path <- Skilift::write_provenance_table(data.table::data.table(value = 1L), root, "case", "curated")
  writeBin(charToRaw("changed"), path)
  cohort <- suppressWarnings(Skilift::Cohort$new(data.table::data.table(pair = "case", jabba_gg = NA_character_)))
  cohort$inputs[, provenance_table_root := root]
  expect_error(Skilift::lift_filtered_events_from_tables(cohort, file.path(root, "data"), cores = 1,
    verbose = TRUE), "SHA-256 mismatch for pair case, stage curated")
})

test_that("the provenance wrapper collects and lifts genuine SNVs without a graph column", {
  root <- tempfile()
  on.exit(unlink(root, recursive = TRUE))
  withr::local_envvar(c(DEFAULT_GENOME = system.file(
    "extdata", "hg19.regularChr.chrom.sizes", package = "gGnome")))
  fixture <- system.file("extdata", "test_data", "oncotable_test_data", "oncokb_snvcn.maf", package = "Skilift")
  cohort <- suppressWarnings(Skilift::Cohort$new(data.table::data.table(
    pair = "case", oncokb_snv = fixture,
    tumor_type = "BRCA", disease = "Breast Cancer", primary_site = "Breast")))
  cohort$type <- "paired"
  expect_false(any(c("events", "balanced_jabba_gg", "jabba_gg") %in% names(cohort$inputs)))
  result <- Skilift::lift_mvp_provenance(cohort, file.path(root, "data"), file.path(root, "provenance"),
    cores = 1, genome_length = 29008, width = 10000,
    gencode = system.file("extdata", "test_data", "test_gencode_v29lift37.rds", package = "Skilift"))
  table <- Skilift::read_provenance_table(file.path(root, "provenance"), "case", "collected")
  expect_gt(nrow(table), 0L)
  events <- jsonlite::fromJSON(file.path(root, "data", "case", "filtered.events.json"))
  expect_gt(nrow(events), 0L)
  expect_true(any(events$vartype == "SNV"))
  metadata <- jsonlite::read_json(file.path(root, "data", "case", "metadata.json"))[[1L]]
  expect_identical(metadata$pair, "case")
  expect_identical(metadata$primary_site, "Breast")
  expect_true(is.character(result$inputs$string_summary))
})
