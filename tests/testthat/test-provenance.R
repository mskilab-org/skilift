test_that("provenance observations retain supplied paths and file facts", {
  source_file <- tempfile("provenance input ")
  writeBin(as.raw(c(1, 2, 3)), source_file)
  on.exit(unlink(source_file))

  observed <- Skilift:::provenance_observation(source_file)
  expect_identical(observed$supplied_path, source_file)
  expect_true(grepl("^/", observed$absolute_path))
  expect_true(file.exists(observed$resolved_path))
  expect_identical(observed$stat_status, "ok")
  expect_equal(observed$size_bytes, 3)
  expect_null(observed$created_at)
  expect_true(grepl("Z$", observed$observed_at))
})

test_that("file bindings omit inline objects and never stringify metadata", {
  inline <- data.frame(gene = "SYNTHETIC_GENE", variant = "SYNTHETIC_VARIANT")
  non_files <- list(
    inline, inline["gene"], list(inline), list(list(gene = "SYNTHETIC_GENE")),
    matrix("SYNTHETIC_GENE"), 0.42, TRUE
  )
  for (value in non_files) {
    expect_length(Skilift:::provenance_input_bindings(list(input = value), c(source = "input")), 0L)
    expect_length(Skilift:::provenance_direct_binding(value, "source"), 0L)
  }

  for (value in list("missing.rds", factor("missing.rds"), list("missing.rds"))) {
    bindings <- Skilift:::provenance_input_bindings(list(input = value), c(source = "input"))
    expect_identical(bindings[[1]]$supplied_path, "missing.rds")
    expect_identical(bindings[[1]]$cohort_column, "input")
  }
  expect_length(Skilift:::provenance_input_bindings(
    list(input = NA_character_, disease = "SYNTHETIC_DIAGNOSIS"), c(source = "input")
  ), 0L)
  expect_length(Skilift:::provenance_direct_binding(c("one.rds", "two.rds"), "source"), 0L)
})

test_that("metadata records Conpair files and only the selected purity and SNV inputs", {
  binding <- Skilift:::provenance_backfill_registry[["metadata.json"]]
  binding$columns["graph"] <- "jabba_gg"
  row <- list(
    jabba_gg = "selected-graph.rds", purple_pp_bestFit = "unused-purple.tsv",
    conpair_contamination = "contamination.txt", conpair_concordance = "concordance.txt",
    oncokb_snv = "selected.maf", somatic_snvs = "unused.vcf", germline_snvs = "unused-germline.vcf"
  )
  selected <- Skilift:::provenance_metadata_columns(row, binding$columns)
  production <- Skilift:::provenance_input_bindings(row, selected)
  backfill <- Skilift:::provenance_backfill_inputs(row, list(inputs = row), binding)
  columns <- function(inputs) vapply(inputs, `[[`, character(1), "cohort_column")
  expected <- c("jabba_gg", "oncokb_snv", "conpair_contamination", "conpair_concordance")
  expect_setequal(columns(production), expected)
  expect_setequal(columns(backfill), expected)

  for (absent in list(NULL, NA_character_)) {
    row["jabba_gg"] <- list(absent)
    row["oncokb_snv"] <- list(absent)
    inputs <- Skilift:::provenance_backfill_inputs(row, list(inputs = row), binding)
    expect_setequal(columns(inputs), c(
      "purple_pp_bestFit", "somatic_snvs", "conpair_contamination", "conpair_concordance"
    ))
  }
})

test_that("legacy list Cohorts backfill BAM and metadata without changing files or inputs", {
  root <- tempfile("provenance-legacy-pairing-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  paths <- file.path(root, c("tumor.bam", "normal.bam", "tumor.wgs", "normal.wgs", "contamination.txt"))
  for (path in c(paths, paste0(paths[1:2], ".bai"))) writeLines("synthetic input", path)

  wide <- data.frame(
    pair = "case-1", tumor_bam = paths[1], normal_bam = paths[2],
    tumor_wgs_metrics = paths[3], normal_wgs_metrics = paths[4], conpair_contamination = paths[5]
  )
  long <- data.frame(
    pair = c("case-1", "case-1"), tumor_type = c("tumor", "normal"),
    bam = paths[1:2], wgs_metrics = paths[3:4], conpair_contamination = paths[5]
  )
  forms <- list(wide, data.table::as.data.table(wide), long, data.table::as.data.table(long))
  for (i in seq_along(forms)) {
    output_dir <- file.path(root, as.character(i))
    case_dir <- file.path(output_dir, "case-1")
    dir.create(case_dir, recursive = TRUE)
    keys <- c("tumor.bam", "normal.bam", "tumor.bam.bai", "normal.bam.bai", "metadata.json")
    outputs <- file.path(case_dir, keys)
    for (path in outputs) writeLines("existing output", path)
    before <- tools::md5sum(outputs)
    cohort <- list(inputs = forms[[i]])
    inputs_before <- serialize(cohort$inputs, NULL)
    cohort_path <- file.path(output_dir, "cohort.rds")
    saveRDS(cohort, cohort_path)

    result <- Skilift::lift_provenance(cohort_path, output_dir)
    expect_identical(result$written, "case-1")
    expect_length(result$failed, 0L)
    expect_identical(tools::md5sum(outputs), before)
    expect_identical(serialize(readRDS(cohort_path)$inputs, NULL), inputs_before)
    document <- jsonlite::fromJSON(file.path(case_dir, "provenance.json"), simplifyVector = FALSE)
    expect_setequal(names(document$outputs), keys)
    for (key in keys[1:4]) {
      expected_path <- file.path(root, key)
      expect_identical(document$outputs[[key]]$inputs[[1]]$supplied_path, expected_path)
    }
    metadata_paths <- vapply(document$outputs[["metadata.json"]]$inputs, `[[`, character(1), "supplied_path")
    expect_setequal(metadata_paths, paths[3:5])

    # In-memory lists use the same path as saved legacy lists and are not mutated.
    result <- Skilift::lift_provenance(cohort, output_dir, force = TRUE)
    expect_identical(result$written, "case-1")
    expect_identical(serialize(cohort$inputs, NULL), inputs_before)
  }
})

test_that("provenance documents merge output entries without changing others", {
  previous <- list(
    schema_version = 1L,
    pair = "case-1",
    scope = "direct-file-inputs",
    extension = list(note = "preserve me"),
    outputs = list(old.json = list(lifter = "old_lifter"))
  )
  update <- list(new.json = list(lifter = "new_lifter"))
  document <- Skilift:::provenance_create_document("case-1", update, previous)

  expect_true("old.json" %in% names(document$outputs))
  expect_true("new.json" %in% names(document$outputs))
  expect_identical(document$extension, previous$extension)
  expect_error(Skilift:::provenance_create_document("../case", update), "Case")
  expect_error(Skilift:::provenance_create_document("case-1", list("../x" = list())), "artifact")
})

test_that("legacy Cohort backfill writes ordinary provenance without changing outputs", {
  output_dir <- tempfile("provenance-backfill-")
  case_dir <- file.path(output_dir, "case-1")
  dir.create(case_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE))

  input_path <- file.path(output_dir, "coverage.rds")
  output_path <- file.path(case_dir, "coverage.arrow")
  saveRDS(list(values = c(1, 2)), input_path)
  writeBin(as.raw(c(4, 5, 6)), output_path)
  output_before <- readBin(output_path, "raw", n = file.info(output_path)$size)

  legacy_cohort <- list(inputs = data.frame(
    pair = "case-1", tumor_coverage = input_path,
    stringsAsFactors = FALSE
  ))
  cohort_path <- file.path(output_dir, "cohort.rds")
  saveRDS(legacy_cohort, cohort_path)

  result <- Skilift::lift_provenance(cohort_path, output_dir)
  expect_identical(result$written, "case-1")
  expect_length(result$skipped, 0L)
  expect_identical(readBin(output_path, "raw", n = file.info(output_path)$size), output_before)

  manifest_path <- file.path(case_dir, "provenance.json")
  document <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  expect_identical(document$pair, "case-1")
  expect_identical(document$scope, "direct-file-inputs")
  expect_true("coverage.arrow" %in% names(document$outputs))
  entry <- document$outputs[["coverage.arrow"]]
  expect_false("capture_mode" %in% names(entry))
  expect_false("input_basis" %in% names(entry))
  expect_identical(entry$lifter, "lift_denoised_coverage")
  expect_identical(entry$inputs[[1]]$cohort_column, "tumor_coverage")
  expect_true(grepl("Z$", entry$run$started_at))

  manifest_before <- readBin(manifest_path, "raw", n = file.info(manifest_path)$size)
  writeBin(as.raw(c(7, 8, 9)), output_path)
  skipped <- Skilift::lift_provenance(cohort_path, output_dir)
  expect_identical(skipped$skipped, "case-1")
  expect_identical(readBin(manifest_path, "raw", n = file.info(manifest_path)$size), manifest_before)

  forced <- Skilift::lift_provenance(cohort_path, output_dir, force = TRUE)
  expect_identical(forced$written, "case-1")
  expect_identical(readBin(output_path, "raw", n = file.info(output_path)$size), as.raw(c(7, 8, 9)))
})

test_that("existing BAM pairing keeps the Case key and tumor/normal paths", {
  inputs <- data.table::data.table(
    pair = c("case-1", "case-1"),
    tumor_type = c("tumor", "normal"),
    bam = c("tumor.bam", "normal.bam")
  )
  paired <- Skilift:::pairify_cohort_inputs(
    inputs, tumor_normal_columns = "bam", sep = "_"
  )

  expect_identical(paired$pair, "case-1")
  expect_identical(paired$tumor_bam, "tumor.bam")
  expect_identical(paired$normal_bam, "normal.bam")
})

test_that("segment-width backfill selection matches the lifter graph check", {
  fields <- Skilift:::segment_width_segstats_fields()
  graph <- structure(
    list(nodes = list(dt = as.data.frame(setNames(rep(list(1L), length(fields)), fields)))),
    class = "gGraph"
  )
  expect_true(Skilift:::segment_width_graph_has_segstats(graph))
  graph$nodes$dt <- graph$nodes$dt[, -1, drop = FALSE]
  expect_false(Skilift:::segment_width_graph_has_segstats(graph))
})

test_that("backfill records only resources selected by saved options", {
  cohort <- structure(list(inputs = data.frame(
    pair = c("masked", "unmasked"),
    tumor_coverage = c("masked.rds", "unmasked.rds"),
    denoised_coverage_apply_mask = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )), class = "Cohort")
  binding <- Skilift:::provenance_backfill_registry[["coverage.arrow"]]
  masked <- Skilift:::provenance_backfill_inputs(cohort$inputs[1, ], cohort, binding)
  unmasked <- Skilift:::provenance_backfill_inputs(cohort$inputs[2, ], cohort, binding)

  expect_true("mask" %in% vapply(masked, `[[`, character(1), "role"))
  expect_false("mask" %in% vapply(unmasked, `[[`, character(1), "role"))
  expect_true("reference_settings" %in% vapply(unmasked, `[[`, character(1), "role"))
})

test_that("production publication commits all staged outputs with one Case manifest", {
  output_dir <- tempfile("provenance-publish-success-")
  case_dir <- file.path(output_dir, "case-1")
  dir.create(case_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE))

  context <- Skilift:::provenance_new_context(output_dir)
  for (key in c("coverage.arrow", "complex.json", "multiplicity.png")) {
    final_path <- file.path(case_dir, key)
    staged_path <- Skilift:::provenance_stage_output(context, final_path)
    expect_true(grepl(paste0("\\.", tools::file_ext(key), "$"), staged_path))
    writeBin(as.raw(c(1, 2, 3)), staged_path)
    Skilift:::provenance_queue_output(
      context, "case-1", key, "smoke", list(), staged_path, final_path
    )
  }

  Skilift:::provenance_finish(context)
  document <- jsonlite::fromJSON(
    file.path(case_dir, "provenance.json"), simplifyVector = FALSE
  )
  expect_setequal(names(document$outputs), c("coverage.arrow", "complex.json", "multiplicity.png"))
  expect_identical(readBin(file.path(case_dir, "coverage.arrow"), "raw", n = 3), as.raw(1:3))
  expect_identical(readBin(file.path(case_dir, "complex.json"), "raw", n = 3), as.raw(1:3))
  expect_identical(readBin(file.path(case_dir, "multiplicity.png"), "raw", n = 3), as.raw(1:3))
})

test_that("production publishes readable records when another artifact record is unreadable", {
  output_dir <- tempfile("provenance-partial-")
  case_dir <- file.path(output_dir, "case-1")
  dir.create(case_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE))
  keys <- c("coverage.arrow", "complex.json")
  previous <- list()
  for (key in keys) {
    writeLines("old output", file.path(case_dir, key))
    previous[[key]] <- list(lifter = "previous_run")
  }
  jsonlite::write_json(
    Skilift:::provenance_create_document("case-1", previous),
    file.path(case_dir, "provenance.json"), auto_unbox = TRUE
  )

  context <- Skilift:::provenance_new_context(output_dir)
  for (key in keys) {
    final <- file.path(case_dir, key)
    staged <- Skilift:::provenance_stage_output(context, final)
    writeLines("new output", staged)
    Skilift:::provenance_queue_output(context, "case-1", key, "current_run", list(), staged, final)
  }
  records <- list.files(context$spool_dir, full.names = TRUE)
  damaged <- records[vapply(records, function(path) identical(readRDS(path)$artifact_key, "complex.json"), logical(1))]
  expect_length(damaged, 1L)
  writeLines("unreadable record", damaged)
  expect_message(suppressWarnings(Skilift:::provenance_finish(context)), "could not read")

  expect_identical(readLines(file.path(case_dir, "coverage.arrow")), "new output")
  expect_identical(readLines(file.path(case_dir, "complex.json")), "old output")
  document <- jsonlite::fromJSON(file.path(case_dir, "provenance.json"), simplifyVector = FALSE)
  expect_identical(document$outputs[["coverage.arrow"]]$lifter, "current_run")
  expect_identical(document$outputs[["complex.json"]]$lifter, "previous_run")
})

test_that("production publication keeps old files when provenance cannot be saved", {
  output_dir <- tempfile("provenance-publish-")
  case_dir <- file.path(output_dir, "case-1")
  dir.create(case_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE))

  output_path <- file.path(case_dir, "coverage.arrow")
  manifest_path <- file.path(case_dir, "provenance.json")
  writeBin(as.raw(c(1, 2, 3)), output_path)
  writeLines("invalid manifest", manifest_path)
  old_output <- readBin(output_path, "raw", n = file.info(output_path)$size)
  old_manifest <- readBin(manifest_path, "raw", n = file.info(manifest_path)$size)

  context <- Skilift:::provenance_new_context(output_dir)
  staged <- Skilift:::provenance_stage_output(context, output_path)
  writeBin(as.raw(c(4, 5, 6)), staged)
  Skilift:::provenance_queue_output(
    context, "case-1", "coverage.arrow", "lift_denoised_coverage",
    inputs = list(), staged_path = staged, final_path = output_path
  )

  expect_message(Skilift:::provenance_finish(context), "provenance failed for Case case-1")
  expect_identical(readBin(output_path, "raw", n = file.info(output_path)$size), old_output)
  expect_identical(readBin(manifest_path, "raw", n = file.info(manifest_path)$size), old_manifest)
})

test_that("backfill skips an existing manifest by default and logs per-Case failures", {
  output_dir <- tempfile("provenance-backfill-failures-")
  dir.create(file.path(output_dir, "case-1"), recursive = TRUE)
  dir.create(file.path(output_dir, "case-2"), recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE))

  writeLines("not-json", file.path(output_dir, "case-1", "provenance.json"))
  writeLines("payload", file.path(output_dir, "case-1", "coverage.arrow"))
  writeLines("payload", file.path(output_dir, "case-2", "coverage.arrow"))
  cohort <- list(inputs = data.frame(
    pair = c("case-1", "case-2"),
    tumor_coverage = c("missing-1.rds", "missing-2.rds"),
    stringsAsFactors = FALSE
  ))

  result <- expect_message(
    Skilift::lift_provenance(cohort, output_dir),
    "already exists"
  )
  expect_identical(result$skipped, "case-1")
  expect_identical(result$written, "case-2")
  expect_length(result$failed, 0L)
})
