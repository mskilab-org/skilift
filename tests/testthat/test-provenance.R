test_that("reader receipts preserve multiple inputs and finalize without changing metrics", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    source1 <- file.path(root, "source1.rds")
    source2 <- file.path(root, "source2.tsv")
    saveRDS(3, source1)
    writeLines("value\n4", source2)
    case_dir <- file.path(root, "case")
    dir.create(case_dir)
    jsonlite::write_json(list(list(pair = "case", tmb = 0)), file.path(case_dir, "metadata.json"), auto_unbox = TRUE)
    worker <- function() {
        .provenance_context <- Skilift:::provenance_context(list(pair = "case", multiplicity = source1))
        a <- Skilift:::readRDS(source1)
        b <- Skilift:::fread(source2)
        Skilift:::provenance_write_json(list(value = a + b$value), file.path(case_dir, "mutations.json"))
    }
    worker()
    Skilift:::provenance_finalize(root)
    first <- jsonlite::read_json(file.path(case_dir, "metadata.json"))
    expect_equal(first[[1]]$tmb, 0)
    artifact <- first[[1]]$provenance$artifacts[[1]]
    expect_equal(artifact$evidence, "captured")
    expect_setequal(vapply(artifact$sources, `[[`, character(1), "path"), c(source1, source2))
    Skilift:::provenance_finalize(root)
    second <- jsonlite::read_json(file.path(case_dir, "metadata.json"))
    expect_length(second[[1]]$provenance$artifacts, 1L)
    worker()
    expect_length(list.files(file.path(case_dir, ".provenance", "history")), 2L)
})

test_that("downstream receipts retain the upstream production snapshot", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    source <- file.path(root, "input.rds")
    intermediate <- file.path(root, "oncotable.rds")
    saveRDS(1, source)
    upstream <- function() {
        .provenance_context <- Skilift:::provenance_context(list(pair = "case"))
        saveRDS(Skilift:::readRDS(source), intermediate)
        Skilift:::provenance_receipt(intermediate)
    }
    downstream <- function() {
        .provenance_context <- Skilift:::provenance_context(list(pair = "case"))
        value <- Skilift:::readRDS(intermediate)
        Skilift:::provenance_write_json(list(value = value), file.path(root, "filtered.events.json"))
    }
    upstream()
    downstream()
    receipt <- jsonlite::read_json(file.path(root, ".provenance", "filtered.events.json.json"))
    expect_equal(receipt$sources[[1]]$production$sources[[1]]$path, source)
})

test_that("metadata step isolates field inputs and preserves false/zero values", {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(0, path)
    worker <- function() {
        .provenance_context <- Skilift:::provenance_context(list(pair = "case"))
        metadata <- data.table::data.table(pair = "case", tmb = NA_real_)
        metadata <- Skilift:::provenance_step(metadata, {
            metadata$tmb <- Skilift:::readRDS(path)
            metadata
        })
        expect_equal(metadata$tmb, 0)
        expect_equal(.provenance_context$fields$tmb$sources[[1]]$path, path)
    }
    worker()
})
test_that("cached cohorts keep exact inputs and remapping requires a matching original", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    old <- file.path(root, "old.rds")
    replacement <- file.path(root, "replacement.rds")
    saveRDS(1, old)
    saveRDS(2, replacement)
    cohort <- Skilift::Cohort$new(data.table::data.table(pair = "P"), merge_tumor_type_db = FALSE)
    cohort$inputs[, allelic_jabba_gg := old]
    cache <- file.path(root, "cohort.rds")
    saveRDS(cohort, cache)
    loaded <- Skilift:::provenance_load_cohort(cache)
    expect_identical(loaded$inputs, cohort$inputs)
    loaded <- Skilift:::provenance_begin_upload(loaded, file.path(root, "output"))
    expect_true(file.exists(loaded$provenance_upload$snapshot$path))
    expect_identical(loaded[1]$provenance_upload$id, loaded$provenance_upload$id)
    remapping <- list(schemaVersion = 1L, caseId = "P", changes = list(list(column = "allelic_jabba_gg", fromPath = "/wrong", toPath = replacement)))
    expect_error(Skilift:::provenance_apply_remapping(loaded, remapping), "no longer matches")
    expect_identical(loaded$inputs$allelic_jabba_gg, old)
    remapping$changes[[1]]$fromPath <- old
    Skilift:::provenance_apply_remapping(loaded, remapping)
    expect_identical(loaded$inputs$allelic_jabba_gg, replacement)
    expect_identical(readRDS(cache)$inputs$allelic_jabba_gg, old)
})
