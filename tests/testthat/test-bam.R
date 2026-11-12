bam_link_fixture <- function(root) {
    work_dir <- file.path(root, "work files")
    staged_dir <- file.path(root, "staged")
    dir.create(work_dir, recursive = TRUE)
    dir.create(staged_dir)

    sources <- file.path(work_dir, c("tumor.bam", "tumor.bam.bai", "normal.bam", "normal.bai"))
    payloads <- c("tumor BAM", "tumor index", "normal BAM", "normal index")
    for (i in seq_along(sources)) writeLines(payloads[i], sources[i])
    stopifnot(
        file.symlink("../work files/tumor.bam", file.path(staged_dir, "tumor-hop")),
        file.symlink("tumor-hop", file.path(staged_dir, "tumor.bam")),
        file.symlink("../work files/tumor.bam.bai", file.path(staged_dir, "tumor.bam.bai")),
        file.symlink(sources[3], file.path(staged_dir, "normal.bam"))
    )
    cohort <- Cohort$new(
        data.table::data.table(
            pair = "T___vs___N",
            tumor_bam = file.path(staged_dir, "tumor.bam"),
            normal_bam = file.path(staged_dir, "normal.bam")
        ),
        merge_tumor_type_db = FALSE
    )
    list(
        cohort = cohort,
        sources = sources,
        payloads = payloads,
        staged_dir = staged_dir,
        output_dir = file.path(root, "gOS data"),
        destinations = file.path(root, "gOS data", "T___vs___N", c("tumor.bam", "tumor.bam.bai", "normal.bam", "normal.bam.bai"))
    )
}

test_that("hardlink migration follows source symlinks for BAMs and indexes", {
    skip_on_os("windows")
    root <- tempfile("bam links ' ")
    on.exit(unlink(root, recursive = TRUE))
    fixture <- bam_link_fixture(root)

    lift_bam(fixture$cohort, fixture$output_dir)
    expect_true(all(file.exists(fixture$destinations)))
    expect_equal(Sys.readlink(fixture$destinations), rep("", 4))

    # Updating the originals distinguishes hardlinks from copies.
    updated <- paste(fixture$payloads, "updated")
    for (i in seq_along(fixture$sources)) writeLines(updated[i], fixture$sources[i])
    unlink(fixture$staged_dir, recursive = TRUE)
    unlink(fixture$sources)
    expect_equal(vapply(fixture$destinations, readLines, character(1), USE.NAMES = FALSE), updated)
})

test_that("symlink migration points directly at absolute original paths", {
    skip_on_os("windows")
    root <- tempfile("bam links ' ")
    on.exit(unlink(root, recursive = TRUE))
    fixture <- bam_link_fixture(root)

    lift_bam(fixture$cohort, fixture$output_dir, symlink = TRUE)
    expect_equal(Sys.readlink(fixture$destinations), normalizePath(fixture$sources))
    unlink(fixture$staged_dir, recursive = TRUE)
    expect_equal(vapply(fixture$destinations, readLines, character(1), USE.NAMES = FALSE), fixture$payloads)
})

test_that("copy migration dereferences sources without sharing their data", {
    skip_on_os("windows")
    root <- tempfile("bam links ' ")
    on.exit(unlink(root, recursive = TRUE))
    fixture <- bam_link_fixture(root)

    lift_bam(fixture$cohort, fixture$output_dir, copy = TRUE)
    for (source in fixture$sources) writeLines("changed original", source)
    expect_equal(Sys.readlink(fixture$destinations), rep("", 4))
    expect_equal(vapply(fixture$destinations, readLines, character(1), USE.NAMES = FALSE), fixture$payloads)
    expect_error(lift_bam(fixture$cohort, fixture$output_dir, copy = TRUE, symlink = TRUE), "Choose either")
})

test_that("overwrite replaces destination symlinks without modifying their targets", {
    skip_on_os("windows")
    root <- tempfile("bam links ' ")
    on.exit(unlink(root, recursive = TRUE))

    for (mode in c("hardlink", "symlink", "copy")) {
        fixture <- bam_link_fixture(file.path(root, mode))
        dir.create(dirname(fixture$destinations[1]), recursive = TRUE)
        sentinel <- file.path(root, paste0(mode, "-sentinel"))
        writeLines("keep target unchanged", sentinel)
        stopifnot(file.symlink("missing-bam", fixture$destinations[1]), file.symlink(sentinel, fixture$destinations[2]))
        args <- list(cohort = fixture$cohort, output_data_dir = fixture$output_dir,
                     symlink = mode == "symlink", copy = mode == "copy")

        invisible(capture.output(do.call(lift_bam, args)))
        expect_equal(Sys.readlink(fixture$destinations[1]), "missing-bam")
        expect_equal(Sys.readlink(fixture$destinations[2]), sentinel)

        do.call(lift_bam, c(args, list(overwrite = TRUE)))
        expect_equal(vapply(fixture$destinations, readLines, character(1), USE.NAMES = FALSE), fixture$payloads)
        expect_equal(readLines(sentinel), "keep target unchanged")
        expect_equal(Sys.readlink(fixture$destinations),
                     if (mode == "symlink") normalizePath(fixture$sources) else rep("", 4))
    }
})

test_that("indexes beside the supplied BAM path take precedence over target indexes", {
    skip_on_os("windows")
    root <- tempfile("bam links ' ")
    on.exit(unlink(root, recursive = TRUE))
    fixture <- bam_link_fixture(root)
    supplied_index <- file.path(fixture$staged_dir, "tumor.bai")
    unlink(file.path(fixture$staged_dir, "tumor.bam.bai"))
    writeLines("published index", supplied_index)

    lift_bam(fixture$cohort, fixture$output_dir, symlink = TRUE)
    expect_equal(Sys.readlink(fixture$destinations[2]), normalizePath(supplied_index))
    expect_equal(readLines(fixture$destinations[2]), "published index")
})
