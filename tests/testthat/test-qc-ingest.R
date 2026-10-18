test_that("QC sources honor overwrite, fill, derived values and key order", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    first <- file.path(root, "first.csv")
    second <- file.path(root, "second.csv")
    writeLines(c("sample,depth,reads", "A,10,100", "B,NA,0", "none,99,99", ",99,99"), first)
    writeLines(c("sample,depth,reads", "A,20,200", "B,30,300"), second)
    fields <- list(annotation_field("depth", transform = as.numeric),
                   annotation_field("reads", transform = as.numeric))
    seed <- qc_source(first, "sample", fields, kind = "csv")
    fill <- qc_source(second, "sample", fields, kind = "csv", precedence = "fill")
    out <- file.path(root, "merged.json")
    merged <- merge_qc_sources(list(seed, fill), out,
        derived = list(ratio = function(record) record$reads / record$depth),
        field_order = c("reads", "depth"))
    expect_identical(names(merged), c("A", "B"))
    expect_identical(names(merged$A), c("reads", "depth", "ratio"))
    expect_equal(merged$A, list(reads = 100, depth = 10, ratio = 10))
    expect_equal(merged$B, list(reads = 0, depth = 30, ratio = 0))
    parsed <- jsonlite::fromJSON(out, simplifyVector = FALSE)
    expect_identical(names(parsed$A), c("reads", "depth", "ratio"))
    overwrite <- qc_source(second, "sample", fields, kind = "csv")
    merged <- merge_qc_sources(list(seed, overwrite), out, field_order = c("depth", "reads"))
    expect_equal(merged$A$depth, 20)
    expect_error(merge_qc_sources(list(seed), out), "field_order")
})

test_that("QC fills retain false, zero and character NA; transforms can omit invalid measurements", {
    path <- tempfile(fileext = ".json")
    out <- tempfile(fileext = ".json")
    on.exit(unlink(c(path, out)))
    jsonlite::write_json(list(id = "sample", zero = 0, flag = FALSE, missing = "bad"), path, auto_unbox = TRUE)
    seed <- qc_source(path, "id", list(
        annotation_field("zero"), annotation_field("flag"),
        annotation_field("missing", transform = function(x) NA_character_),
        annotation_field("absent", transform = function(x) NULL)
    ))
    fill <- qc_source(path, "id", list(
        annotation_field("zero", transform = function(x) 10),
        annotation_field("flag", transform = function(x) TRUE),
        annotation_field("missing", transform = function(x) "present")
    ), precedence = "fill")
    merged <- merge_qc_sources(list(seed, fill), out, field_order = c("zero", "flag", "missing"))
    expect_equal(merged$sample$zero, 0)
    expect_identical(merged$sample$flag, FALSE)
    expect_identical(merged$sample$missing, NA_character_)
    expect_false("absent" %in% names(merged$sample))
})

test_that("one QC source resolves duplicate samples before fill precedence", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    first <- file.path(root, "first.csv")
    last <- file.path(root, "last.csv")
    writeLines(c("sample,depth", "A,1", "A,2"), first)
    writeLines(c("sample,depth", "A,3"), last)
    source <- qc_source(c(first, last), "sample",
        list(annotation_field("depth", transform = as.numeric)), kind = "csv", precedence = "fill")
    merged <- merge_qc_sources(list(source), file.path(root, "out.json"), field_order = "depth")
    expect_equal(merged$A$depth, 3)
})

test_that("Picard reader strips comments and the plugin hook selects PAIR", {
    path <- tempfile()
    out <- tempfile()
    on.exit(unlink(c(path, out)))
    writeLines(c("# Picard", "", "CATEGORY\tTOTAL_READS", "FIRST_OF_PAIR\t50", "PAIR\t100", "# end"), path)
    table <- read_picard_metrics(path)
    expect_equal(table$TOTAL_READS, c(50L, 100L))
    source <- qc_source(path, "CATEGORY", list(annotation_field("TOTAL_READS", as = "total_reads")),
        kind = "picard", locate = function(row, key) if (key == "PAIR") row else NULL,
        key_transform = function(key) "DNA")
    merged <- merge_qc_sources(list(source), out, field_order = "total_reads")
    expect_equal(merged$DNA$total_reads, 100)
    expect_identical(names(merged), "DNA")
})

test_that("merged QC is ingested before flags and records field-level provenance", {
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path))
    jsonlite::write_json(list(DNA = list(pair = "not-published", mean_depth = 600,
        fraction_of_reads_aligned = 0.98)), path, auto_unbox = TRUE)
    flags <- list(list("FAIL", "mean_depth", `<`, 500, "Mean target depth", "500X"))
    worker <- function() {
        .provenance_context <- Skilift:::provenance_context(list(pair = "TM", qc_merged = path))
        metadata <- suppressWarnings(create_metadata(pair = "TM", cohort_type = "heme",
            seqnames_genome_width_or_genome_length = 29008, qc_flags_config = flags,
            qc_merged = path, qc_sample_key = "DNA"))
        expect_equal(metadata$tumor_median_coverage, 600)
        expect_false("pair" %in% names(metadata$coverage_qc[[1]]))
        expect_equal(metadata$qcMetrics[[1]][[1]]$code, "PASS")
        expect_equal(.provenance_context$fields$coverage_qc$sources[[1]]$path, path)
        expect_equal(.provenance_context$fields$tumor_median_coverage$sources[[1]]$path, path)
        parsed <- jsonlite::fromJSON(jsonlite::toJSON(metadata, auto_unbox = TRUE), simplifyVector = FALSE)[[1]]
        expect_equal(parsed$coverage_qc$mean_depth, 600)
        qc <- parsed$qcMetrics
        if (is.list(qc) && length(qc) == 1L && is.list(qc[[1]]) &&
            length(qc[[1]]) > 0L && is.list(qc[[1]][[1]]) &&
            all(c("code", "title") %in% names(qc[[1]][[1]]))) qc <- qc[[1]]
        expect_identical(names(qc[[1]]), c("code", "title"))
    }
    worker()
    before <- data.table::data.table(pair = "TM", tumor_median_coverage = 12)
    expect_identical(add_coverage_metrics_from_qc_file(before, NULL, "DNA"), before)
    expect_identical(add_coverage_metrics_from_qc_file(before, path, "absent"), before)
    expect_identical(Skilift:::qc_sample_key_from_row(list(pair = "TM",
        ir_analysis_id = "TM-DNA_v3_TM-RNA_RNA_v2")), "TM-DNA")
    expect_identical(Skilift:::qc_sample_key_from_row(list(pair = "TM", qc_sample_key = "explicit")), "explicit")
    expect_identical(Skilift:::qc_sample_key_from_row(list(pair = "ordinary")), "ordinary")
})

test_that("lift_metadata passes merged QC and IR sample identity to the first write", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    path <- file.path(root, "qc_merged.json")
    jsonlite::write_json(list(DNA = list(mean_depth = 700)), path, auto_unbox = TRUE)
    cohort <- Skilift::Cohort$new(data.table::data.table(pair = "TM"), merge_tumor_type_db = FALSE)
    cohort$type <- "heme"
    cohort$inputs$qc_merged <- path
    cohort$inputs$ir_analysis_id <- "DNA_v2_OTHER_RNA_v1"
    flags <- list(list("FAIL", "mean_depth", `<`, 500, "Mean target depth", "500X"))
    data.table::set(cohort$inputs, j = "qc_flags", value = list(list(flags)))
    cohort <- suppressWarnings(lift_metadata(cohort, file.path(root, "data"),
        genome_length = 29008, do_lift_datafiles_json = FALSE))
    output <- jsonlite::fromJSON(file.path(root, "data", "TM", "metadata.json"), simplifyVector = FALSE)[[1]]
    expect_equal(output$coverage_qc$mean_depth, 700)
    expect_equal(output$tumor_median_coverage, 700)
    expect_false("pair" %in% names(output$coverage_qc))
    expect_equal(output$qcMetrics, list(list(code = "PASS", title = "Mean target depth (700) >= 500X")))
})

test_that("Myeloseqer plugin preserves API BAM PDF precedence and numeric parsing", {
    pipeline <- Sys.getenv("MYELOSEQER_PIPELINE_ROOT")
    skip_if(!nzchar(pipeline), "Set MYELOSEQER_PIPELINE_ROOT to exercise the pipeline QC plugin")
    plugin <- new.env(parent = globalenv())
    sys.source(file.path(pipeline, "lib", "mys_qc_plugin.R"), envir = plugin)
    sys.source(file.path(pipeline, "myeloseqer_qc_flags.R"), envir = plugin)
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    qc_json <- file.path(root, "qc_metrics.json")
    bam_csv <- file.path(root, "bam_stats.csv")
    report_csv <- file.path(root, "report_coverage.csv")
    coverage_csv <- file.path(root, "coverage_from_qc_pdf.csv")
    alignment <- file.path(root, "DNA_v3_OTHER_RNA_v1", "qc", "alignment_summary_metrics")
    dir.create(dirname(alignment), recursive = TRUE)
    writeLines(c("# Picard", "CATEGORY\tTOTAL_READS\tPF_READS_ALIGNED\tPF_ALIGNED_BASES\tMEAN_READ_LENGTH",
        "PAIR\t1000\t900\t100000\t100"), alignment)
    jsonlite::write_json(list(samples = list(DNA = "DNA_v3"), qc_metrics = list(DNA_v3 = list(
        "Mean Coverage Depth(fold)" = " 700.5 ", "Total number of Reads" = " 800 ",
        "Total number of Bases(Mbp)" = "not numeric", "Total number of Bases(AQ20)(Mbp)" = "0.096",
        "Mean Read Length(AQ20)" = "", "Mean Read Length(AQ30)" = " 115 ",
        "Coverage within Target Region" = "100"
    ))), qc_json, auto_unbox = TRUE)
    writeLines(c("sample,fraction_of_reads_aligned,mapped_reads", "DNA,0.98,980"), bam_csv)
    writeLines(c("sample,mean_depth,uniformity,on_target", "DNA,100,80%,85%"), report_csv)
    writeLines(c("sample,mean_depth,uniformity,on_target", "DNA,650, 91.2% , 97.5% "), coverage_csv)
    merge <- function() merge_qc_sources(
        plugin$myeloseqer_qc_sources(qc_json, c(report_csv, coverage_csv), bam_csv, alignment),
        file.path(root, "merged.json"), derived = plugin$MYELOSEQER_QC_DERIVED,
        field_order = plugin$MYELOSEQER_QC_FIELD_ORDER)
    merged <- merge()
    expect_equal(merged$DNA, list(mean_depth = 700.5, total_reads = 800, pf_aligned_bases = 96000,
        mean_read_length = 115, fraction_of_reads_aligned = 0.98, pf_reads_aligned = 980,
        uniformity = 0.912, on_target = 0.975, percent_aligned = 1))
    expect_identical(names(merged$DNA), plugin$MYELOSEQER_QC_FIELD_ORDER)
    expect_identical(plugin$myeloseqer_qc_paths(root, "json"), qc_json)
    metadata <- suppressWarnings(create_metadata(pair = "TM", cohort_type = "heme",
        seqnames_genome_width_or_genome_length = 29008,
        qc_merged = file.path(root, "merged.json"), qc_sample_key = "DNA",
        qc_flags_config = plugin$myeloseqer_qc_flag_thresholds))
    output <- jsonlite::fromJSON(jsonlite::toJSON(metadata, auto_unbox = TRUE), simplifyVector = FALSE)[[1]]
    expect_length(output$qcMetrics, 6L)
    expect_identical(vapply(output$qcMetrics, `[[`, character(1), "code"),
        c("PASS", "PASS", "PASS", "PASS", "WARN", "PASS"))
    expect_identical(vapply(output$qcMetrics, `[[`, character(1), "title"), c(
        "Fraction of reads aligned (0.98) >= 90%",
        "Percent on target (0.975) >= 80%",
        "Mean target depth (700) >= 500X",
        "Uniformity of coverage (0.912) >= 85%",
        "Total reads (800) < 500k",
        "Mean read length (115) >= 100 bp"
    ))
    writeLines(c("sample,fraction_of_reads_aligned,mapped_reads", "DNA,1.02,NA"), bam_csv)
    expect_message(merged <- merge(), "\\[bam_stats warn\\] DNA: fraction_of_reads_aligned=1.02 is > 1; ignoring")
    expect_equal(merged$DNA$fraction_of_reads_aligned, 0.9)
    expect_equal(merged$DNA$pf_reads_aligned, 900)
    no_api <- merge_qc_sources(plugin$myeloseqer_qc_sources(coverage_csv = c(report_csv, coverage_csv)),
        file.path(root, "no_api.json"), derived = plugin$MYELOSEQER_QC_DERIVED,
        field_order = plugin$MYELOSEQER_QC_FIELD_ORDER)
    expect_equal(no_api$DNA$mean_depth, 650)
    expect_false("percent_aligned" %in% names(no_api$DNA))
})

test_that("datafiles retains legacy core omission and rounding without simplifying provenance", {
    root <- tempfile()
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    records <- list(
        list(pair = "P1", tumor_type = "NA", disease = "NA", explicit_null = NULL,
             patient_id = "00123", value = 1.23456789, zero = 0, flag = FALSE,
             coverage_qc = list(mean_depth = 701.234567),
             qcMetrics = list(list(code = "PASS", title = "single flag")),
             tumor_details = "NA annotation"),
        list(pair = "P2", tumor_type = "AML", disease = "NA", explicit_null = NULL,
             patient_id = "00987", value = 9.87654321, zero = 0, flag = FALSE,
             coverage_qc = list(mean_depth = 999.876543),
             qcMetrics = list(list(code = "WARN", title = "another flag")),
             tumor_details = "NA annotation")
    )
    provenance <- list(schemaVersion = 1L, caseId = "P1", text = "NA",
        precision = 0.123456789123, empty = list(), absent = NULL,
        artifacts = list(list(path = "one", sources = list(list(path = "input", sha256 = "NA"))),
                         list(path = "two", values = list(TRUE, 2.5, "NA"))))
    legacy <- vector("list", length(records))
    input_paths <- character(length(records))
    for (i in seq_along(records)) {
        pair_dir <- file.path(root, records[[i]]$pair)
        dir.create(pair_dir)
        input_paths[i] <- file.path(pair_dir, "metadata.json")
        core <- jsonlite::toJSON(list(records[[i]]), auto_unbox = TRUE, null = "null", digits = NA)
        legacy[[i]] <- jsonlite::unbox(jsonlite::fromJSON(core))
        record <- records[[i]]
        if (i == 1L) record$provenance <- provenance
        jsonlite::write_json(list(record), input_paths[i], auto_unbox = TRUE, null = "null", digits = NA)
    }
    expected <- jsonlite::fromJSON(jsonlite::toJSON(legacy, auto_unbox = TRUE,
        pretty = TRUE, null = "null"), simplifyVector = FALSE)
    before <- lapply(input_paths, readLines)
    lift_datafiles_json(root, cores = 1L)
    output <- jsonlite::fromJSON(file.path(root, "datafiles.json"), simplifyVector = FALSE)
    expect_length(output, 2L)
    expect_identical(output[[1]]$provenance, provenance)
    core_output <- output
    core_output[[1]]$provenance <- NULL
    expect_identical(core_output, expected)
    expect_false(any(c("tumor_type", "disease", "explicit_null") %in% names(output[[1]])))
    expect_identical(output[[2]]$tumor_type, "AML")
    expect_identical(output[[1]]$patient_id, "00123")
    expect_identical(output[[1]]$tumor_details, "NA annotation")
    expect_identical(output[[1]]$zero, 0L)
    expect_identical(output[[1]]$flag, FALSE)
    expect_equal(output[[1]]$value, 1.2346)
    expect_identical(lapply(input_paths, readLines), before)
    expect_true(file.exists(file.path(root, "datafiles.arrow")))
})
