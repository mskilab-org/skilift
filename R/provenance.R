# Production-time provenance is scoped to one sample worker. Readers record
# inputs before consuming them; writers publish receipts only after success.
# Worker-local environments and per-artifact atomic receipts avoid concurrent
# writes to metadata.json. The finalizer folds receipts into the existing row.

provenance_load_cohort <- function(path) {
    cohort <- Skilift::Cohort$new(base::readRDS(path))
    cohort$provenance_cached_cohort <- list(path = path, fingerprint = provenance_fingerprint(path),
                                          sha256 = digest::digest(file = path, algo = "sha256"))
    cohort
}

provenance_apply_remapping <- function(cohort, remapping) {
    stopifnot(remapping$schemaVersion == 1L, length(remapping$changes) > 0L)
    index <- which(cohort$inputs$pair == remapping$caseId)
    if (length(index) != 1L) stop("Remapping must identify exactly one cached case")
    columns <- vapply(remapping$changes, function(x) x$column, character(1))
    if (anyDuplicated(columns)) stop("Only one replacement is allowed per input column")
    # Validate the entire proposal before applying any changes.
    for (change in remapping$changes) {
        if (!change$column %in% names(cohort$inputs) ||
            !identical(cohort$inputs[[change$column]][[index]], change$fromPath)) {
            stop("Recorded input no longer matches remapping: ", change$column)
        }
        if (!is.character(change$toPath) || length(change$toPath) != 1L || !file.exists(change$toPath)) {
            stop("Replacement is unavailable: ", change$column)
        }
    }
    for (change in remapping$changes) data.table::set(cohort$inputs, i = index, j = change$column, value = change$toPath)
    cohort$provenance_remapping <- remapping
    invisible(cohort)
}

provenance_begin_upload <- function(cohort, output_data_dir, arguments = list()) {
    # Older serialized R6 objects predate these declared fields.
    if (!"provenance_upload" %in% names(cohort)) cohort <- Skilift::Cohort$new(cohort)
    root <- file.path(output_data_dir, ".uploads")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    directory <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-"), tmpdir = root)
    dir.create(directory)
    snapshot <- file.path(directory, "cohort.rds")
    base::saveRDS(cohort, snapshot)
    rows <- lapply(seq_len(nrow(cohort$inputs)), function(i) as.list(cohort$inputs[i]))
    provenance_atomic_json(rows, file.path(directory, "inputs.json"))
    record <- list(id = basename(directory), startedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
                   cachedCohort = cohort$provenance_cached_cohort,
                   remapping = cohort$provenance_remapping,
                   snapshot = list(path = snapshot, sha256 = digest::digest(file = snapshot, algo = "sha256")),
                   reference = cohort$reference_name, study = cohort$type,
                   arguments = arguments,
                   producer = list(package = "Skilift", version = as.character(utils::packageVersion("Skilift"))),
                   options = options()[grep("^skilift", names(options()))])
    provenance_atomic_json(record, file.path(directory, "upload.json"))
    cohort$provenance_upload <- record
    cohort
}

provenance_context <- function(row, cohort = NULL) {
    context <- new.env(parent = emptyenv())
    context$row <- as.list(row)
    context$sources <- list()
    context$fields <- list()
    context$startedAt <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
    context$reference <- if (!is.null(cohort)) cohort$reference_name else NULL
    context$inputTable <- if (!is.null(cohort)) cohort$provenance_input_path else NULL
    context$columnMapping <- if (!is.null(cohort)) cohort$provenance_column_mapping else NULL
    context$upload <- if (!is.null(cohort)) cohort$provenance_upload else NULL
    context$selection <- list()
    context$parameters <- list()
    context
}

provenance_current <- function() {
    for (frame in rev(sys.frames())) {
        if (exists(".provenance_context", envir = frame, inherits = FALSE)) {
            return(get(".provenance_context", envir = frame, inherits = FALSE))
        }
    }
    NULL
}

provenance_fingerprint <- function(path) {
    info <- file.info(path)
    if (is.na(info$size[1])) return(list(state = "missing"))
    list(state = "present", resolvedPath = normalizePath(path, mustWork = FALSE),
         size = unname(info$size[1]), mtime = as.numeric(info$mtime[1]))
}

provenance_note <- function(path) {
    context <- provenance_current()
    if (is.null(context) || !is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path)) {
        return(invisible(NULL))
    }
    supplied <- path.expand(path)
    if (!startsWith(supplied, "/")) supplied <- file.path(getwd(), supplied)
    role <- names(context$row)[vapply(context$row, function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) && identical(value, path)
    }, logical(1))]
    key <- supplied
    if (is.null(context$sources[[key]])) {
        fingerprint <- provenance_fingerprint(supplied)
        context$sources[[key]] <- list(
            id = supplied, path = supplied,
            role = if (length(role)) paste(role, collapse = ", ") else "resource",
            tablePath = context$inputTable,
            originalColumn = if (length(role)) paste(unlist(context$columnMapping[role]), collapse = ", ") else NULL,
            cohortColumns = unname(as.list(role)),
            fingerprintAtUse = fingerprint
        )
        # Carry the producer snapshot of transformed inputs (e.g. oncotable)
        # only when it still describes the file actually being consumed.
        receipt <- file.path(dirname(supplied), ".provenance", paste0(basename(supplied), ".json"))
        if (file.exists(receipt)) {
            upstream <- jsonlite::read_json(receipt, simplifyVector = FALSE)
            used <- upstream$outputFingerprint
            if (identical(used$resolvedPath, fingerprint$resolvedPath) &&
                identical(as.numeric(used$size), fingerprint$size) &&
                isTRUE(abs(used$mtime - fingerprint$mtime) < 0.00001)) {
                context$sources[[key]]$production <- upstream
            }
        }
    }
    invisible(NULL)
}

provenance_literal <- function(field, value, origin = "cohort-input") {
    context <- provenance_current()
    if (!is.null(context)) {
        context$fields[[field]] <- list(sources = list(), evidence = "literal",
                                       value = value, origin = origin)
    }
    invisible(NULL)
}

# These local bindings preserve the reader interfaces used by skilift. Calls
# outside a sample lifter behave identically and do not create receipts.
readRDS <- function(file, ...) {
    provenance_note(file)
    base::readRDS(file, ...)
}

fread <- function(input = "", ...) {
    provenance_note(input)
    data.table::fread(input = input, ...)
}

readLines <- function(con = stdin(), ...) {
    provenance_note(con)
    base::readLines(con, ...)
}

readVcf <- function(file, ...) {
    provenance_note(file)
    VariantAnnotation::readVcf(file, ...)
}

provenance_step <- function(metadata, expression, metadata_inputs = character()) {
    context <- provenance_current()
    if (is.null(context)) return(force(expression))
    previous <- context$sources
    context$sources <- list()
    on.exit({ context$sources <- utils::modifyList(previous, context$sources) }, add = TRUE)
    before <- data.table::copy(metadata)
    dependencies <- context$fields[intersect(metadata_inputs, names(context$fields))]
    for (dependency in dependencies) {
        for (source in dependency$sources) context$sources[[source$path]] <- source
    }
    result <- force(expression)
    changed <- names(result)[vapply(names(result), function(field) {
        !identical(before[[field]], result[[field]])
    }, logical(1))]
    for (field in changed) {
        context$fields[[field]] <- list(sources = unname(context$sources),
                                      evidence = if (length(context$sources)) "captured" else "literal",
                                      value = result[[field]], origin = "metadata-transform")
        if (length(metadata_inputs)) {
            context$fields[[field]]$metadataInputs <- as.list(before[, intersect(metadata_inputs, names(before)), with = FALSE])
        }
    }
    result
}

provenance_atomic_json <- function(value, path, ...) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(pattern = ".provenance-", tmpdir = dirname(path))
    on.exit(unlink(temporary), add = TRUE)
    jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA, ...)
    if (!file.rename(temporary, path)) stop("Could not publish provenance: ", path)
}

provenance_publish_receipt <- function(record, path) {
    history <- file.path(dirname(path), "history")
    dir.create(history, recursive = TRUE, showWarnings = FALSE)
    archive <- tempfile(pattern = "production-", tmpdir = history, fileext = ".json")
    record$recordId <- basename(archive)
    provenance_atomic_json(record, archive)
    provenance_atomic_json(record, path)
}

provenance_receipt <- function(path) {
    context <- provenance_current()
    if (is.null(context) || !file.exists(path)) return(invisible(NULL))
    row <- context$row
    params <- utils::modifyList(row[intersect(names(row), Skilift::config_parameter_names)], context$parameters)
    record <- list(
        schemaVersion = 1L, path = basename(path), evidence = "captured",
        sourceScope = "reader-session", caseId = as.character(row$pair),
        producedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        startedAt = context$startedAt, reference = context$reference,
        producer = list(package = "Skilift", version = as.character(utils::packageVersion("Skilift"))),
        inputTable = context$inputTable,
        upload = context$upload, selection = context$selection,
        parameters = params, sources = unname(context$sources),
        outputFingerprint = provenance_fingerprint(path)
    )
    receipt_dir <- file.path(dirname(path), ".provenance")
    if (basename(path) == "metadata.json") {
        # Snapshot values separately from file mtimes, which finalization changes.
        values <- jsonlite::read_json(path, simplifyVector = FALSE)[[1]]
        values$provenance <- NULL
        record$outputValues <- values
    }
    provenance_publish_receipt(record, file.path(receipt_dir, paste0(basename(path), ".json")))
    if (basename(path) == "metadata.json") {
        for (field in names(context$fields)) {
            field_record <- utils::modifyList(record, context$fields[[field]], keep.null = TRUE)
            # modifyList merges unnamed lists; explicitly replace the input list.
            field_record$sources <- context$fields[[field]]$sources
            field_record$value <- context$fields[[field]]$value
            field_record$outputValues <- NULL
            field_record$outputValue <- values[[field]]
            field_record$origin <- context$fields[[field]]$origin
            field_record$pointer <- paste0("/0/", gsub("/", "~1", gsub("~", "~0", field, fixed = TRUE), fixed = TRUE))
            field_record$sourceScope <- "metadata-reader"
            provenance_publish_receipt(field_record, file.path(receipt_dir, paste0("field-", field, ".json")))
        }
    }
    invisible(NULL)
}

provenance_write_json <- function(x, path, ...) {
    if (basename(path) == "metadata.json" && file.exists(path)) {
        previous <- jsonlite::read_json(path, simplifyVector = FALSE)
        if (length(previous) && !is.null(previous[[1]]$provenance) && inherits(x, "data.frame")) {
            x$provenance <- list(previous[[1]]$provenance)
        }
    }
    jsonlite::write_json(x, path, ...)
    provenance_receipt(path)
    invisible(NULL)
}

ggsave <- function(filename, plot = ggplot2::last_plot(), ..., path = NULL) {
    result <- ggplot2::ggsave(filename = filename, plot = plot, ..., path = path)
    provenance_receipt(if (is.null(path)) filename else file.path(path, filename))
    invisible(result)
}

provenance_finalize <- function(output_data_dir) {
    for (case_dir in list.dirs(output_data_dir, recursive = FALSE, full.names = TRUE)) {
        path <- file.path(case_dir, "metadata.json")
        receipts <- list.files(file.path(case_dir, ".provenance"), pattern = "\\.json$", full.names = TRUE)
        if (!file.exists(path) || !length(receipts)) next
        metadata <- jsonlite::read_json(path, simplifyVector = FALSE)
        if (!length(metadata)) next
        provenance <- metadata[[1]]$provenance
        if (is.null(provenance)) provenance <- list(schemaVersion = 1L, caseId = basename(case_dir), artifacts = list())
        records <- lapply(receipts, jsonlite::read_json, simplifyVector = FALSE)
        key <- function(x) paste(x$path, if (is.null(x$pointer)) "" else x$pointer, sep = "#")
        old <- provenance$artifacts
        # A full metadata write supersedes all previous field records.
        fresh_metadata <- Filter(function(x) identical(x$path, "metadata.json") && is.null(x$pointer), records)
        if (length(fresh_metadata)) {
            produced <- fresh_metadata[[1]]$producedAt
            records <- Filter(function(x) !identical(x$path, "metadata.json") || identical(x$producedAt, produced), records)
            old <- Filter(function(x) !identical(x$path, "metadata.json"), old)
        }
        keep <- !vapply(old, key, character(1)) %in% vapply(records, key, character(1))
        provenance$artifacts <- c(old[keep], records)
        provenance$recordedAt <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
        metadata[[1]]$provenance <- provenance
        provenance_atomic_json(metadata, path)
    }
    invisible(NULL)
}
