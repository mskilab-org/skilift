#' Describe a source of sample-keyed QC measurements
#'
#' Sources are merged in list order. Paths may be a vector; duplicate sample
#' records within one source use the last record, before applying precedence.
#' A field transform returning NULL omits an unavailable measurement; NA is an
#' explicit missing value. This lets plugins retain a prior usable measurement.
#' @param path Input path or paths.
#' @param sample_key JSON name path resolving to a sample id, or a table column.
#' @param fields List of annotation_field specifications.
#' @param kind JSON, CSV, or Picard metrics.
#' @param locate Optional function(object, key) returning a JSON record, or
#'   function(row, key) selecting/transforming a tabular record; NULL skips it.
#' @param precedence Whether to overwrite or only fill NULL/numeric NA values.
#' @param key_transform Optional function(key) normalizing the output sample id.
#' @return A QC source specification.
#' @export
qc_source <- function(path, sample_key, fields, kind = c("json", "csv", "picard"),
                      locate = NULL, precedence = c("overwrite", "fill"),
                      key_transform = NULL) {
    if (is.character(fields)) fields <- lapply(fields, annotation_field)
    structure(list(path = path, sample_key = sample_key, fields = fields,
                   kind = match.arg(kind), locate = locate,
                   precedence = match.arg(precedence), key_transform = key_transform),
              class = "qc_source")
}

#' Read a Picard metrics table, ignoring comments and empty lines
#'
#' Uses the reader from process_qc_metrics without changing that legacy reader.
#' @param path Path to a Picard metrics file.
#' @return A data.table of metrics.
#' @export
read_picard_metrics <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path)) {
        stop("Invalid Picard metrics path")
    }
    provenance_note(path)
    fcon <- file(path, "r")
    on.exit(close(fcon))
    txt <- character(0)
    l <- readLines(fcon, 1)
    nr <- NROW(l)
    while (nr > 0) {
        is_comment <- startsWith(l, "#")
        is_empty <- all(nr == 0 || !nzchar(l))
        if (!is_comment && !is_empty) txt <- c(txt, l)
        l <- readLines(fcon, 1)
        nr <- NROW(l)
    }
    data.table::setDT(data.table::fread(text = txt, header = TRUE, sep = "\t"))
}

qc_is_missing <- function(x) {
    is.null(x) || (is.numeric(x) && length(x) == 1L && is.na(x))
}

qc_source_records <- function(source) {
    records <- list()
    add_record <- function(record, key) {
        if (!is.null(source$key_transform)) key <- source$key_transform(key)
        if (is.null(record) || length(key) != 1L || is.na(key) || !nzchar(key) || key == "none") {
            return(invisible(NULL))
        }
        fields <- list()
        for (field in source$fields) {
            value <- record[[field$source]]
            if (!is.null(field$transform)) value <- field$transform(value)
            if (!is.null(value) && length(value)) fields[field$as] <- list(value)
        }
        records[[as.character(key)]] <<- fields
        invisible(NULL)
    }
    for (path in source$path) {
        provenance_note(path)
        if (source$kind == "json") {
            object <- jsonlite::fromJSON(path, simplifyVector = FALSE)
            key <- object
            for (name in source$sample_key) {
                if (is.null(key)) break
                key <- key[[name]]
            }
            if (is.null(key)) {
                message(sprintf("[qc_json skip] %s: no %s field", path,
                                paste(source$sample_key, collapse = ".")))
                next
            }
            record <- if (is.null(source$locate)) object else source$locate(object, key)
            add_record(record, key)
        } else {
            table <- if (source$kind == "picard") read_picard_metrics(path) else
                utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
            if (!(source$sample_key %in% names(table))) {
                message(sprintf("[qc skip] %s: no '%s' column", path, source$sample_key))
                next
            }
            for (i in seq_len(nrow(table))) {
                key <- trimws(as.character(table[[source$sample_key]][i]))
                if (is.na(key) || key == "" || key == "none") next
                record <- as.list(table[i, , drop = FALSE])
                if (!is.null(source$locate)) record <- source$locate(record, key)
                add_record(record, key)
            }
        }
    }
    records
}

#' Merge QC sources into an ordered, sample-keyed JSON file
#' @param sources Ordered list of qc_source specifications.
#' @param out_path Output JSON path.
#' @param derived Named functions of the post-merge record; NULL omits a field.
#' @param field_order Required output field order; other keys follow in encounter order.
#' @return The merged records, invisibly.
#' @export
merge_qc_sources <- function(sources, out_path, derived = list(), field_order = NULL) {
    if (is.null(field_order)) stop("field_order is required for stable QC output")
    records <- list()
    for (source in sources) {
        incoming <- qc_source_records(source)
        for (key in names(incoming)) {
            record <- records[[key]]
            if (is.null(record)) record <- list()
            for (name in names(incoming[[key]])) {
                if (source$precedence == "overwrite" || qc_is_missing(record[[name]])) {
                    record[name] <- incoming[[key]][name]
                }
            }
            records[[key]] <- record
        }
    }
    for (key in names(records)) {
        record <- records[[key]]
        for (name in names(derived)) {
            value <- derived[[name]](record)
            if (!is.null(value) && length(value)) record[name] <- list(value)
        }
        records[[key]] <- record[c(intersect(field_order, names(record)),
                                  setdiff(names(record), field_order))]
    }
    jsonlite::write_json(records, out_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    invisible(records)
}

#' Add premerged QC before metadata QC flag calculation
#' @param metadata A metadata row.
#' @param qc_merged Path to the sample-keyed JSON emitted by merge_qc_sources.
#' @param sample_key Key identifying the sample record.
#' @param median_coverage_field Record field used for tumor_median_coverage.
#' @return Updated metadata, or unchanged metadata when no record matches.
#' @export
add_coverage_metrics_from_qc_file <- function(metadata, qc_merged, sample_key,
                                             median_coverage_field = "mean_depth") {
    if (is.null(qc_merged) || !length(qc_merged) || is.na(qc_merged) || !nzchar(qc_merged)) {
        return(metadata)
    }
    provenance_note(qc_merged)
    records <- jsonlite::fromJSON(qc_merged, simplifyVector = FALSE)
    if (is.null(sample_key) || length(sample_key) != 1L || is.na(sample_key)) return(metadata)
    record <- records[[sample_key]]
    if (is.null(record)) return(metadata)
    record$pair <- NULL
    # set's outer list selects the column; the inner list is its single row.
    # Unlike `$<- list(record)`, this does not unwrap a one-field record.
    if (data.table::is.data.table(metadata)) {
        data.table::set(metadata, j = "coverage_qc", value = list(list(as.list(record))))
    } else {
        metadata$coverage_qc <- list(as.list(record))
    }
    if (!is.null(record[[median_coverage_field]])) {
        metadata$tumor_median_coverage <- record[[median_coverage_field]]
    }
    metadata
}

# The generic lifter accepts explicit keys. The IR fallback retains the legacy
# DNA-id convention for cohorts whose displayed pair is a shortened TM id.
qc_sample_key_from_row <- function(row) {
    key <- row[["qc_sample_key"]]
    if (length(key) == 1L && !is.na(key) && nzchar(key)) return(as.character(key))
    candidate <- row[["ir_analysis_id"]]
    if (length(candidate) != 1L || is.na(candidate) || !nzchar(candidate)) return(row[["pair"]])
    match <- regmatches(candidate, regexec("^(.+?)_v[0-9]+_.+_RNA_v[0-9]+$", candidate, perl = TRUE))[[1]]
    if (length(match) >= 2L) match[2] else row[["pair"]]
}
