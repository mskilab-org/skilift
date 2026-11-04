# Simple per-Case provenance helpers. Provenance is kept outside data payloads.

provenance_timestamp <- function(time = Sys.time()) {
  format(time, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC", usetz = FALSE)
}

provenance_os_account <- function() {
  account <- tryCatch(unname(Sys.info()[["user"]]), error = function(e) NULL)
  if (is.null(account) || !length(account) || is.na(account) || !nzchar(account)) {
    account <- Sys.getenv("USER", unset = "")
  }
  if (!nzchar(account)) account <- Sys.getenv("USERNAME", unset = "")
  if (!nzchar(account)) return(NULL)
  account
}

provenance_run_context <- function(context = NULL) {
  if (!is.null(context)) return(context)
  run <- list(started_at = provenance_timestamp())
  operator <- provenance_os_account()
  if (!is.null(operator)) run$operator <- operator
  run
}

provenance_new_context <- function(output_data_dir, enabled = TRUE, run = NULL) {
  if (!isTRUE(enabled)) return(NULL)
  if (missing(output_data_dir) || length(output_data_dir) != 1L || is.na(output_data_dir)) {
    stop("output_data_dir must be a single path when provenance is enabled", call. = FALSE)
  }
  if (!dir.exists(output_data_dir) && !dir.create(output_data_dir, recursive = TRUE)) {
    stop("Could not create output directory for provenance: ", output_data_dir, call. = FALSE)
  }
  context <- new.env(parent = emptyenv())
  context$output_data_dir <- normalizePath(output_data_dir, winslash = "/", mustWork = TRUE)
  context$run <- provenance_run_context(run)
  context$spool_dir <- tempfile("skilift-provenance-")
  if (!dir.create(context$spool_dir, recursive = TRUE)) {
    stop("Could not create provenance staging directory", call. = FALSE)
  }
  context$closed <- FALSE
  context
}

provenance_acquire_context <- function(output_data_dir, enabled = TRUE, context = NULL) {
  owns <- is.null(context) && isTRUE(enabled)
  if (owns) context <- provenance_new_context(output_data_dir, enabled = TRUE)
  list(context = context, owns = owns)
}

provenance_safe_pair <- function(pair) {
  if (length(pair) != 1L || is.na(pair) || !nzchar(pair) ||
      pair %in% c(".", "..") || grepl("[/\\\\]", pair) ||
      grepl("^[A-Za-z]:", pair)) {
    stop("Invalid Case identifier for provenance", call. = FALSE)
  }
  pair
}

provenance_safe_artifact_key <- function(key) {
  if (length(key) != 1L || is.na(key) || !nzchar(key) ||
      grepl("^[/\\\\]", key) || grepl("^[A-Za-z]:", key) ||
      grepl("\\\\", key) || any(strsplit(key, "/", fixed = TRUE)[[1]] %in% c("", ".", "..")) ||
      identical(key, "provenance.json")) {
    stop("Invalid Case-relative provenance artifact key", call. = FALSE)
  }
  key
}

provenance_is_uri <- function(path) {
  grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)
}

provenance_safe_path <- function(path) {
  # Do not persist URL user-info, query strings or fragments that may contain
  # credentials or signed tokens.
  if (!provenance_is_uri(path)) return(path)
  path <- sub("([?#]).*$", "", path)
  sub("(://)[^/@]*@", "\\1", path)
}

provenance_absolute_path <- function(path) {
  path <- path.expand(path)
  if (provenance_is_uri(path)) return(path)
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (is_absolute) path else file.path(getwd(), path)
}

provenance_format_time <- function(value) {
  if (length(value) != 1L || is.na(value)) return(NULL)
  provenance_timestamp(as.POSIXct(value, tz = "UTC"))
}

provenance_observation <- function(path, output = FALSE) {
  observed_at <- provenance_timestamp()
  supplied <- provenance_safe_path(as.character(path)[1])
  absolute <- provenance_absolute_path(supplied)
  uri <- provenance_is_uri(supplied)
  file_path <- if (uri) NA_character_ else path.expand(supplied)
  info <- if (!uri) tryCatch(file.info(file_path), error = function(e) NULL) else NULL
  is_file <- !is.null(info) && nrow(info) == 1L && !is.na(info$isdir) && !info$isdir
  status <- if (uri) "unsupported" else if (is_file) "ok" else "unavailable"
  size <- if (is_file && !is.na(info$size)) unname(as.numeric(info$size)) else NULL
  modified <- if (is_file && !is.na(info$mtime)) provenance_format_time(info$mtime) else NULL
  resolved <- NULL
  if (is_file) {
    resolved <- tryCatch(
      normalizePath(file_path, winslash = "/", mustWork = TRUE),
      error = function(e) NULL
    )
  }

  if (isTRUE(output)) {
    return(list(
      observed_at = observed_at,
      size_bytes = size,
      modified_at = modified
    ))
  }
  list(
    supplied_path = supplied,
    absolute_path = absolute,
    resolved_path = resolved,
    observed_at = observed_at,
    stat_status = status,
    size_bytes = size,
    modified_at = modified,
    created_at = NULL
  )
}

provenance_file_value <- function(value) {
  # A selected file column may be a factor or a single-cell list column.
  # Never flatten tables/objects or stringify clinical values into paths.
  if (is.list(value) && !is.data.frame(value) && is.null(dim(value)) &&
      is.null(names(value)) && length(value) == 1L) value <- value[[1L]]
  if (is.factor(value)) value <- as.character(value)
  if (!is.character(value) || !is.null(dim(value)) || length(value) != 1L) {
    return(character())
  }
  if (is.na(value) || !nzchar(trimws(value))) return(character())
  value
}

provenance_direct_binding <- function(path, role) {
  path <- provenance_file_value(path)
  if (!length(path)) return(list())
  observation <- provenance_observation(path)
  observation$role <- role
  observation["cohort_column"] <- list(NULL)
  list(observation)
}

provenance_input_bindings <- function(row, columns) {
  if (is.null(columns) || !length(columns)) return(list())
  roles <- names(columns)
  columns <- as.character(columns)
  if (is.null(roles)) roles <- columns
  inputs <- list()
  for (i in seq_along(columns)) {
    column <- columns[[i]]
    if (is.null(row[[column]])) next
    values <- provenance_file_value(row[[column]])
    for (path in values) {
      observation <- provenance_observation(path)
      observation$role <- roles[[i]]
      observation$cohort_column <- column
      inputs[[length(inputs) + 1L]] <- observation
    }
  }
  inputs
}

provenance_entry <- function(lifter, run, inputs, output_observation) {
  entry <- list(lifter = lifter, run = run,
                output_observation = output_observation,
                inputs = inputs)
  entry
}

provenance_read_document <- function(path, pair) {
  if (!file.exists(path)) return(NULL)
  document <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) stop("Cannot read existing provenance document ", path,
                             ": ", conditionMessage(e), call. = FALSE)
  )
  if (!is.list(document) || !identical(document$schema_version, 1L) &&
      !identical(document$schema_version, 1)) {
    stop("Unsupported or malformed provenance document: ", path, call. = FALSE)
  }
  if (!identical(as.character(document$pair), pair) || !is.list(document$outputs)) {
    stop("Existing provenance document has invalid Case identity or outputs: ", path,
         call. = FALSE)
  }
  document
}

provenance_create_document <- function(pair, updates, previous = NULL) {
  pair <- provenance_safe_pair(as.character(pair)[1])
  if (is.null(previous)) {
    previous <- list(
      schema_version = 1L,
      pair = pair,
      scope = "direct-file-inputs",
      outputs = list()
    )
  }
  if (!is.list(previous$outputs)) {
    stop("Existing provenance outputs must be a JSON object", call. = FALSE)
  }
  for (key in names(updates)) {
    provenance_safe_artifact_key(key)
    previous$outputs[[key]] <- updates[[key]]
  }
  previous$schema_version <- 1L
  previous$pair <- pair
  if (is.null(previous$scope)) previous$scope <- "direct-file-inputs"
  previous
}

provenance_write_json_temp <- function(document, directory) {
  path <- tempfile(".provenance-stage-", tmpdir = directory, fileext = ".json")
  tryCatch(
    jsonlite::write_json(
      document, path, auto_unbox = TRUE, pretty = TRUE,
      null = "null", na = "null", digits = NA
    ),
    error = function(e) {
      unlink(path)
      stop(e)
    }
  )
  if (!file.exists(path)) stop("Provenance JSON was not written", call. = FALSE)
  path
}

provenance_log_failure <- function(pair, error) {
  message("Skilift provenance failed for Case ", pair, ": ", conditionMessage(error))
  invisible(NULL)
}

provenance_stage_output <- function(context, final_path) {
  if (is.null(context)) return(final_path)
  if (isTRUE(context$closed)) stop("Provenance context is already closed", call. = FALSE)
  extension <- tools::file_ext(basename(final_path))
  tempfile(
    ".skilift-output-stage-",
    tmpdir = dirname(final_path),
    fileext = if (nzchar(extension)) paste0(".", extension) else ""
  )
}

provenance_queue_output <- function(context, pair, artifact_key, lifter,
                                    inputs, staged_path, final_path) {
  if (is.null(context)) return(invisible(NULL))
  pair <- provenance_safe_pair(pair)
  artifact_key <- provenance_safe_artifact_key(artifact_key)
  if (!file.exists(staged_path) || dir.exists(staged_path)) {
    stop("Staged output is missing: ", staged_path, call. = FALSE)
  }
  record <- list(
    pair = pair,
    artifact_key = artifact_key,
    case_dir = file.path(context$output_data_dir, pair),
    staged_path = staged_path,
    final_path = final_path,
    entry = provenance_entry(
      lifter = lifter,
      run = context$run,
      inputs = inputs,
      output_observation = provenance_observation(staged_path, output = TRUE)
    )
  )
  queue_path <- tempfile("record-", tmpdir = context$spool_dir, fileext = ".rds")
  saveRDS(record, queue_path)
  invisible(NULL)
}

provenance_restore_files <- function(installed, backups) {
  for (path in rev(installed)) {
    if (file.exists(path) || dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  }
  for (path in names(backups)) {
    backup <- backups[[path]]
    if (file.exists(backup) || dir.exists(backup)) {
      if (file.exists(path) || dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
      file.rename(backup, path)
    }
  }
}

provenance_publish_case <- function(pair, records) {
  pair <- provenance_safe_pair(pair)
  case_dir <- unique(vapply(records, `[[`, character(1), "case_dir"))
  if (length(case_dir) != 1L || !dir.exists(case_dir)) {
    stop("Invalid output directory for Case ", pair, call. = FALSE)
  }
  keys <- vapply(records, `[[`, character(1), "artifact_key")
  paths <- vapply(records, `[[`, character(1), "final_path")
  if (anyDuplicated(keys) || anyDuplicated(paths)) {
    stop("More than one staged replacement targets an artifact for Case ", pair,
         call. = FALSE)
  }
  for (record in records) {
    if (!file.exists(record$staged_path) || dir.exists(record$staged_path)) {
      stop("Staged output is missing: ", record$staged_path, call. = FALSE)
    }
  }

  manifest <- file.path(case_dir, "provenance.json")
  previous <- provenance_read_document(manifest, pair)
  updates <- setNames(lapply(records, `[[`, "entry"), keys)
  document <- provenance_create_document(pair, updates, previous)
  manifest_stage <- provenance_write_json_temp(document, case_dir)
  on.exit(unlink(manifest_stage), add = TRUE)

  backups <- list()
  installed <- character()
  staged_paths <- c(
    vapply(records, `[[`, character(1), "staged_path"),
    manifest_stage
  )
  final_paths <- c(paths, manifest)
  ok <- FALSE
  tryCatch({
    for (i in seq_along(final_paths)) {
      final <- final_paths[[i]]
      staged <- staged_paths[[i]]
      if (file.exists(final) || dir.exists(final)) {
        backup <- tempfile(".skilift-backup-", tmpdir = dirname(final))
        if (!file.rename(final, backup)) {
          stop("Could not preserve previous file before replacing: ", final,
               call. = FALSE)
        }
        backups[[final]] <- backup
      }
      if (!file.rename(staged, final)) {
        stop("Could not publish staged file: ", final, call. = FALSE)
      }
      installed <- c(installed, final)
    }
    ok <- TRUE
  }, error = function(e) {
    provenance_restore_files(installed, backups)
    stop(e)
  })

  if (ok) {
    if (length(backups)) unlink(unlist(backups, use.names = FALSE), force = TRUE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

provenance_finish <- function(context) {
  if (is.null(context) || isTRUE(context$closed)) return(invisible(NULL))
  context$closed <- TRUE
  queued <- list.files(context$spool_dir, pattern = "\\.rds$", full.names = TRUE)
  records <- lapply(queued, function(path) tryCatch(readRDS(path), error = function(e) e))
  valid <- !vapply(records, inherits, logical(1), "error")
  if (any(!valid)) {
    message("Skilift provenance: could not read one or more staged Case records")
  }
  records <- records[valid]
  if (length(records)) {
    pairs <- unique(vapply(records, `[[`, character(1), "pair"))
    for (pair in pairs) {
      case_records <- records[vapply(records, function(x) identical(x$pair, pair), logical(1))]
      tryCatch(
        provenance_publish_case(pair, case_records),
        error = function(e) {
          for (record in case_records) {
            if (file.exists(record$staged_path)) unlink(record$staged_path)
          }
          provenance_log_failure(pair, e)
        }
      )
    }
  }
  unlink(context$spool_dir, recursive = TRUE, force = TRUE)
  invisible(NULL)
}

# Explicit output/input bindings used by legacy-output backfill. Extend this
# registry only when the lifter's actual selection and output key are understood.
provenance_backfill_registry <- list(
  "coverage.arrow" = list(
    lifter = "lift_denoised_coverage",
    columns = c(coverage = "tumor_coverage"),
    resources = c(reference_settings = ".default_settings"),
    conditional_resources = list(
      mask = list(column = "denoised_coverage_apply_mask", path = ".mask")
    )
  ),
  "hetsnps.arrow" = list(
    lifter = "lift_hetsnps",
    columns = c(het_pileups = "het_pileups"),
    resources = c(
      reference_settings = ".default_settings", scatterplot_mask = ".mask"
    )
  ),
  "tumor.bam" = list(
    lifter = "lift_bam",
    columns = c(tumor_bam = "tumor_bam")
  ),
  "tumor.bam.bai" = list(
    lifter = "lift_bam",
    columns = c(tumor_bai = "tumor_bam_bai")
  ),
  "normal.bam" = list(
    lifter = "lift_bam",
    columns = c(normal_bam = "normal_bam")
  ),
  "normal.bam.bai" = list(
    lifter = "lift_bam",
    columns = c(normal_bai = "normal_bam_bai")
  ),
  "mutations.json" = list(
    lifter = "lift_multiplicity",
    columns = c(multiplicity = "multiplicity", oncokb_snv = "oncokb_snv")
  ),
  "germline_mutations.json" = list(
    lifter = "lift_multiplicity",
    columns = c(multiplicity = "germline_multiplicity")
  ),
  "sage.qc.json" = list(
    lifter = "lift_variant_qc",
    columns = c(
      oncokb_snv = "oncokb_snv", multiplicity = "multiplicity",
      sbs_post_prob_signatures = "sbs_post_prob_signatures",
      indel_post_prob_signatures = "indel_post_prob_signatures"
    )
  ),
  "mutation_catalog.json" = list(
    lifter = "lift_signatures",
    columns = c(sbs = "matrix_sbs_signatures")
  ),
  "sbs_decomposed_prob.json" = list(
    lifter = "lift_signatures",
    columns = c(sbs = "decomposed_sbs_signatures")
  ),
  "id_mutation_catalog.json" = list(
    lifter = "lift_signatures",
    columns = c(indel = "matrix_indel_signatures")
  ),
  "id_decomposed_prob.json" = list(
    lifter = "lift_signatures",
    columns = c(indel = "decomposed_indel_signatures")
  ),
  "somatic_altered_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "multiplicity"),
    resources = c(mask = ".mask")
  ),
  "somatic_total_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "multiplicity"),
    resources = c(mask = ".mask")
  ),
  "germline_altered_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "germline_multiplicity"),
    resources = c(mask = ".mask")
  ),
  "germline_total_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "germline_multiplicity"),
    resources = c(mask = ".mask")
  ),
  "hetsnps_major_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "hetsnps_multiplicity"),
    resources = c(mask = ".mask")
  ),
  "hetsnps_minor_hist.json" = list(
    lifter = "lift_multiplicity_fits",
    columns = c(multiplicity = "hetsnps_multiplicity"),
    resources = c(mask = ".mask")
  ),
  "multiplicity.png" = list(
    lifter = "lift_2d_purity_ploidy_plot",
    columns = c(
      somatic_multiplicity = "multiplicity",
      hetsnps_multiplicity = "hetsnps_multiplicity"
    ),
    resources = c(mask = ".mask")
  ),
  "hetsnps_major_minor.png" = list(
    lifter = "lift_allelic_pp_fit",
    columns = c(graph = ".default_jabba", het_pileups = "het_pileups")
  ),
  "hetsnps_major_minor.json" = list(
    lifter = "lift_allelic_pp_fit",
    columns = c(graph = ".default_jabba", het_pileups = "het_pileups")
  ),
  "ppfit.json" = list(
    lifter = "lift_segment_width_distribution",
    columns = c(graph = ".default_jabba")
  ),
  "complex.json" = list(
    lifter = "lift_copy_number_graph",
    columns = c(graph = ".default_jabba"),
    resources = c(settings = ".default_settings")
  ),
  "allelic.json" = list(
    lifter = "lift_copy_number_graph",
    columns = c(graph = "allelic_jabba_gg"),
    resources = c(settings = ".default_settings")
  ),
  "filtered.events.json" = list(
    lifter = "lift_filtered_events",
    columns = c(oncotable = "oncotable", graph = ".default_jabba"),
    resources = c(gene_locations = ".gene_locations", heme_database = ".heme_database")
  ),
  "highlights.json" = list(
    lifter = "lift_filtered_events",
    columns = c(
      oncotable = "oncotable", graph = ".default_jabba", itdseek = "itdseek_rds"
    ),
    resources = c(heme_database = ".heme_database", duncavage_database = ".duncavage_database")
  ),
  "metadata.json" = list(
    lifter = "lift_metadata",
    columns = c(
      graph = ".default_jabba", events = "events",
      oncokb_snv = "oncokb_snv", somatic_snvs = "somatic_snvs",
      tumor_coverage = "tumor_coverage",
      estimate_library_complexity = "estimate_library_complexity",
      alignment_summary_metrics = "alignment_summary_metrics",
      insert_size_metrics = "insert_size_metrics", tumor_wgs_metrics = "tumor_wgs_metrics",
      normal_wgs_metrics = "normal_wgs_metrics", het_pileups = "het_pileups",
      decomposed_sbs_signatures = "decomposed_sbs_signatures",
      decomposed_indel_signatures = "decomposed_indel_signatures",
      matrix_sbs_signatures = "matrix_sbs_signatures",
      matrix_indel_signatures = "matrix_indel_signatures",
      activities_sbs_signatures = "activities_sbs_signatures",
      activities_indel_signatures = "activities_indel_signatures",
      hrdetect = "hrdetect", onenesstwoness = "onenesstwoness",
      msisensorpro = "msisensorpro", purple_pp_bestFit = "purple_pp_bestFit",
      conpair_contamination = "conpair_contamination", conpair_concordance = "conpair_concordance",
      purple_qc = "purple_qc"
    )
  )
)

provenance_backfill_columns <- function(cohort, columns) {
  roles <- names(columns)
  columns <- as.character(columns)
  names(columns) <- roles
  default_jabba <- columns == ".default_jabba"
  if (any(default_jabba)) {
    columns[default_jabba] <- tryCatch(
      Skilift::DEFAULT_JABBA(object = cohort),
      error = function(e) "jabba_gg"
    )
  }
  columns
}

provenance_metadata_columns <- function(row, columns) {
  graph_column <- unname(columns["graph"])
  graph <- if (length(graph_column) && !is.na(graph_column)) row[[graph_column]] else NULL
  # Match add_purity_ploidy: Purple is read only when no graph is selected.
  if (!is.null(graph) && !any(is.na(graph))) {
    columns <- columns[columns != "purple_pp_bestFit"]
  }
  if (all(c("oncokb_snv", "somatic_snvs") %in% columns)) {
    oncokb <- row[["oncokb_snv"]]
    oncokb_selected <- !is.null(oncokb) && length(oncokb) > 0L &&
      !all(is.na(oncokb))
    if (oncokb_selected) {
      columns <- columns[columns != "somatic_snvs"]
    } else {
      columns <- columns[columns != "oncokb_snv"]
    }
  }
  columns
}

provenance_backfill_inputs <- function(row, cohort, binding) {
  columns <- provenance_backfill_columns(cohort, binding$columns)
  if (identical(binding$lifter, "lift_metadata")) {
    columns <- provenance_metadata_columns(row, columns)
  }
  inputs <- provenance_input_bindings(row, columns)
  if (identical(binding$lifter, "lift_segment_width_distribution")) {
    graph_column <- unname(columns["graph"])
    graph <- if (length(graph_column) && !is.na(graph_column)) row[[graph_column]] else NULL
    ggraph <- tryCatch(Skilift:::process_jabba(graph), error = function(e) NULL)
    if (!is.null(ggraph) && inherits(ggraph, "gGraph") &&
        !Skilift:::segment_width_graph_has_segstats(ggraph)) {
      inputs <- c(inputs, provenance_input_bindings(row, c(coverage = "tumor_coverage")))
    }
  }
  if (identical(binding$lifter, "lift_hetsnps")) {
    filter_mask <- row[["hetsnps_mask"]]
    if (is.null(filter_mask) || !length(filter_mask) || is.na(filter_mask)) {
      filter_mask_path <- system.file("extdata", "data", "maskA_re.rds", package = "Skilift")
      inputs <- c(inputs, provenance_direct_binding(filter_mask_path, "filter_mask"))
    } else {
      inputs <- c(inputs, provenance_input_bindings(row, c(filter_mask = "hetsnps_mask")))
    }
  }
  resources <- binding$resources
  if (length(resources)) {
    paths <- as.character(resources)
    names(paths) <- names(resources)
    paths[paths == ".default_settings"] <- Skilift:::default_settings_path
    paths[paths == ".gene_locations"] <- system.file(
      "extdata", "data", "gene_locations.rds", package = "Skilift"
    )
    paths[paths == ".mask"] <- system.file(
      "extdata", "data", "maskA_re.rds", package = "Skilift"
    )
    paths[paths == ".heme_database"] <- if (identical(cohort$type, "heme")) {
      Skilift:::HEMEDB()
    } else ""
    paths[paths == ".duncavage_database"] <- if (identical(cohort$type, "heme")) {
      Skilift:::DUNCAVAGEDB()
    } else ""
    for (role in names(paths)) {
      inputs <- c(inputs, provenance_direct_binding(paths[[role]], role))
    }
  }
  conditional_resources <- binding$conditional_resources
  if (length(conditional_resources)) {
    for (role in names(conditional_resources)) {
      resource <- conditional_resources[[role]]
      value <- row[[resource$column]]
      enabled <- if (is.null(value) || !length(value)) FALSE else {
        tryCatch(isTRUE(as.logical(value[[1L]])), error = function(e) FALSE)
      }
      if (enabled) {
        resource_path <- resource$path
        if (identical(resource_path, ".mask")) {
          resource_path <- system.file("extdata", "data", "maskA_re.rds", package = "Skilift")
        }
        inputs <- c(inputs, provenance_direct_binding(resource_path, role))
      }
    }
  }
  inputs
}

provenance_write_document_atomic <- function(path, document) {
  directory <- dirname(path)
  temporary <- provenance_write_json_temp(document, directory)
  backup <- NULL
  installed <- FALSE
  tryCatch({
    if (file.exists(path)) {
      backup <- tempfile(".skilift-backup-", tmpdir = directory)
      if (!file.rename(path, backup)) stop("Could not preserve existing manifest", call. = FALSE)
    }
    if (!file.rename(temporary, path)) stop("Could not publish provenance manifest", call. = FALSE)
    installed <- TRUE
  }, error = function(e) {
    if (installed && file.exists(path)) unlink(path)
    if (!is.null(backup) && file.exists(backup)) file.rename(backup, path)
    if (file.exists(temporary)) unlink(temporary)
    stop(e)
  })
  if (!is.null(backup) && file.exists(backup)) unlink(backup)
  invisible(path)
}

#' Generate per-Case provenance for existing Skilift outputs
#'
#' This function backfills provenance from a saved Cohort and an existing output
#' directory. It never reruns a lifter or modifies data outputs. Only output
#' keys in Skilift's documented provenance registry are considered.
#'
#' @details Currently registered keys include coverage, hetsnp, BAM/index,
#'   copy-number graph, segment-width JSON, filtered-event, Heme highlights,
#'   metadata, multiplicity histograms, selected purity/ploidy plots, signatures
#'   and variant-QC files. The segment-width entry inspects the saved graph to
#'   decide whether coverage is also a direct input; it does not rerun the lifter.
#'   Other supported lifters and custom output names may not be covered.
#'
#' @param cohort A Cohort object or path to a serialized Cohort (including legacy
#'   `cohort.rds` files).
#' @param output_data_dir Directory containing the Case output directories.
#' @param force Replace existing Case provenance documents when TRUE. Defaults
#'   to FALSE, which preserves and skips them.
#' @param cores Number of Cases to process at once. Defaults to 1.
#' @return Invisibly, a list with `written`, `skipped`, and `failed` Case IDs.
#' @export
lift_provenance <- function(cohort, output_data_dir, force = FALSE, cores = 1) {
  if (length(force) != 1L || is.na(force) || !is.logical(force)) {
    stop("force must be TRUE or FALSE", call. = FALSE)
  }
  if (length(cores) != 1L || is.na(cores) || cores < 1) {
    stop("cores must be a positive integer", call. = FALSE)
  }
  if (is.character(cohort) && length(cohort) == 1L) {
    if (!file.exists(cohort)) stop("Cohort RDS does not exist: ", cohort, call. = FALSE)
    cohort <- readRDS(cohort)
  }
  if (is.null(cohort$inputs) || !"pair" %in% names(cohort$inputs)) {
    stop("cohort must contain an inputs table with a pair column", call. = FALSE)
  }
  if (!dir.exists(output_data_dir)) {
    stop("output_data_dir does not exist: ", output_data_dir, call. = FALSE)
  }
  output_data_dir <- normalizePath(output_data_dir, winslash = "/", mustWork = TRUE)
  inputs <- cohort$inputs
  # Legacy saved Cohorts can be plain lists. Adapt only the pairing input;
  # leave the saved Cohort, table and existing pairing rules unchanged.
  pairing_source <- if (R6::is.R6(cohort)) cohort else {
    data.table::as.data.table(data.table::copy(inputs))
  }
  metadata_inputs <- tryCatch(
    pairify_cohort_inputs(
      pairing_source, tumor_normal_columns = c("wgs_metrics"),
      keep_remaining = TRUE, sep_cast = "_"
    ),
    error = function(e) NULL
  )
  bam_inputs <- tryCatch({
    paired <- pairify_cohort_inputs(
      pairing_source, tumor_normal_columns = c("bam"), sep = "_"
    )
    for (sample in c("tumor", "normal")) {
      bam_column <- paste0(sample, "_bam")
      index_column <- paste0(bam_column, "_bai")
      if (bam_column %in% names(paired)) {
        bam <- as.character(paired[[bam_column]])
        adjacent <- paste0(bam, ".bai")
        sans_extension <- paste0(tools::file_path_sans_ext(bam), ".bai")
        selected <- ifelse(
          file.exists(adjacent), adjacent,
          ifelse(file.exists(sans_extension), sans_extension, NA_character_)
        )
        paired[[index_column]] <- selected
      }
    }
    paired
  }, error = function(e) NULL)
  pairs <- unique(c(
    as.character(inputs[["pair"]]),
    if (!is.null(metadata_inputs)) as.character(metadata_inputs[["pair"]]),
    if (!is.null(bam_inputs)) as.character(bam_inputs[["pair"]])
  ))
  if (anyNA(pairs) || any(!nzchar(pairs))) stop("Cohort contains an invalid Case pair", call. = FALSE)
  run <- provenance_run_context()
  written <- skipped <- failed <- character()

  process_case <- function(pair) {
    tryCatch({
      pair <- provenance_safe_pair(pair)
      case_dir <- file.path(output_data_dir, pair)
      if (!dir.exists(case_dir)) return(list(status = "skipped", pair = pair))
      manifest <- file.path(case_dir, "provenance.json")
      if (file.exists(manifest) && !force) {
        message("Skilift provenance backfill skipped Case ", pair,
                ": provenance.json already exists")
        return(list(status = "skipped", pair = pair))
      }
      previous <- provenance_read_document(manifest, pair)
      updates <- list()
      for (key in names(provenance_backfill_registry)) {
        provenance_safe_artifact_key(key)
        output_path <- file.path(case_dir, key)
        if (!file.exists(output_path) || dir.exists(output_path)) next
        binding <- provenance_backfill_registry[[key]]
        source <- if (identical(binding$lifter, "lift_metadata")) {
          metadata_inputs
        } else if (identical(binding$lifter, "lift_bam")) {
          bam_inputs
        } else {
          inputs
        }
        if (is.null(source)) {
          stop("Could not apply saved-Cohort input selection for ", binding$lifter, call. = FALSE)
        }
        row_index <- match(pair, as.character(source[["pair"]]))
        if (is.na(row_index)) next
        row <- if (inherits(source, "data.table")) {
          source[row_index]
        } else {
          source[row_index, , drop = FALSE]
        }
        input_bindings <- provenance_backfill_inputs(row, cohort, binding)
        output_observation <- provenance_observation(output_path, output = TRUE)
        updates[[key]] <- provenance_entry(
          lifter = binding$lifter,
          run = run,
          inputs = input_bindings,
          output_observation = output_observation
        )
      }
      if (!length(updates)) return(list(status = "skipped", pair = pair))
      document <- provenance_create_document(pair, updates, previous)
      provenance_write_document_atomic(manifest, document)
      list(status = "written", pair = pair)
    }, error = function(e) {
      provenance_log_failure(pair, e)
      list(status = "failed", pair = pair)
    })
  }

  results <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(pairs, process_case, mc.cores = as.integer(cores),
                       mc.preschedule = TRUE)
  } else {
    lapply(pairs, process_case)
  }
  for (result in results) {
    if (identical(result$status, "written")) written <- c(written, result$pair)
    else if (identical(result$status, "failed")) failed <- c(failed, result$pair)
    else skipped <- c(skipped, result$pair)
  }
  invisible(list(written = written, skipped = skipped, failed = failed))
}
