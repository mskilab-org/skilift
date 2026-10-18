#' Resolve the conventional provenance-table artifact paths
#'
#' @param root Directory containing per-pair provenance tables.
#' @param pair Case identifier.
#' @param stage Either "collected" or "curated".
#' @return A list with table and manifest paths. The writer may choose RDS instead.
#' @export
provenance_table_paths <- function(root, pair, stage) {
  stage <- match.arg(stage, c("collected", "curated"))
  list(table = file.path(root, pair, paste0(stage, ".parquet")),
       manifest = file.path(root, pair, "manifest.json"))
}

# Store declarative fields and executable transforms without serializing closures.
provenance_table_field_parameters <- function(fields) {
  lapply(fields, function(field) {
    out <- unclass(field)
    if (is.function(out$transform)) out$transform <- paste(deparse(out$transform), collapse = "\n")
    out
  })
}

#' Write an immutable-stage provenance table and its production receipt
#'
#' Collected and curated outputs occupy separate files; rerunning a stage only
#' replaces that stage's manifest entry. Atomic JSON publication preserves the
#' other stage. List columns, or options(skilift_provenance_format = "rds"), use
#' RDS; otherwise parquet preserves quotes, newlines, and column types.
#'
#' @inheritParams provenance_table_paths
#' @param table A data.frame or data.table to cache.
#' @param inputs Named source paths, or records containing a path member.
#' @param parameters Parameters used to produce this stage.
#' @return The actual table path, invisibly.
#' @export
write_provenance_table <- function(table, root, pair, stage, inputs = list(),
                                   parameters = list()) {
  paths <- provenance_table_paths(root, pair, stage)
  dir.create(dirname(paths$table), recursive = TRUE, showWarnings = FALSE)
  paths$table <- file.path(normalizePath(dirname(paths$table)), basename(paths$table))
  paths$manifest <- file.path(dirname(paths$table), "manifest.json")
  manifest <- if (file.exists(paths$manifest)) {
    jsonlite::read_json(paths$manifest, simplifyVector = FALSE)
  } else {
    list(schemaVersion = 1L, caseId = pair, stages = list())
  }

  .provenance_context <- provenance_current()
  if (is.null(.provenance_context)) {
    .provenance_context <- provenance_context(c(list(pair = pair), as.list(inputs)))
  }
  .provenance_context$parameters <- utils::modifyList(
    .provenance_context$parameters, c(list(provenance_stage = stage), parameters))

  input_records <- list()
  for (role in names(inputs)) {
    path <- inputs[[role]]
    if (is.list(path)) path <- path$path
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !file.exists(path) || dir.exists(path)) next
    provenance_note(path)
    input_records[[role]] <- list(path = normalizePath(path),
      sha256 = digest::digest(file = path, algo = "sha256"))
  }
  # Reader wrappers also capture references and any pluggable annotation files.
  for (source in .provenance_context$sources) {
    path <- source$path
    recorded <- vapply(input_records, function(x) x$path, character(1))
    if (!file.exists(path) || dir.exists(path) || normalizePath(path) %in% recorded) next
    role <- source$role
    if (is.null(role) || !nzchar(role) || role %in% names(input_records)) role <- path
    input_records[[role]] <- list(path = normalizePath(path),
      sha256 = digest::digest(file = path, algo = "sha256"))
  }

  format <- if (any(vapply(table, is.list, logical(1))) ||
                identical(getOption("skilift_provenance_format"), "rds")) "rds" else "parquet"
  if (format == "rds") {
    paths$table <- sub("\\.parquet$", ".rds", paths$table)
    saveRDS(table, paths$table)
  } else {
    arrow::write_parquet(table, paths$table)
  }
  manifest$stages[[stage]] <- list(path = paths$table, format = format,
    rows = nrow(table), columns = unname(as.list(names(table))),
    sha256 = digest::digest(file = paths$table, algo = "sha256"),
    producedAt = base::format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    inputs = input_records, parameters = parameters)
  provenance_atomic_json(manifest, paths$manifest)
  provenance_receipt(paths$table)
  invisible(paths$table)
}

# Resolve without reading the table, so cohort lifters can skip absent artifacts.
provenance_table_artifact <- function(root, pair, stage, manifest = NULL) {
  paths <- provenance_table_paths(root, pair, stage)
  if (is.null(manifest) && file.exists(paths$manifest)) manifest <- paths$manifest
  if (is.character(manifest)) manifest <- jsonlite::read_json(manifest, simplifyVector = FALSE)
  if (!is.null(manifest)) {
    record <- manifest$stages[[stage]]
    if (is.null(record)) stop("No provenance table for pair ", pair, ", stage ", stage)
    if (!record$format %in% c("parquet", "rds")) {
      stop("Unsupported provenance table format for pair ", pair, ", stage ", stage, ": ", record$format)
    }
    return(record)
  }
  path <- paths$table
  format <- "parquet"
  rds <- sub("\\.parquet$", ".rds", path)
  if (!file.exists(path) && file.exists(rds)) {
    path <- rds
    format <- "rds"
  }
  list(path = path, format = format)
}

#' Read a provenance table, rejecting a cache whose recorded SHA-256 has changed
#'
#' @inheritParams provenance_table_paths
#' @param manifest Optional manifest path or parsed manifest object.
#' @return A data.table; the source receipt is attached to the active lineage.
#' @export
read_provenance_table <- function(root, pair, stage = "curated", manifest = NULL) {
  record <- provenance_table_artifact(root, pair, stage, manifest)
  path <- record$path
  if (!file.exists(path)) stop("Missing provenance table for pair ", pair, ", stage ", stage, ": ", path)
  if (!is.null(record$sha256)) {
    actual <- digest::digest(file = path, algo = "sha256")
    if (!identical(actual, record$sha256)) {
      stop(structure(list(message = paste0("Provenance table SHA-256 mismatch for pair ",
        pair, ", stage ", stage, ": recorded ", record$sha256, "; actual ", actual),
        call = NULL), class = c("provenance_table_hash_mismatch", "error", "condition")))
    }
  }
  provenance_note(path)
  table <- if (record$format == "rds") base::readRDS(path) else arrow::read_parquet(path)
  data.table::as.data.table(table)
}

#' Collect a case's annotated variants into a provenance table
#'
#' This is the additive table path; the legacy oncotable collector is unchanged.
#' @inheritParams oncotable
#' @param sources Extra allele-keyed annotation sources.
#' @param fields Declarative annotation fields, including post-collection transforms.
#' @return A data.table with the original collector ordering and column projection.
#' @export
create_provenance_table <- function(
    pair,
    somatic_variant_annotations = NULL,
    fusions = NULL,
    jabba_gg = NULL,
    karyograph = NULL,
    events = NULL,
    signature_counts = NULL,
    multiplicity = NULL,
    oncokb_snv = NULL,
    oncokb_cna = NULL,
    oncokb_fusions = NULL,
    oncokb_intragenic_deletions = NULL,
    oncokb_intragenic_duplications = NULL,
    gencode,
    cytoband,
    verbose = TRUE,
    amp.thresh = 4,
    filter = "PASS",
    del.thresh = 0.5,
    somatic_variant_annotations_fields = Skilift:::echtvar_vcf_fields,
    sources = list(),
    fields = list()
    ) {
  out <- data.table()

  ## if ("type" %in% names(mcols(gencode))) {
  ##     pge <- gencode %Q% (type == "gene" & gene_type == "protein_coding")
  ## } else {
  ##   pge <- gencode %Q% (gene_type == "protein_coding")
  ## }

  
  if ("type" %in% names(mcols(gencode))) {
    ## pge <- gencode %Q% (type == "gene" & gene_type == "protein_coding")
    pge = (
      gencode
      %Q% (type == "gene")
      %Q% (order(
        GenomeInfoDb::rankSeqlevels(as.character(seqnames)),
        -rank(width))
      ) %Q% (!duplicated(gene_name))
    )
  } else {
    pge <- gencode %Q% (gene_type == "protein_coding")
  }




  ## collect gene fusions
  # prefer fusions from oncokb
  if (!is.na(oncokb_fusions) && !is.null(oncokb_fusions) && file.exists(oncokb_fusions)) {
      ## FIXME: the logic to grab coordinates doesn't make sense for fusions. Should take into account the actual junction + gene name..
    out <- rbind(
      out,
      collect_oncokb_fusions(oncokb_fusions, pge, cytoband, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_fusions(fusions, pge, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  if ( (
      !any(is.na(oncokb_intragenic_deletions))
    && !is.null(oncokb_intragenic_deletions)
    && all(file.exists(oncokb_intragenic_deletions))
  ) || (
    !any(is.na(oncokb_intragenic_duplications))
    && !is.null(oncokb_intragenic_duplications)
    && all(file.exists(oncokb_intragenic_duplications))    
  ) ) {
      ## FIXME: the logic to grab coordinates doesn't make sense for fusions. Should take into account the actual junction + gene name..
    out <- rbind(
      out,
      collect_oncokb_intragenic(
        oncokb_intragenic_deletions = oncokb_intragenic_deletions,
        oncokb_intragenic_duplications = oncokb_intragenic_duplications,
        pge,
        cytoband,
        verbose
      ),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_fusions(fusions, pge, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect complex events
  out <- rbind(
    out,
    collect_complex_events(events, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect copy number
  # prefer oncokb cna
  if (!is.na(oncokb_cna) && !is.null(oncokb_cna) && file.exists(oncokb_cna)) {
    out <- rbind(
      out,
      collect_oncokb_cna(oncokb_cna, jabba_gg, pge, amp.thresh, del.thresh, karyograph, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_copy_number_jabba(jabba_gg, pge, amp.thresh, del.thresh, verbose, karyograph),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect gene mutations
  # prefer oncokb snv
  if (!is.na(oncokb_snv) && !is.null(oncokb_snv) && file.exists(oncokb_snv)) {
    out <- rbind(
      out,
      collect_oncokb_with_sources(oncokb_snv, multiplicity, sources = sources, verbose = verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_mutations(somatic_variant_annotations, jabba_gg, filter, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  is_vcf_present = (
    !is.null(somatic_variant_annotations) 
    && !any(is.na(somatic_variant_annotations)) 
    && all(file.exists(somatic_variant_annotations)) 
  )

  ## is_echtvar_vcf_present = is_vcf_present && all(grepl("echtvar", somatic_variant_annotations))

  ## if (is_echtvar_vcf_present) {
  ##   ## Added redundancy just in case the fields that get fed in don't contain the canonical nf-gos echtvar annotated fields
  ##   somatic_variant_annotations_fields = c(somatic_variant_annotations_fields, Skilift:::echtvar_vcf_fields)
  ##   somatic_variant_annotations_fields = unique(somatic_variant_annotations_fields)
  ## }

  are_small_variants_populated = (
    NROW(out$type) > 0
    && any(out$track %in% "variants")
  )

  are_somatic_variant_annotations_fields_in_vcf = FALSE
  if (is_vcf_present) {
    provenance_note(somatic_variant_annotations)
    vcfheader = VariantAnnotation::scanVcfHeader(somatic_variant_annotations)
    fields_info = rownames(VariantAnnotation::info(vcfheader))
    fields_geno = rownames(VariantAnnotation::geno(vcfheader))
    fields_to_pass = unique(
      c(
        intersect(fields_info, Skilift:::echtvar_vcf_fields), ## redundancy to ensure these are added by default
        intersect(fields_geno, Skilift:::echtvar_vcf_fields), ## redundancy to ensure these are added by default
        intersect(fields_info, somatic_variant_annotations_fields),
        intersect(fields_geno, somatic_variant_annotations_fields)
      )
    )
    for (field in setdiff(somatic_variant_annotations_fields, fields_to_pass)) {
      message(paste("field", field, "not found in vcf"))
    }
    somatic_variant_annotations_fields = fields_to_pass
    are_somatic_variant_annotations_fields_in_vcf = NROW(somatic_variant_annotations_fields) > 0
  }

  if (is_vcf_present && are_small_variants_populated && are_somatic_variant_annotations_fields_in_vcf) {
    ## merge echtvar
    out <- Skilift:::merge_annotations(
      onco = out, 
      annotated_vcf = somatic_variant_annotations, 
      fields = somatic_variant_annotations_fields
    )
  }

  for (field in fields) {
    if (!is.null(field$transform) && field$as %in% names(out)) {
      data.table::set(out, j = field$as, value = field$transform(out[[field$as]]))
    }
  }

  ## add gene locations
  gene_locations <- readRDS(system.file("extdata", "data", "gene_locations.rds", package = "Skilift"))

  # Merge gene locations and create location string
  if ("gene" %in% names(out)) {
    out[gene_locations,
      gene_location := paste0(i.seqnames, ":", i.start, "-", i.end),
      on = c(gene = "gene_name")
    ]
  }

  out$id <- pair

  # remove all rows for which data was not passed
  out <- out[!is.na(type)]

  # coerce all empty strings to NA without affecting columns of levels
  out <- out[, lapply(.SD, function(x) {
    if (is.factor(x)) x <- as.character(x) # Convert factors to characters
    ifelse(x == "", NA_character_, x) # Replace empty strings with NA
  })]

  if (verbose) message("done processing sample")
  return(out)
}

#' Collect and cache provenance tables for a cohort
#'
#' Gencode and cytoband are loaded once before workers fork. Successful workers
#' publish collected-stage files, manifests, and lineage receipts.
#' @inheritParams create_oncotable
#' @param root Directory containing per-pair provenance tables.
#' @param sources_for Optional function(row) returning extra annotation sources.
#' @param fields Declarative fields and post-collection transforms.
#' @return The cohort with provenance_table_root and provenance_table inputs.
#' @export
lift_provenance_tables <- function(
    cohort,
    root,
    amp_thresh_multiplier = 1.5,
    gencode = Skilift::get_default_gencode(),
    cytoband = system.file("extdata", "data", "cytoband.rds", package = "Skilift"),
    cores = 1,
    somatic_variant_annotations_fields = Skilift:::echtvar_vcf_fields,
    sources_for = NULL,
    fields = list()) {
  if (!inherits(cohort, "Cohort")) {
    stop("Input must be a Cohort object")
  }

  if (nrow(cohort$inputs) == 0) {
    stop("No samples found in the cohort")
  }

  if (Sys.which("bcftools") == "") {
    stop("bcftools is not available on the system PATH. Try `module load htslib` first or install it.")
  } else {
    message("bcftools is available.")
  }

  if (gencode == "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds") {
    message("using default gencode: /gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds")
  }

  .provenance_context <- provenance_context(list(), cohort)
  gencode <- process_gencode(gencode)
  cytoband <- process_cytoband(cytoband)
  reference_sources <- .provenance_context$sources

  if (amp_thresh_multiplier == 1.5) {
    message("using default amp_thres_multiplier: 1.5")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(root)) {
    dir.create(root, recursive = TRUE)
  }

  # Create a copy of the cohort to modify
  updated_cohort <- cohort

  # Add provenance_table column if it doesn't exist
  if (!"provenance_table" %in% names(updated_cohort$inputs)) {
    updated_cohort$inputs[, provenance_table := NA_character_]
  }

  updated_cohort$inputs[, provenance_table_root := root]

  jabba_column = Skilift::DEFAULT_JABBA(object = cohort)

  futile.logger::flog.threshold("ERROR")
  iterfun = function(i) {
    main = function() {
      row <- cohort$inputs[i]
      .provenance_context <- provenance_context(row, cohort)
      .provenance_context$sources <- reference_sources
      .provenance_context$selection <- list(graphColumn = jabba_column)
      .provenance_context$parameters <- list(amp.thresh = amp_thresh_multiplier,
        del.thresh = 0.5, filter = "PASS", somatic_variant_annotations_fields = somatic_variant_annotations_fields)

      # Create output directory for this pair
      pair_root <- file.path(root, row$pair)
      if (!dir.exists(pair_root)) {
        dir.create(pair_root, recursive = TRUE)
      }

      # Get ploidy from jabba output, default to 2 if missing
      ploidy <- 2 # Default ploidy
      jabba_path = row[[jabba_column]]
      if (!is.null(jabba_path) && all(file.exists(jabba_path))) {
        ploidy_ggraph <- tryCatch(
          {
            process_jabba(row[[jabba_column]])
          },
          error = function(e) {
            msg <- sprintf("Error reading JaBbA file for %s: %s. Using default ploidy of 2.", row$pair, e$message)
            warning(msg)
          }
        )

        if (!is.null(ploidy_ggraph)) {
          ploidy <- ifelse(
            !is.null(ploidy_ggraph$meta$ploidy),
            ploidy_ggraph$meta$ploidy,
            ploidy_ggraph$ploidy
          )
        }
      }

      # amp_thresh <- amp_thresh_multiplier * ploidy
      amp_thresh = amp_thresh_multiplier
      message(paste("Processing", row$pair, "using amp.thresh of", amp_thresh))


      # Run provenance_table for this pair
      provenance_table_result <- create_provenance_table(
        pair = row[["pair"]],
        somatic_variant_annotations = row[["somatic_variant_annotations"]],
        fusions = row[["fusions"]], 
        jabba_gg = row[[jabba_column]],
        karyograph = row[["karyograph"]],
        events = row[["events"]],
        signature_counts = row[["signature_counts"]],
        multiplicity = row[["multiplicity"]],
        oncokb_snv = row[["oncokb_snv"]],
        oncokb_cna = row[["oncokb_cna"]],
        oncokb_fusions = row[["oncokb_fusions"]],
        oncokb_intragenic_deletions = row[["oncokb_intragenic_deletions"]],
        oncokb_intragenic_duplications = row[["oncokb_intragenic_duplications"]],
        gencode = gencode,
        cytoband = cytoband,
        verbose = TRUE,
        amp.thresh = amp_thresh,
        filter = "PASS",
        del.thresh = 0.5,
        somatic_variant_annotations_fields = somatic_variant_annotations_fields,
        sources = if (is.null(sources_for)) list() else sources_for(row),
        fields = fields
      )

      if (!is.null(provenance_table_result)) {
        # Save successful results
        path <- write_provenance_table(
          provenance_table_result, root, row$pair, "collected",
          inputs = as.list(row)[intersect(names(row), c(
            "somatic_variant_annotations", "fusions", jabba_column, "karyograph",
            "events", "signature_counts", "multiplicity", "oncokb_snv", "oncokb_cna",
            "oncokb_fusions", "oncokb_intragenic_deletions", "oncokb_intragenic_duplications"))],
          parameters = list(amp.thresh = amp_thresh, del.thresh = 0.5,
            filter = "PASS", somatic_variant_annotations_fields = somatic_variant_annotations_fields,
            fields = provenance_table_field_parameters(fields))
        )
        return(list(index = i, path = path))
      }
    }
    tryCatchLog(
      main(),
      error = function(e) {
        print(sprintf("Unexpected error processing %s: %s", cohort$inputs[i]$pair, e$message))
        return(NULL)
      }
    )
  }
  results <- mclapply(seq_len(nrow(cohort$inputs)), iterfun, mc.cores = cores, mc.preschedule = FALSE)


  # Update provenance_table paths in the cohort
  results <- Filter(Negate(is.null), results)
  for (result in results) {
    updated_cohort$inputs[result$index, provenance_table := result$path][]
  }
  message(sprintf("\nProcessing complete - results written to %s", root))

  return(updated_cohort)
}
