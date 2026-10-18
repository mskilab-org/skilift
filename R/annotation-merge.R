#' Declare a source column and its table/JSON projection
#' @export
annotation_field <- function(source, as = source, json = NULL, title = NULL,
                             format = NULL, transform = NULL) {
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    stop("source must be a length-one character column name")
  }
  structure(list(source = source, as = as, json = json, title = title,
                 format = format, transform = transform), class = "annotation_field")
}

annotation_fields <- function(fields) {
  if (is.character(fields)) lapply(fields, annotation_field) else fields
}

#' Construct an allele key using the existing MAF anchor-base convention
#' @export
allele_key <- function(x, chromosome = "Chromosome", position = "vcf_pos",
                       ref = "vcf_ref", alt = "vcf_alt", sep = " ") {
  if (!all(c(chromosome, position) %in% names(x))) {
    stop("allele_key requires columns: ", paste(c(chromosome, position), collapse = ", "))
  }
  reference <- x[[ref]]
  alternate <- x[[alt]]
  if (is.null(reference) || is.null(alternate)) {
    maf_columns <- c("Variant_Type", "Allele", "Reference_Allele", "flanking_bps", "Tumor_Seq_Allele2")
    if (!all(maf_columns %in% names(x))) {
      stop("allele_key requires VCF REF/ALT columns or MAF columns: ", paste(maf_columns, collapse = ", "))
    }
    reconstructed <- maf_reconstruct_refalt(data.table::as.data.table(x))
    if (is.null(reference)) reference <- reconstructed$vcf_ref_reconstructed
    if (is.null(alternate)) alternate <- reconstructed$vcf_alt_reconstructed
  }
  paste(x[[chromosome]], x[[position]], reference, alternate, sep = sep)
}

read_annotation_source <- function(source) {
  if (is.character(source) && length(source) == 1L) {
    if (grepl("\\.rds$", source, ignore.case = TRUE)) readRDS(source) else fread(source)
  } else source
}

#' Merge projected annotations with deterministic duplicate-key semantics
#' @export
merge_annotations_by_allele <- function(target, source, fields,
                                        target_key = NULL, source_key = NULL,
                                        fallback = c("none", "oncokb_multiplicity"),
                                        overwrite = FALSE,
                                        direction = c("target_from_source", "source_into_matched_target")) {
  fallback <- match.arg(fallback)
  direction <- match.arg(direction)
  target <- data.table::copy(data.table::as.data.table(read_annotation_source(target)))
  source <- read_annotation_source(source)
  fields <- annotation_fields(fields)
  if (fallback == "oncokb_multiplicity") {
    if (direction != "target_from_source") stop("multiplicity fallback requires target_from_source direction")
    # Preserve all six legacy geometric matching passes; only projection is generic.
    matched <- merge_oncokb_multiplicity(data.table::copy(target), source,
      overwrite = TRUE, cols.keep = unique(vapply(fields, `[[`, character(1), "source")))
    target$multiplicity_id_match <- matched$multiplicity_id_match
    source <- matched
    rows <- seq_len(nrow(target))
    idx <- rows
  } else {
    source <- data.table::as.data.table(source)
    if (is.null(target_key)) target_key <- allele_key(target)
    if (is.null(source_key)) source_key <- allele_key(source)
    if (length(target_key) != nrow(target) || length(source_key) != nrow(source)) {
      stop("annotation key lengths must match their tables")
    }
    if (direction == "source_into_matched_target") {
      rows <- match(source_key, target_key)
      idx <- which(!is.na(rows))
      rows <- rows[idx]
      # Assignment into the first matching target: last source row wins.
      keep <- !duplicated(rows, fromLast = TRUE)
      rows <- rows[keep]
      idx <- idx[keep]
    } else {
      rows <- seq_len(nrow(target))
      idx <- match(target_key, source_key)
    }
  }
  for (f in fields) {
    if (!f$source %in% names(source)) {
      # The legacy collector derives absent fields (notably VAF) from counts.
      if (fallback == "oncokb_multiplicity") next
      message("Annotation source column missing: ", f$source)
      values <- rep(NA_character_, length(idx))
    } else {
      values <- source[[f$source]][idx]
    }
    if (!is.null(f$transform)) values <- f$transform(values)
    if (overwrite && direction == "target_from_source") {
      data.table::set(target, j = f$as, value = values)
    } else {
      column <- target[[f$as]]
      if (is.null(column)) column <- values[rep(NA_integer_, nrow(target))]
      use <- if (overwrite) seq_along(rows) else which(is.na(column[rows]))
      column[rows[use]] <- values[use]
      data.table::set(target, j = f$as, value = column)
    }
  }
  target
}

#' Build filtered-events JSON columns from the same projection specification
#' @export
annotation_json_columns <- function(fields, base = default_echtvar_json_fields) {
  mapped <- lapply(Filter(function(f) !is.null(f$json), annotation_fields(fields)), function(f) {
    col <- if (is.null(f$format)) list() else template_filtered_events_formatted_number
    col$id <- f$json
    col$title <- if (is.null(f$title)) f$json else f$title
    col$dataIndex <- f$as
    if (!is.null(f$format)) col$rendererProps$format <- f$format
    col
  })
  c(base, mapped)
}

MULTIPLICITY_FIELDS <- lapply(c(
  "annotation", "ref", "alt", "ref_denoised", "alt_denoised", "normal.ref", "normal.alt",
  "variant.g", "major.count", "minor.count", "major_snv_copies", "minor_snv_copies",
  "total_snv_copies", "total_copies", "VAF", "cn", "altered_copies"
), annotation_field)

#' Collect OncoKB annotations with optional allele-keyed source projections
#' @export
collect_oncokb_with_sources <- function(oncokb_maf, multiplicity = NA_character_, sources = list(), verbose = TRUE) {
  
  empty_output_oncokb = data.table(type = NA, source = "oncokb_maf")
  if (missing(oncokb_maf) || is.null(oncokb_maf) || !file.exists(oncokb_maf)) {
    if (verbose) message("OncoKB MAF file is missing or does not exist.")
    return(empty_output_oncokb)
  }

  # snpeff_ontology = readRDS(system.file("extdata", "data", "snpeff_ontology.rds", package = "Skilift"))
  oncokb <- fread(oncokb_maf, colClasses = c("Chromosome" = "character"))
  is_multiplicity_present = is.character(multiplicity) && NROW(multiplicity) == 1 && !is_loosely_na(multiplicity) && file.exists(multiplicity)
  is_oncokb_populated = NROW(oncokb) > 0

  is_multiplicity_populated = FALSE
  if (is_multiplicity_present) {
    multiplicity <- readRDS(multiplicity)
    is_multiplicity_populated = NROW(multiplicity) > 0
  }

  # if (is_oncokb_populated && !is_multiplicity_populated) {
  #   stop("Something's off - oncokb is populated with variants, but not multiplicity.")
  # }
  
  if (is_oncokb_populated && is_multiplicity_populated) {
    oncokb = oncokb[
      order(
        GenomeInfoDb::rankSeqlevels(Chromosome),
        vcf_pos
      )
    ]
    oncokb$Chromosome = ifelse(oncokb$Chromosome == "", "na", oncokb$Chromosome)
    lst = Skilift:::maf_reconstruct_refalt(oncokb)
    oncokb$vcf_ref = lst$vcf_ref_reconstructed
    oncokb$vcf_alt = lst$vcf_alt_reconstructed
    order_ix_multiplicity = with(data = multiplicity, expr = {
      order(
        GenomeInfoDb::rankSeqlevels(as.character(seqnames)),
        start
      )
    })
    multiplicity = multiplicity[order_ix_multiplicity]
    oncokb <- merge_annotations_by_allele(oncokb, multiplicity, fields = MULTIPLICITY_FIELDS,
      fallback = "oncokb_multiplicity", overwrite = TRUE)
    oncokb = Skilift:::annotate_multihit(oncokb)
  }

  all_source_as_names <- character()
  for (s in sources) {
    s$fields <- annotation_fields(s$fields)
    oncokb <- merge_annotations_by_allele(oncokb, s$table, s$fields,
      source_key = s$key, overwrite = isTRUE(s$overwrite))
    all_source_as_names <- union(all_source_as_names, vapply(s$fields, `[[`, character(1), "as"))
  }

  for (col in c("segment_cn", "segment_cn_low", "segment_cn_high")) {
    if (!col %in% names(oncokb)) {
      oncokb[[col]] <- rep_len(NA, NROW(oncokb))
    }
  }

  # if (is_oncokb_populated && ! is_multiplicity_populated) {
  #   # oncokb$major.count = NA_integer_
  #   # oncokb$minor.count = NA_integer_
  #   # oncokb$major_snv_copies = NA_real_
  #   # oncokb$minor_snv_copies = NA_real_
  #   # oncokb$altered_copies = NA_real_
  #   oncokb$segment_cn = NA_real_
  #   oncokb$is_multi_hit_per_gene = NA
  # }

  if (is_oncokb_populated) {
    lst = maf_reconstruct_refalt(oncokb)
    vcf_ref = lst$vcf_ref_reconstructed
    vcf_alt = lst$vcf_alt_reconstructed
    # mat_vcf_id = stringr::str_split_fixed(oncokb$vcf_id, "_|/", n = 3)
    # vcf_ref = mat_vcf_id[,2]
    # vcf_alt = mat_vcf_id[,3]
    oncokb$vcf_ref = vcf_ref
    oncokb$vcf_alt = vcf_alt
    # rm(mat_vcf_id)
    oncokb[["variant.g"]] = paste0(
      oncokb$vcf_ref, ">", oncokb$vcf_alt
    )

      ## oncokb$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
  
    oncokb$short <- dplyr::case_when(
      grepl("frameshift", oncokb$Consequence) & grepl("fs$", oncokb$HGVSp) ~ "trunc",
      grepl("stop", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "trunc",
      grepl("lost", oncokb$Consequence) ~ "trunc",
      grepl("missense", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "missense",
      grepl("splice", oncokb$Consequence) ~ "splice",
      grepl("inframe_(insertion|deletion)", oncokb$Consequence) ~ "inframe",
      grepl("upstream", oncokb$Consequence) ~ "upstream",
      grepl("(5|3)_prime_UTR_variant", oncokb$Consequence) ~ "UTR",
      TRUE ~ NA_character_
    )

    oncokb <- parse_oncokb_tier(
      oncokb,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )

    # coerce_column_tuples <- list(
    #   c("altered_copies", "total_copies")
    # )
    # current_nms <- names(oncokb)
    # for (col in coerce_column_tuples) {
    #   # if (col[1] %in% current_nms && !col[2] %in% current_nms) {
    #   current_nms[current_nms == col[1]] <- col[2]
    #   # }
    # }
    # names(oncokb) <- current_nms
    return(
      oncokb[, {
        ENV = environment()
        main = function() {
          if (!exists("major.count", where = ENV)) {
            major.count = NA_real_
          }
          if (!exists("minor.count", where = ENV)) {
            minor.count = NA_real_
          }
          if (!exists("major_snv_copies", where = ENV)) {
            major_snv_copies = NA_real_
          }
          if (!exists("minor_snv_copies", where = ENV)) {
            minor_snv_copies = NA_real_
          }
          if (!exists("altered_copies", where = ENV)) {
            altered_copies = NA_real_
          }
          if (!exists("segment_cn", where = ENV)) {
            segment_cn = NA_real_
          }
          if (
              !exists("ref", mode = "character", where = ENV)
              && !exists("ref", mode = "numeric", where = ENV)
          ) {
            ref = NA_real_
          }
          if (
              !exists("alt", mode = "character", where = ENV)
              && !exists("alt", mode = "numeric", where = ENV)
          ) {
            alt = NA_real_
          }
          if (
            !exists("VAF", mode = "character", where = ENV)
            && !exists("VAF", mode = "numeric", where = ENV)
          ) {
            VAF = alt / (ref + alt)
          }
          if (!exists("is_multi_hit_per_gene", mode = "character", where = ENV)) {
            is_multi_hit_per_gene = NA
          }
          out = list(
            gene = Hugo_Symbol,
            gene_summary = GENE_SUMMARY,
            role = Role,
            variant.g = paste(Chromosome, ":", Start_Position, "-", End_Position, " ", variant.g, sep = ""),
            variant.c = HGVSc,
            variant.p = HGVSp,
            annotation = Consequence,
            type = short,
            tier = tier,
            tier_description = tier_factor,
            variant_summary = VARIANT_SUMMARY,
            therapeutics = tx_string,
            resistances = rx_string,
            diagnoses = dx_string,
            prognoses = px_string,
            distance = NA_integer_,
            effect = MUTATION_EFFECT,
            effect_description = MUTATION_EFFECT_DESCRIPTION,
            major_count = major.count,
            minor_count = minor.count,
            major_snv_copies = major_snv_copies,
            minor_snv_copies = minor_snv_copies,
            # total_copies, ## total_copies gets converted to estimated_altered_copies downstream
            altered_copies = altered_copies,
            segment_cn = segment_cn,
            ref = ref,
            alt = alt,
            VAF = VAF,
            vartype = "SNV",
            track = "variants",
            source = "oncokb_maf",
            is_multi_hit_per_gene = is_multi_hit_per_gene,
            Chromosome = Chromosome,
            vcf_pos = vcf_pos,
            Tumor_Seq_Allele2 = Tumor_Seq_Allele2,
            vcf_ref = vcf_ref,
            vcf_alt = vcf_alt,
            ## transcript_id: prefer MAF Transcript_ID, fall back to Feature.
            ## Both are populated by the upstream OncoKB annotator from the
            ## snpEff-derived transcript (e.g. "NM_022552.4").
            transcript_id = {
              tx_id_col <- if (exists("Transcript_ID", where = ENV)) {
                as.character(Transcript_ID)
              } else {
                rep(NA_character_, length(Hugo_Symbol))
              }
              feature_col <- if (exists("Feature", where = ENV)) {
                as.character(Feature)
              } else {
                rep(NA_character_, length(Hugo_Symbol))
              }
              ifelse(!is.na(tx_id_col) & nzchar(tx_id_col), tx_id_col, feature_col)
            }
          )
          extra_names <- setdiff(intersect(all_source_as_names, names(oncokb)), names(out))
          out <- c(out, as.list(oncokb[, extra_names, with = FALSE]))
          return(out)
        }
        main()
      }]
    )
  }
  return(empty_output_oncokb)
}
