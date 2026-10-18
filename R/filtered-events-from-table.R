#' Create filtered events from a collected or curated provenance table
#'
#' The legacy filtered-events projection is retained verbatim after reading the
#' hash-validated cache, including its heme selection and fusion coordinate logic.
#' @inheritParams create_filtered_events
#' @param provenance_table_root Directory containing per-pair provenance tables.
#' @param provenance_stage Stage to consume; defaults to curated.
#' @export
create_filtered_events_from_table <- function(
    pair,
    provenance_table_root,
    provenance_stage = "curated",
    jabba_gg,
    out_file,
    return_table = FALSE,
    cohort_type = "paired",
    additional_columns = NULL,
    fusion_coord_padding = 250
    ) {
  ot <- read_provenance_table(provenance_table_root, pair, provenance_stage)

  possible_drivers = empty_oncotable = structure(list(gene = character(0), gene_summary = character(0), 
    role = character(0), value = integer(0), vartype = character(0), 
    type = character(0), tier = integer(0), tier_description = character(0), 
    variant_summary = character(0), therapeutics = character(0), 
    resistances = logical(0), diagnoses = character(0), prognoses = character(0), 
    effect = character(0), effect_description = character(0), 
    fusion_genes = character(0), footprint = character(0), fusion_gene_coords = character(0), 
    track = character(0), source = character(0), variant.g = character(0), 
    variant.c = character(0), variant.p = character(0), annotation = character(0), 
    distance = logical(0), major_count = numeric(0), minor_count = numeric(0), 
    major_snv_copies = numeric(0), minor_snv_copies = numeric(0), 
    altered_copies = numeric(0), segment_cn = integer(0), ref = integer(0), 
    alt = integer(0), VAF = numeric(0), gene_location = character(0), 
    id = character(0)), row.names = c(NA, 0L), class = c("data.table", 
                                                         "data.frame"))

    oncotable_col_to_filtered_events_col <- c(
        "id" = "id",
        "gene" = "gene",
        "fusion_genes" = "fusion_genes",
        "fusion_gene_coords" = "fusion_gene_coords",
        "value" = "fusion_cn",
        "vartype" = "vartype",
        "type" = "type",
        "variant.g" = "Variant_g",
        "variant.p" = "Variant",
        "altered_copies" = "estimated_altered_copies",
        "segment_cn" = "segment_cn",
        "ref" = "ref",
        "alt" = "alt",
        "VAF" = "VAF",
        "gene_location" = "Genome_Location",
        "tier" = "Tier",
        "role" = "role",
        "gene_summary" = "gene_summary",
        "variant_summary" = "variant_summary",
        "effect" = "effect",
        "effect_description" = "effect_description",
        "therapeutics" = "therapeutics",
        "resistances" = "resistances",
        "diagnoses" = "diagnoses",
        "prognoses" = "prognoses",
        "is_multi_hit_per_gene" = "is_multi_hit_per_gene",
        "echtvar_fields" = "tertiary_field"
    )

    if (!is.null(additional_columns)) {
      is_list = inherits(additional_columns, "list")
      is_json_schema = FALSE
      if (is_list) {
        has_json_schema_fields = sapply(
          X = additional_columns, 
          FUN = function(x) {
            all(c("id", "title", "dataIndex") %in% names(x))
          }
        )
        is_json_schema = all(has_json_schema_fields)
      }
      if (is_json_schema) {
        IX = seq_len(NROW(additional_columns))
        lst = gGnome::transpose(lapply(
          X = IX, 
          FUN = function(ii) {
            x = additional_columns[ii]
            name_in_json = names(x)
            y = x[[1]]
            if (is.null(name_in_json)) {
              name_in_json = y[["dataIndex"]]
            }
            added_colnms = y[["dataIndex"]]
            viewType = y[["viewType"]]
            list(
              added_colnms = added_colnms,
              name_in_json = name_in_json,
              viewType = viewType
            )
          }
        ), ffun = c)
        added_colnms = lst$added_colnms
        name_in_json_vector = gsub("[[:space:]]+", "_", lst$name_in_json, perl = TRUE)
        viewType = lst$viewType
        # dt_schema = data.table::data.table(
        #   added_colnms,
        #   name_in_json = name_in_json_vector,
        #   viewType
        # )
      } else {
        added_colnms = names(additional_columns)
        if (is.null(added_colnms)) {
          added_colnms = additional_columns
        }
        added_colnms = ifelse(
          is.na(added_colnms) | added_colnms == "", 
          additional_columns, 
          added_colnms
        )
      }
      for (j in seq_len(NROW(additional_columns))) {
        if (is_json_schema) {
          name_in_json = name_in_json_vector[j]
        } else if (is_list) {
          name_in_json = additional_columns[[j]]
        } else {
          name_in_json = additional_columns[j]
        }
        col = added_colnms[j]
        oncotable_col_to_filtered_events_col[col] <- name_in_json
      }
    }


  res.final = Skilift:::change_names(possible_drivers, oncotable_col_to_filtered_events_col)

  if (NROW(ot) > 0) {
      # add a fusion_gene_coords column of NAs if no fusions
      if (!"fusion_gene_coords" %in% colnames(ot)) {
        ot[, fusion_genes := NA]
        ot[, fusion_gene_coords := NA]
      }
      # snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation, perl = TRUE)]
      snvs <- ot[ot$vartype == "SNV" & !is.na(ot$type)]
      if ("tier" %in% colnames(ot)) {
        snvs <- snvs[order(is.na(tier)), ]
      }
      if (NROW(snvs) > 0) {
        snvs[, is_unique_p := !is.na(variant.p) & !duplicated(cbind(gene, variant.p))]
        snvs[, is_unique_g := !duplicated(cbind(gene, variant.g))]
        # FIXME: variant.c should not be absent..
        ## remove_variant_c = FALSE
        ## if (is.null(snvs$variant.c)) {
        ##   snvs$variant.c = 1:NROW(snvs)
        ##   remove_variant_c = TRUE
        ## }
        snvs[, is_unique_c := !duplicated(cbind(gene, variant.c))]
        snvs <-  snvs[is_unique_p | (is_unique_g & is_unique_c)]
        snvs$is_unique_p = NULL
        snvs$is_unique_g = NULL
        snvs$is_unique_c = NULL
        ## if (remove_variant_c) snvs$variant.c = NA
        variant_concatenated = snvs[, paste(variant.p, "/", variant.c, sep = " ")]
        variant_concatenated = gsub(" / NA", "",
          gsub("NA / ", "",
            variant_concatenated,
            perl = TRUE
          ), perl = TRUE
        )
        snvs$variant.p = variant_concatenated
        ## snvs$variant.c = NULL
        rm("variant_concatenated")

        snvs = Skilift:::parse_echtvar_oncotable(snvs)
        
      }

      homdels <- ot[ot$type == "homdel",][, vartype := "HOMDEL"][, type := "SCNA"]
      amps <- ot[ot$type == "amp",][, vartype := "AMP"][, type := "SCNA"]
      fusions <- ot[ot$type == "fusion",] ## fusion vartype is either fusion or outframe_fusion
      complex_events = ot[ot$source == "complex",][, vartype := type][, type := "Complex SV"]
      fp = get0("footprint", as.environment(as.list(complex_events)), ifnotfound = rep_len(NA_character_, NROW(complex_events)))
      complex_events$gene_location = fp
    
      possible_drivers <- rbind(snvs, homdels, amps, fusions, complex_events, fill = TRUE)
      
  }


  filtered_events_columns <- names(possible_drivers)[names(possible_drivers) %in% names(oncotable_col_to_filtered_events_col)]

  res <- possible_drivers[, ..filtered_events_columns]
  intersected_columns <- intersect(filtered_events_columns, names(res))
  setnames(res, old = intersected_columns, new = oncotable_col_to_filtered_events_col[intersected_columns])

  dedup_field = c("gene", "vartype", "Variant")
  dedup_df = list()
  for (f in dedup_field) {
      dedup_df[[f]] = res[[f]]
  }
  dedup_df = as.data.frame(dedup_df)
  ix_duplicated = which(base::duplicated.data.frame(dedup_df))
  if (NROW(ix_duplicated) > 0) {
      res = res[-ix_duplicated]
  }
  ## res <- res %>% unique(., by = c("gene", "vartype", "Variant"))
  if (nrow(res) > 0) {
    res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
    res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
    res.mut <- res[vartype == "SNV"] ## FIXME: This is all small mutations, not just SNV
    if (nrow(res.mut) > 0) {
      ## FIXME: Doesn't do anything? maybe parse echtvar further here..

      NULL
    }
    res.fus = res[type == "fusion",] ## need to deal each class explicitly
    if (NROW(res.fus) > 0) {
      res.fus$fusion_gene_coords = gsub("(?<=,)[+-]|[+-](?=,)|[+-]$", "", res.fus$fusion_gene_coords, perl = TRUE)
      res.fus$gene = res.fus$fusion_genes
      is_intragenic = ! grepl("::", res.fus$fusion_genes)
      is_intragenic_del = is_intragenic & grepl("Del", res.fus$Variant)
      is_intragenic_ptd = is_intragenic & grepl("PTD", res.fus$Variant)
      is_noframe = res.fus$vartype == "fusion"
      is_inframe = res.fus$vartype == "inframe_fusion"
      is_outframe = res.fus$vartype == "outframe_fusion"
      ## If variant.p isn't found
      ## i.e. OncoKB isn't updated
      ## then only in frame or out of frame
      ## will be labeled
      variant.p.fus = base::get0(
        "Variant",
        as.environment(as.list(res.fus)),
        ifnotfound = rep_len("", NROW(res.fus))
      )
      fus_frame_label = data.table::fcase(
        is_inframe, "In-Frame",
        is_outframe, "Out-of-Frame",
        is_noframe, "",
        default = NA_character_
      )
      fus_intragenic_label = data.table::fcase(
        is_intragenic_del, "",
        is_intragenic_ptd, "",
        !is_intragenic, "Fusion",
        default = NA_character_
      )
      fus_frame_label = trimws(paste(fus_frame_label, fus_intragenic_label))
        
      # fus_frame_label
      if (any(is.na(fus_frame_label))) stop("A fusion was not labeled as inframe, outframe, or noframe")
      variant_label = paste(
          fus_frame_label,
          variant.p.fus
      )
      variant_label = trimws(variant_label)
      variant_label = gsub("[[:space:]]{2,}", "", variant_label)
      res.fus$Variant = variant_label
      ## Fudging around the coordinates a bit because of the frontend keying behavior
      ## gOS frontend expects that the genome location is unique..
      ## This may be violated by fusions that use the same breakpoint
      ## But may annotate a different event due to transcript ambiguity.
      coords = res.fus$fusion_gene_coords
      ix_present = ( function(x) { which(! ( is.na(x) | !nzchar(x) )) } )(coords)
      grl_coords = gUtils::parse.grl(coords[ix_present])
      if (!is.numeric(fusion_coord_padding) || length(fusion_coord_padding) != 1 ||
          is.na(fusion_coord_padding) || fusion_coord_padding < 0) {
        stop("fusion_coord_padding must be a single non-negative number")
      }
      grl_coords_red = GenomicRanges::reduce(grl_coords + fusion_coord_padding, ignore.strand = TRUE)
      padded_coords = gUtils::grl.string(grl_coords_red)
      ## padded_coords = GenomicRanges::reduce(gUtils::parse.grl(res.fus$fusion_gene_coords) + 250, ignore.strand = TRUE) %>% gUtils::grl.string()
      coord_tbl = data.table::data.table(fusion_gene_coords = padded_coords)[, .(listid = .I, .GRP, iix = seq_len(.N)), by = fusion_gene_coords][order(rank(listid))]
      
      ## coord_tbl = data.table::data.table(fusion_gene_coords = res.fus$fusion_gene_coords)[, .(listid = .I, .GRP, iix = seq_len(.N)), by = fusion_gene_coords][order(rank(listid))]
      coord_tbl[, is_duplicated := iix > 1]
      ## dlfus = Skilift:::dunlist(strsplit(res.fus$fusion_gene_coords, ","))
      dlfus = Skilift:::dunlist(unname(strsplit(padded_coords, ",")))
      dlfus = data.table::merge.data.table(dlfus, coord_tbl, by = "listid")
      matfus = stringr::str_split_fixed(dlfus$V1, pattern = ":|-", 3)
      dlfus$seqnames = matfus[,1]
      dlfus$start = as.integer(matfus[,2])
      dlfus$end = as.integer(matfus[,3])
      ## dlfus[, start_padded := start - 250L][]
      ## dlfus[, end_padded := end + 250L][]
      set.seed(10)
      ## dlfus[is_duplicated == TRUE, start_padded := round(jitter(start_padded, amount = 5)) - round(runif(1) * 5)][]
      ## dlfus[is_duplicated == TRUE, end_padded := round(jitter(end_padded, amount = 5)) + round(runif(1) * 5)]
      dlfus$start_padded = dlfus$start
      dlfus$end_padded = dlfus$end
      dlfus[is_duplicated == TRUE, start_padded := round(jitter(start_padded, amount = 5)) - round(runif(1) * 5)][]
      dlfus[is_duplicated == TRUE, end_padded := round(jitter(end_padded, amount = 5)) + round(runif(1) * 5)]      
      dlfus[, padded_fus_gene_coords := paste(
        seqnames,
        ":",
        start_padded, "-", end_padded, sep = ""
        ), by = listid][]
      var_gene_coords = character(NROW(res.fus))
      var_gene_coords[ix_present] = dlfus[, paste(padded_fus_gene_coords, collapse = ","), by = listid]$V1
      ## res.fus$fusion_gene_coords = dlfus[, paste(padded_fus_gene_coords, collapse = ","), by = listid]$V1
      res.fus$fusion_gene_coords = var_gene_coords
      dlfus[, padded_genome_location_coords := paste(
        seqnames,
        ":",
        start_padded, "-", seqnames, ":", end_padded, sep = ""
      ), by = listid][]
      ## res.fus$Genome_Location = dlfus[, paste(padded_genome_location_coords, collapse = "|"), by = listid]$V1
      res.fus$Genome_Location = res.fus$fusion_gene_coords

      res.fus$estimated_altered_copies = res.fus$fusion_cn

      ## Filtering Fusions based on COSMIC and Duncavage:
      ## If the exact match is found -
      ## remove any reverse reciprocal fusion
      fusions_curated = Skilift:::get_curated_fusions()
      fusions_recurrent = fusions_curated$fusions_recurrent

      blacklist = fusions_curated$blacklist
      fusions_recurrent = fusions_recurrent[!fusions_recurrent %in% blacklist]

      fus_lst = data.table::tstrsplit(fusions_recurrent, "::")
      forwards = paste(fus_lst[[1]], "::", fus_lst[[2]], sep = "")
      reverse = paste(fus_lst[[2]], "::", fus_lst[[1]], sep = "")

      fgenes = gsub("@.*", "", res.fus$fusion_genes)
      fgenes = gsub("\\([0-9]+\\)", "", fgenes, perl = TRUE)

      is_in_forward = fgenes %in% forwards
      is_in_reverse = fgenes %in% reverse
      ## Some are present in both directions,
      ## but are not well-curated. So leave them in.
      is_both_f_and_r = fgenes %in% forwards[forwards %in% reverse]
      is_valid = (
          is_in_forward & !is_in_reverse
      ) | is_both_f_and_r
      ## Other is anything that's completely outside of the list
      ## in either direction.
      is_other = ! fgenes %in% c(forwards, reverse)

      res.fus = base::subset(res.fus, subset = is_valid | is_other)
      
    }
    res.cn.dt = res.cn <- res[(
      type == "SCNA" ## FIXME: redundant logic for now
      | vartype %in% c("AMP", "HOMDEL")
    )]
    if (nrow(res.cn) > 0) {
      # jab <- readRDS(jabba_gg)
      # res.cn.gr <- GRanges(res.cn)
      # res.cn.dt
      # res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"))
      # res.cn.dt <- as.data.table(res.cn.gr)
      # FIXME: value is actually min_cn is renamed fusion_cn!!
      # Also, we already picked out the corresponding segment cn to use (min_cn)
      res.cn.dt$min_cn = res.cn.dt$fusion_cn 
      res.cn.dt[, estimated_altered_copies := abs(min_cn - 2)]
      res.cn.dt[, segment_cn := min_cn]
      res.cn.dt[, Variant := vartype]
      res.cn.dt$fusion_cn = NA_real_

      # To fix very small CNAs that can appear due to forcing in reciprocal
      # junctions with very small gaps.
      ## res.cn.dt = res.cn.dt[!res.cn.dt$node_width_jab < 1000,]
      # remove redundant columns since already added to Variant
      
      res.cn.dt[, c("min_cn", "cn", "cn.high", "cn.low", "width", "strand") := NULL]
    }
    res.complex = res[type == "Complex SV"]
    has_complex = NROW(res.complex) > 0
    ## initial pass to determine if anything is worth showing (non decoy chrom)
    if (has_complex) {
      loc_complex = res.complex[["Genome_Location"]]
      loc_complex = gsub(";", ",", loc_complex)
      grl_footprint = gUtils::parse.grl(loc_complex)
      snm = GenomeInfoDb::seqnames(grl_footprint)
      stdchrom = GenomeInfoDb::standardChromosomes(snm)
      is_snm_in_stdchrom = S4Vectors::`%in%`(snm, stdchrom)
      grl_footprint_subset = grl_footprint[is_snm_in_stdchrom]
      numelem = S4Vectors::elementNROWS(grl_footprint_subset)
      is_all_grl_in_stdchrom = numelem > 0
      ## has_complex = any(is_all_grl_in_stdchrom)
      grl_footprint = grl_footprint_subset[is_all_grl_in_stdchrom]
      res.complex = res.complex[is_all_grl_in_stdchrom]
      has_complex = NROW(res.complex) > 0
    }
    if (has_complex) {
      
      grlfp_1mb = GenomicRanges::reduce(grl_footprint + 0.5e6)
      grlfp_5mb = GenomicRanges::reduce(grl_footprint + 2.5e6)
      grlfp_10mb = GenomicRanges::reduce(grl_footprint + 5e6)
      
      get_arbitrary_loc = function(grl, return_string = TRUE) {
        grlfp = gUtils::grl.unlist(grl)
        grlfp = data.table::setDT(BiocGenerics::as.data.frame(grlfp))[]
        grlfp_subset = grlfp[, {
          byfun = function() {
            nr = NROW(.SD)
            has_three_or_more = nr >= 3
            if (!has_three_or_more) {
              return(.SD)
            }
            set.seed(10)
            gr_rank = rank(-width, ties.method = "random")
            out = .SD[order(gr_rank)][1:3]
            return(out)
          }
          byfun()
        }, keyby = grl.ix]
        grlfp_subset$strand = rep_len("*", NROW(grlfp_subset))
        grfp_subset = GenomicRanges::trim(dt2gr(grlfp_subset))
        out = split(grfp_subset, grlfp_subset$grl.ix)
        if (return_string) {
          loc_complex = gUtils::grl.string(out)
          return(loc_complex)
        }
        return(out)
      }
      
      loc_complex_1mb = get_arbitrary_loc(grlfp_1mb)
      loc_complex_5mb = get_arbitrary_loc(grlfp_5mb)
      loc_complex_10mb = get_arbitrary_loc(grlfp_10mb)

      enr_orig = S4Vectors::elementNROWS(grl_footprint)
      enr_1mb = S4Vectors::elementNROWS(grlfp_1mb)
      wid_1mb = sum(BiocGenerics::width(grlfp_1mb))
      enr_5mb = S4Vectors::elementNROWS(grlfp_5mb)
      wid_5mb = sum(BiocGenerics::width(grlfp_5mb))
      enr_10mb = S4Vectors::elementNROWS(grlfp_10mb)
      wid_10mb = sum(BiocGenerics::width(grlfp_10mb))
      use_1mb = (
        (enr_1mb == enr_5mb) & (enr_5mb == enr_10mb) 
      )
      use_5mb = (
        (enr_1mb > enr_5mb) & (enr_5mb == enr_10mb)
      )
      use_10mb = (
        (enr_5mb > enr_10mb)
      )
      loc_complex = data.table::fcase(
        use_1mb, loc_complex_1mb,
        use_5mb, loc_complex_5mb,
        use_10mb, loc_complex_10mb,
        default = NA_character_
      )
      if (anyNA(loc_complex)) {
        stop("Error determining complex event location")
      }
      # enrs = c(enr_1mb, enr_5mb, enr_10mb)
      # delta_numwindows = c(0, diff(enrs))
      # is_biggest_change = delta_numwindows == min(delta_numwindows, na.rm = TRUE)
      # choices = list(
      #   loc_complex_1mb,
      #   loc_complex_5mb,
      #   loc_complex_10mb
      # )
      
      # loc_complex = choices[[which(is_biggest_change & enrs > 1)[1]]]
      # pick_first_anyways = NROW(loc_complex) == 0
      # if (pick_first_anyways) {
      #   loc_complex = loc_complex_1mb
      # }
      

      res.complex[["Genome_Location"]] = loc_complex
      res.complex$Variant = res.complex$vartype
    }
    res.final <- rbind(res.mut, res.cn.dt, res.fus, res.complex, fill = TRUE)
    res.final[, sample := pair]
    if (identical(cohort_type, "heme")) {
      res.final <- select_heme_events(res.final)
    }
    res.final$type = tools::toTitleCase(res.final$type)
    ### FIXME: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    res.final <- res.final[res.final$seqnames != "Y"]
    ### FIXME ^^^: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    
    
    ## FIXME: Decide whether to either propagate return_table to lift_mvp in lift_heme somehow, or change this
    ## conditional logic to be based on cohort_type. Can also be left alone.
    ## Note that return_table for debugging purposes as well.

  }

  provenance_write_json(
    res.final, 
    out_file, 
    pretty = TRUE, 
    auto_unbox = TRUE
  )

  if (return_table) {
      return(res.final)
  }
}

#' Lift filtered events, highlights, and summaries from provenance tables
#'
#' @inheritParams lift_filtered_events
#' @param provenance_stage Stage to consume; defaults to curated.
#' @return The cohort with string_summary populated for successful cases.
#' @export
lift_filtered_events_from_tables <- function(
  cohort, 
  output_data_dir, 
  cores = 1, 
  return_table = TRUE,
  additional_columns = Skilift:::default_echtvar_json_fields,
  show_all_tags = FALSE,
  show_all_small_muts_tags = FALSE,
  show_all_cna_tags = FALSE,
  show_all_svs_tags = FALSE,
  verbose = FALSE,
  fusion_coord_padding = 250,
  provenance_stage = "curated"
) {
  if (!inherits(cohort, "Cohort")) {
      stop("Input must be a Cohort object")
  }
  
  if (!dir.exists(output_data_dir)) {
      dir.create(output_data_dir, recursive = TRUE)
  }
  
  # Validate required columns exist
  jabba_column = Skilift::DEFAULT_JABBA(object = cohort)
  required_cols <- c("pair", "provenance_table_root")
  missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
  if (length(missing_cols) > 0) {
    print(
      paste(
        "Missing required columns in cohort:",
        paste(missing_cols, collapse = ", "),
        collapse = ""
      )
    )
  }
    
  cohort_type = cohort$type
  # Create a copy of the cohort to modify
  # Add summary column if it doesn't exist
  if (!"string_summary" %in% names(cohort$inputs)) {
    cohort$inputs[, string_summary := ""]
  }

  # Process each sample in parallel
  futile.logger::flog.threshold("ERROR")
  iterate_fun = function(i) {
      main = function() {
      row <- cohort$inputs[i,]
      .provenance_context <- provenance_context(row, cohort)
      pair_dir <- file.path(output_data_dir, row$pair)
      
      if (!dir.exists(pair_dir)) {
          dir.create(pair_dir, recursive = TRUE)
      }
      
      out_file <- file.path(pair_dir, "filtered.events.json")
      highlights_out_file <- file.path(pair_dir, "highlights.json")

      out = NULL
      string_summary = ""
      artifact <- provenance_table_artifact(row$provenance_table_root, row$pair, provenance_stage)
      if (!file.exists(artifact$path)) {
        return(
          list(
            pair = row$pair,
	    strin_summary = ""
	  )
        )
      }
      out <- create_filtered_events_from_table(
        pair = row$pair,
        provenance_table_root = row$provenance_table_root,
        provenance_stage = provenance_stage,
        jabba_gg = row[[jabba_column]],
        out_file = out_file,
        return_table = return_table,
        cohort_type = cohort_type,
        additional_columns = additional_columns,
        fusion_coord_padding = fusion_coord_padding
      )
      if (identical(cohort_type, "heme")) {
        create_heme_highlights(
          events_tbl = out,
          jabba_gg = row[[jabba_column]],
          out_file = highlights_out_file,
          tumor_type = row$tumor_type,
          cohorttuple = row
        )
      }
      string_summary = create_summary(
        events_tbl = out,
        cohort_type = cohort_type,
        cohorttuple = row,
        show_all = show_all_tags,
        show_all_small_muts_tags = show_all_small_muts_tags,
        show_all_cna_tags = show_all_cna_tags,
        show_all_svs_tags = show_all_svs_tags
      )
      return(
        list(
          pair = row$pair,
          string_summary = string_summary
        )
      )
    }
    tryCatchLog({
      if (verbose) {
        out = main()
      } else {
        Skilift::shutup(suppressMessages({out = main()}))
      }
      return(out)
    }, error = function(e) {
      if (inherits(e, "provenance_table_hash_mismatch")) stop(e)
      print(sprintf("Error processing %s: %s", cohort$inputs[i]$pair, e$message))
      NULL
    })
  }
  results = mclapply(seq_len(nrow(cohort$inputs)), iterate_fun, mc.cores = cores, mc.preschedule = TRUE)
  failures <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failures)) stop(paste(as.character(results[failures]), collapse = "\n"))
  results = results[!vapply(results, is.null, FALSE)]
  if (NROW(results) > 0) {
    # Transpose the list
    results = do.call(Map, c(f = c, results))
    setkey(cohort, pair)
    cohort$inputs[results$pair, string_summary := results$string_summary][]
  }
  return(cohort)
}
