#' @name initialize_metadata_columns
#' @title Initialize Metadata Columns
#' @description Creates a data.table with all expected columns initialized to NA
#' @param pair The sample pair identifier
#' @return A data.table with all possible columns initialized to NA
initialize_metadata_columns <- function(pair) {
    if (is.null(pair)) {
        stop("pair must be a character string")
    }
    if (!is.character(pair)) {
        stop("pair must be a character string")
    }
    if (nchar(pair) == 0) {
        stop("pair cannot be empty")
    }
    
    dt <- data.table(
        # Basic metadata
        pair = pair,
        tumor_type = NA_character_,
        disease = NA_character_,
        primary_site = NA_character_,
        inferred_sex = NA_character_,
        
        # Coverage QC
        coverage_qc = list(NULL),
        
        # Variant counts
        snv_count = NA_integer_,
        snv_count_normal_vaf_greater0 = NA_integer_,
        
        # Coverage parameters
        cov_slope = NA_real_,
        cov_intercept = NA_real_,
        
        # Signature lists
        deconstructsigs_sbs_fraction = list(NULL),
        sigprofiler_indel_fraction = list(NULL),
        sigprofiler_sbs_count = list(NULL),
        sigprofiler_sbs_fraction = list(NULL),
        
        # HRD scores
        hrd_score = NA_real_,
        hrd = list(NULL),
        b1_2 = NA_real_,
        b1 = NA_real_,
        b2 = NA_real_,
        tmb = NA_real_  # Add this line
    )
    return(dt)
}

#' @name add_basic_metadata
#' @title Add Basic Metadata
#' @description
#' Adds basic metadata information such as tumor type, disease, and primary site.
#'
#' @param metadata A data.table containing metadata.
#' @param tumor_type The type of tumor.
#' @param disease The disease associated with the sample.
#' @param primary_site The primary site of the tumor.
#' @return Updated metadata with basic information added.
add_basic_metadata <- function(
    metadata,
    input_tumor_type,
    input_disease,
    input_primary_site
) {
    # Input validation
    if (is.null(metadata)) {
        stop("metadata cannot be NULL")
    }
    if (!is.data.table(metadata)) {
        stop("metadata must be a data.table")
    }
    
    # Validate tumor_type if provided
    if (!is.null(input_tumor_type)) {
        if (!is.character(input_tumor_type)) {
            stop("tumor_type must be NULL or a character string")
        }
        metadata[, tumor_type := input_tumor_type]
    }
    
    # Validate disease if provided
    if (!is.null(input_disease)) {
        if (!is.character(input_disease)) {
            stop("disease must be NULL or a character string")
        }
        metadata[, disease := input_disease]
    }
    
    # Validate primary_site if provided
    if (!is.null(input_primary_site)) {
        if (!is.character(input_primary_site)) {
            stop("primary_site must be NULL or a character string")
        }
        metadata[, primary_site := input_primary_site]
    }
    
    return(metadata)
}

#' @name add_sex_information
#' @title Add Sex Information
#' @description
#' Adds inferred sex information to the metadata based on various inputs.
#'
#' @param metadata A data.table containing metadata.
#' @param input_inferred_sex The inferred sex of the sample.
#' @param jabba_gg A file path to the jabba_gg data.
#' @param tumor_coverage Coverage data for the tumor.
#' @return Updated metadata with inferred sex information added.
add_sex_information <- function(
    metadata,
    input_inferred_sex = NULL,
    jabba_gg = NULL,
    tumor_coverage = NULL 
) {
    # Direct sex specification
    if (!is.null(input_inferred_sex)) {
        if (!input_inferred_sex %in% c("male", "female")) {
            warning("inferred sex must be one of `male` or `female`")
        }
        metadata[, inferred_sex := input_inferred_sex]
        return(metadata)
    }
    
    # Infer from jabba_gg if available
    if (!is.null(jabba_gg)) {
        gg <- readRDS(jabba_gg)
        ncn.x <- gg$nodes$dt[
            (seqnames == "X" | seqnames == "chrX"),
            weighted.mean(cn, w = end - start + 1, na.rm = TRUE)
        ]
        metadata[, inferred_sex := ifelse(ncn.x < 1.4, "male", "female")]
        return(metadata)
    }
    
    # Infer from tumor coverage if available
    if (!is.null(tumor_coverage)) {
        t_cov <- gGnome:::readCov(tumor_coverage)
        ncn.x <- unique(gr2dt(t_cov)[, foreground.chr := mean(foreground), by = seqnames][, .(seqnames, foreground.chr)])[seqnames %in% c("X", "chrX")]$foreground.chr
        if (length(ncn.x) == 0) {
            warning("Could not extract X chromosome coverage from tumor coverage data")
            return(metadata)
        }
        metadata[, inferred_sex := ifelse(ncn.x < 0.7, "male", "female")]
        return(metadata)
    }
    
    return(metadata)
}

#' @name dlrs
#' @title derivative log ratio spread
#' @description
#'
#' function to get the dlrs for coverage data
#' used in meta_data_json
#'
#' @param x foreground from dryclean
#' @return dlrs
#' @export
#' @author Joel Rosiene
dlrs <- function(x) {
    nx <- length(x)
    if (nx < 3) {
        stop("Vector length>2 needed for computation")
    }
    tmp <- embed(x, 2)
    diffs <- tmp[, 2] - tmp[, 1]
    dlrs <- IQR(diffs, na.rm = TRUE) / (sqrt(2) * 1.34)
    return(dlrs)
}

#' @name extract_metrics
#' @title Extract Metrics from QC File
#' @description Helper function to extract specified metrics from a QC file
#' @param qc_data Data.table containing QC metrics
#' @param metrics Named vector of column names to extract
#' @param pair Sample pair identifier
#' @return Data.table with extracted metrics
extract_metrics <- function(qc_data, metrics, pair) {
    # Verify all metrics exist in data
    missing_cols <- setdiff(metrics, names(qc_data))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing columns in QC data: %s", paste(missing_cols, collapse = ", ")))
    }
    
    # Select and rename columns
    result <- qc_data[, metrics, with = FALSE]
    setnames(result, old = metrics, new = names(metrics))
    
    # Add pair identifier
    result[, pair := pair]
    
    return(result)
}

#' @name process_qc_metrics
#' @title Process Quality Control Metrics
#' @description Processes quality control metrics from multiple QC files and combines them
#' @param estimate_library_complexity Path to the estimate_library_complexity_metrics file
#' @param alignment_summary_metrics Path to the alignment_summary_metrics file
#' @param insert_size_metrics Path to the insert_size_metrics file
#' @param wgs_metrics Path to the wgs_metrics file
#' @param pair Sample pair identifier
#' @return A list containing processed QC metrics
process_qc_metrics <- function(
    estimate_library_complexity,
    alignment_summary_metrics,
    insert_size_metrics,
    wgs_metrics,
    pair
) {
    # Define metric mappings for each file type
    complexity_metrics_cols <- c(
        read_pairs_examined = "READ_PAIRS_EXAMINED",
        read_pair_duplicates = "READ_PAIR_DUPLICATES",
        read_pair_optical_duplicates = "READ_PAIR_OPTICAL_DUPLICATES",
        percent_duplication = "PERCENT_DUPLICATION"
    )
    
    alignment_metrics_cols <- c(
        total_reads = "TOTAL_READS",
        pf_reads_aligned = "PF_READS_ALIGNED",
        pf_aligned_bases = "PF_ALIGNED_BASES",
        mean_read_length = "MEAN_READ_LENGTH"
    )
    
    insert_metrics_cols <- c(
        median_insert_size = "MEDIAN_INSERT_SIZE"
    )
    
    wgs_metrics_cols <- c(
        median_coverage = "MEAN_COVERAGE",
        pct_30x = "PCT_30X",
        pct_50x = "PCT_50X"
    )
    
    # Read and extract metrics from each file
    complexity_data <- extract_metrics(
        fread(estimate_library_complexity),
        complexity_metrics_cols,
        pair
    )
    
    alignment_data <- extract_metrics(
        fread(alignment_summary_metrics)[CATEGORY=="PAIR"],
        alignment_metrics_cols,
        pair
    )
    
    insert_data <- extract_metrics(
        fread(insert_size_metrics),
        insert_metrics_cols,
        pair
    )
    
    wgs_data <- extract_metrics(
        fread(wgs_metrics),
        wgs_metrics_cols,
        pair
    )
    
    # Merge all metrics on pair
    qc_metrics <- Reduce(function(x, y) merge(x, y, by = "pair"), 
                        list(complexity_data, alignment_data, insert_data, wgs_data))
    
    # Calculate additional metrics
    qc_metrics[, `:=`(
        m_reads = total_reads / 1e6,
        m_reads_mapped = pf_reads_aligned / 1e6,
        percent_optical_duplication = read_pair_optical_duplicates / read_pairs_examined,
        percent_aligned = pf_aligned_bases / (total_reads * mean_read_length),
        percent_optical_dups_of_dups = read_pair_optical_duplicates / read_pair_duplicates
    )]
    
    return(as.list(qc_metrics))
}

#' @name add_coverage_metrics
#' @title Add Coverage Metrics
#' @description
#' Adds coverage metrics to the metadata based on tumor coverage and QC metrics.
#'
#' @param metadata A data.table containing metadata.
#' @param tumor_coverage Coverage data for the tumor.
#' @param foreground_col_name Name of the column in the coverage data to use for foreground coverage
#' @param estimate_library_complexity Path to library complexity metrics file
#' @param alignment_summary_metrics Path to alignment summary metrics file
#' @param insert_size_metrics Path to insert size metrics file
#' @param wgs_metrics Path to WGS metrics file
#' @return Updated metadata with coverage metrics added.
add_coverage_metrics <- function(
    metadata,
    tumor_coverage = NULL,
    foreground_col_name = "foreground",
    estimate_library_complexity = NULL,
    alignment_summary_metrics = NULL,
    insert_size_metrics = NULL,
    wgs_metrics = NULL
) {
    coverage_variance <- NULL
    if (!is.null(tumor_coverage)) {
        foreground = as.data.table(readRDS(tumor_coverage))[[foreground_col_name]]
        cov_var <- dlrs(foreground) /100
        coverage_variance <- list(coverage_variance = cov_var)
        metadata$coverage_qc <- list(list(coverage_variance))
    }
    
    processed_metrics <- NULL
    if (!is.null(estimate_library_complexity) && 
        !is.null(alignment_summary_metrics) && 
        !is.null(insert_size_metrics) && 
        !is.null(wgs_metrics)) {
        
        processed_metrics <- process_qc_metrics(
            estimate_library_complexity,
            alignment_summary_metrics,
            insert_size_metrics,
            wgs_metrics,
            metadata$pair
        )
        
        metadata$coverage_qc <- list(as.list(c(processed_metrics, coverage_variance)))
    }

    return(metadata)
}

#' @name vcf_count
#' @title vcf_count
#' @description
#' takes in a vcf and returns a data.table with total count and count of variants with normal vaf greater than 0
#'
#' @param vcf_path Path to the VCF file
#' @param genome Reference genome name (e.g., "hg19", "hg38")
#' @return data.table
#' @export
#' @author Shihab Dider, Tanubrata Dey
vcf_count <- function(
    vcf_path,
    genome) {
    if (!file.exists(vcf_path)) {
        stop("VCF file does not exist.")
    }
    vcf <- readVcf(vcf_path, genome)

    # Filter for PASS variants
    pass_variants <- rowRanges(vcf)$FILTER == "PASS"
    vcf <- vcf[pass_variants, ]

    snv_count <- length(vcf)

    return(
        data.table(
            category = c("snv_count", "snv_count_normal_vaf_greater0"),
            counts = c(snv_count, NA)
        )
    )
}
    


#' @name add_variant_counts
#' @title Add Variant Counts
#' @description
#' Adds variant counts to the metadata based on somatic and germline SNVs.
#'
#' @param metadata A data.table containing metadata.
#' @param somatic_snvs Path to somatic SNV VCF file
#' @param genome The genome reference used.
#' @return Updated metadata with variant counts added.
add_variant_counts <- function(
    metadata,
    somatic_snvs = NULL,
    genome = "hg19"
) {
    if (!is.null(somatic_snvs)) {
        # snv_counts_dt <- sage_count(somatic_snvs, genome = genome)
        snv_counts_dt <- vcf_count(somatic_snvs, genome = genome)
        
        metadata$snv_count <- snv_counts_dt[category == "snv_count", ]$counts
        ## count up filtered variants in VCF 

        metadata$snv_count_normal_vaf_greater0 <- snv_counts_dt[category == "snv_count_normal_vaf_greater0", ]$counts
    }
    return(metadata)
}

#' @name add_sv_counts
#' @title Add Structural Variant Counts
#' @description Adds junction, loose end, and total SV counts from JaBbA graph
#' @param meta_dt A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @return Updated metadata with SV counts
add_sv_counts <- function(metadata, jabba_gg = NULL) {
    if (is.null(jabba_gg)) {
        return(metadata)
    }
    
    gg <- readRDS(jabba_gg)
    
    # Count junctions
    if (nrow(gg$junctions$dt) > 0) {
        metadata$junction_count <- nrow(gg$junctions$dt[type != "REF", ])
    } else {
        warning("No junctions found in JaBbA graph, skipping SV counts...")
        metadata$junction_count <- 0
    }
    
    # Count loose ends
    if (length(gg$loose) > 0) {
        metadata$loose_count <- nrow(as.data.table(gg$loose)[terminal == FALSE, ])
    } else {
        warning("No loose ends found in JaBbA graph, skipping loose counts...")
        metadata$loose_count <- 0
    }
    
    # Calculate total SV count
    metadata[, sv_count := (junction_count + (loose_count / 2))]
    
    return(metadata)
}

#' @name add_purity_ploidy
#' @title Add Purity and Ploidy Metrics
#' @description Adds purity and ploidy information from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @return Updated metadata with purity and ploidy information
add_purity_ploidy <- function(metadata, jabba_gg = NULL) {
    if (is.null(jabba_gg)) {
        return(metadata)
    }
    
    gg <- readRDS(jabba_gg)
    metadata$purity <- gg$meta$purity
    metadata$ploidy <- gg$meta$ploidy
    metadata$beta <- gg$meta$ploidy
    metadata$gamma <- 1 - gg$meta$purity
    
    return(metadata)
}

#' @name add_loh
#' @title Add Loss of Heterozygosity Metrics
#' @description Adds LOH metrics from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param seqnames_loh Vector of sequence names to consider for LOH
#' @return Updated metadata with LOH metrics
add_loh <- function(
    metadata,
    jabba_gg = NULL,
    seqnames_loh = c(1:22)
) {
    if (is.null(jabba_gg) || is.null(seqnames_loh)) {
        return(metadata)
    }
    
    gg <- readRDS(jabba_gg)
    nodes.dt <- gg$nodes$dt
    nodes.dt[, seqnames := gsub("chr", "", seqnames)]
    nodes.dt <- nodes.dt[seqnames %in% seqnames_loh]
    totalseglen <- nodes.dt$width %>% sum()
    
    if ("cn.low" %in% names(nodes.dt)) {
        lohsegs <- nodes.dt[cn.low == 0, ] %>%
            .[cn.high > 0] %>%
            .$width %>%
            sum()
        loh_frc <- lohsegs / totalseglen
        metadata[, `:=`(
            loh_fraction = loh_frc,
            loh_seglen = lohsegs,
            loh_total_genome = totalseglen
        )]
    } else {
        metadata$loh_fraction <- "Not Allelic Jabba"
    }
    
    return(metadata)
}

#' @name add_genome_length
#' @title Add Total Genome Length
#' @description Adds total genome length from JaBbA graph
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param seqnames_genome_width_or_genome_length Vector of sequence names to include in genome length or genome length as a numeric
#' @return Updated metadata with genome length information
add_genome_length <- function(
    metadata,
    jabba_gg = NULL,
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
) {
    if (is.numeric(seqnames_genome_width_or_genome_length) && NROW(seqnames_genome_width_or_genome_length) == 1) {
        metadata$total_genome_length <- seqnames_genome_width_or_genome_length
        return(metadata)
    }

    if (is.null(jabba_gg) || is.null(seqnames_genome_width_or_genome_length)) {
        return(metadata)
    }
    
    gg <- readRDS(jabba_gg)
    nodes.gr <- gg$nodes$gr
    seqlengths.dt <- suppressWarnings(
        as.data.table(
            seqinfo(nodes.gr),
            keep.rownames = "seqnames"
        )
    )
    seqlengths.dt[, seqnames := gsub("chr", "", seqnames)]
    seqlengths.dt <- seqlengths.dt[seqnames %in% seqnames_genome_width_or_genome_length, ]
    metadata$total_genome_length <- sum(seqlengths.dt$seqlengths)
    
    return(metadata)
}

#' @name add_sv_types
#' @title Add Structural Variant Types
#' @description Adds counts of different SV types from JaBbA graph and complex events
#' @param metadata A data.table containing metadata
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param complex Path to complex events RDS file
#' @return Updated metadata with SV type counts
add_sv_types <- function(metadata, jabba_gg = NULL, complex = NULL) {
    sv_like_types_count <- NULL
    if (!is.null(jabba_gg)) {
        gg <- readRDS(jabba_gg)
        if (!is.null(gg$edges$dt) && nrow(gg$edges$dt) > 0) {
            sv_like_types_count <- table(gg$edges$dt[type == "ALT" & class != "REF"]$class)
            metadata$sv_types_count <- list(list(list(sv_like_types_count)))
        } else {
            warning("No edges found in JaBbA graph, skipping SV-like types count...")
        }
    }
    
    qrp_counts <- data.table(qrpmin = 0, qrpmix = 0, qrppos = 0)
    
    sv_types_count <- NULL
    if (!is.null(complex)) {
        complex_data <- readRDS(complex)
        if (!is.null(complex_data$meta$events$type)) {
            sv_types_count <- table(complex_data$meta$events$type)
            
            # Update QRP counts
            qrp_types <- c("qrpmin", "qrpmix", "qrppos")
            for (qrp_type in qrp_types) {
                if (qrp_type %in% names(sv_types_count)) {
                    qrp_counts[[qrp_type]] <- sv_types_count[[qrp_type]]
                }
            }
            
        }
    } else {
        warning("Complex events and JaBbA graph not found as inputs, skipping SV types count...")
    }
    
    metadata$sv_types_count <- list(as.list(c(sv_types_count, sv_like_types_count)))
    metadata$hrd <- list(as.list(qrp_counts))
    return(metadata)
}

#' @name add_coverage_parameters
#' @title Add Coverage Parameters
#' @description Adds coverage parameters to the metadata based on tumor coverage
#' @param metadata A data.table containing metadata
#' @param tumor_coverage Coverage data for the tumor
#' @return Updated metadata with coverage parameters added
add_coverage_parameters <- function(metadata, tumor_coverage) {
    if (!is.null(tumor_coverage)) {
        if (is.null(metadata$purity) || is.null(metadata$ploidy)) {
            warning("Purity and ploidy not found in metadata, cov_slope and cov_intercept will not be calculated")
        }
        rel2abs.cov <- skitools::rel2abs(readRDS(tumor_coverage),
            field = "foreground",
            purity = metadata$purity,
            ploidy = metadata$ploidy,
            return.params = TRUE
        )
        metadata$cov_slope <- rel2abs.cov[1] %>% unname()
        metadata$cov_intercept <- rel2abs.cov[2] %>% unname()
    }
    return(metadata)
}

#' @name add_het_pileups_parameters
#' @title Add Heterozygous Pileups Parameters
#' @description Adds slope and intercept parameters calculated from heterozygous pileups
#' @param metadata A data.table containing metadata
#' @param het_pileups_wgs Path to heterozygous pileups WGS data
#' @return Updated metadata with heterozygous pileups parameters
add_het_pileups_parameters <- function(metadata, het_pileups) {
    if (!is.null(het_pileups)) {
        hets.read <- grab.hets(het_pileups) %>% gr2dt()
        hets.read <- dt2gr(hets.read[, .(count = sum(count)), by = c("seqnames", "start", "end")])
        if (is.null(metadata$purity) || is.null(metadata$ploidy)) {
            warning("Purity and ploidy not found in metadata, het_slope and het_intercept will not be calculated")
        }
        rel2abs.hets <- skitools::rel2abs(hets.read,
            field = "count",
            purity = metadata$purity,
            ploidy = metadata$ploidy,
            return.params = TRUE
        )
        metadata$hets_slope <- rel2abs.hets[1] %>% unname()
        metadata$hets_intercept <- rel2abs.hets[2] %>% unname()
    }
    return(metadata)
}

#' @name add_tmb
#' @title Add Tumor Mutation Burden
#' @description
#' Calculates and adds tumor mutation burden (TMB) to the metadata based on SNV count and genome length
#'
#' @param metadata A data.table containing metadata
#' @param somatic_snvs Path to somatic SNV VCF file
#' @param jabba_gg Path to JaBbA graph RDS file
#' @param genome The genome reference used
#' @param seqnames_genome_width_or_genome_length Sequence names for genome width calculation
#' @return Updated metadata with TMB value added
add_tmb <- function(
    metadata,
    somatic_snvs = NULL,
    jabba_gg = NULL,
    genome = "hg19",
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
) {
    meta_dt <- copy(metadata)
    # First ensure we have snv_count
    if (is.null(meta_dt$snv_count) || is.na(meta_dt$snv_count)) {
        meta_dt <- add_variant_counts(meta_dt, somatic_snvs, genome)
    }
    
    # Then ensure we have total_genome_length
    if (is.null(meta_dt$total_genome_length) || is.na(meta_dt$total_genome_length)) {
        meta_dt <- add_genome_length(meta_dt, jabba_gg, seqnames_genome_width_or_genome_length)
    }
    
    # Calculate TMB if we have both required values
    is_tmb_computable <- !is.na(meta_dt$snv_count) && !is.na(meta_dt$total_genome_length) && !is.null(meta_dt$snv_count) && !is.null(meta_dt$total_genome_length)
    if (is_tmb_computable) {
        tmb_vals = meta_dt$snv_count / (as.numeric(meta_dt$total_genome_length) / 1e6)
        metadata[, tmb := tmb_vals]
        metadata[, tmb := round(tmb, digits = 3)]
    } else {
        warning("Cannot calculate TMB without both snv_count and total_genome_length")
    }
    
    return(metadata)
}

#' @name compute_signature_averages
#' @title Compute Signature Averages
#' @description
#' Computes signature averages from various signature files
#'
#' @param sig_file Path to signature file
#' @param sample Sample or pair name needed to subset to the correct sample
#' @param is_indel Boolean indicating if processing indel signatures
#' @param is_deconstruct_sigs Boolean indicating if using deconstructSigs format
#' @return List of signature averages
compute_signature_averages <- function(
    sig_file,
    is_indel = FALSE,
    is_deconstruct_sigs = FALSE
) {
    if (is_deconstruct_sigs) {
        sig_data <- readRDS(sig_file)
        weights <- as.data.table(sig_data$weights)
        weights <- t(weights)
        sigs.dt <- as.data.table(weights, keep.rownames = "Signature") %>% 
            setnames(., c("Signature", "weights"))
        sigs.dt[, Signature := gsub("Signature.", "SBS", Signature)]
        sigs.vect <- sigs.dt$weights
        names(sigs.vect) <- sigs.dt$Signature
        return(sigs.vect)
    } else {
        sig.dt <- fread(sig_file)
        if (nrow(sig.dt) == 0) {
            stop("No signatures found in file")
        }
        pair <- sig.dt$Samples[1]
        sig.dt[, pair := pair]
        sig.dt[, Samples := NULL]
        sig.dt_avg <- copy(sig.dt)
        sig.dt_avg[, pair := NULL]
        row_sum <- sig.dt_avg %>% rowSums()
        sigs.dt <- melt.data.table(sig.dt_avg, measure.vars = names(sig.dt_avg)) %>% 
            setnames(., c("signature", "value"))
        sigs.dt[, avg_value := value / row_sum]
        sigs.vect <- sigs.dt$avg_value
        names(sigs.vect) <- sigs.dt$signature
        sigs_counts.vect <- sigs.dt$value
        names(sigs_counts.vect) <- sigs.dt$signature
        
        if (is_indel) {
            return(list(
                indel_fraction = sigs.vect,
                indel_count = sigs_counts.vect
            ))
        } else {
            return(list(
                sbs_fraction = sigs.vect,
                sbs_count = sigs_counts.vect
            ))
        }
    }
}

#' @name add_signatures
#' @title Add Signatures
#' @description
#' Adds mutational signatures to the metadata based on various signature matrices.
#'
#' @param metadata A data.table containing metadata.
#' @param signatures_pair_name The name of the signature pair.
#' @param activities_sbs_signatures Activities of SBS signatures.
#' @param activities_indel_signatures Activities of indel signatures.
#' @param deconstructsigs_sbs_signatures DeconstructSigs SBS signatures.
#' @return Updated metadata with mutational signatures added.
add_signatures <- function(
    metadata,
    activities_sbs_signatures,
    activities_indel_signatures,
    deconstructsigs_sbs_signatures
) {
    # Initialize list columns with proper length
    n <- nrow(metadata)
    metadata[, `:=`(
        deconstructsigs_sbs_fraction = list(list()),
        deletionInsertion = list(list()),
        sigprofiler_indel_fraction = list(list()),
        sigprofiler_indel_count = list(list()),
        sigprofiler_sbs_fraction = list(list()),
        sigprofiler_sbs_count = list(list()),
        signatures = list(list())
    )]

    if (!is.null(activities_sbs_signatures)) {
        signatures <- compute_signature_averages(
            sig_file = activities_sbs_signatures,
            is_indel = FALSE,
            is_deconstruct_sigs = FALSE
        )

        metadata$sigprofiler_sbs_fraction <- list(as.list(signatures[["sbs_fraction"]]))
        metadata$sigprofiler_sbs_count <- list(as.list(signatures[["sbs_count"]]))
        metadata$signatures <- list(as.list(signatures[["sbs_fraction"]]))
    }

    if (!is.null(activities_indel_signatures)) {
        deletionInsertion <- compute_signature_averages(
            sig_file = activities_indel_signatures,
            is_indel = TRUE,
            is_deconstruct_sigs = FALSE
        )
        metadata$deletionInsertion <- list(as.list(deletionInsertion[["indel_fraction"]]))
        metadata$sigprofiler_indel_fraction <- list(as.list(deletionInsertion[["indel_fraction"]]))
        metadata$sigprofiler_indel_count <- list(as.list(deletionInsertion[["indel_count"]]))
    }

    if (!is.null(deconstructsigs_sbs_signatures)) {
        signatures <- compute_signature_averages(
            sig_file = deconstructsigs_sbs_signatures,
            is_indel = FALSE,
            is_deconstruct_sigs = TRUE
        )
        metadata$deconstructsigs_sbs_fraction <- list(as.list(signatures))
    }

    return(metadata)
}

#' @name add_hrd_scores
#' @title Add HRD Scores
#' @description
#' Adds homologous recombination deficiency (HRD) scores to the metadata.
#'
#' @param metadata A data.table containing metadata.
#' @param hrdetect HRDetect scores.
#' @param onenesstwoness Oneness and twoness scores.
#' @return Updated metadata with HRD scores added.
add_hrd_scores <- function(metadata, hrdetect, onenesstwoness) {
    if (!is.null(hrdetect)) {
        hrd <- readRDS(hrdetect)
        hrd_values <- data.table(
            dels_mh = hrd$indels_classification_table$del.mh.count,
            rs3 = hrd$exposures_rearr["RefSigR3", ],
            rs5 = hrd$exposures_rearr["RefSigR5", ],
            sbs3 = hrd$exposures_subs["SBS3", ],
            sbs8 = hrd$exposures_subs["SBS8", ],
            del_rep = hrd$indels_classification_table$del.rep.count
        )
        hrd_score <- hrd$hrdetect_output[1, "Probability"]
        metadata$hrd_score <- hrd_score

        if (is.null(metadata$hrd[[1]]) || is.na(metadata$hrd[[1]])) {
            metadata$hrd <- list(as.list(hrd_values))
        } else {
            metadata$hrd <- list(c(metadata$hrd[[1]], as.list(hrd_values)))
        }
    } else {
        warning("HRDetect scores not found, skipping HRD scores...")
    }

    if (!is.null(onenesstwoness)) {
        onetwo <- readRDS(onenesstwoness)
        metadata$b1_2 <- onetwo$ot_scores[, "SUM12"]
        metadata$b1 <- onetwo$ot_scores[, "BRCA1"]
        metadata$b2 <- onetwo$ot_scores[, "BRCA2"]
    } else {
        warning("Oneness and twoness scores not found, skipping...")
    }
    
    return(metadata)
}

#' @title Create Metadata for a Sample
#' @description
#' Creates a comprehensive metadata object for a single sample pair by aggregating various data inputs.
#'
#' @param pair The sample pair identifier.
#' @param tumor_type The type of tumor.
#' @param disease The disease associated with the sample.
#' @param primary_site The primary site of the tumor.
#' @param inferred_sex The inferred sex of the sample.
#' @param jabba_gg A file path to the jabba_gg data.
#' @param events Structural variant events.
#' @param somatic_snvs Somatic single nucleotide variants.
#' @param germline_snvs Germline single nucleotide variants.
#' @param tumor_coverage Coverage data for the tumor.
#' @param estimate_library_complexity Path to the estimate_library_complexity_metrics file.
#' @param alignment_summary_metrics Path to the alignment_summary_metrics file.
#' @param insert_size_metrics Path to the insert_size_metrics file.
#' @param wgs_metrics Path to the wgs_metrics file.
#' @param het_pileups Heterozygous pileups data.
#' @param signatures_pair_name The name of the signature pair.
#' @param matrix_indel_signatures Matrix of indel signatures.
#' @param matrix_sbs_signatures Matrix of SBS signatures.
#' @param activities_sbs_signatures Activities of SBS signatures.
#' @param hrdetect HRDetect scores.
#' @param onenesstwoness Oneness and twoness scores.
#' @param genome The genome reference used.
#' @param seqnames_loh Sequence names for loss of heterozygosity.
#' @param seqnames_genome_width_or_genome_length Sequence names and genome width in list or genome length as a numeric
#' @return A data.table containing the metadata for a single sample pair.
#' @export
create_metadata <- function(
    pair,
    tumor_type = NULL,
    disease = NULL, 
    primary_site = NULL,
    inferred_sex = NULL,
    jabba_gg = NULL,
    events = NULL,
    somatic_snvs = NULL,
    germline_snvs = NULL,
    tumor_coverage = NULL,
    foreground_col_name = "foreground",
    estimate_library_complexity = NULL,
    alignment_summary_metrics = NULL,
    insert_size_metrics = NULL,
    wgs_metrics = NULL,
    het_pileups = NULL,
    activities_indel_signatures = NULL,
    deconstructsigs_sbs_signatures = NULL,
    activities_sbs_signatures = NULL,
    hrdetect = NULL,
    onenesstwoness = NULL,
    genome = "hg19",
    seqnames_loh = c(1:22),
    seqnames_genome_width_or_genome_length = c(1:22, "X", "Y")
) {
    # Initialize metadata with all possible columns
    metadata <- initialize_metadata_columns(pair)
    # change NA to NULL
    fix_entries = c("tumor_type", "disease", "primary_site", "inferred_sex", "jabba_gg", "events", "somatic_snvs", "germline_snvs", "tumor_coverage", "estimate_library_complexity", "alignment_summary_metrics", "insert_size_metrics", "wgs_metrics", "het_pileups", "activities_indel_signatures", "deconstructsigs_sbs_signatures", "activities_sbs_signatures", "hrdetect", "onenesstwoness")
    for (x in fix_entries) {
        if (!exists(x) || is.null(get(x)) || is.na(get(x))) {
            assign(x, NULL)
        }
    }    
    # Add each component sequentially
    metadata <- add_basic_metadata(metadata, tumor_type, disease, primary_site)
    metadata <- add_sex_information(metadata, inferred_sex, jabba_gg, tumor_coverage)
    # Add coverage metrics
    metadata <- add_coverage_metrics(
        metadata,
        tumor_coverage,
        foreground_col_name,
        estimate_library_complexity,
        alignment_summary_metrics,
        insert_size_metrics,
        wgs_metrics
    )
    metadata <- add_variant_counts(metadata, somatic_snvs, genome)
    
    # New SV-related function calls
    metadata <- add_sv_counts(metadata, jabba_gg)
    metadata <- add_purity_ploidy(metadata, jabba_gg)
    metadata <- add_loh(metadata, jabba_gg, seqnames_loh)
    metadata <- add_genome_length(metadata, jabba_gg, seqnames_genome_width_or_genome_length)
    metadata <- add_sv_types(metadata, jabba_gg, events)
    metadata <- add_coverage_parameters(metadata, tumor_coverage)
    metadata <- add_het_pileups_parameters(metadata, het_pileups)
    
    # Add TMB calculation
    metadata <- add_tmb(metadata, somatic_snvs, jabba_gg, genome, seqnames_genome_width_or_genome_length)
    
    metadata <- add_signatures(
        metadata,
        activities_sbs_signatures,
        activities_indel_signatures,
        deconstructsigs_sbs_signatures
    )
    
    # Add HRD scores
    metadata <- add_hrd_scores(metadata, hrdetect, onenesstwoness)
    
    return(metadata)
}

#' @name lift_metadata
#' @title lift_metadata
#' @description
#' Create metadata files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param genome_length Genome length for the samples (for targeted panels or WES data)
#' @return None
#' @export
lift_metadata <- function(cohort, output_data_dir, cores = 1, genome_length = NULL) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Define all possible columns
    all_cols <- c(
        "pair", "tumor_type", "disease", "primary_site", "inferred_sex",
        "jabba_gg", "events", "somatic_snvs", "germline_snvs", "tumor_coverage",
        "estimate_library_complexity", "alignment_summary_metrics",
        "insert_size_metrics", "wgs_metrics", "het_pileups",
        "activities_sbs_signatures", "activities_indel_signatures",
        "hrdetect", "onenesstwoness"
    )
    
    # Check for required column
    if (!"pair" %in% names(cohort$inputs)) {
        stop("Missing required column 'pair' in cohort")
    }
    
    # Warn about missing optional columns
    missing_cols <- all_cols[!all_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        warning("Missing optional columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "metadata.json")
        
        tryCatch({
            # Create metadata object
            metadata <- create_metadata(
                pair = row$pair,
                tumor_type = row$tumor_type,
                disease = row$disease,
                primary_site = row$primary_site,
                inferred_sex = row$inferred_sex,
                jabba_gg = row$jabba_gg,
                events = row$events,
                somatic_snvs = row$somatic_snvs,
                germline_snvs = row$germline_snvs,
                tumor_coverage = row$tumor_coverage,
                estimate_library_complexity = row$estimate_library_complexity,
                alignment_summary_metrics = row$alignment_summary_metrics,
                insert_size_metrics = row$insert_size_metrics,
                wgs_metrics = row$wgs_metrics,
                het_pileups = row$het_pileups,
                activities_sbs_signatures = row$activities_sbs_signatures,
                activities_indel_signatures = row$activities_indel_signatures,
                hrdetect = row$hrdetect,
                onenesstwoness = row$onenesstwoness,
                seqnames_genome_width_or_genome_length = genome_length # if genome length is provided
            )

            if (is.null(metadata)) {
                warning(sprintf("No metadata generated for %s", row$pair))
                return()
            }
            
            # Write to JSON
            jsonlite::write_json(
                metadata,
                out_file,
                auto_unbox = TRUE,
                pretty = TRUE,
                null = "null"
            )
            
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s\nCallstack: %s", row$pair, e$message, deparse(e$call)))
        })
    }, mc.cores = cores, mc.preschedule = FALSE)
    
    invisible(NULL)
}
