#' @name create_variant_qc
#' @title create_variant_qc
#' @description
#'
#' Function to create quality control metrics for somatic variants
#'
#' @param somatic_snvs Path to sage VCF file to be used (paired T-N/ Tumor only)
#' @param reference_name Reference genome name (e.g., "hg19", "hg38")
#' @return data.table containing variant QC metrics
#' @export
#' @author Shihab Dider, Tanubrata Dey, Stanley Clarke
create_variant_qc <- function(
    somatic_snvs,
    reference_name = "hg19") {
    vcf <- readVcf(somatic_snvs, reference_name)

    # Filter for PASS variants
    pass_variants <- rowRanges(vcf)$FILTER == "PASS"
    vcf <- vcf[pass_variants, ]

    # Extract necessary information from VCF object
    chrom <- as.character(seqnames(rowRanges(vcf)))
    pos <- start(rowRanges(vcf))
    ref <- as.character(VariantAnnotation::ref(vcf))
    alt <- as.character(unlist(VariantAnnotation::alt(vcf)))
    filter <- as.character(rowRanges(vcf)$FILTER)
    qual <- as.numeric(rowRanges(vcf)$QUAL)

    # Extract depth and allele count information from the genotype (geno) slot
    geno_data <- VariantAnnotation::geno(vcf)

    # check if tumor-only or paired analysis 
    is_tumor_only <- length(colnames(geno_data$DP)) == 1
    if (is_tumor_only) {
        tumor <- colnames(geno_data$DP)[1]
    } else {
        normal <- colnames(geno_data$DP)[1]
        tumor <- colnames(geno_data$DP)[2]
    }

    T_DP <- as.numeric(geno_data$DP[, tumor])
    alt_count_T <- sapply(geno_data$AD[, tumor], function(x) as.numeric(x[2])) # Extract the second element for alternate allele depth
    T_ABQ <- as.numeric(geno_data$ABQ[, tumor])
    VAF_T <- as.numeric(geno_data$AF[, tumor])

    # Check if normal sample data exists
    if (is_tumor_only) {
        sq <- data.table(
            chromosome = chrom,
            position = pos,
            reference = ref,
            alternate = alt,
            filter = filter,
            mapping_quality = qual,
            tumor_depth = T_DP,
            tumor_alt_counts = alt_count_T,
            tumor_abq = T_ABQ,
            tumor_vaf = VAF_T
        )

        consider_numeric <- c("tumor_depth", "tumor_alt_counts", "tumor_abq", "tumor_vaf")
    } else {
        N_DP <- as.numeric(geno_data$DP[, normal])
        alt_count_N <- sapply(geno_data$AD[, normal], function(x) as.numeric(x[2])) # Extract the second element for alternate allele depth
        VAF_N <- as.numeric(geno_data$AF[, normal])

        sq <- data.table(
            chromosome = chrom,
            position = pos,
            reference = ref,
            alternate = alt,
            filter = filter,
            mapping_quality = qual,
            tumor_depth = T_DP,
            normal_depth = N_DP,
            normal_alt_counts = alt_count_N,
            tumor_alt_counts = alt_count_T,
            tumor_abq = T_ABQ,
            tumor_vaf = VAF_T,
            normal_vaf = VAF_N
        )

        consider_numeric <- c("tumor_depth", "normal_depth", "normal_alt_counts", "tumor_alt_counts", "tumor_abq", "tumor_vaf", "normal_vaf")
    }

    sq[, (consider_numeric) := lapply(.SD, as.numeric), .SDcols = consider_numeric]
    sq[, id := .I]

    return(sq)
}
#' @name lift_variant_qc
#' @title lift_variant_qc
#' @description
#' Creates variant QC metrics for all samples in a cohort and writes them to JSON files
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_variant_qc <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    if (!"somatic_snvs" %in% names(cohort$inputs)) {
        stop("Missing required column 'somatic_snvs' in cohort")
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "sage.qc.json")
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            # Generate QC metrics
            qc_data <- create_variant_qc(
                somatic_snvs = row$somatic_snvs,
                reference_name = cohort$reference_name
            )
            
            # Write to JSON
            jsonlite::write_json(qc_data, out_file, pretty = TRUE)
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}

.infer_type <- function(x) {
    if (is.character(x) || is.factor(x) || is.logical(x)) {
        return("enum")
    }
    if (is.numeric(x)) {
        if (all(na.omit(x) == floor(na.omit(x)))) {
            return("int")
        } else {
            return("float")
        }
    }
    return(NA_character_)
}

#' @name create_variant_qc_schema
#' @title create_variant_qc_schema
#' @description
#'
#' Creates a schema for variant QC metrics from a data.table.
#' The schema maps field names to their inferred types (enum, int, float).
#'
#' @param qc_dt data.table containing variant QC metrics (e.g., from create_variant_qc)
#' @return A named list representing the schema.
#' @export
create_variant_qc_schema <- function(qc_dt) {
    # Columns to include in schema, based on user example
    schema_cols <- c("filter", "mapping_quality", "tumor_depth", "normal_depth",
                     "normal_alt_counts", "tumor_alt_counts", "tumor_abq",
                     "tumor_vaf", "normal_vaf")

    # Filter for columns that exist in the data.table
    schema_cols_exist <- intersect(schema_cols, names(qc_dt))

    if (length(schema_cols_exist) == 0) {
        return(list())
    }

    # Infer types for existing columns
    schema <- sapply(qc_dt[, ..schema_cols_exist], .infer_type)

    # Filter out any columns where type could not be inferred
    schema <- schema[!is.na(schema)]

    return(as.list(schema))
}

#' @name add_variant_qc_schema_to_dataset
#' @title add_variant_qc_schema_to_dataset
#' @description
#'
#' Generates a variant QC schema and adds it to a specified dataset in a JSON file.
#'
#' @param qc_object A data.table of variant QC metrics or a path to a variant QC JSON file.
#' @param datasets_json_path Path to the datasets.json file.
#' @param dataset_id The ID of the dataset to update.
#' @return None
#' @export
add_variant_qc_schema_to_dataset <- function(qc_object, datasets_json_path, dataset_id) {
    if (is.character(qc_object) && length(qc_object) == 1 && file.exists(qc_object)) {
        qc_dt <- as.data.table(jsonlite::fromJSON(qc_object))
    } else if (is.data.table(qc_object)) {
        qc_dt <- qc_object
    } else {
        stop("qc_object must be a data.table or a valid file path to the sage.qc JSON file.")
    }

    schema <- create_variant_qc_schema(qc_dt)

    if (!file.exists(datasets_json_path)) {
        stop("datasets.json file not found at: ", datasets_json_path)
    }

    datasets <- jsonlite::read_json(datasets_json_path, simplifyVector = FALSE)

    dataset_indices <- which(sapply(datasets, function(d) d$id == dataset_id))

    if (length(dataset_indices) == 0) {
        stop(paste("Dataset with id", dataset_id, "not found in", datasets_json_path))
    }

    for (idx in dataset_indices) {
        datasets[[idx]][["variant_qc_dropdown_schema"]] <- schema
    }

    jsonlite::write_json(datasets, datasets_json_path, pretty = TRUE, auto_unbox = TRUE)

    invisible(NULL)
}
