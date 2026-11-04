#' @name lift_bam
#' @title lift_bam
#' @description
#' Lift bam file to output data directory.
#' Currently hardlinking.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param provenance Enable per-Case provenance recording (default: TRUE).
#' @param .provenance_context Internal shared per-run provenance collector.
#' @return None
#' @export
lift_bam <- function(cohort, output_data_dir, cores = 1, overwrite = FALSE, copy = FALSE,
                     provenance = TRUE, .provenance_context = NULL) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    transfer_method = "ln -n"
    # force_flag = if (overwrite) "f" else ""
    do_ln_and_add_f_flag = !all(copy) && all(overwrite)
    do_copy = all(copy)
    if (do_ln_and_add_f_flag) {
        transfer_method = paste(transfer_method, "-f")
    }
    if (all(copy)) {
        transfer_method = "cp"
    }

    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Check for required column
    if (!"pair" %in% names(cohort$inputs)) {
        stop("Missing required column 'pair' in cohort")
    }
    
    # Checks if cohort needs to be paired
    # Ensures that "tumor_bam" and "normal_bam" columns are present
    lift_inputs = pairify_cohort_inputs(cohort, tumor_normal_columns = c("bam"), sep = "_")

    if (!any(c("tumor_bam", "normal_bam") %in% names(lift_inputs))) {
        message("No tumor_bam or normal_bam columns present in cohort inputs")
        return(NULL)
    }
    if (!is.null(lift_inputs$tumor_bam)) {
      bambai_path = paste(lift_inputs$tumor_bam, ".bai", sep = "")
      bai_path = paste(tools::file_path_sans_ext(lift_inputs$tumor_bam), ".bai", sep = "")
      lift_inputs$tumor_bam_bai = dplyr::case_when(
                                           file.exists(bambai_path) ~ bambai_path,
                                           file.exists(bai_path) ~ bai_path,
                                           TRUE ~ NA_character_
                                         )
    }
    if (!is.null(lift_inputs$normal_bam)) {
      bambai_path = paste(lift_inputs$normal_bam, ".bai", sep = "")
      bai_path = paste(tools::file_path_sans_ext(lift_inputs$normal_bam), ".bai", sep = "")
      lift_inputs$normal_bam_bai = dplyr::case_when(
                                           file.exists(bambai_path) ~ bambai_path,
                                           file.exists(bai_path) ~ bai_path,
                                           TRUE ~ NA_character_
                                         )
    }

    provenance_call <- provenance_acquire_context(
        output_data_dir, enabled = provenance, context = .provenance_context
    )
    provenance_context <- provenance_call$context
    if (provenance_call$owns) on.exit(provenance_finish(provenance_context), add = TRUE)

    # Process each sample in parallel
    upload_bam = function(i) {
        row <- lift_inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        if (!dir.exists(pair_dir)) dir.create(pair_dir, recursive = TRUE)

        transfer_one <- function(source, destination, artifact_key, source_column,
                                 input_role) {
            if (is.null(source) || !length(source) || anyNA(source) || !nzchar(as.character(source)[1])) {
                return(invisible(NULL))
            }
            if (file.exists(destination) && !all(overwrite) && !all(copy)) {
                return(invisible(NULL))
            }
            staged <- provenance_stage_output(provenance_context, destination)
            command <- paste(transfer_method, shQuote(as.character(source)[1]), shQuote(staged))
            status <- system(command)
            if (status != 0L || !file.exists(staged) || dir.exists(staged)) {
                if (!is.null(provenance_context)) {
                    if (file.exists(staged)) unlink(staged)
                }
                return(invisible(NULL))
            }
            if (!is.null(provenance_context)) {
                inputs <- provenance_input_bindings(
                    row, stats::setNames(source_column, input_role)
                )
                provenance_queue_output(
                    provenance_context,
                    pair = as.character(row$pair),
                    artifact_key = artifact_key,
                    lifter = "lift_bam",
                    inputs = inputs,
                    staged_path = staged,
                    final_path = destination
                )
            }
            invisible(NULL)
        }

        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            if (!is.null(row$tumor_bam) && !any(is.na(row$tumor_bam))) {
                transfer_one(row$tumor_bam, file.path(pair_dir, "tumor.bam"),
                             "tumor.bam", "tumor_bam", "tumor_bam")
                transfer_one(row$tumor_bam_bai, file.path(pair_dir, "tumor.bam.bai"),
                             "tumor.bam.bai", "tumor_bam_bai", "tumor_bai")
            }
            if (!is.null(row$normal_bam) && !any(is.na(row$normal_bam))) {
                transfer_one(row$normal_bam, file.path(pair_dir, "normal.bam"),
                             "normal.bam", "normal_bam", "normal_bam")
                transfer_one(row$normal_bam_bai, file.path(pair_dir, "normal.bam.bai"),
                             "normal.bam.bai", "normal_bam_bai", "normal_bai")
            }
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }
    mclapply(seq_len(nrow(lift_inputs)), upload_bam, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
