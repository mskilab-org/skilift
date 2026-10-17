#' @name lift_bam
#' @title lift_bam
#' @description
#' Lift BAM files and indexes to the gOS output data directory. Source symlinks
#' (including relative links and chains) are resolved to their original files.
#' Hardlinks are used by default; symbolic links and copies are optional.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param overwrite Replace existing destination files or symlinks (default: FALSE).
#'   Replacement removes the destination link, not the file it points to.
#' @param copy Copy original file contents instead of linking (default: FALSE).
#'   Cannot be combined with \code{symlink = TRUE}.
#' @param symlink Create absolute symbolic links directly to the original BAM
#'   and index files instead of hardlinks (default: FALSE).
#' @details
#' Outputs are \code{<output_data_dir>/<pair>/tumor.bam}, \code{normal.bam},
#' and their \code{.bam.bai} indexes. Index lookup prefers \code{.bam.bai},
#' then \code{.bai}, beside the supplied BAM path; if neither exists, it checks
#' beside the resolved BAM. Each index is resolved independently.
#'
#' Hardlinks require the original and destination to be on the same filesystem
#' and survive removal of the original path. Symbolic links can cross filesystems
#' but require the original paths to remain available to gOS. Neither link mode
#' copies BAM contents. Use \code{overwrite = TRUE} to replace previously broken
#' destination links. Filesystem operations use GNU \code{ln} and \code{cp}.
#' @examples
#' \dontrun{
#' lift_bam(cohort, output_data_dir, overwrite = TRUE) # Hardlink original files
#' lift_bam(cohort, output_data_dir, symlink = TRUE, overwrite = TRUE)
#' }
#' @return None
#' @export
lift_bam <- function(cohort, output_data_dir, cores = 1, overwrite = FALSE, copy = FALSE, symlink = FALSE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }

    if (all(copy) && all(symlink)) {
        stop("Choose either copy = TRUE or symlink = TRUE, not both")
    }

    transfer_command = if (all(copy)) "cp" else "ln"
    transfer_args = if (all(copy)) {
        c("-L", "-T", if (all(overwrite)) "--remove-destination" else "--no-clobber")
    } else {
        c("-n", "-T", if (all(symlink)) "-s" else "-L", if (all(overwrite)) "-f")
    }

    transfer_file = function(source, destination) {
        status = system2(transfer_command, c(transfer_args, "--", shQuote(source), shQuote(destination)))
        if (status != 0L) {
            stop(sprintf("Failed to transfer %s to %s (exit status %s)", source, destination, status))
        }
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

    # Process each sample in parallel
    
    upload_bam = function(i) {
        row <- lift_inputs[i,]
        .provenance_context <- provenance_context(row, cohort)
        .provenance_context$parameters <- list(copy = copy, symlink = symlink, overwrite = overwrite)
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            for (sample_type in c("tumor", "normal")) {
                bam_path = row[[paste0(sample_type, "_bam")]]
                if (is.null(bam_path) || anyNA(bam_path)) next

                # Resolve relative symlinks and chains before linking into the output.
                bam_source = normalizePath(bam_path, mustWork = TRUE)
                index_candidates = c(
                    paste0(bam_path, ".bai"),
                    paste0(tools::file_path_sans_ext(bam_path), ".bai"),
                    paste0(bam_source, ".bai"),
                    paste0(tools::file_path_sans_ext(bam_source), ".bai")
                )
                bai_path = index_candidates[file.exists(index_candidates)][1L]
                if (is.na(bai_path)) {
                    stop(sprintf("No BAM index found for %s", bam_path))
                }
                bai_source = normalizePath(bai_path, mustWork = TRUE)

                for (kind in c("bam", "bai")) {
                    source <- if (kind == "bam") bam_source else bai_source
                    destination <- file.path(pair_dir, paste0(sample_type, if (kind == "bam") ".bam" else ".bam.bai"))
                    link <- Sys.readlink(destination)
                    existed <- file.exists(destination) || (!is.na(link) && nzchar(link))
                    .provenance_context$sources <- list()
                    provenance_note(if (kind == "bam") bam_path else bai_path)
                    transfer_file(source, destination)
                    if (!existed || all(overwrite)) provenance_receipt(destination)
                }
            }
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }
    mclapply(seq_len(nrow(lift_inputs)), upload_bam, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
