# Test-only readers and independent checks for the committed synthetic corpus.

fixture_read_table <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "", colClasses = "character")
}

synthetic_case_inputs <- function(root, profile) {
  stopifnot(profile %in% c("ngs", "myeloseq"))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  x <- fixture_read_table(file.path(root, profile, "inputs.tsv"))
  for (name in setdiff(names(x), "pair")) {
    present <- !is.na(x[[name]])
    x[[name]][present] <- file.path(root, profile, x[[name]][present])
  }
  x
}

# gGraph/gWalk objects contain executable R6 environments. Never clone clinical
# objects to manufacture fixtures. Construct fresh graphs when gGnome is available.
materialize_synthetic_graphs <- function(root, output) {
  if (!requireNamespace("gGnome", quietly = TRUE) || !requireNamespace("data.table", quietly = TRUE)) {
    stop("Graph fixtures require gGnome and data.table; no fake gGraph RDS is substituted.")
  }
  cases <- fixture_read_table(file.path(root, "ngs/cases.tsv"))
  inputs <- synthetic_case_inputs(root, "ngs")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  paths <- setNames(character(nrow(cases)), cases$pair)
  for (i in seq_len(nrow(cases))) {
    nodes <- readRDS(inputs$segments_cbs[i])
    nodes$cn <- 2 * 2 ^ nodes$seg.mean
    edges <- data.table::fread(file.path(root, "ngs", cases$case_dir[i], "jabba/mock.edges.tsv"))
    graph <- gGnome::gG(nodes = nodes, edges = edges,
                       meta = list(purity = 0.75, ploidy = 2.2, gr.colorfield = "cn"))
    path <- file.path(output, cases$pair[i], "jabba.simple.gg.rds")
    if (file.exists(path)) stop("Refusing to overwrite graph fixture: ", path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(graph, path, version = 2)
    paths[i] <- normalizePath(path, winslash = "/")
  }
  paths
}

fixture_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

fixture_check_text <- function(text) {
  # Generic source-path and clinical identifier patterns, not a copied blacklist
  # of real patient identifiers. Generation from constants is the primary safeguard.
  fixture_assert(!any(grepl("/gpfs/|/Users/|/home/|TM[0-9]{2}-[0-9]|WG-[0-9]{2}-[0-9]", text)),
                 "Unexpected private path or non-synthetic sample identifier")
}

validate_synthetic_cohorts <- function(root) {
  for (package in c("GenomicRanges", "Rsamtools", "jsonlite")) {
    fixture_assert(requireNamespace(package, quietly = TRUE), paste("Missing fixture dependency:", package))
  }
  fixture_assert(dir.exists(root), "Missing synthetic fixture directory")
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  fixture_assert(length(files) > 0L, "Empty fixture corpus")
  fixture_assert(sum(file.info(files)$size) < 256 * 1024, "Synthetic corpus must remain below 256 KiB")
  fixture_assert(!any(grepl("\\.(bam|bai|cram|fastq|fq)(\\.|$)", files, ignore.case = TRUE)), "Raw sequencing files are forbidden")
  fixture_assert(!any(nzchar(Sys.readlink(files))), "Fixtures must not depend on symlinks")
  provenance <- fixture_read_table(file.path(root, "provenance.tsv"))
  fixture_assert(!anyDuplicated(provenance$path), "Duplicate provenance paths")
  for (profile in c("ngs", "myeloseq")) {
    cases <- fixture_read_table(file.path(root, profile, "cases.tsv"))
    raw_inputs <- fixture_read_table(file.path(root, profile, "inputs.tsv"))
    inputs <- synthetic_case_inputs(root, profile)
    fixture_assert(nrow(cases) == 3L && !anyDuplicated(cases$pair), "Each profile needs three distinct Cases")
    fixture_assert(identical(cases$pair, inputs$pair), "Case and input identities must agree")
    fixture_assert(all(grepl("^MOCK-", cases$pair)), "Case identifiers must be synthetic")
    fixture_assert(cases$patient_id[1] == cases$patient_id[2] && is.na(cases$patient_id[3]), "Repeated assay and unknown patient examples are required")
    fixture_assert(cases$assay_date[1] != cases$assay_date[2], "Repeated assays must have distinct dates")
    fixture_assert(!any(cases$pair %in% cases$patient_id), "Pair ID must not be a patient ID")
    late <- fixture_read_table(file.path(root, "updates", profile, "patient-links.tsv"))
    fixture_assert(identical(late$pair, cases$pair[3]) && !is.na(late$patient_id), "Late linkage must retain the original pair ID")
    for (i in seq_len(nrow(cases))) {
      for (role in setdiff(names(inputs), "pair")) {
        path <- inputs[[role]][i]
        if (is.na(path)) next
        fixture_assert(file.exists(path) && file.info(path)$size > 0, paste("Missing or zero-byte input:", path))
        rel <- raw_inputs[[role]][i]
        fixture_assert(!grepl("^/|(^|/)\\.\\.(/|$)", rel), "Manifest paths must be portable and contained")
        record <- provenance[provenance$profile == profile & provenance$pair == cases$pair[i] & provenance$role == role, ]
        fixture_assert(nrow(record) == 1L && record$path == file.path(profile, rel), "Input provenance must match its file")
        fixture_assert(record$frontend_updated_at == sprintf("2000-%02d-03T00:00:00Z", i), "Unexpected initial ingestion date")
        fixture_assert(if (i == 3L) is.na(record$source_created_at) else record$source_created_at == sprintf("2000-%02d-02T00:00:00Z", i), "Unknown creation date must remain unknown")
      }
      maf <- fixture_read_table(inputs$oncokb_snv[i])
      gr <- readRDS(inputs$multiplicity[i])
      expected <- c(3L, 2L, 0L)[i]
      fixture_assert(inherits(gr, "GRanges") && isTRUE(methods::validObject(gr)), "Multiplicity must be a valid GRanges")
      fixture_assert(nrow(maf) == expected && length(gr) == expected, "Unexpected variant count")
      fixture_assert(ncol(maf) == if (profile == "ngs") 154L else 152L, "MAF should retain the observed profile's column set")
      fixture_assert(identical(as.character(GenomicRanges::seqnames(gr)), maf$Chromosome), "MAF/multiplicity chromosomes differ")
      fixture_assert(identical(as.integer(GenomicRanges::start(gr)), as.integer(maf$vcf_pos)), "MAF/multiplicity coordinates differ")
      fixture_assert(identical(gr$REF, maf$Reference_Allele) && identical(gr$ALT, maf$Tumor_Seq_Allele2), "MAF/multiplicity alleles differ")
      fixture_assert(isTRUE(all.equal(gr$VAF, gr$alt / (gr$ref + gr$alt))), "VAF must match read counts")
      fixture_assert(all(maf$Tumor_Sample_Barcode == cases$tumor_sample[i]), "Wrong tumor sample")
      breakends <- unlist(Rsamtools::scanTabix(inputs$structural_variants[i]), use.names = FALSE)
      fixture_assert(length(breakends) == c(0L, 2L, 0L)[i], "Unexpected breakend count")
      if (length(breakends)) {
        fields <- strsplit(breakends, "\t")
        ids <- vapply(fields, `[`, character(1), 3L)
        mates <- sub(".*MATEID=", "", vapply(fields, `[`, character(1), 8L))
        fixture_assert(identical(sort(ids), sort(mates)), "Breakends must have reciprocal mates")
      }
      if (!is.na(cases$normal_sample[i]) && nzchar(cases$normal_sample[i])) {
        fixture_assert(all(maf$Matched_Norm_Sample_Barcode == cases$normal_sample[i]), "Wrong normal sample")
      }
      if (profile == "ngs") {
        cov <- readRDS(inputs$tumor_coverage[i])
        fixture_assert(inherits(cov, "GRanges") && length(cov) == 4L, "NGS coverage should contain four bins")
        fixture_assert(all(c("foreground", "background", "input.read.counts") %in% names(S4Vectors::mcols(cov))), "Missing coverage columns")
        if (expected > 0L) {
          hits <- GenomicRanges::findOverlaps(gr, cov)
          fixture_assert(length(hits) == expected && identical(gr$cn, cov$foreground[S4Vectors::subjectHits(hits)]), "Variant CN must agree with its coverage bin")
        }
        segments <- readRDS(inputs$segments_cbs[i])
        fixture_assert(identical(names(S4Vectors::mcols(segments)), c("ID", "num.mark", "seg.mean")), "CBS columns must match the producer schema")
        fixture_assert(file.exists(paste0(inputs$somatic_snvs[i], ".tbi")), "Missing VCF index")
        records <- unlist(Rsamtools::scanTabix(inputs$somatic_snvs[i]), use.names = FALSE)
        fixture_assert(length(records) == c(2L, 2L, 0L)[i], "BGZF/index or PASS selection is invalid")
        fixture_assert(nrow(fixture_read_table(inputs$oncokb_fusions[i])) == 0L, "Empty fusion output must retain its header")
      } else {
        fixture_assert(all(is.na(gr$segment_cn)) && all(is.na(gr$altered_copies)), "Panel CN must remain unknown, not fabricated")
        fixture_assert(is.na(inputs$tumor_coverage[i]), "Panel coverage GRanges are deliberately absent")
        qc <- jsonlite::read_json(inputs$qc_metrics[i])
        fixture_assert(identical(names(qc), cases$tumor_sample[i]), "QC must use a synthetic sample key")
        fixture_assert(qc[[1]]$mean_depth == 200, "Unexpected synthetic QC depth")
        reported <- fixture_read_table(file.path(dirname(inputs$oncokb_snv[i]), "merged_oncokb.maf___REPORTED"))
        fixture_assert(nrow(reported) == min(1L, expected), "Reported subset must remain distinct from the full table")
      }
    }
    replacement <- fixture_read_table(file.path(root, "updates", profile, "replacement.tsv"))
    fixture_assert(replacement$pair == cases$pair[1] && replacement$role == "oncokb_snv", "Replacement must target one input of the same Case")
    original <- fixture_read_table(inputs$oncokb_snv[1])
    revised <- fixture_read_table(file.path(root, replacement$path))
    fixture_assert(revised$VARIANT_SUMMARY[1] != original$VARIANT_SUMMARY[1], "Replacement must change an annotation")
    retained <- setdiff(names(original), c("VARIANT_SUMMARY", "LEVEL_1"))
    fixture_assert(identical(original[retained], revised[retained]), "Annotation replacement must preserve variant identity and measurements")
    fixture_assert(replacement$frontend_updated_at == "2000-02-02T00:00:00Z", "Replacement must have new ingestion provenance")
  }
  for (path in files) {
    if (grepl("\\.tbi$", path)) next
    if (grepl("\\.rds$", path)) {
      x <- readRDS(path)
      fixture_assert(inherits(x, "GRanges") && isTRUE(methods::validObject(x)), "Only generated GRanges RDS files belong in the base corpus")
      fixture_check_text(c(names(x), as.character(GenomicRanges::seqnames(x)), unlist(as.list(S4Vectors::mcols(x)), use.names = FALSE)))
    } else {
      con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
      text <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
      fixture_check_text(text)
      if (grepl("\\.vcf(\\.gz)?$", path)) {
        fixture_assert(text[1] == "##fileformat=VCFv4.2" && sum(startsWith(text, "#CHROM\t")) == 1L, "Invalid VCF header")
        header <- strsplit(text[startsWith(text, "#CHROM\t")], "\t")[[1]]
        records <- text[!startsWith(text, "#")]
        fixture_assert(all(lengths(strsplit(records, "\t")) == length(header)), "VCF column count differs from header")
      }
    }
  }
  invisible(TRUE)
}

check_synthetic_parallel_reads <- function(root) {
  paths <- unlist(lapply(c("ngs", "myeloseq"), function(profile) synthetic_case_inputs(root, profile)$multiplicity), use.names = FALSE)
  inspect <- function(path) {
    suppressPackageStartupMessages(requireNamespace("GenomicRanges", quietly = TRUE))
    x <- readRDS(path)
    list(rows = length(x), ref = sum(x$ref), alt = sum(x$alt))
  }
  serial <- lapply(paths, inspect)
  cluster <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cluster))
  fixture_assert(identical(serial, parallel::parLapply(cluster, paths, inspect)), "Parallel fixture reads differ from serial reads")
  invisible(TRUE)
}
