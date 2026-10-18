# Original suite retained verbatim in tests/deprecated/test-filtered-events.R.
# Active tests use temporary, schema-complete fixtures and assert consumer output.
library(testthat)

setup({
  ot_test_paths <<- list(
    oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/oncotable.rds', package='Skilift'),
    no_oncokb_unit_oncotable = system.file('extdata/test_data/oncotable_test_data/no_oncokb_unit_oncotable.rds', package='Skilift'),
    unit_oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/unit_oncotable.rds', package='Skilift'),
    annotated_bcf = system.file('extdata/test_data/oncotable_test_data/annotated.bcf', package='Skilift'),
    unit_annotated_bcf = system.file('extdata/test_data/oncotable_test_data/unit_annotated.bcf', package='Skilift'),
    jabba_simple_gg = system.file('extdata/test_data/oncotable_test_data/jabba.simple.gg.rds', package='Skilift'),
    complex = system.file('extdata/test_data/oncotable_test_data/complex.rds', package='Skilift'),
    fusions = system.file('extdata/test_data/oncotable_test_data/fusions.rds', package='Skilift'),
    karyograph = system.file('extdata/test_data/oncotable_test_data/karyograph.rds', package='Skilift'),
    # oncokb_maf = system.file('extdata/test_data/oncotable_test_data/oncokb.maf', package='Skilift'),
    oncokb_snvcn_maf = system.file('extdata/test_data/oncotable_test_data/oncokb_snvcn.maf', package='Skilift'),
    oncokb_cna = system.file('extdata/test_data/oncotable_test_data/oncokb_cna.txt', package='Skilift'),
    oncokb_fusions = system.file('extdata/test_data/oncotable_test_data/oncokb_fusions.tsv', package='Skilift')
  )

  # gencode <<- process_gencode('/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds')
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  gencode <<- process_gencode(test_rds_path)

  cytoband_path = system.file("extdata/data/cytoband.rds", package = "Skilift")
  cytoband <<- readRDS(cytoband_path)
})

legacy_event_fixture <- function(root) {
  table <- data.table::copy(base::readRDS(ot_test_paths$unit_oncotable))
  # The old cache predates numeric copy counts and per-event complex footprints.
  numeric_columns <- intersect(c("value", "tier", "major_count", "minor_count",
    "major_snv_copies", "minor_snv_copies", "total_copies", "segment_cn", "ref", "alt", "VAF"), names(table))
  for (column in numeric_columns) data.table::set(table, j = column, value = as.numeric(table[[column]]))
  graph_events <- base::readRDS(ot_test_paths$complex)$meta$events
  footprints <- graph_events[, .(footprint = paste(footprint, collapse = ",")), by = type]
  complex_rows <- which(table$source == "complex")
  matched <- match(table$type[complex_rows], footprints$type)
  stopifnot(!anyNA(matched))
  table$footprint <- NA_character_
  table$footprint[complex_rows] <- footprints$footprint[matched]
  # Missing optional evidence has the same typed NA representation as collectors.
  table$altered_copies <- NA_real_
  table$is_multi_hit_per_gene <- NA
  path <- file.path(root, "legacy_events.rds")
  saveRDS(table, path)
  path
}

legacy_event_environment <- function(envir = parent.frame()) {
  withr::local_envvar(c(DEFAULT_GENOME = system.file(
    "extdata", "hg19.regularChr.chrom.sizes", package = "gGnome")), .local_envir = envir)
  # complex.rds is an events-only graph without purity/ploidy; the paired graph
  # stores these values and must be selected for copy-number calculations.
  withr::local_options(list(skilift_jabba_columns = "jabba_gg"), .local_envir = envir)
}

# Finite, shared seqinfo keeps graph and annotation overlaps on the same genome.
cna_test_ranges <- function(start, end, type, transcript_id = NA_character_,
                            gene_name = "GENE1", strand = "*") {
  GenomicRanges::GRanges(
    "1", IRanges::IRanges(start, end), strand = strand,
    seqinfo = GenomeInfoDb::Seqinfo("1", seqlengths = 10000),
    gene_name = gene_name, gene_id = gene_name, type = type,
    transcript_id = transcript_id
  )
}

cna_test_graph <- function(start, end, cn) {
  nodes <- GenomicRanges::GRanges(
    "1", IRanges::IRanges(start, end), strand = "+",
    seqinfo = GenomeInfoDb::Seqinfo("1", seqlengths = 10000),
    cn = cn, type = "REF", col = "#FF0000"
  )
  gGnome::gG(
    nodes = nodes,
    edges = data.table::data.table(
      n1 = seq_len(length(nodes) - 1L), n2 = seq_len(length(nodes) - 1L) + 1L,
      n1.side = 1L, n2.side = 0L
    ),
    meta = list(purity = 0.8, ploidy = 2)
  )
}

cna_test_transcript <- function(id, exon_start, exon_end,
                                cds_start = numeric(), cds_end = numeric()) {
  cna_test_ranges(
    c(min(exon_start), exon_start, cds_start),
    c(max(exon_end), exon_end, cds_end),
    c("transcript", rep("exon", length(exon_start)), rep("CDS", length(cds_start))),
    transcript_id = id
  )
}

cna_test_oncokb <- function(path, alteration, value) {
  template <- data.table::fread(system.file(
    "extdata", "test_data", "oncotable_test_data", "oncokb_cna.txt", package = "Skilift"
  ))
  symbol <- if (alteration == "Amplification") "CDK6" else "CDKN2A"
  row <- template[HUGO_SYMBOL == symbol & ALTERATION == alteration][1]
  row[, `:=`(HUGO_SYMBOL = "GENE1", ALTERATION = alteration, min_cn = value)]
  data.table::fwrite(row, path, sep = "\t")
}

test_that("process_gencode handles .rds input", {
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  result <- process_gencode(test_rds_path)
  expect_true(is(result, "GRanges"))
})

test_that("process_cytoband handles NULL input", {
  expect_warning(process_cytoband(NULL), "Using default cytoband")
})

test_that("process_cytoband handles .rds input", {
  cytoband_path <- system.file("extdata/data/cytoband.rds", package = "Skilift")
  result <- process_cytoband(cytoband_path)
  expect_true(is(result, "GRanges"))
  expect_true(all(c("band", "stain", "chrom_name", "chromband") %in% names(mcols(result))))
})

test_that("process_cytoband handles text file input", {
  # Create a temporary cytoband text file
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  
  result <- process_cytoband(temp_file)
  expect_true(is(result, "GRanges"))
  expect_true(all(c("band", "stain", "chrom_name", "chromband") %in% names(mcols(result))))
  expect_equal(length(result), 3)
  
  unlink(temp_file)
})

test_that("process_cytoband handles coarse parameter", {
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  result_fine <- process_cytoband(temp_file, coarse = FALSE)
  result_coarse <- process_cytoband(temp_file, coarse = TRUE)
  
  # Check that coarse bands have no decimal points
  expect_true(any(grepl("\\.", result_fine$chromband)))
  expect_false(any(grepl("\\.", result_coarse$chromband)))
})

test_that("process_cytoband handles zero-based coordinates", {
  # Create a temporary cytoband text file with zero-based coordinates
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  
  result <- process_cytoband(temp_file)
  expect_true(min(start(result)) > 0)  # Should be 1-based
  expect_equal(start(result)[1], 1)    # First start position should be 1
  
  unlink(temp_file)
})
test_that("collect_gene_fusions handles valid input", {
  fusions_path <- ot_test_paths$fusions
  result_fusions <- collect_gene_fusions(fusions_path, gencode, verbose = FALSE)
  expect_true(nrow(result_fusions) > 0)
  expect_true(all(c("gene", "vartype", "fusion_genes", "track", "type", "source", "fusion_gene_coords") %in% colnames(result_fusions)))
})

test_that("complex event collection preserves genomic evidence and excludes simple events", {
  result <- collect_complex_events(ot_test_paths$complex, verbose = FALSE)
  original <- base::readRDS(ot_test_paths$complex)$meta$events
  expected <- original[!type %in% c("del", "dup", "invdup", "tra", "inv")]
  expect_identical(result$footprint, expected$footprint)
  expect_true(all(nzchar(result$footprint)))
  expect_false(any(result$type %in% c("del", "dup", "invdup", "tra", "inv")))
})


# not testing this function as signatures needs to be updated to use sigprofiler
# test_that("collect_signatures handles valid input", {
#   signature_counts_path <- system.file('extdata/test_data/oncotable_test_data/signature_counts.txt', package='Skilift')
#   result_signatures <- collect_signatures(signature_counts_path, verbose = FALSE)
#   expect_true(nrow(result_signatures) > 0)
#   expect_true(all(c("value", "type", "track", "source") %in% colnames(result_signatures)))
# })

test_that("collect_gene_mutations handles valid input", {
  annotated_bcf_path <- ot_test_paths$unit_annotated_bcf
  result_mutations <- collect_gene_mutations(annotated_bcf_path, ot_test_paths$jabba_simple_gg, filter = 'PASS', verbose = FALSE)
  expect_true(nrow(result_mutations) > 0)
  expect_true(all(c("value", "type", "track", "source") %in% colnames(result_mutations)))
})

test_that("parse_oncokb_tier correctly assigns tiers", {
  # Create test data
  test_oncokb <- data.table(
    LEVEL_1 = c("drug1,drug2", NA, NA, NA),
    LEVEL_2 = c("drug3,drug4", "drug1,drug2", NA, NA),
    LEVEL_R1 = c(NA, NA, "LevelR1", NA),
    LEVEL_Dx1 = c(NA, NA, NA, NA),
    LEVEL_Px1 = c(NA, NA, NA, NA),
    ONCOGENIC = c("Oncogenic", "Likely Oncogenic", "Unknown", "Likely Neutral")
  )
  
  result <- parse_oncokb_tier(test_oncokb)
  
  # Check tier assignments
  expect_equal(result$tier, c(1, 1, 1, 3))
  expect_equal(result$tier_factor, factor(c("Clinically Actionable", "Clinically Actionable", 
                                          "Clinically Actionable", "VUS"),
                                        levels = c("Clinically Actionable", "Clinically Significant", "VUS")))
  
  # Check string concatenation
  expect_equal(result$tx_string, c("drug1,drug2,drug3,drug4", "drug1,drug2", NA, NA))
  expect_equal(result$rx_string, c(NA, NA, "LevelR1", NA))
  expect_true(all(is.na(result$dx_string)))
  expect_true(all(is.na(result$px_string)))
})

test_that("collect_oncokb_fusions handles missing file", {
  result <- collect_oncokb_fusions(NULL, gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
  
  result <- collect_oncokb_fusions("nonexistent.tsv", gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
})

test_that("collect_oncokb_fusions handles valid input", {
  result <- collect_oncokb_fusions(
    ot_test_paths$oncokb_fusions,
    gencode,
    cytoband,
    verbose = FALSE
  )
  
  # Check basic structure
  expect_true(is.data.table(result))
  expect_true(nrow(result) > 0)
  
  # Check required columns exist
  expected_cols <- c("gene", "gene_summary", "role", "value", "vartype", "type",
                    "tier", "tier_description", "variant_summary", "therapeutics",
                    "resistances", "diagnoses", "prognoses", "effect",
                    "effect_description", "fusion_genes", "fusion_gene_coords",
                    "track", "source")
  expect_true(all(expected_cols %in% names(result)))
  
  # Check specific values
  expect_equal(result$track[1], "variants")
  expect_equal(result$source[1], "oncokb_fusions")
  expect_true(all(result$tier %in% 1:3))
  expect_true(all(result$tier_description %in% c("Clinically Actionable", "Clinically Significant", "VUS")))
  source <- data.table::fread(ot_test_paths$oncokb_fusions)[silent == FALSE]
  fusion_key <- gsub("-", "::", source$FUSION, fixed = TRUE)
  frame <- source$in.frame[match(result$fusion_genes, fusion_key)]
  expect_true(all(result$vartype[which(frame)] == "inframe_fusion"))
  expect_true(all(result$vartype[which(!frame)] == "outframe_fusion"))
  
  # Check fusion gene coordinates format
  expect_true(all(grepl("^[0-9XY]+:[0-9]+-[0-9]+.,", result$fusion_gene_coords)))
})

test_that("collect_oncokb_fusions handles empty input file", {
  # Create temporary empty TSV file
  temp_file <- tempfile(fileext = ".tsv")
  empty_dt <- data.table(
    FUSION = character(),
    GENE_SUMMARY = character(),
    Role = character(),
    LEVEL_1 = character(),
    LEVEL_2 = character(),
    LEVEL_R1 = character(),
    LEVEL_Dx1 = character(),
    LEVEL_Px1 = character(),
    ONCOGENIC = character(),
    silent = logical(),
    in.frame = logical()
  )
  fwrite(empty_dt, temp_file)
  
  result <- collect_oncokb_fusions(temp_file, gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
  
  unlink(temp_file)
})

test_that("collect_oncokb_fusions handles malformed gene names", {
  # Create temporary TSV file with invalid gene names
  temp_file <- tempfile(fileext = ".tsv")
  malformed_dt <- data.table(
    FUSION = c("INVALID1-INVALID2", "NONEXISTENT1-NONEXISTENT2"),
    Hugo_Symbol = c("INVALID1-INVALID2", "NONEXISTENT1-NONEXISTENT2"),
    GENE_SUMMARY = c("Test summary 1", "Test summary 2"),
    VARIANT_SUMMARY = c("Test summary 1", "Test summary 2"),
    min_cn = c(2, 2),
    Role = c("Oncogene", "TSG"),
    LEVEL_1 = c(NA, NA),
    LEVEL_2 = c(NA, NA),
    LEVEL_R1 = c(NA, NA),
    LEVEL_Dx1 = c(NA, NA),
    LEVEL_Px1 = c(NA, NA),
    ONCOGENIC = c("Oncogenic", "Likely Oncogenic"),
    silent = c(FALSE, FALSE),
    in.frame = c(TRUE, FALSE),
    MUTATION_EFFECT = c("Missense", "Truncating"),
    MUTATION_EFFECT_DESCRIPTION = c("Test effect 1", "Test effect 2")
  )
  fwrite(malformed_dt, temp_file)
  
  result <- collect_oncokb_fusions(
    temp_file,
    gencode,
    cytoband,
    verbose = FALSE
  )
  # Should still return a valid data.table with NA coordinates for invalid genes
  expect_true(is.data.table(result))
  expect_true(all(is.na(result$fusion_gene_coords)))
  
  unlink(temp_file)
})

test_that("collect_oncokb_cna handles missing file", {
  result <- collect_oncokb_cna(NULL, verbose = FALSE)
  expect_equal(result$type, NA)
  expect_equal(result$source, "oncokb_cna")
})

test_that("transcript CNA coverage selects CDS or exons per transcript", {
  segments <- data.table::data.table(
    seqnames = "1", start = 1L, end = 100L, gene_name = "GENE1", type = "amp"
  )
  gencode <- cna_test_ranges(
    c(1, 1, 151, 1, 1, 1), c(100, 200, 200, 200, 200, 20),
    c("CDS", "exon", "CDS", "exon", "exon", "exon"),
    c("TX_CDS", "TX_CDS", "TX_CDS_MISS", "TX_CDS_MISS", "TX_EXON", "TX_SHORT")
  )
  detail <- get_transcript_cna_coverage(segments, gencode)
  expect_identical(detail$transcript_id, c("TX_CDS", "TX_CDS_MISS", "TX_EXON", "TX_SHORT"))
  expect_equal(detail$coverage_basis, c("CDS", "CDS", "exon", "exon"))
  expect_equal(detail$feature_width, c(100, 50, 200, 20))
  expect_equal(detail$overlapped_width, c(100, 0, 100, 20))
  expect_equal(detail$coverage_fraction, c(1, 0, 0.5, 1))
})

test_that("transcript CNA coverage unions intervals without cross-gene leakage", {
  segments <- data.table::data.table(
    seqnames = "1", start = c(1L, 51L, 1L), end = c(75L, 100L, 150L),
    gene_name = c("GENE1", "GENE1", "GENE2"), type = "amp", strand = "+"
  )
  gencode <- cna_test_ranges(
    c(1, 1, 51, 101), c(100, 100, 150, 150), "CDS", "TX1",
    gene_name = c("GENE1", "GENE1", "GENE1", "GENE2"), strand = "-"
  )
  detail <- get_transcript_cna_coverage(segments, gencode)
  expect_equal(detail$gene_name, c("GENE1", "GENE2"))
  expect_equal(detail$transcript_id, c("TX1", "TX1"))
  expect_equal(detail$feature_width, c(150, 50))
  expect_equal(detail$overlapped_width, c(100, 50))
  expect_equal(detail$coverage_fraction, c(2 / 3, 1))
  expect_equal(get_transcript_cna_coverage(segments[c(3, 1, 2)], gencode[c(4, 2, 3, 1)]), detail)
})

test_that("transcript CNA coverage handles unevaluable annotations", {
  segments <- data.table::data.table(
    seqnames = "1", start = 1L, end = 100L, gene_name = "GENE1", type = "amp"
  )
  empty <- get_transcript_cna_coverage(segments[0], NULL)
  expect_equal(nrow(empty), 0L)
  expect_identical(vapply(empty, typeof, ""), c(
    gene_name = "character", type = "character", transcript_id = "character",
    coverage_basis = "character", feature_width = "double",
    overlapped_width = "double", coverage_fraction = "double"
  ))
  unnamed <- data.table::copy(segments)[, type := NA_character_]
  expect_equal(get_transcript_cna_coverage(unnamed, NULL), empty)
  gencode <- cna_test_ranges(1, 200, "transcript", "TX1")
  expect_equal(get_transcript_cna_coverage(segments, gencode), empty)
  gencode <- cna_test_ranges(151, 200, "CDS", "TX1")
  expect_equal(get_transcript_cna_coverage(segments, gencode)$coverage_fraction, 0)
  # A zero-width CDS union is unevaluable, not a reason to fall back to exons.
  gencode <- cna_test_ranges(c(101, 1), c(100, 200), c("CDS", "exon"), "TX1")
  expect_equal(get_transcript_cna_coverage(segments, gencode), empty)

  gencode <- cna_test_ranges(1, 100, "CDS", "TX1")
  gencode$transcript_id <- NULL
  expect_error(get_transcript_cna_coverage(segments, gencode),
    "^gencode must contain gene_name, type, and transcript_id metadata$")
  for (bad_id in c("", NA_character_)) {
    gencode <- cna_test_ranges(1, 100, "CDS", bad_id)
    expect_error(get_transcript_cna_coverage(segments, gencode),
      "^CDS/exon records must have non-empty gene_name and transcript_id$")
  }
  for (bad_gene in c("", NA_character_)) {
    gencode <- cna_test_ranges(1, 100, "exon", "TX1", gene_name = bad_gene)
    expect_error(get_transcript_cna_coverage(segments, gencode),
      "^CDS/exon records must have non-empty gene_name and transcript_id$")
  }
  gencode <- cna_test_ranges(1, 100, "CDS", NA_character_, gene_name = "OTHER")
  expect_equal(get_transcript_cna_coverage(segments, gencode), empty)

  gencode <- cna_test_ranges(1, 1010, "gene")
  graph <- cna_test_graph(c(1, 1001), c(1000, 1010), c(10, 2))
  calls <- get_gene_ampdels_from_jabba(graph, gencode, nseg = NULL)
  amp <- calls[gene_name == "GENE1" & type == "amp"]
  expect_equal(nrow(amp), 1L)
  expect_identical(amp$max_transcript_fraction, 0)
  expect_equal(amp$transcript_coverage[[1]], empty[, .(
    transcript_id, coverage_basis, feature_width, overlapped_width, coverage_fraction
  )])
})

test_that("transcript CNA coverage gates AMP and HOMDEL by any transcript", {
  oncokb_path <- tempfile(fileext = ".tsv")
  withr::defer(unlink(oncokb_path))
  check_case <- function(graph, gene, features, alteration, event_type, fraction, passes) {
    gencode <- c(gene, features)
    value <- if (event_type == "amp") 10 else 0
    cna_test_oncokb(oncokb_path, alteration, value)
    pge <- gencode[gencode$type == "gene"]
    calls <- get_gene_ampdels_from_jabba(graph, pge, nseg = NULL)
    event <- calls[gene_name == "GENE1" & type == event_type]
    expect_equal(nrow(event), 1L)
    expect_equal(event$max_transcript_fraction, fraction)
    evidence <- event$transcript_coverage[[1]]
    if (nrow(evidence) > 0L) {
      expect_equal(max(evidence$coverage_fraction), event$max_transcript_fraction)
    } else {
      expect_equal(event$max_transcript_fraction, 0)
    }
    result <- collect_oncokb_cna(
      oncokb_path, graph, pge, amp.thresh = 4, del.thresh = 0.5, karyograph = NULL
    )
    retained <- result[result$gene %in% "GENE1" & result$type %in% event_type]
    expect_equal(nrow(retained), as.integer(passes))
    if (passes) {
      expect_equal(retained$value, value)
      expect_equal(retained$source, "oncokb_cna")
    }
  }

  amp_graph <- cna_test_graph(c(1, 1001), c(1000, 1010), c(10, 2))
  amp_gene <- cna_test_ranges(1, 1010, "gene")
  amp_fail <- cna_test_transcript("TX_FAIL", 1, 1010, 1001, 1010)
  amp_pass <- cna_test_transcript("TX_PASS", 1, 1010, 910, 1009)
  check_case(amp_graph, amp_gene, cna_test_transcript("TX90", 1, 1010, 911, 1010),
    "Amplification", "amp", 0.90, FALSE)
  check_case(amp_graph, amp_gene, amp_pass, "Amplification", "amp", 0.91, TRUE)
  check_case(amp_graph, amp_gene, amp_fail, "Amplification", "amp", 0, FALSE)
  check_case(amp_graph, amp_gene, c(amp_fail, amp_pass), "Amplification", "amp", 0.91, TRUE)
  check_case(amp_graph, amp_gene, cna_test_transcript("TX_EXON", 910, 1009),
    "Amplification", "amp", 0.91, TRUE)
  check_case(amp_graph, amp_gene, amp_gene[0], "Amplification", "amp", 0, FALSE)

  homdel_graph <- cna_test_graph(c(1, 11), c(10, 1000), c(0, 2))
  homdel_gene <- cna_test_ranges(1, 1000, "gene")
  homdel_fail <- cna_test_transcript("TX_FAIL", 1, 110, 11, 110)
  check_case(homdel_graph, homdel_gene, homdel_fail, "Deletion", "homdel", 0, FALSE)
  check_case(homdel_graph, homdel_gene, cna_test_transcript("TX_ONE", 1, 110, 10, 109),
    "Deletion", "homdel", 0.01, TRUE)
  check_case(homdel_graph, homdel_gene,
    c(homdel_fail, cna_test_transcript("TX_PASS", 1, 10, 1, 10)),
    "Deletion", "homdel", 1, TRUE)
  check_case(homdel_graph, homdel_gene, cna_test_transcript("TX_EXON", 1, 10),
    "Deletion", "homdel", 1, TRUE)
  check_case(homdel_graph, homdel_gene, homdel_gene[0], "Deletion", "homdel", 0, FALSE)
})

test_that("transcript CNA coverage preserves gene-body nomination and raw calls", {
  graph <- cna_test_graph(c(1, 101), c(100, 1000), c(10, 2))
  gencode <- c(cna_test_ranges(1, 1000, "gene"), cna_test_transcript("TX1", 1, 100, 1, 100))
  pge <- gencode[gencode$type == "gene"]
  calls <- get_gene_ampdels_from_jabba(graph, pge, nseg = NULL)
  expect_equal(nrow(calls[gene_name == "GENE1" & type == "amp"]), 0L)
  oncokb_path <- tempfile(fileext = ".tsv")
  withr::defer(unlink(oncokb_path))
  cna_test_oncokb(oncokb_path, "Amplification", 10)
  result <- collect_oncokb_cna(
    oncokb_path, graph, pge, amp.thresh = 4, del.thresh = 0.5, karyograph = NULL
  )
  expect_equal(nrow(result[result$gene %in% "GENE1" & result$type %in% "amp"]), 0L)

  graph <- cna_test_graph(c(1, 1001), c(1000, 1010), c(10, 2))
  gencode <- c(cna_test_ranges(1, 1010, "gene"), cna_test_transcript("TX1", 1, 1010, 1001, 1010))
  pge <- gencode[gencode$type == "gene"]
  calls <- get_gene_ampdels_from_jabba(graph, pge, nseg = NULL)
  expect_equal(nrow(calls), 1L)
  amp <- calls[gene_name == "GENE1" & type == "amp"]
  expect_equal(nrow(amp), 1L)
  expect_equal(amp$min_cn, 10)
  expect_equal(amp$min_quantile_cn, 10)
  expect_equal(amp$max_cn, 10)
  expect_equal(amp$total_node_width, 1000)
  evidence <- amp$transcript_coverage[[1]]
  expect_equal(evidence$coverage_basis, "CDS")
  expect_equal(evidence$feature_width, 10)
  expect_equal(evidence$overlapped_width, 0)
  expect_equal(evidence$coverage_fraction, 0)
  expect_equal(amp$max_transcript_fraction, max(evidence$coverage_fraction))
  graph_path <- tempfile(fileext = ".rds")
  withr::defer(unlink(graph_path))
  saveRDS(graph, graph_path)
  raw <- collect_copy_number_jabba(
    graph_path, pge, amp.thresh = 4, del.thresh = 0.5, karyograph = NULL
  )
  raw_amp <- raw[gene == "GENE1" & type == "amp"]
  expect_equal(nrow(raw_amp), 1L)
  expect_equal(raw_amp$value, 10)
})

test_that("create_oncotable persists each case's annotated variants and source evidence", {
  legacy_event_environment()
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  test_inputs <- data.table::data.table(
    pair = c("397089", "397090"),
    somatic_variant_annotations = ot_test_paths$unit_annotated_bcf,
    oncokb_fusions = ot_test_paths$oncokb_fusions,
    jabba_gg = ot_test_paths$jabba_simple_gg,
    events = ot_test_paths$complex,
    oncokb_snv = ot_test_paths$oncokb_snvcn_maf,
    oncokb_cna = ot_test_paths$oncokb_cna
  )
  test_cohort <- suppressWarnings(Cohort$new(test_inputs))
  updated <- create_oncotable(test_cohort, outdir = temp_dir, cores = 1,
    gencode = system.file("extdata", "test_data", "test_gencode_v29lift37.rds", package = "Skilift"),
    cytoband = system.file("extdata", "data", "cytoband.rds", package = "Skilift"))
  expected_snv <- Skilift:::collect_oncokb(ot_test_paths$oncokb_snvcn_maf, verbose = FALSE)
  expected_snv <- expected_snv[!is.na(type)]
  selected_fields <- c("Chromosome", "vcf_pos", "vcf_ref", "vcf_alt", "gene", "type")
  expected_rows <- as.data.frame(expected_snv[, ..selected_fields])
  expected_rows[] <- lapply(expected_rows, function(x) { x[x == ""] <- NA; x })
  for (pair in test_inputs$pair) {
    path <- updated$inputs$oncotable[match(pair, updated$inputs$pair)]
    expect_true(file.exists(path))
    table <- base::readRDS(path)
    expect_true(all(table$id == pair))
    selected_rows <- as.data.frame(table[source == "oncokb_maf", ..selected_fields])
    expect_equal(selected_rows, expected_rows)
    expect_true(any(table$source == "oncokb_fusions"))
    expect_true(all(nzchar(table$footprint[table$source == "complex"])))
    expect_true(file.exists(file.path(temp_dir, pair, "oncotable.txt")))
    expect_true(file.exists(file.path(temp_dir, pair, ".provenance", "oncotable.rds.json")))
  }
})

test_that("create_filtered_events writes complete gOS events from all variant classes", {
  legacy_event_environment()
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  fixture <- legacy_event_fixture(temp_dir)
  out_file <- file.path(temp_dir, "filtered.events.json")
  result <- create_filtered_events("test_sample", fixture, ot_test_paths$jabba_simple_gg,
    out_file, return_table = TRUE)
  expect_gt(nrow(result), 0L)
  expect_true(all(c("SNV", "AMP", "HOMDEL") %in% result$vartype))
  expect_setequal(result$vartype[result$type == "Complex SV"], c("chromothripsis", "rigma", "pyrgo", "bfb"))
  expect_true(any(result$type == "Fusion"))
  expect_true(all(result$sample == "test_sample"))
  expect_true(all(!is.na(result$type)))
  expect_true(all(!is.na(result$gene[result$type != "Complex SV"])))
  expect_true(all(nzchar(result$Genome_Location[result$type == "Complex SV"])))
  json <- jsonlite::fromJSON(out_file)
  expect_equal(nrow(json), nrow(result))
  expect_identical(json$Variant, result$Variant)
  expect_identical(json$Genome_Location, result$Genome_Location)
  expect_equal(json$Tier, result$Tier)
  no_return_file <- file.path(temp_dir, "no_return.json")
  create_filtered_events("test_sample", fixture, ot_test_paths$jabba_simple_gg,
    no_return_file, return_table = FALSE)
  expect_identical(jsonlite::read_json(no_return_file), jsonlite::read_json(out_file))
})

test_that("lift_filtered_events writes case JSON and summaries into new and existing directories", {
  legacy_event_environment()
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  fixture <- legacy_event_fixture(temp_dir)
  inputs <- data.table::data.table(pair = "397089", oncotable = fixture,
    jabba_gg = ot_test_paths$jabba_simple_gg, tumor_type = "BRCA", disease = "Breast Cancer")
  cohort <- suppressWarnings(Cohort$new(inputs))
  cohort$type <- "paired"
  outputs <- file.path(temp_dir, c("new", "existing"))
  dir.create(outputs[2L])
  for (output in outputs) {
    updated <- lift_filtered_events(cohort, output, cores = 1)
    path <- file.path(output, "397089", "filtered.events.json")
    events <- jsonlite::fromJSON(path)
    expect_gt(nrow(events), 0L)
    expect_true(all(events$sample == "397089"))
    expect_true(any(events$vartype == "SNV"))
    expect_true(any(events$type == "Complex SV"))
    expect_true(nzchar(updated$inputs$string_summary[1L]))
    expect_match(updated$inputs$string_summary[1L], "Complex SV", fixed = TRUE)
    receipt <- jsonlite::read_json(file.path(output, "397089", ".provenance", "filtered.events.json.json"))
    expect_true(normalizePath(fixture) %in% vapply(receipt$sources, `[[`, character(1), "path"))
  }
  expect_identical(jsonlite::read_json(file.path(outputs[1L], "397089", "filtered.events.json")),
    jsonlite::read_json(file.path(outputs[2L], "397089", "filtered.events.json")))
})

