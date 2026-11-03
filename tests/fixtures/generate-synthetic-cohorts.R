# Synthetic inputs only: this generator never opens clinical/source data.
# Run from any directory: Rscript <this-file> <new-empty-output-directory>

fixture_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing fixture dependencies: ", paste(missing, collapse = ", "))
}

fixture_variants <- function(n) {
  data.frame(
    chromosome = c("1", "1", "2"), position = c(1101L, 11101L, 10101L),
    reference = c("A", "C", "G"), alternate = c("G", "T", "A"),
    gene = c("MOCK1", "MOCK2", "MOCK3"),
    ref = c(30L, 32L, 38L), alt = c(10L, 8L, 2L),
    filter = c("PASS", "PASS", "LowQual"), stringsAsFactors = FALSE
  )[seq_len(n), , drop = FALSE]
}

fixture_ranges <- function(chromosome, start, width = 1L) {
  # Seqinfo moved out of GenomeInfoDb in recent Bioconductor releases.
  seqinfo_package <- if (requireNamespace("GenomeInfoDb", quietly = TRUE)) "GenomeInfoDb" else "Seqinfo"
  Seqinfo <- getExportedValue(seqinfo_package, "Seqinfo")
  GenomicRanges::GRanges(
    seqnames = chromosome, ranges = IRanges::IRanges(start = start, width = width),
    seqinfo = Seqinfo(c("1", "2"), c(249250621L, 243199373L), genome = "hg19")
  )
}

fixture_maf <- function(v, tumor, normal, profile) {
  # Column names are format metadata observed in the two source profiles.
  common <- strsplit(paste(
    "Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand",
    "Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2",
    "dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status",
    "Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_File Sequencer",
    "Tumor_Sample_UUID Matched_Norm_Sample_UUID HGVSc HGVSp HGVSp_Short Transcript_ID Exon_Number",
    "t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count all_effects Allele Gene Feature",
    "Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons",
    "Existing_variation ALLELE_NUM DISTANCE STRAND_VEP SYMBOL SYMBOL_SOURCE HGNC_ID BIOTYPE",
    "CANONICAL CCDS ENSP SWISSPROT TREMBL UNIPARC RefSeq SIFT PolyPhen EXON INTRON DOMAINS",
    "AF AFR_AF AMR_AF ASN_AF EAS_AF EUR_AF SAS_AF AA_AF EA_AF CLIN_SIG SOMATIC PUBMED",
    "MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE IMPACT PICK VARIANT_CLASS TSL",
    "HGVS_OFFSET PHENO MINIMISED GENE_PHENO FILTER flanking_bps vcf_id vcf_qual",
    "gnomADe_AF gnomADe_AFR_AF gnomADe_AMR_AF gnomADe_ASJ_AF gnomADe_EAS_AF gnomADe_FIN_AF",
    "gnomADe_NFE_AF gnomADe_OTH_AF gnomADe_SAS_AF vcf_pos"
  ), " ")[[1]]
  annotation <- strsplit(paste(
    "ANNOTATED GENE_IN_ONCOKB VARIANT_IN_ONCOKB MUTATION_EFFECT MUTATION_EFFECT_CITATIONS ONCOGENIC",
    "LEVEL_1 LEVEL_2 LEVEL_3A LEVEL_3B LEVEL_4 LEVEL_R1 LEVEL_R2 HIGHEST_LEVEL",
    "HIGHEST_SENSITIVE_LEVEL HIGHEST_RESISTANCE_LEVEL TX_CITATIONS LEVEL_Dx1 LEVEL_Dx2 LEVEL_Dx3",
    "HIGHEST_DX_LEVEL DX_CITATIONS LEVEL_Px1 LEVEL_Px2 LEVEL_Px3 HIGHEST_PX_LEVEL PX_CITATIONS",
    "GENE_SUMMARY VARIANT_SUMMARY TUMOR_TYPE_SUMMARY DIAGNOSTIC_SUMMARY PROGNOSTIC_SUMMARY",
    "MUTATION_EFFECT_DESCRIPTION ROW_ID Role"
  ), " ")[[1]]
  extra <- if (profile == "ngs") c("IS-A-HOTSPOT", "IS-A-3D-HOTSPOT") else c("REPORTED", "FUNC", "OLD_ID")
  columns <- c(common, extra, annotation, if (profile == "ngs") c("segment_cn", "segment_cn_low", "segment_cn_high"))
  n <- nrow(v)
  x <- as.data.frame(setNames(rep(list(rep("", n)), length(columns)), columns), check.names = FALSE)
  values <- list(
    Hugo_Symbol = v$gene, Center = "SYNTHETIC", NCBI_Build = "GRCh37",
    Chromosome = v$chromosome, Start_Position = v$position, End_Position = v$position,
    Strand = "+", Variant_Classification = "Missense_Mutation", Variant_Type = "SNP",
    Reference_Allele = v$reference, Tumor_Seq_Allele1 = v$reference, Tumor_Seq_Allele2 = v$alternate,
    Tumor_Sample_Barcode = tumor, Matched_Norm_Sample_Barcode = normal,
    Mutation_Status = "Somatic", Sequence_Source = "SYNTHETIC",
    HGVSc = paste0("c.1", v$reference, ">", v$alternate),
    HGVSp = c("p.Lys1Glu", "p.Pro1Ser", "p.Ala1Thr")[seq_len(n)],
    HGVSp_Short = c("p.K1E", "p.P1S", "p.A1T")[seq_len(n)],
    Transcript_ID = "MOCK_TRANSCRIPT", t_depth = v$ref + v$alt,
    t_ref_count = v$ref, t_alt_count = v$alt, Consequence = "missense_variant",
    FILTER = v$filter, vcf_pos = v$position,
    vcf_id = paste0(v$chromosome, ":", v$position, "_", v$reference, "/", v$alternate),
    flanking_bps = paste0("T", v$reference, "T"), vcf_qual = 60,
    ONCOGENIC = rep(c("Oncogenic", "Likely Oncogenic", "Unknown"), length.out = n),
    LEVEL_1 = rep(c("MOCK_DRUG", "", ""), length.out = n),
    GENE_SUMMARY = "Synthetic gene description; not clinical evidence.",
    VARIANT_SUMMARY = "Synthetic variant annotation.",
    MUTATION_EFFECT = "Unknown", MUTATION_EFFECT_DESCRIPTION = "Synthetic effect.",
    Role = "Oncogene", ROW_ID = seq_len(n)
  )
  if (profile == "myeloseq") values$REPORTED <- seq_len(n) == 1L
  if (profile == "ngs") values$segment_cn <- rep(c(2, 6, 2), length.out = n)
  for (name in names(values)) x[[name]] <- rep(values[[name]], length.out = n)
  x
}

fixture_multiplicity <- function(v, profile) {
  n <- nrow(v)
  gr <- fixture_ranges(v$chromosome, v$position)
  common <- list(REF = v$reference, ALT = v$alternate, ref = v$ref, alt = v$alt,
                 VAF = v$alt / (v$ref + v$alt), gene = v$gene,
                 variant.g = rep(paste0(v$chromosome, ":", v$position, "_", v$reference, "/", v$alternate), length.out = n))
  if (profile == "myeloseq") {
    # Observed panel estimates have unknown CN (logical all-NA columns), not zero CN.
    common$ref <- as.numeric(v$ref)
    common$alt <- as.numeric(v$alt)
    for (name in c("major.count", "minor.count", "altered_copies", "major_snv_copies",
                   "minor_snv_copies", "segment_cn", "segment_cn_low", "segment_cn_high")) common[[name]] <- rep(NA, n)
    common$vkey <- paste(v$chromosome, v$position, v$reference, v$alternate, sep = "_")
    common$IX <- seq_len(n)
  } else {
    for (name in c("paramRangeID", "allele", "gene_id", "feature_id", "rank", "cdna_pos",
                   "cds_pos", "protein_pos", "distance")) common[[name]] <- rep(NA_character_, n)
    common$QUAL <- rep(60, n)
    common$FILTER <- v$filter
    common$annotation <- rep("missense_variant", n)
    common$impact <- rep("MODERATE", n)
    common$feature_type <- rep("transcript", n)
    common$transcript_type <- rep("protein_coding", n)
    common$variant.c <- rep(paste0("c.1", v$reference, ">", v$alternate), length.out = n)
    common$variant.p <- c("p.Lys1Glu", "p.Pro1Ser", "p.Ala1Thr")[seq_len(n)]
    common$normal.ref <- rep(40L, n)
    common$normal.alt <- rep(0L, n)
    common$gt <- rep("0/1", n)
    common$cn <- rep(c(2, 6, 2), length.out = n)
    common$ref_denoised <- as.numeric(v$ref)
    common$alt_denoised <- as.numeric(v$alt)
    common$major.count <- as.numeric(pmax(v$ref, v$alt))
    common$minor.count <- as.numeric(pmin(v$ref, v$alt))
    for (name in c("total_snv_copies", "major_snv_copies", "altered_copies")) common[[name]] <- rep(1, n)
    common$minor_snv_copies <- rep(0, n)
  }
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(common)
  gr
}

fixture_write_vcf <- function(v, path, tumor, normal = "", compressed = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  header <- c(
    "##fileformat=VCFv4.2", "##source=SkiliftSyntheticFixtures", "##reference=GRCh37",
    '##contig=<ID=1,length=249250621>', '##contig=<ID=2,length=243199373>',
    '##FILTER=<ID=LowQual,Description="Synthetic low quality example">',
    '##INFO=<ID=ANN,Number=.,Type=String,Description="Synthetic SnpEff annotation">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", tumor,
            if (nzchar(normal)) normal), collapse = "\t")
  )
  records <- vapply(seq_len(nrow(v)), function(i) {
    annotation <- paste(c(v$alternate[i], "missense_variant", "MODERATE", v$gene[i], v$gene[i],
                         "transcript", "MOCK_TRANSCRIPT", "protein_coding", "1/1",
                         paste0("c.1", v$reference[i], ">", v$alternate[i]),
                         c("p.Lys1Glu", "p.Pro1Ser", "p.Ala1Thr")[i], "1/30", "1/30", "1/10", "", ""), collapse = "|")
    paste(c(v$chromosome[i], v$position[i], paste0("MOCK_VARIANT_", i), v$reference[i], v$alternate[i],
            60, v$filter[i], paste0("ANN=", annotation), "GT:AD:DP",
            paste0("0/1:", v$ref[i], ",", v$alt[i], ":", v$ref[i] + v$alt[i]),
            if (nzchar(normal)) "0/0:40,0:40"), collapse = "\t")
  }, character(1))
  plain <- if (compressed) tempfile(fileext = ".vcf") else path
  on.exit(if (compressed) unlink(plain), add = TRUE)
  writeLines(c(header, records), plain)
  if (compressed) {
    Rsamtools::bgzip(plain, dest = path, overwrite = FALSE)
    Rsamtools::indexTabix(path, format = "vcf")
  }
  invisible(path)
}

fixture_write_breakends <- function(path, present) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c("##fileformat=VCFv4.2", "##source=SkiliftSyntheticFixtures",
    '##contig=<ID=1,length=249250621>', '##contig=<ID=2,length=243199373>',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate breakend ID">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  if (present) lines <- c(lines,
    "1\t1250\tMOCK_BND_1\tA\tA]2:10250]\t60\tPASS\tSVTYPE=BND;MATEID=MOCK_BND_2",
    "2\t10250\tMOCK_BND_2\tC\tC]1:1250]\t60\tPASS\tSVTYPE=BND;MATEID=MOCK_BND_1")
  plain <- tempfile(fileext = ".vcf")
  on.exit(unlink(plain))
  writeLines(lines, plain)
  Rsamtools::bgzip(plain, dest = path, overwrite = FALSE)
  Rsamtools::indexTabix(path, format = "vcf")
  invisible(path)
}

generate_synthetic_cohorts <- function(output) {
  fixture_require(c("GenomicRanges", "IRanges", "S4Vectors", "Rsamtools", "jsonlite"))
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Output must be a new or empty directory; existing fixtures are never overwritten.")
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  table <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  }
  rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(x, path, version = 2)
  }
  json <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  }
  all_cases <- list()
  provenance <- list()
  for (profile in c("ngs", "myeloseq")) {
    cases <- inputs <- list()
    for (i in 1:3) {
      label <- sprintf("MOCK-%s-%03d", if (profile == "ngs") "NGS" else "MYELO", i)
      tumor <- paste0(label, if (profile == "ngs") "-T" else "-DNA")
      normal <- if (profile == "ngs" && i < 3) paste0(label, "-N") else ""
      pair <- if (profile == "ngs") {
        if (nzchar(normal)) paste0(tumor, "___vs___", normal) else tumor
      } else paste0(tumor, "_v1_", label, "-RNA_v1")
      cohort_type <- if (profile == "myeloseq") "heme" else if (nzchar(normal)) "paired" else "tumor_only"
      case_dir <- if (profile == "ngs") file.path(cohort_type, "results", pair) else label
      processed <- if (profile == "ngs") case_dir else file.path(case_dir, pair)
      root <- file.path(output, profile)
      case <- data.frame(pair = pair, patient_id = if (i == 3) "" else paste0("MOCK-PATIENT-", profile),
                         assay = if (profile == "ngs") "WGS" else "MYELOSEQ",
                         assay_date = sprintf("2000-%02d-01", i), cohort_type = cohort_type,
                         tumor_sample = tumor, normal_sample = normal, case_dir = case_dir)
      cases[[i]] <- case
      v <- fixture_variants(c(3L, 2L, 0L)[i])
      maf <- fixture_maf(v, tumor, normal, profile)
      input <- list(pair = pair, oncokb_snv = file.path(processed, if (profile == "ngs") "oncokb/merged_oncokb.maf" else "merged_oncokb.maf"),
                    multiplicity = file.path(processed, if (profile == "ngs") "snv_multiplicity/est_snv_cn_somatic.rds" else "est_snv_cn.rds"),
                    somatic_snvs = "", somatic_variant_annotations = "", tumor_coverage = "", segments_cbs = "",
                    oncokb_cna = "", oncokb_fusions = "", structural_variants = "", qc_metrics = "")
      table(maf, file.path(root, input$oncokb_snv))
      rds(fixture_multiplicity(v, profile), file.path(root, input$multiplicity))
      if (profile == "ngs") {
        input$somatic_snvs <- file.path(processed, "sage/somatic", paste0(pair, ".sage.pass_filtered.vcf.gz"))
        fixture_write_vcf(v[v$filter == "PASS", ], file.path(root, input$somatic_snvs), tumor, normal, TRUE)
        input$somatic_variant_annotations <- file.path(processed, "snpeff/somatic", paste0(pair, "_somatic.ann.vcf"))
        fixture_write_vcf(v, file.path(root, input$somatic_variant_annotations), tumor, normal)
        input$tumor_coverage <- file.path(processed, "dryclean/tumor/drycleaned.cov.rds")
        cov <- fixture_ranges(c("1", "1", "1", "2"), c(1L, 10001L, 20001L, 10001L), 10000L)
        for (name in c("background.log", "foreground.log", "input.read.counts", "signal", "center.all",
                       "median.chr", "foreground", "background", "log.reads")) S4Vectors::mcols(cov)[[name]] <- c(2, 6, 2, 2)
        rds(cov, file.path(root, input$tumor_coverage))
        input$segments_cbs <- file.path(processed, "cbs/seg.rds")
        seg <- fixture_ranges(rep("1", 3), c(1L, 10001L, 20001L), 10000L)
        seg$ID <- tumor
        seg$num.mark <- c(10, 10, 10)
        seg$seg.mean <- log2(c(2, 6, 2) / 2)
        rds(seg, file.path(root, input$segments_cbs))
        # Constructor ingredients, NOT a list masquerading as a serialized gGraph.
        table(data.frame(n1 = c(1L, 2L, 2L), n2 = c(2L, 3L, 2L), n1.side = 1L, n2.side = 0L,
                         type = c("REF", "REF", "ALT"), cn = c(2, 2, 4)),
              file.path(root, processed, "jabba/mock.edges.tsv"))
        input$oncokb_cna <- file.path(processed, "oncokb/merged_oncokb_cna.tsv")
        cna <- data.frame(HUGO_SYMBOL = "MOCK2", SAMPLE_ID = tumor, CANCER_TYPE = "SYNTHETIC",
                          ALTERATION = "AMPLIFICATION", ONCOGENIC = "Oncogenic", LEVEL_1 = "", LEVEL_2 = "",
                          LEVEL_R1 = "", LEVEL_Dx1 = "", LEVEL_Px1 = "", GENE_SUMMARY = "Synthetic gene.",
                          VARIANT_SUMMARY = "Synthetic amplification.", MUTATION_EFFECT = "Unknown",
                          MUTATION_EFFECT_DESCRIPTION = "Synthetic effect.", min_cn = 6, Role = "Oncogene")
        table(cna, file.path(root, input$oncokb_cna))
        input$oncokb_fusions <- file.path(processed, "oncokb/merged_oncokb_fusions.tsv")
        # Header-only fusion output is an observed producer shape, not a missing file.
        fusion <- as.data.frame(setNames(rep(list(character()), 12),
          c("silent", "in.frame", "geneA", "exonA", "aminoA", "geneB", "exonB", "aminoB", "min_cn", "mincn", "bp5p", "bp3p")))
        table(fusion, file.path(root, input$oncokb_fusions))
      } else {
        input$somatic_snvs <- file.path(processed, "SmallVariants.filtered.vcf")
        fixture_write_vcf(v[v$filter == "PASS", ], file.path(root, input$somatic_snvs), tumor)
        input$somatic_variant_annotations <- file.path(processed, "SmallVariants.vcf")
        fixture_write_vcf(v, file.path(root, input$somatic_variant_annotations), tumor)
        table(maf[seq_len(min(1L, nrow(maf))), , drop = FALSE], file.path(root, processed, "merged_oncokb.maf___REPORTED"))
        fixture_write_vcf(v[seq_len(min(1L, nrow(v))), ], file.path(root, processed, ".reported_set/effective_filtered.vcf"), tumor)
        table(data.frame(tsv_row = seq_len(nrow(v)), gene = v$gene, selected = seq_len(nrow(v)) == 1L,
                         decision = ifelse(seq_len(nrow(v)) == 1L, "selected", "excluded")),
              file.path(root, processed, ".reported_set/selection_manifest.tsv"))
        input$qc_metrics <- file.path(case_dir, "qc_merged.json")
        json(setNames(list(list(mean_depth = 200, total_reads = 10000L, pf_aligned_bases = 900000L,
                               mean_read_length = 100, fraction_of_reads_aligned = 0.9, pf_reads_aligned = 9000L,
                               uniformity = 0.95, on_target = 0.9, percent_aligned = 90)), tumor), file.path(root, input$qc_metrics))
      }
      input$structural_variants <- file.path(processed, if (profile == "ngs")
        paste0("sv_chimera_filter/tumor/", tumor, ".gridss.filtered.ffpe_filtered.vcf.gz") else "fusions_filtered.vcf.gz")
      fixture_write_breakends(file.path(root, input$structural_variants), present = i == 2L)
      inputs[[i]] <- as.data.frame(input)
      for (role in setdiff(names(input), "pair")) {
        if (!nzchar(input[[role]])) next
        provenance[[length(provenance) + 1L]] <- data.frame(profile = profile, pair = pair, role = role,
          path = file.path(profile, input[[role]]), source_created_at = if (i == 3) "" else sprintf("2000-%02d-02T00:00:00Z", i),
          frontend_updated_at = sprintf("2000-%02d-03T00:00:00Z", i))
      }
      if (i == 1) {
        revised <- maf
        revised$VARIANT_SUMMARY[1] <- "Synthetic revised annotation."
        revised$LEVEL_1[1] <- "MOCK_DRUG_V2"
        update_path <- file.path("updates", profile, "merged_oncokb.maf")
        table(revised, file.path(output, update_path))
        table(data.frame(pair = pair, role = "oncokb_snv", path = update_path,
                         source_created_at = "2000-02-01T00:00:00Z", frontend_updated_at = "2000-02-02T00:00:00Z"),
              file.path(output, "updates", profile, "replacement.tsv"))
      }
    }
    cases <- do.call(rbind, cases)
    table(cases, file.path(root, "cases.tsv"))
    table(do.call(rbind, inputs), file.path(root, "inputs.tsv"))
    table(data.frame(pair = cases$pair[3], patient_id = paste0("MOCK-PATIENT-LATE-", profile)),
          file.path(output, "updates", profile, "patient-links.tsv"))
    all_cases[[profile]] <- cases
  }
  table(do.call(rbind, provenance), file.path(output, "provenance.tsv"))
  genes <- fixture_ranges(c("1", "1", "2"), c(1001L, 11001L, 10001L), 500L)
  genes$gene_name <- genes$gene_id <- c("MOCK1", "MOCK2", "MOCK3")
  genes$type <- "gene"
  rds(genes, file.path(output, "reference/genes.rds"))
  # Uniform, explicitly synthetic source mtimes; consumers should use provenance.tsv for dates.
  files <- list.files(output, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  Sys.setFileTime(files, as.POSIXct("2000-01-01", tz = "UTC"))
  invisible(all_cases)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Usage: Rscript generate-synthetic-cohorts.R <new-empty-output-directory>")
  generate_synthetic_cohorts(args[1])
}
