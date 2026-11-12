test_that("OncoKB frameshifts survive oncotable creation across HGVS formats", {
  cases <- data.frame(
    protein = c(
      "p.Phe1104fs", "p.F1104fs", "legacyfs",
      "p.Phe1104LeufsTer27", "p.F1104Lfs*27",
      "p.Trp288CysfsTer12", "p.W288Cfs*12",
      "p.Arg97ProfsTer?", "p.R97Pfs*?", "p.(Arg97ProfsTer23)",
      "p.Ser1586Ter", "p.S1586*", "p.(Ser1586Ter)",
      "p.Val600Glu", "p.Trp24Ter", "p.Met1?", "p.Ter10Glnext*5",
      "p.?", "p.Glu10del", "", "", "",
      "", NA_character_, "p.?", "p.Val600Glu", "p.F1104Lfs*27",
      "p.F1104Lfs*27", "p.Ser1586Ter",
      "p.Phe1104fs", "p.Phe1104LeufsTer27", "p.Phe1104LeufsTer27",
      "p.Ser1586Ter", "p.Phe1104LeufsTer27", "p.Phe1104LeufsTer27"
    ),
    consequence = c(
      rep("frameshift_variant", 13),
      "missense_variant", "stop_gained", "start_lost", "stop_lost",
      "splice_acceptor_variant", "inframe_deletion", "upstream_gene_variant",
      "5_prime_UTR_variant", "synonymous_variant",
      rep("frameshift_variant", 4), "missense_variant", NA_character_, "synonymous_variant",
      "frameshift_variant&splice_region_variant", "frameshift_variant&splice_region_variant",
      "frameshift_variant&missense_variant", "frameshift_variant&stop_gained",
      "frameshift_variant&inframe_deletion", "frameshift_variant&5_prime_UTR_variant"
    ),
    expected = c(
      rep("trunc", 13), "missense", "trunc", "trunc", "trunc",
      "splice", "inframe", "upstream", "UTR", NA_character_,
      rep(NA_character_, 4), "missense", NA_character_, NA_character_,
      "trunc", "splice", "missense", "trunc", "inframe", "UTR"
    ),
    stringsAsFactors = FALSE
  )
  n <- nrow(cases)
  maf <- data.table::data.table(
    Hugo_Symbol = "TET2", Chromosome = "4", Start_Position = seq_len(n),
    End_Position = seq_len(n), vcf_pos = seq_len(n),
    Reference_Allele = "C", Allele = "G", Tumor_Seq_Allele2 = "G",
    flanking_bps = "ACA", Variant_Type = "SNP",
    HGVSc = paste0("c.", seq_len(n), "C>G"), HGVSp = cases$protein,
    Consequence = cases$consequence, Transcript_ID = "NM_001127208.3",
    GENE_SUMMARY = "", Role = "Tumor Suppressor", VARIANT_SUMMARY = "",
    MUTATION_EFFECT = "Unknown", MUTATION_EFFECT_DESCRIPTION = "",
    ONCOGENIC = rep(c("Oncogenic", "Unknown"), length.out = n),
    LEVEL_1 = rep(c("drug", "", ""), length.out = n),
    LEVEL_2 = "", LEVEL_R1 = "", LEVEL_Dx1 = "", LEVEL_Px1 = ""
  )
  # Indels are real MAF INS/DEL alleles, not SNP stand-ins. Classification must
  # still require the consequence/protein evidence rather than allele length.
  indels <- seq_len(13)
  maf[indels, `:=`(Variant_Type = "INS", Reference_Allele = "-", Allele = "GG",
                  Tumor_Seq_Allele2 = "GG")]
  maf[11:13, `:=`(Variant_Type = "DEL", Reference_Allele = "C", Allele = "-",
                  Tumor_Seq_Allele2 = "-")]
  path <- tempfile(fileext = ".maf")
  on.exit(unlink(path), add = TRUE)
  data.table::fwrite(maf, path, sep = "\t", na = "NA")

  collected <- Skilift:::collect_oncokb(path, verbose = FALSE)
  expect_identical(collected$type, cases$expected)
  expect_identical(collected$variant.p, cases$protein)
  expect_identical(collected$variant.c, maf$HGVSc)
  expected_tier <- ifelse(nzchar(maf$LEVEL_1), 1L,
                          ifelse(maf$ONCOGENIC == "Oncogenic", 2L, 3L))
  expect_identical(collected$tier, expected_tier)

  genes <- GenomicRanges::GRanges(
    "4", IRanges::IRanges(1, 1000), gene_name = "TET2",
    type = "gene", gene_type = "protein_coding"
  )
  result <- Skilift::oncotable(
    pair = "frameshift-regression", oncokb_snv = path,
    oncokb_fusions = NA_character_, oncokb_cna = NA_character_,
    fusions = NA_character_, jabba_gg = NA_character_, events = NA_character_,
    gencode = genes, cytoband = NULL, verbose = FALSE
  )
  retained <- !is.na(cases$expected)
  expect_identical(result$variant.c, maf$HGVSc[retained])
  expect_identical(result$type, cases$expected[retained])
  expect_equal(result$tier, expected_tier[retained])
})
