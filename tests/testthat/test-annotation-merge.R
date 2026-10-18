testthat::test_that("MAF allele keys preserve anchoring for substitutions and indels", {
  maf <- data.table::data.table(
    Chromosome = c("1", "2", "3", "X"), vcf_pos = c(10L, 20L, 30L, 40L),
    Variant_Type = c("SNP", "INS", "DEL", "MNP"),
    Allele = c("T", "AG", "A", "TG"), Reference_Allele = c("C", "-", "AC", "CC"),
    flanking_bps = c("[C]", "[A]", "[AC]", "[CC]"),
    Tumor_Seq_Allele2 = c("T", "G", "-", "TG")
  )
  testthat::expect_identical(allele_key(maf), c("1 10 C T", "2 20 A AG", "3 30 AC A", "X 40 CC TG"))
  rebuilt <- Skilift:::maf_reconstruct_refalt(maf)
  maf$vcf_ref <- rebuilt$vcf_ref_reconstructed
  maf$vcf_alt <- rebuilt$vcf_alt_reconstructed
  testthat::expect_identical(allele_key(maf), paste(maf$Chromosome, maf$vcf_pos, maf$vcf_ref, maf$vcf_alt))
})

testthat::test_that("annotation merges preserve fill, overwrite and duplicate direction", {
  target <- data.table::data.table(vkey = c("a", "b", "b", "c"), value = c(1L, NA_integer_, 3L, 4L))
  source <- data.table::data.table(vkey = c("a", "b", "b"), field = c(2.5, 5.5, 9.5))
  fields <- list(annotation_field("field", as = "value"))
  filled <- merge_annotations_by_allele(target, source, fields, target$vkey, source$vkey)
  testthat::expect_equal(filled$value, c(1, 5.5, 3, 4))
  replaced <- merge_annotations_by_allele(target, source, fields, target$vkey, source$vkey, overwrite = TRUE)
  testthat::expect_equal(replaced$value, c(2.5, 5.5, 5.5, NA_real_))
  assigned <- merge_annotations_by_allele(target, source, fields, target$vkey, source$vkey,
    overwrite = TRUE, direction = "source_into_matched_target")
  testthat::expect_equal(assigned$value, c(2.5, 9.5, 3, 4))
  testthat::expect_identical(target$value, c(1L, NA_integer_, 3L, 4L))
})

testthat::test_that("missing fields create named NA columns and transforms are vectorised", {
  target <- data.table::data.table(vkey = c("a", "b"))
  source <- data.table::data.table(vkey = c("a", "b"), raw = c(0.1, 0.2))
  testthat::expect_message(result <- merge_annotations_by_allele(target, source,
    list(annotation_field("absent", as = "optional"), annotation_field("raw", as = "pct", transform = function(x) x * 100)),
    target$vkey, source$vkey), "absent")
  testthat::expect_identical(result$optional, c(NA_character_, NA_character_))
  testthat::expect_equal(result$pct, c(10, 20))
})

testthat::test_that("JSON projection renames the event key without losing its table source", {
  spec <- annotation_json_columns(list(annotation_field("gnomad_211_af", json = "gnomad_af", title = "gnomAD_AF", format = ".2%")), base = list())[[1]]
  testthat::expect_identical(spec$id, "gnomad_af")
  testthat::expect_identical(spec$dataIndex, "gnomad_211_af")
  testthat::expect_identical(spec$rendererProps$format, ".2%")
})

testthat::test_that("multiplicity fallback preserves absent derived-field semantics", {
  target <- data.table::data.table(Chromosome = "1", Start_Position = 10L,
    End_Position = 10L, vcf_pos = 10L, Allele = "T", Reference_Allele = "A",
    Variant_Type = "SNP", Hugo_Symbol = "GENE", vcf_ref = "A", vcf_alt = "T")
  source <- GenomicRanges::GRanges("1", IRanges::IRanges(10L, 10L),
    REF = "A", ALT = "T", gene = "GENE", ref = 8L, alt = 2L)
  merged <- merge_annotations_by_allele(target, source,
    Skilift:::MULTIPLICITY_FIELDS, fallback = "oncokb_multiplicity", overwrite = TRUE)
  testthat::expect_identical(merged$ref, 8L)
  testthat::expect_identical(merged$alt, 2L)
  testthat::expect_false("VAF" %in% names(merged))
  testthat::expect_equal(merged$alt / (merged$ref + merged$alt), 0.2)
})
