# Skilift
Faciliates easier loading and manipulation of genomic data for use with [PGV](https://github.com/mskilab-org/pgv) and[Case Report](https://github.com/mskilab-org/case-report).

# Installation
- To install directly from github: `devtools::install_github("mskilab-org/skilift")`
- To install from source:
    1. `git clone https://github.com/mskilab-org/skilift.git`
    2. `devtools::load_all("path/to/clone")`

# Usage with Case Reports
The basic workflow is as follows:
1. Load the cohort data into a Cohort object
2. Use lifter methods to generate the Case Report plot files

## Cohort
The Cohort class is a representation of a cohort of cases. It can be instantiated from a data table (a pairs table) or the path to the output directory of our nextflow case-reports pipeline (https://github.com/mskilab-org/nf-casereports/). It also takes a reference name (e.g 'hg38') to associate with the cohort (default reference is "hg19").

```r
# load a pairs table
cohort <- Cohort$new(pairs_table, reference_name = 'hg19')

# load a case-reports output directory
cohort <- Cohort$new(output_dir, reference_name = 'hg19')
```

The Cohort does a prefix search on the column names of the pairs table to extract the columns containing paths to data files relevant to Case Reports. It maps these columns to a set of well-defined names:

```r
default_col_mapping <- list(
    pair = c("pair", "patient_id", "pair_id", "sample"),
    tumor_type = c("tumor_type"),
    tumor_details = c("tumor_details"),
    disease = c("disease"),
    primary_site = c("primary_site"),
    inferred_sex = c("inferred_sex"),
    structural_variants = c("structural_variants", "gridss_somatic", "gridss_sv", "svaba_sv", "sv", "svs"),
    tumor_coverage = c("tumor_coverage", "dryclean_tumor", "tumor_dryclean_cov"),
    somatic_snvs = c("somatic_snvs", "sage_somatic_vcf", "strelka_somatic_vcf", "strelka2_somatic_vcf", "somatic_snv", "snv_vcf", "somatic_snv_vcf"),
    germline_snvs = c("germline_snvs", "sage_germline_vcf", "germline_snv", "germline_snv_vcf"),
    het_pileups = c("het_pileups", "hets", "sites_txt", "hets_sites"),
    somatic_snv_cn = c("somatic_snv_cn", "multiplicity"),
    germline_snv_cn = c("germline_snv_cn", "germline_multiplicity"),
    somatic_variant_annotations = c("somatic_variant_annotations", "annotated_bcf"),
    germline_variant_annotations = c("germline_variant_annotations", "annotated_vcf_germline"),
    oncokb_snv = c("oncokb_snv", "oncokb_maf", "maf"),
    oncokb_cna = c("oncokb_cna", "cna"),
    oncokb_fusions = c("oncokb_fusions", "oncokb_fusion", "fusion_maf"),
    jabba_gg = c("jabba_gg", "jabba_simple", "jabba_rds", "jabba_simple_gg"),
    karyograph = c("karyograph"),
    balanced_jabba_gg = c("balanced_jabba_gg", "non_integer_balance", "balanced_gg"),
    events = c("events", "complex"),
    fusions = c("fusions"),
    allelic_jabba_gg = c("allelic_jabba_gg", "lp_phased_balance", "allelic_gg"),
    activities_sbs_signatures = c("activities_sbs_signatures", "sbs_activities"),
    matrix_sbs_signatures = c("matrix_sbs_signatures", "sbs_matrix"),
    decomposed_sbs_signatures = c("decomposed_sbs_signatures", "sbs_decomposed"),
    activities_indel_signatures = c("activities_indel_signatures", "indel_activities"),
    matrix_indel_signatures = c("matrix_indel_signatures", "indel_matrix"),
    decomposed_indel_signatures = c("decomposed_indel_signatures", "indel_decomposed"),
    hrdetect = c("hrdetect", "hrd"),
    oncotable = c("oncotable"),
    estimate_library_complexity = c("estimate_library_complexity", "library_complexity_metrics", "est_lib_complex"),
    alignment_summary_metrics = c("alignment_summary_metrics", "alignment_metrics"),
    insert_size_metrics = c("insert_size_metrics", "insert_metrics"),
    wgs_metrics = c("wgs_metrics", "wgs_stats")
)
```

If the column names in the pairs table match (either exactly or by prefix) any of those in the above mapping, the Cohort object will automatically map the columns to the corresponding names. Priority is established by order (e.g if there are columns "pair" and "pair_id", the "pair" field of the cohort will be mapped to "pair" of the pairs table, not "pair_id").

In case the column names in the pairs table do not match those in default column mapping (or you want to change the default priority of the column names), you can pass a custom column mapping to the Cohort object:

```r
custom_col_mapping <- list(
    pair = "pid"
)

cohort <- Cohort$new(pairs_table, reference_name = 'hg19', col_mapping = custom_col_mapping)
```

Once the Cohort object is instantiated, you can access the data:

```r
# get the data table
cohort_table <- cohort$inputs
```

You can also validate the data (e.g check for missing values, invalid paths, etc) by calling the `validate_inputs` method on the Cohort object:

```r
# validate the data
cohort$validate_inputs()
```

## Lifters

THe lifter methods are used to generate the plot files for Case Reports. Each lifter takes a cohort object as input and generates the plot files (in parallel) for all samples in the cohort to a specified `output_data_dir`. Some lifters take additional parameters. 

The lifter methods are:

- `lift_copy_number_graph`: Create copy number gGraphs as JSON files. Use `is_allelic` to generate allelic copy number graph.
- `lift_denoised_coverage`: Create coverage plot as Apache Arrow files from GRanges.
- `lift_hetsnps`: Create hetSNP plots as Apace Arrow files from GRanges.
- `lift_filtered_events`: Create filtered events as JSON files from oncotable.
- `lift_multiplicity`: Create multiplicity plots as JSON files from SNVplicity (https://github.com/mskilab-org/SNVplicity). Set `is_germline = TRUE` to generate germline multiplicity plots.
- `lift_segment_width_distribution`: Create segment width distribution plots as JSON files. 
- `lift_signatures`: Create SBS and Indel signatures as JSON files from SigProfilerAssignment (https://github.com/AlexandrovLab/SigProfilerAssignment?tab=readme-ov-file)
- `lift_variant_qc`: Create variant QC plots as JSON files from Strelka2 or Sage VCF files.
- `lift_metadata`: Create metadata JSON files from assorted data (included qc metrics, etc).

E.g

```r
# instantiate a cohort object
pairs_table <- readRDS("./pairs.rds")
cohort <- Cohort$new(pairs_table, reference_name = 'hg19')

output_data_dir <- "~/public_html/case_reports_cohort_data"

# create copy number gGraphs
lift_copy_number_graph(cohort, output_data_dir = output_data_dir, is_allelic = FALSE, cores = 1)
lift_copy_number_graph(cohort, output_data_dir = output_data_dir, is_allelic = TRUE, cores = 1) # allelic

# create scatterplots
lift_denoised_coverage(cohort, output_data_dir = output_data_dir, cores = 1) # coverage
lift_hetsnps(cohort, output_data_dir = output_data_dir, cores = 1) #hetSNPs

# create filtered events
oncotable_dir <- "~/path/to/oncotable_dir"
cohort_w_oncotable <- create_oncotable(cohort, outdir = oncotable_dir, cores = 1)
lift_filtered_events(cohort_w_oncotable, output_data_dir = output_data_dir, cores = 1)

# create multiplicity plots
lift_multiplicity(cohort, output_data_dir = output_data_dir, is_germline = FALSE, cores = 1)
lift_multiplicity(cohort, output_data_dir = output_data_dir, is_germline = TRUE, cores = 1) # germline

# create segment width distribution plots
lift_segment_width_distribution(cohort, output_data_dir = output_data_dir, cores = 1)

# create signatures
lift_signatures(cohort, output_data_dir = output_data_dir, cores = 1)

# create variant qc plots
lift_variant_qc(cohort, output_data_dir = output_data_dir, cores = 1)

# create metadata
lift_metadata(cohort, output_data_dir = output_data_dir, cores = 1)
# if targeted panel or exome data, set genome_length to the size of the target region
lift_metadata(cohort, output_data_dir = output_data_dir, cores = 1, genome_length = 3e6)
```

Note that cohort only needs to be instantiated once, and can be reused for all
lifter methods. The lifter will automatically use the appropriate columns in
the cohort for that plot type. The lifter methods can be run in parallel, and
the number of cores can be specified by setting the `cores` parameter. By
default, the lifter methods will use 1 core.

## Production provenance

Skilift records **production receipts** for instrumented lifters automatically;
this does not require `use_provenance_table = TRUE`. Receipts describe inputs
observed by readers during production, not inputs inferred later from filenames.
The optional [table cache](#cached-provenance-tables-opt-in) adds separate
collected/curated artifacts and SHA-256 validation.

### Interaction points

| Interaction | Provenance behavior |
| --- | --- |
| `Cohort$new(...)` | Retains matched input-column names in `provenance_column_mapping`; construction from a pipeline output directory also records the input CSV path in `provenance_input_path`. An in-memory table has no input-file path to record. |
| `lift_all(...)` | Starts an upload record under `<output_data_dir>/.uploads/<id>/`: `cohort.rds` snapshots the cohort before lifting, `inputs.json` holds its input rows, and `upload.json` records the snapshot hash, reference, study, arguments, Skilift version and `skilift*` options. The returned cohort exposes `provenance_upload`. |
| Individual lifters | Instrumented readers capture sources; successful writers publish a receipt beside the output. Calling a lifter directly does not start an upload snapshot. Standalone readers outside a lifter context do not create receipts. |
| `lift_metadata(...)` | Captures field-level sources and literal overrides as well as the whole metadata artifact. By default it also calls `lift_datafiles_json()`; set `do_lift_datafiles_json = FALSE` to defer aggregation. |
| `lift_datafiles_json(output_data_dir)` | Folds current receipts into each existing case's `metadata.json`, then preserves the provenance object in cohort `datafiles.json`. Call this after additional standalone lifters to refresh the published view. `lift_all()` does this after its final outputs. |
| `write_provenance_table()` / `read_provenance_table()` | Write stage manifests and receipts, or hash-check and read a stage. When consumed inside an active lifter, the cached table joins its lineage. See [cache interrogation](#interrogate-and-curate-a-cached-table). |

For a manually composed workflow, finish with aggregation:

```r
lift_metadata(cohort, output_data_dir = output_data_dir, cores = 1,
              do_lift_datafiles_json = FALSE)
# Run any other required lifters before publishing the combined metadata.
lift_datafiles_json(output_data_dir = output_data_dir, cores = 1)
```

### Where evidence lives and what it means

- `<output_data_dir>/<pair>/.provenance/<output-basename>.json` is the
  latest artifact receipt. For example, `filtered.events.json` has receipt
  `.provenance/filtered.events.json.json`.
- `.provenance/field-<field>.json` records a metadata field, with a JSON Pointer
  such as `/0/tumor_median_coverage`. `.provenance/history/production-*.json`
  retains each published receipt; `recordId` identifies that historical record.
  This is receipt history, **not a backup of old output bytes**.
- After aggregation, `metadata.json` is still a one-row array. Its first row's
  `provenance` has `schemaVersion`, `caseId`, `recordedAt` and `artifacts`.
  `artifacts` includes both file receipts and field receipts; distinguish them
  by `path` and optional `pointer`. Aggregation requires an existing metadata
  row; it does not create one for receipt-only case directories.
- Receipts include production time, producer version, reference, parameters,
  selection details where supplied (for example `graphColumn`), upload context,
  sources and an `outputFingerprint`. Each source includes its supplied absolute
  `path`, role/cohort columns, original input column and table path when known,
  and `fingerprintAtUse` (`resolvedPath`, byte `size`, numeric `mtime`, `state`).
- A source's nested `production` is the upstream receipt captured when it was
  read. It is attached only if the upstream output's resolved path, size and
  modification time still match. Follow `sources[].production.sources[]` to
  trace transformed inputs back through earlier producers.
- `evidence = "captured"` means reader-observed sources. Artifact
  `sourceScope = "reader-session"` can include all reads accumulated in that
  worker, not a proven minimal dependency set. Field receipts use
  `sourceScope = "metadata-reader"`; `metadataInputs` records declared derived
  metadata dependencies. `evidence = "literal"` records a value without a file
  source, with origins such as `cohort-input`, `metadata-transform` or
  `added-field-override`. It does not mean missing evidence of a file read.

Fingerprints are **not content hashes**: equal size/mtime does not establish
byte identity. Receipts describe production-time state, not current validity;
aggregation does not revalidate files or rerun stale outputs. Old/uninstrumented
outputs can lack receipts, and an absent field receipt is not proof that a value
was manually supplied. Metadata receipts also retain `outputValues` (or
`outputValue` for a field), because aggregation itself changes `metadata.json`.
Keep hidden receipt/upload directories when retaining provenance, and review
paths, input rows and snapshots before sharing a report outside its intended
audience.

### Interrogate a report in R

These examples use an existing lifted case. Keep JSON as nested lists rather
than simplifying heterogeneous receipts into data frames:

```r
output_data_dir <- "build/data"
pair <- "CASE_ID"                         # a pair present in this output
case_dir <- file.path(output_data_dir, pair)
metadata <- jsonlite::read_json(file.path(case_dir, "metadata.json"),
                                simplifyVector = FALSE)[[1]]
records <- metadata$provenance$artifacts

# Which inputs and settings produced the metadata file?
artifact <- Filter(function(x) identical(x$path, "metadata.json") &&
                     is.null(x$pointer), records)
lapply(artifact, function(x) x[c("producedAt", "producer", "reference",
                                "selection", "parameters", "sources")])

# Why does this coverage value appear? Distinguish a file read from an override.
coverage <- Filter(function(x) identical(x$path, "metadata.json") &&
                     identical(x$pointer, "/0/tumor_median_coverage"), records)
lapply(coverage, function(x) x[c("outputValue", "evidence", "origin",
                                "sources", "metadataInputs")])
```

An empty result means no matching receipt was published. To inspect a receipt
before aggregation, read its `.provenance/<output-basename>.json` directly.
For a cohort-wide query, read `datafiles.json` with `simplifyVector = FALSE`
and select the row by `pair`; the same `provenance$artifacts` structure is there.

Trace upstream producers recursively, including sources of intermediate tables:

```r
source_tree <- function(record) {
  lapply(record$sources, function(source) {
    list(path = source$path, role = source$role,
         originalColumn = source$originalColumn,
         fingerprintAtUse = source$fingerprintAtUse,
         upstream = if (is.null(source$production)) NULL else
           source_tree(source$production))
  })
}
lapply(records, source_tree)
```

Recover the exact upload input row from a receipt, rather than a later cohort:

```r
uploads <- Filter(function(x) !is.null(x$upload$snapshot$path), records)
if (length(uploads)) {
  upload <- uploads[[1]]$upload
  snapshot <- upload$snapshot$path
  stopifnot(identical(digest::digest(file = snapshot, algo = "sha256"),
                      upload$snapshot$sha256))
  original_cohort <- readRDS(snapshot)
  original_cohort$inputs[original_cohort$inputs$pair == pair, ]
  upload[c("id", "startedAt", "arguments", "options", "cachedCohort", "remapping")]
}
```

A standalone lift may have no upload. Cases lifted across multiple uploads can
have different upload IDs on different artifact receipts; inspect the receipt
for the output of interest rather than assuming one upload owns the whole case.

Inspect previous productions of one output (optionally add a `pointer` filter
to select a metadata field):

```r
history_paths <- list.files(file.path(case_dir, ".provenance", "history"),
                            pattern = "\\.json$", full.names = TRUE)
history <- lapply(history_paths, jsonlite::read_json, simplifyVector = FALSE)
history <- Filter(function(x) identical(x$path, "metadata.json") &&
                    is.null(x$pointer), history)
lapply(history, function(x) x[c("recordId", "producedAt", "outputValues", "sources")])
```

### Cached-cohort replay and explicit remapping (internal integration hooks)

`Skilift:::provenance_load_cohort(path)` loads a saved cohort without
rediscovering inputs and records its path, fingerprint and SHA-256 in
`provenance_cached_cohort`. `Skilift:::provenance_apply_remapping(cohort, proposal)`
validates a complete replacement proposal before changing the in-memory inputs:

```r
replay <- Skilift:::provenance_load_cohort("cached-cohort.rds")
proposal <- list(schemaVersion = 1L, caseId = "CASE_ID", changes = list(
  list(column = "tumor_coverage", fromPath = "/original/coverage.rds",
       toPath = "/replacement/coverage.rds")
))
Skilift:::provenance_apply_remapping(replay, proposal)
# A subsequent lift_all(replay, ...) records the proposal in its upload.
```

The case must match exactly one row, each column can appear only once, each
`fromPath` must exactly match the recorded input, and every `toPath` must exist.
The saved cohort is not rewritten. Replacement content/reference compatibility
is not validated by this hook. These `:::` helpers are internal, not exported
stable APIs; there is no automatic path search or remapping when a source moves.

## Cached provenance tables (opt-in)

The provenance-table path is additive. Existing `oncotable()`,
`create_oncotable()`, `create_filtered_events()` and `lift_filtered_events()`
remain available with their existing behavior. `lift_all()` defaults to
`use_provenance_table = FALSE`; opt in explicitly:

```r
cohort <- lift_all(cohort, output_data_dir = "build/data", cores = 1,
                   use_provenance_table = TRUE)
```

This chooses `lift_mvp_provenance()`, which collects tables and consumes
the **collected** stage directly. It does not invent a curation step.
Its cache root is `<nextflow_results_path>/provenance` when available,
otherwise `<output_data_dir>/provenance`. For explicit cache location or
application-specific collection/curation boundaries, use the APIs below.
Annotation-only cohorts do not need a JaBbA graph; graph-dependent tracks
retain their own input requirements.

### Allele-keyed annotation projection

Declare an extra field once, then reuse it for table projection and the
filtered-events `additional_columns` interface:

```r
fields <- list(annotation_field(
  source = "review_score", as = "review_score",
  json = "review_confidence", title = "Review confidence", format = ".2%"
))

target <- data.table::data.table(
  Chromosome = "1", vcf_pos = 100L, vcf_ref = "A", vcf_alt = "T"
)
source <- data.table::data.table(
  Chromosome = "1", vcf_pos = 100L, vcf_ref = "A", vcf_alt = "T",
  review_score = 0.95
)
annotated <- merge_annotations_by_allele(target, source, fields)
json_columns <- annotation_json_columns(fields)
```

`source` names the incoming column; `as` names its table column; `json`
names the output field/column ID. The JSON spec's `dataIndex` reads `as`.
Optional `transform` is applied to projected values; `title` and `format`
control presentation. Fields without `json` remain table-only.
`annotation_json_columns()` includes the default echtvar columns unless
its `base` argument is overridden.

- `allele_key()` defaults to chromosome, **1-based VCF position**, REF and
  ALT. When REF/ALT are absent it uses the existing MAF reconstruction
  (including indel anchor-base conventions); callers must still supply
  chromosome/position and the required MAF fields. It does not reinterpret
  arbitrary coordinates or genome builds.
- Default matching uses the first source row for a duplicate key.
  `overwrite = FALSE` fills only missing target values; TRUE replaces them.
  Missing source fields emit a message and yield missing projected values.
  Explicit `target_key` / `source_key` support non-coordinate keys, such as
  fusion identifiers. `direction = "source_into_matched_target"` assigns
  into the first matching target row, with the last source assignment winning.
- `fallback = "oncokb_multiplicity"` delegates pairing to the existing
  `merge_oncokb_multiplicity()` geometry cascade before projecting fields:
  exact allele, gene/ALT overlap, adjusted indels and the existing gene/
  interval fallbacks are preserved, not replaced by a new matching rule.
  Unlike default projection, fallback mode leaves source-absent fields
  absent, preserving collector derivation such as VAF from allele counts
  when the multiplicity input does not supply a VAF column.

`collect_oncokb_with_sources()` accepts extra sources as
`list(table = <table or path>, fields = <field list>, key = NULL,
overwrite = FALSE)`. RDS paths use `readRDS`; other tabular paths use
`fread`. `lift_provenance_tables()` accepts `sources_for = function(row)`
to supply these per case, and `fields` for post-collection transforms.
For example, with one allele-keyed review TSV per pair:

```r
cache_root <- "build/provenance"
cohort <- lift_provenance_tables(
  cohort, root = cache_root, cores = 1,
  sources_for = function(row) list(list(
    table = file.path("review", paste0(row$pair, ".tsv")),
    fields = fields, overwrite = TRUE
  ))
)
cohort <- lift_filtered_events_from_tables(
  cohort, output_data_dir = "build/data", cores = 1,
  provenance_stage = "collected", additional_columns = json_columns
)
```

### Two-stage cache and integrity

The explicit stage names are `collected` and `curated`. For split workflows,
read with `read_provenance_table(root, pair, stage = "collected")`, apply
the application's curation, then write that result with
`write_provenance_table(table, root, pair, stage = "curated", inputs = ...,
parameters = ...)`. The two stages occupy **separate files**; writing
curated preserves the collected artifact and manifest entry. Rerunning a
stage replaces that stage, not an append-only history. Do not run competing
writers for the same pair/stage.

`create_provenance_table()` collects a single case; the cohort producer
`lift_provenance_tables()` adds `provenance_table_root` and `provenance_table`
to cohort inputs. Both `create_filtered_events_from_table()` and
`lift_filtered_events_from_tables()` default to consuming **curated**;
set `provenance_stage = "collected"` explicitly when no curation is needed.

Each pair directory contains stage tables and `manifest.json`. The manifest
records format/path, row count, columns, SHA-256, production time, input
hashes and parameters. Normal tables use parquet, preserving column types
and interpretation strings with embedded quotes/newlines. Tables containing
list columns use RDS automatically; `options(skilift_provenance_format =
"rds")` explicitly selects RDS. The manifest, not a guessed extension,
determines reader dispatch. `provenance_table_paths()` returns conventional
parquet paths; `write_provenance_table()` invisibly returns the actual path.

`read_provenance_table()` returns a data.table and verifies a manifest's
recorded SHA-256. Altered bytes raise `provenance_table_hash_mismatch` with
pair, stage and both hashes. Restore/rebuild the intended cache rather than
editing its digest. Without a manifest the reader can resolve conventional
filenames, but there is no recorded hash to check. Input receipts join the
active provenance lineage. Hash validation is **not** automatic invalidation
when upstream inputs or curation parameters change; rerun the owning stage.

### Interrogate and curate a cached table

For an existing collected cache, inspect the manifest before choosing a stage.
This example assumes the collected table has a `gene` column and illustrates
an application-defined gene selection, not a built-in clinical curation rule:

```r
cache_root <- "build/provenance"            # use the actual cache root
pair <- "CASE_ID"
manifest_path <- provenance_table_paths(cache_root, pair, "collected")$manifest
manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
names(manifest$stages)
manifest$stages$collected[c("path", "format", "rows", "columns", "sha256",
                            "producedAt", "inputs", "parameters")]

collected <- read_provenance_table(cache_root, pair, stage = "collected")
genes_to_report <- c("TP53", "KRAS")
collected[gene %in% genes_to_report]

# Persist an explicit curation boundary and its dependency on the collected file.
curated_path <- write_provenance_table(
  collected[gene %in% genes_to_report], cache_root, pair, stage = "curated",
  inputs = list(collected = manifest$stages$collected$path),
  parameters = list(genes_to_report = genes_to_report)
)
curated <- read_provenance_table(cache_root, pair, stage = "curated")
manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
manifest$stages$curated$inputs$collected   # recorded parent path and SHA-256
manifest$stages$curated$parameters
```

The manifest records file-level inputs and parameters, not per-cell annotation
lineage. Passing an in-memory table cannot recover its original file sources;
supply existing source paths through `inputs` when writing a manual stage.
Reading a manifest alone does not verify its table; `read_provenance_table()`
performs that check. Ordinary stage receipts and their history are also written
under the cache pair's `.provenance/` directory, separate from report outputs.

### Generic QC premerge and metadata ingest

QC sources are merged before metadata lift, not patched into emitted JSON.
A plugin declares its source keys, normalization and units with `qc_source()`
and `annotation_field()`; the generic API handles JSON, CSV and Picard input:

```r
qc_fields <- list(annotation_field("depth", as = "mean_depth",
                                   transform = as.numeric))
sources <- list(
  qc_source("primary_qc.csv", sample_key = "sample", fields = qc_fields,
            kind = "csv", precedence = "overwrite"),
  qc_source("fallback_qc.csv", sample_key = "sample", fields = qc_fields,
            kind = "csv", precedence = "fill")
)
merge_qc_sources(sources, out_path = "build/qc_merged.json",
                 field_order = "mean_depth")

cohort$inputs$qc_merged <- "build/qc_merged.json"
cohort$inputs$qc_sample_key <- cohort$inputs$pair
lift_metadata(cohort, output_data_dir = "build/data", cores = 1)
```

`sample_key` is a CSV/Picard column or a JSON name path resolving to the
sample ID. Optional `locate(object, key)` selects/transforms the record
(return NULL to skip); `key_transform(key)` normalizes the output ID.
Tabular input is read as character, so declare numeric/unit transforms.
Blank/`none` sample keys are skipped. Within a source's path list the last
record for a key wins; sources are then applied in list order with
overwrite/fill precedence. `merge_qc_sources()` requires `field_order` for
stable key order and accepts named `derived` functions of the merged record.

`add_coverage_metrics_from_qc_file(metadata, qc_merged, sample_key)` ingests
the record as `coverage_qc` (without a `pair` member), and uses `mean_depth`
for `tumor_median_coverage` by default. `lift_metadata()` consumes cohort
`qc_merged` / `qc_sample_key` columns before QC flags are evaluated; omitting
merged QC retains the legacy path. Set an explicit key for generic cohorts.
The Myeloseqer pipeline's `lib/mys_qc_plugin.R` is an application plugin
for IR/API, BAM and PDF precedence, not a universal metric mapping.

# Usage with PGV
The Skilift class manages a database of patient metadata and genomic plot data. It provides an interface to load, update, query, and manipulate the datafiles.json of PGV programmatically. This object is only necessary if you want to push data to PGV.

The key methods allow converting to/from JSON for data storage, validating the data, adding/removing plots, and generating plot JSON for visualization. The metadata and plots are stored as data tables for easy manipulation.

The conversion methods can also be used in isolation without having to instantiate a Skilift object. See tests (plot generation methods section) for runnable code. This is particularly useful for generating JSON files for use with Case Reports.

## Methods
The following methods are for the Skilift class:
- `initialize(datafiles_json_path, datadir, settings)`: Initialize a new Skilift object

- `update_datafiles_json()`: Update JSON data files on disk 

- `to_datatable(filter)`: Convert to a single data table, apply optional filter

- `add_plots(new_plots_dt)`: Add new plots to PGV

- `remove_plots(plots_to_remove_dt)`: Remove plots

- `list_higlass_tilesets(endpoint, username, password)`: List bigwig tilesets on higlass

- `validate()`: Validate metadata and plots

Methods for converting to JSON:
plot_metadata is a data table containing plot metadata and should have the following columns:
    - `patient.id`: Patient identifier (e.g '0124')
    - `source`: output datafile json filename (e.g complex.json)
    - `x`: filepath to the RDS object with the raw data (e.g ~/gg.rds)
    - `overwrite`: Whether to overwrite existing plot files
    - `ref`: Reference genome (e.g 'hg38')

datadir is the parent directory of the sample directory where the raw data file is stored 
e.g if the data file is stored in `data/0124/gg.rds` then the datadir is `data/`

settings is the path to the settings file (by default it will use the one included with the package). It's required for parsing seqlengths

## Fields

- `metadata`: Data table containing patient metadata 

- `plots`: Data table containing plot metadata

- `datadir`: Path to data directory

- `settings`: Path to settings file

- `higlass_metadata`: List containing endpoint, username, and password to Higlass server

### Constructor
```r
pgvdb <- Skilift$new(datafiles_json_path, datadir, settings)
```

Create a new Skilift object by passing the path to the datafiles.json, data directory, and settings file.

### update_datafiles_json
```r
pgvdb$update_datafiles_json()  
```

Update the datafiles JSON on disk with the current metadata and plots.

### to_datatable
```r
dt <- pgvdb$to_datatable(filter = c("patient.id", "HCC"))
```

Convert metadata and plots to a single data table, applying an optional filter.

### add_plots

```r
pgvdb$add_plots(plots_to_add, cores=2)  
```

Add new plots by passing a `data.table` containing required columns:

- `patient.id`: Patient identifier
- `x`: can be
    - list(list(server="", uuid="")): a list of the [server, uuid] 
    - list(GRanges), list(gWalk), list(gGraph): a list containing an object
    - Filepath to an RDS object
- `visible`: Whether plot is visible

If required columns are missing, an empty table will be returned with the
required column names. The `type` column (indiciating the type of plot:
scatterplot, genome, walk, bigwig, etc) is derived from `x` if not supplied by
the user. For `GRanges`, the `type` and `field` (the name of the score column)
must be specified by the user (either bigwig or scatterplot). The `source`
column (indicating the name of the plot file inside the pgv data directory) is
derived from the `type`, if not supplied by the user. 

Unless an `overwrite` column is set, it will not overwrite existing plot files,
instead it will just increment the filename of the new plot file (i.e
coverage.json -> coverage2.json). For GRanges containing bigwig data, the
GRange will be converted into a bigwig (with a unique name) and automaticaly
uploaded to the default mskilab higlass server (see above for changing to a
different server). If overwrite is set to TRUE, the converted bigwig file will
be deleted after uploading.

The `cores` parameter determines how many cores to use for parallel execution.
By default, it will not run in parallel (i.e use 1 core).

### remove_plots

```r
pgvdb$remove_plots(plots_to_remove, delete = FALSE)
```

Remove plots by passing a `data.table` with either:

- `patient.id`: Remove all plots for a patient 
- `patient.id` and `source`: Remove specific plots

Or alternatively with `patient.id`, `server`, and `uuid`.

The `delete` flag determines whether to also delete the source data files.

Passing just `patient.id` will remove all plots for that patient.

### validate

```r
pgvdb$validate()
```

Validate metadata and plot data, removing invalid entries. Is automatically
called when adding/removing plots. Useful if you make manual changes to the
pgvdb plots or metadata.


### init_pgv

```r 
pgvdb$init_pgv(pgv_dir, build=FALSE)
```

Initialize a pgv instance loaded with the data in pgvdb at `pgv_dir`. This
method will clone/pull the ![pgv repo](https://github.com/mskilab-org/pgv) and
create symlinks in the pgv data directory that point to your pgvdb data. 

If the `build` flag is set to `TRUE` it will build pgv instead of launching a
local instance (useful when running on a remote server or hpc). In that case,
you should set `pgv_dir` to be a directory inside whichever directory is served
by your remote server (e.g `public_html`).

### list_higlass_tilesets

```r 
pgvdb$list_higlass_tilesets(endpoint, username, password)
```

Return all bigwig tileset info on the higlass server as a data.table.

```r 
endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/" # dev endpoint
tilesets <- pgvdb$list_higlass_tilesets(
    endpoint,
    username = "username_here",
    password = "password_here"
)
```



### upload_to_higlass

```r 
pgvdb$upload_to_higlass(endpoint, datafile, filetype, datatype, coordSystem, name, username, password)
```

Upload a file to the higlass server. Will also add the file to the current
pgvdb instance. Note that you will need to upload a `chromSizes.tsv` file
first, before uploading other files. 

```r 
endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/" # dev endpoint
pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "chromSizes.tsv", package = "Skilift"),
    filetype = "chromsizes-tsv",
    datatype = "chromsizes",
    coordSystem = "hg38",
    name = "hg38",
    username = "username_here",
    password = "password_here"
)


pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "higlass_test_bigwig.bw", package = "Skilift"),
    name = "test_bigwig",
    filetype = "bigwig",
    datatype = "vector",
    coordSystem = "hg38",
    username = "username_here",
    password = "password_here"
)
```

### delete_from_higlass

```r 
 pgvdb$delete_from_higlass(endpoint, uuid, username, password)
```

Delete a file from the higlass server. Will also remove the plot from the pgvdb
instance. Tilesets can only be deleted by their uuid.
