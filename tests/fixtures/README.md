# Synthetic cohort fixtures

## Existing data review

The existing corpus at `/Users/user1/projects/skilift/inst/extdata/test_data`
occupies about 194 MB. It mainly contains gOS exports, Case Reports data,
oncotables, and isolated input examples. It has no complete Myeloseq input
cohort. Existing fixtures are unchanged.

Some existing tests require files or paths that are not available in this
checkout, including `cytoband.rds` and `unit_annotated.bcf`. The existing
`/Users/user1/projects/skilift/tests/run-all-tests.R` also names `PGVdb` rather
than `Skilift`. Those pre-existing issues are outside this fixture-only change.

## New corpus

Location: `/Users/user1/projects/skilift/inst/extdata/test_data/synthetic_cohorts`.

The generator uses constants only. Downloaded clinical files were inspected
for layout, headers, object classes, and column types; no clinical records,
serialized graph environments, identifiers, measured values, or paths were
copied. `MOCK*` names, positions, read counts, annotations, and dates are invented.
Chromosome names and lengths are public GRCh37 reference metadata. Fictional
genes must not be used for biological interpretation or clinical decisions.

There are three Cases per profile:

| Profile | Layout and contents |
| --- | --- |
| NGS | Two paired Cases and one tumor-only Case; native results-directory layout; 154-column OncoKB MAF; multiplicity, coverage and CBS GRanges; small-variant VCFs; BGZF/indexed SV VCFs; reduced CNA table; header-only fusion tables; fresh graph construction ingredients. |
| Myeloseq | Three nested DNA/RNA-named Cases; 152-column OncoKB MAF; panel multiplicity GRanges with unknown CN; full, filtered and reported small-variant subsets; BGZF/indexed fusion VCFs; sample-keyed QC JSON. No RNA reads or BAMs. |

In each profile, Case 1 has three SNVs, Case 2 has two SNVs and one reciprocal
breakend pair, and Case 3 has no SNVs or breakends. Empty results retain their
headers/classes. Cases 1 and 2 share a patient but have different pair IDs and
assay dates. Case 3 initially has no patient association.

All paths in the manifests are relative. Input paths are relative to their
profile directory; provenance and replacement paths are relative to the corpus
root. Missing optional inputs and unknown creation dates are blank, not invented.
The separate update examples link a patient later or replace only the first
Case's annotation file. The original inputs remain selected in the initial
manifest. There is no newest-file selection.

These manifests are test scaffolding, **not a new production domain schema**.
Small VCF, CNA, QC and selection-manifest examples retain relevant format
features, not every upstream annotation field. Annotation VCFs are plain text
instead of the source BCF so they can be inspected without additional tooling.
The panel's reported subset is an explicit input, not a proposed Skilift filter.

## Check or regenerate

Commands below assume the checkout location above; adjust it for another clone.
Required R packages: GenomicRanges (and its dependencies), Rsamtools, jsonlite.
No downloads or installations are performed.

```sh
Rscript /Users/user1/projects/skilift/tests/fixtures/check-synthetic-cohorts.R
Rscript /Users/user1/projects/skilift/tests/fixtures/generate-synthetic-cohorts.R /tmp/skilift-synthetic-fresh
```

The generator refuses a non-empty destination. Check regeneration in a fresh
directory before deliberately replacing the committed corpus. The check compares
committed and regenerated contents, then requires two fresh builds to have
identical bytes. RDS/BGZF bytes are not a cross-version binary contract.

The focused check verifies schemas, row counts, cross-file variant identity,
read counts/VAF, CN/coverage agreement, tabix reads, reciprocal breakends,
optional inputs, empty results, patient links, replacement provenance, generic
privacy guards, reproducible generation, and equivalent reads on two PSOCK
workers. It also injects an incorrect VAF into a temporary copy and requires
validation to fail. This is **not** a parallel production-conversion benchmark.

## Use in tests

```r
source("/Users/user1/projects/skilift/tests/fixtures/synthetic-cohorts.R")
root <- "/Users/user1/projects/skilift/inst/extdata/test_data/synthetic_cohorts"
inputs <- synthetic_case_inputs(root, "ngs") # absolute paths, no filesystem discovery
```

Patient metadata stays separate from the legacy input table. When using
`Cohort$new()`, map `pair` explicitly; do not let the legacy `patient_id` alias
choose Case identity. Split paired/tumor-only rows before running conversions
that use a single cohort-wide type.

Real gGraph objects can be constructed from the invented CBS nodes and edge
tables, but are not checked in as fake RDS stand-ins:

```r
graph_paths <- materialize_synthetic_graphs(root, tempfile("skilift-graphs-"))
inputs$jabba_gg <- unname(graph_paths[inputs$pair])
```

The optional test requires gGnome and data.table. gGnome was unavailable when
this corpus was added, so graph construction and full `lift_*()` conversions
have not been verified. No clinical graph or gWalk object is copied as a
shortcut. Native complex-event/gWalk, fusion-RDS, karyograph, signature and
full QC-report fixtures remain future coverage.

The testthat entry is
`/Users/user1/projects/skilift/tests/testthat/test-synthetic-cohorts.R`.
Its legacy Cohort and OncoKB collection smoke tests require a loadable Skilift
package. The standalone check does not require Skilift, testthat, devtools,
gUtils, or gGnome.
