# Simple per-Case provenance

Status: **Draft for review; no implementation approved by this document.**

Repository: `/Users/user1/projects/skilift`.
Baseline: `refactor` at `d4ada9e8bcf76c5f34766ee689a7387e11a75425`.
The separate `provenance` branch is not a dependency, implementation template,
or migration target. Do not merge or copy its tracking framework.

## 1. Purpose

For each published gOS output, answer:

1. Who ran Skilift?
2. When did that run start?
3. Which immediate file inputs were selected for this output's conversion?
4. What paths and filesystem details were observed for those inputs at that time?

Export one small JSON document per Case. A frontend component can read its
output's entry without opening the input files or parsing the associated data.
The same document must also be usable from R or a CLI.

This is **file-level provenance**, not a claim that every value or variant in an
output came from every listed file. Suggested tooltip heading: **Inputs supplied
to this output's conversion**.

## 2. Proposed decisions

- Keep the existing Cohort, case-level table, `Cohort$new()` and `lift_*()` APIs.
- Use the existing normalized `pair` value as Case identity. Do not require a
  patient ID or introduce a new Case class.
- Add one ordinary package source file:
  `/Users/user1/projects/skilift/R/provenance.R`.
- Keep the pure record-building helpers, filesystem helpers and provenance lifter
  in that file, following the existing `create_*()` / `lift_*()` convention.
- Store only two new run-level facts: `operator` (the OS account, omitted if
  unavailable) and `started_at`. Do not accept a caller-supplied operator override.
- Derive file observations from the selected paths, but **save the observations
  with the output's provenance**. Do not recompute old provenance from today's
  filesystem or today's Cohort.
- Store provenance outside the data outputs, in one `provenance.json` per Case.
  For example: `/data/demo/MOCK-001/provenance.json`.
- Record immediate inputs only. Do not inspect another file's provenance, follow
  upstream manifests, or reconstruct how an input was made.
- Do not hash input contents by default. Paths, realpath, size and modification
  time are cheap observations, not cryptographic content identities.
- Keep only the provenance of current published outputs. No run archive,
  publication history, database, receipt directory or recursive lineage graph.
- Support **backfill** from an existing `cohort.rds` plus the directory where
  outputs were already published. Backfill generates `provenance.json` without
  rerunning the scientific conversions or modifying their data files.

The operator and start time are the only new run-level facts. The document also
needs structural information (Case, output, input roles) and saved filesystem
observations. Those cannot be omitted if the result is to describe a particular
published output rather than the current input selection.

## 3. Scope and limits

### Included

- Per-output mappings to direct file arguments selected by the lifter.
- Direct resource files selected through explicit arguments or package defaults,
  when a conversion uses them. A reference file read by Skilift is a direct input,
  not upstream provenance.
- Original supplied paths, absolute supplied paths and resolved paths.
- Size, modification time and an observation timestamp.
- Creation/birth time only when genuinely available; otherwise JSON `null`.
- Partial reruns: replace entries only for outputs actually replaced by the run.
- Backfill from a saved `cohort.rds` and an existing output directory, without
  rerunning conversions or modifying data outputs.
- A common representation for JSON, Arrow, images and linked/copied files.
- New, file-backed custom fields through the same explicit input-binding helper.

### Not included

- Earlier pipelines, raw sequencing ancestry, PMKB matching history, or automatic
  reading of provenance stored beside or inside input files.
- Row-level, variant-level or attribute-level attribution.
- Reproducing a complete analysis from this document. File paths and stat results
  do not capture all software, parameters, references or immutable bytes needed
  for reproducibility.
- Copying clinical values, whole Cohort rows, command lines, environment variables
  or arbitrary configuration into the document.
- Changing filtering, interpretation, input precedence or file discovery rules.
- A general executor, scheduler, event protocol, new world model or new cache.
- Recovering the original publication time or exact input selections during
  backfill. Backfill uses the saved Cohort's paths and current filesystem state;
  it does not scan directories to reconstruct dynamically selected files. The
  provenance timestamp means the latest provenance-generation run, not a verified
  historical output-publication time.
- Automatic download, remote object listing or credential-based source access.
- General rollback of past outputs. Only best-effort rollback within the current
  Case publication is required. The feature must log failures and must not falsely
  claim that a failed commit left all previous data intact.

An output may also depend on supplied values and settings. This first version
records **file inputs**, not a complete provenance record for those values.
An empty input list means no direct file inputs were recorded for that conversion;
it must not be presented as proof that the output has no external origin.

For a normal conversion, the entry records selected paths and observations from
that run. For backfill, use the same schema and labels, without a backfill marker
or qualification. The OS account and timestamp identify the latest provenance
generation (the backfill operation); paths are taken from the saved Cohort and
stat observations are made now. Do not claim a verified historical publication
time or exact historical input use. This product convention treats the latest
provenance record as the relevant run, including for backfill.

## 4. Corrections to the two-facts-plus-stat idea

### 4.1 Capture observations at production time

Filesystem facts are derived, but they are not permanent. A path can be replaced,
a symlink can change target, and a file can disappear after conversion.

Observe selected inputs just before the conversion uses them, and retain that
snapshot when the output succeeds. Export the saved observations. A later
inspection may compare current filesystem facts with the snapshot, but must not
replace the snapshot or describe those new facts as the original provenance.

The input files are expected to remain stable while a conversion runs. Cheap
before/after observations can flag detectable changes; they cannot prove stable
bytes. Strong content identity would need a separate requirement.

### 4.2 Filesystem access is not pure

Calling `file.info()`, resolving symlinks, reading the clock, discovering the OS
user and reading/writing JSON all depend on external state. They belong in small
I/O helpers.

Selecting bindings from explicit values, assembling records, validating records
and merging records can be pure. Tests should pass fixed times, user names and
file observations to those functions.

Do not parse `ls` output or launch a shell process per file. Prefer the existing R
filesystem APIs, with platform-specific birth-time support only if justified.

### 4.3 The Cohort does not define the complete output-to-input mapping

A Cohort contains paths for many outputs, alongside non-path values. Some inputs
come from function arguments or package resources. Some paths are selected by
fallback rules. Therefore:

- Do not attach every Cohort path to every output.
- Do not infer a path solely because a character value happens to name a file.
- Use the actual column/argument selection made by the relevant lifter.
- Record a fallback source instead of a higher-priority source that was not selected.
- Keep input bindings next to the code that owns the selection. Pass the resolved
  choices into the shared provenance helpers.

A declared file argument is not proof that every row or attribute used its
contents. This feature does not instrument every read or establish cell-level
causality. Known alternative branches must still be represented accurately.

### 4.4 Creation time and file ownership have limited meaning

- POSIX `ctime` is metadata-change time, **not creation time**.
- Base R's `file.info()` does not provide a portable filesystem birth time.
- Even a genuine birth time describes the local filesystem object, not necessarily
  when its data was originally generated. A copy can have a new birth time.
- The file owner is not necessarily its creator, and is not the Skilift operator.
- Do not infer a download/retrieval date from mtime, ctime or birth time.

For the initial portable implementation, `created_at: null` is acceptable on every
platform. Do not add a dependency merely to fill that field. The implementation
must never substitute another date under the creation-time label.

## 5. Existing code and integration boundaries

These findings come from the current `refactor` checkout, not the other branch.

| Existing owner | Relevant behavior | Consequence for provenance |
| --- | --- | --- |
| `/Users/user1/projects/skilift/R/Cohort.R` | Maps inputs to conventional columns; `DEFAULT_JABBA()` selects a graph column using options and filesystem checks. | Use the selected column from the lifter. Do not independently repeat that decision inside a supposedly pure helper. |
| `/Users/user1/projects/skilift/R/lift-wrappers.R` | Has `required_columns`, `lift_all()` and `lift_mvp()`. | Reuse orchestration conventions, but do not treat `required_columns` as a complete dependency specification. It is primarily an execution gate and includes non-file values. |
| `/Users/user1/projects/skilift/R/filtered-events.R` | `lift_filtered_events()` passes an oncotable and selected graph to the converter, then may write highlights and a summary in the same error boundary. Heme filtering reads an additional reference. | Record immediate table/graph/resource inputs, not the files that originally made the oncotable. Account for each successfully written file separately. A later highlights failure must not erase evidence of a successful events write. |
| `/Users/user1/projects/skilift/R/metadata.R` | `lift_metadata()` prefers `oncokb_snv` over `somatic_snvs`, selects a graph, pairifies some inputs and can locate a Purple QC file. Other producers can locate reference files. | Bind the effective inputs after selection, including directly selected resources. A final scan of the original Cohort is insufficient. Preserve existing value and fallback behavior. |
| `/Users/user1/projects/skilift/R/scatterplot.R` | `lift_denoised_coverage()` selects `tumor_coverage`, writes Arrow and returns `invisible(NULL)`. | A useful small integration test. Preserve its public return value; reporting must happen internally. Account for any additional resources used by the conversion. |
| `/Users/user1/projects/skilift/R/metadata.R` | `lift_datafiles_json()` combines Case metadata. | Do not embed this provenance document into Case metadata or enlarge the catalog with its contents. |

### Why a new file alone is not enough

The shared logic belongs in the new source file. However, small integration calls
are needed at existing selection and successful-write points. Otherwise a
standalone exporter cannot know which files a particular output received, whether
an old file was skipped, or whether a replacement succeeded.

Do not solve that gap by overriding `readRDS()`, `fread()` or other readers globally,
walking call stacks, or copying the other branch's tracking machinery.

### Output identity

Use the output's path relative to its Case directory as its artifact key. This is
a data identifier, not an input filesystem path. It remains useful when a dataset
is moved or served under another URL. For example, `filtered.events.json` is the
key for the file `/data/demo/MOCK-001/filtered.events.json` in the example dataset.
Nested outputs need their complete Case-relative key; a basename alone can collide.

Record every actual emitted file, not just one entry per lifter. One lifter can
produce several files or a variable number of signature/plot outputs. The writer
supplies the key; the exporter must not guess filenames or recursively inventory
the output directory to infer success.

## 6. Run context

A run context is an ordinary R list with two values:

| Field | Meaning |
| --- | --- |
| `operator` | OS-reported execution account, omitted if unavailable. Do not substitute Git author, file owner or scheduler identity. |
| `started_at` | UTC timestamp captured at the start of the provenance-generation run, before Case workers start. |

An OS account may be a service account. This is an attribution label, not verified
human identity. Do not substitute the Git author, file owner or scheduler job
owner. If the OS account cannot be read, omit `operator` without warning.

`lift_all()` should create one context and pass it to its nested lifters. A
standalone lifter creates its own context unless the caller supplies one to share
across a sequence of calls. Do not capture this information in `Cohort$new()`:
construction and conversion may happen on different days or under different users.

Each output entry carries the provenance-generation context for its latest
record. There is no single document-level operator/date that overwrites the
identity of retained outputs. A partial rerun can legitimately leave different
run contexts in one Case's manifest.

Suggested compatibility surface:

- Append an optional `provenance` argument to integrated lifters and `lift_all()`;
  its default is `TRUE`. A caller may explicitly set `FALSE` to disable provenance
  I/O and recover legacy behavior.
- Current default-on production hooks cover coverage, hetsnps, BAM/index,
  copy-number graph, filtered events, Heme highlights, metadata, multiplicity,
  signatures, variant QC, segment-width output, multiplicity histograms, and the
  allelic and 2D purity/ploidy plots. `lift_coverage_jabba_cn()` and
  `lift_purple_sunrise_plot()` remain uninstrumented; custom output names are not
  in the backfill registry. The `ppfit.json` backfill entry uses the saved graph
  structure to determine whether coverage is also a direct input. Do not interpret
  a missing output entry as a conversion failure or proof of no input.
- Resolve the OS account and start time once for the run, before workers start.
  Tests may inject a fixed context internally; callers do not supply an operator
  override.
- Backfill is an explicit mode of the public provenance function and accepts an
  existing Cohort or path to a serialized `cohort.rds`, plus the directory where
  outputs were published. It does not rerun output lifters.

No new required constructor argument, global option, global mutable recorder or
persistent run registry is needed.

## 7. File observations

For each selected direct file input, record:

| Field | Meaning |
| --- | --- |
| `role` | Input role in this conversion, such as `oncotable` or `graph`. |
| `cohort_column` | Effective normalized Cohort column, or `null` for a direct argument/resource. |
| `supplied_path` | Path exactly as supplied to the selection/conversion boundary. |
| `absolute_path` | Supplied path made absolute using the conversion's working directory; it must not discard the supplied symlink spelling. |
| `resolved_path` | Realpath observed then, or `null` when resolution is unavailable. |
| `observed_at` | UTC time of the filesystem observation, not export time. |
| `stat_status` | `ok`, `unavailable`, or `unsupported`. |
| `size_bytes` | Observed target size, or `null`. |
| `modified_at` | Observed target modification time in UTC, or `null`. |
| `created_at` | Genuine target birth time when available, otherwise `null`. |

The target's stat values are used consistently for symlinks. This version does not
need the complete symlink chain, link ownership, permissions or access time.
Preserve distinct supplied paths/roles even when they resolve to one target;
underlying stat observations may be shared internally.

A failed stat is not automatically proof that a file is missing. Permission or
filesystem errors can look similar through base APIs. Use `unavailable` when the
reason cannot be established. Unknown values must stay `null`, never fabricated
zero sizes or dates. A stat failure does not by itself decide whether the data
conversion can proceed.

Missing optional bindings (`NULL`, `NA`, blank values) are not paths and are omitted.
Inspect only designated file-input columns and explicit file arguments. Accept
scalar character/factor paths, including single-cell list columns. Omit inline
tables, objects, numeric/logical values and multi-valued cells; never flatten or
convert their contents into path strings. Metadata records include Conpair files,
exclude the unused germline-SNV argument, and include Purple best-fit only when
the graph branch is not selected.
An explicitly selected path whose stat fails remains recorded as unavailable if
the output is produced. A zero-byte input is a present file, not a missing value.

This first version performs no network requests. An unsupported URI must not be
passed to local path normalization as if it were a filesystem path. If such a
binding is accepted by a converter, preserve a safe URI and report unavailable
stat details; implementing remote metadata retrieval is a separate change.

## 8. JSON contract

One document per Case, for example `/data/demo/MOCK-001/provenance.json`.
This example uses invented paths and values. Its source list is abbreviated, not
a declaration that these are the only inputs of the existing events converter.

```json
{
  "schema_version": 1,
  "pair": "MOCK-001",
  "scope": "direct-file-inputs",
  "outputs": {
    "filtered.events.json": {
      "lifter": "lift_filtered_events",
      "run": {
        "operator": "demo-operator",
        "started_at": "2026-01-02T09:00:00Z"
      },
      "output_observation": {
        "observed_at": "2026-01-02T09:00:04Z",
        "size_bytes": 4096,
        "modified_at": "2026-01-02T09:00:04Z"
      },
      "inputs": [
        {
          "role": "oncotable",
          "cohort_column": "oncotable",
          "supplied_path": "/work/demo/MOCK-001/current.rds",
          "absolute_path": "/work/demo/MOCK-001/current.rds",
          "resolved_path": "/work/demo/MOCK-001/oncotable-v1.rds",
          "observed_at": "2026-01-02T09:00:01Z",
          "stat_status": "ok",
          "size_bytes": 8192,
          "modified_at": "2026-01-01T16:00:00Z",
          "created_at": null
        }
      ]
    }
  }
}
```

Contract rules:

- `outputs` is an object keyed by Case-relative artifact identifiers; `inputs` is
  always an array, including zero or one input. Empty objects must stay objects.
- The output observation is captured after a successful write, or at backfill
  time for an existing file. Its size and mtime can detect some later changes;
  they do not establish content equality.
- `run.started_at` identifies when this provenance record was generated. For
  backfill it is not a verified historical output-publication time.
- The JSON does not label records as backfilled. Backfill uses the same shape and
  presents current saved-Cohort paths as the latest provenance, per product policy.
- The lifter name is supplied by code, not an extra operator-maintained fact.
- All generated timestamps use ISO 8601 UTC with `Z`; retain available precision.
- JSON `null` represents unknown facts. Do not write `NA`, `NaN` or guessed dates.
- Source roles and output keys have deterministic ordering for stable inspection.
- No hashes, upstream records, row locators, attempt history or full data values
  are required by schema version 1.
- Validate Case identity and artifact keys. Reject traversal or keys outside the
  configured Case directory.
- Reject unsupported schema versions or malformed existing documents rather than
  silently replacing them. Preserve unrelated output entries and extension fields
  within a supported schema when merging an update.

## 9. Proposed functions and ownership

All common logic initially lives in
`/Users/user1/projects/skilift/R/provenance.R`. Names below are proposed, not APIs
that exist today.

### Pure helpers

- An internal binding helper takes a Case row, explicit selected column names and
  resolved resource arguments. It returns named source paths. It does not call
  `DEFAULT_JABBA()`, scan directories, stat files or inspect the clock.
- An internal record builder takes the artifact key, lifter name, run context,
  selected input observations and output observation.
- `create_provenance(pair, records, previous = NULL)` validates and assembles a
  Case document. It replaces only the supplied artifact entries and preserves
  unrelated entries. It performs no I/O.

### I/O helpers

- Resolve the run context from the OS account and clock; permit fixed context
  injection internally for tests without exposing an operator override.
- Observe selected paths through filesystem APIs. Accept injected observations or
  an observation function for tests. Keep cache lifetime explicit.
- Read an existing Case document when present.
- Write JSON to a temporary sibling file, then rename it into place. Do not write
  directly over the only valid copy. Clean up the temporary file on failure.
- Stage and publish a Case's changed data files and manifest with recoverable
  backups, preserving unrelated entries and restoring the previous Case state on
  reported commit failures where possible.

### Public lifter

Current public backfill entry point:

```r
lift_provenance(cohort, output_data_dir, force = FALSE, cores = 1)
```

`cohort` may be an existing Cohort or a path to a serialized `cohort.rds`.
`output_data_dir` identifies where outputs were previously published. The
function loads the serialized object with `readRDS()` and accepts the legacy
Cohort/table representation already supported by Skilift; it does not require
new provenance fields in the saved RDS. It identifies existing supported Case
outputs using an explicit output/input registry, applies those lifters' input
selection rules to the saved Cohort, observes the selected paths and existing
output files, and writes per-Case `provenance.json`. It does not rerun scientific
conversion or modify data outputs. Entries use the same JSON fields as production
records, with no backfill label. The OS account and `started_at` describe this
provenance generation. Input paths come from the saved Cohort, and file stats are
observed at backfill time. This is the accepted product meaning of the latest
provenance; it does not prove the historical operator, publication timestamp or
exact source use. Do not guess a historical time from file mtimes.

Production lifters use an internal per-run context to collect staged records and
publish the Case's data and manifest together. This context is not a public
standalone `publish` mode.

If a Case's `provenance.json` already exists, backfill leaves that Case unchanged
and reports it as skipped by default. `force = TRUE` explicitly permits replacing
that Case's existing provenance document. This force behavior applies to backfill;
normal production updates continue to merge successful output entries and preserve
unrelated entries.

The Cohort supplies Case membership and input paths. For captured production
records, it is not used to replace source paths. For backfill, its current inputs
are used to populate the same output/input mapping as a normal record. The
function groups updates by Case, calls the pure assembler and publishes documents.
Existing lifters may use the same export path after collecting successful
production records while preserving their current public return values.
Re-exporting a saved production record must not change its original run or
observation times.

## 10. Production and publication behavior

### Normal flow

1. Resolve the run context in the parent process.
2. Let the lifter make its existing input selections.
3. Capture the corresponding file observations before conversion.
4. Run the existing conversion without changing its scientific behavior.
5. For each successfully written artifact, capture its output observation and
   build a provenance record.
6. Merge and export records for successful artifacts. Preserve other entries.

Do not infer success from output-file existence, a changed mtime, or a non-`NULL`
return alone. An output may already exist, and a successful lifter may return
`invisible(NULL)`. Integrate at the actual write/replace boundary.

### Important cases

| Event | Required behavior |
| --- | --- |
| New successful output | Add an entry with its captured sources and current run. |
| Valid empty output | Record it as a successful output. Empty data is not absent provenance. |
| Successful replacement | Replace that output's entry, including its source observations and run context. |
| Output skipped or reused without rewriting | Keep the previous entry unchanged. Do not attribute the old file to the new operator/run. |
| Failed attempt; previous output known to be untouched | Retain or restore the previous entry. Do not insert the attempted new source paths. |
| Failure while replacing selected outputs and their manifest | Restore the previous data files and manifest where possible; log the Case failure and continue with other Cases. Do not knowingly leave mismatched output/provenance versions published. |
| One artifact succeeds before another conversion or queued record fails | Publish the successful artifact and its record. Leave the failed artifact's prior output/entry unchanged, or absent when none existed. Log the failure; it must not discard other available artifacts. |
| Explicit output removal | Remove only that entry when removal is confirmed. Do not treat a missing input as an instruction to delete output provenance. |
| Existing output with no captured provenance, normal mode | Leave it unrecorded; the UI reports provenance unavailable. Do not infer production provenance from current paths. |
| Existing output with no provenance, explicit backfill mode | Add an entry using the saved Cohort's input selection and current stat observations, with the ordinary schema and no backfill marker. The product convention treats this as the latest provenance. |
| Existing `provenance.json` during backfill | Skip that Case and preserve the file by default; log the skip. Replace only when `force = TRUE`. |
| A pipeline patches an output after Skilift | The saved record describes Skilift's production, not the later edit. A detected mismatch is a warning, not a reason to fabricate new sources. |

### Two files are not one filesystem transaction

Atomic replacement of the provenance JSON alone does **not** atomically replace
both the data and its provenance. Current lifters also write some data outputs
directly. This limitation must be addressed explicitly, not hidden by the word
"atomic".

When provenance is enabled, stage all changed outputs for a Case and the complete
updated provenance document before replacing published files. Keep recoverable
prior outputs and the old manifest until the new Case publication completes.
A conversion or individual queue-record failure does not require a complete Case:
publish the other successfully staged and recorded artifacts so the frontend can
display what is available. Do not replace the failed artifact or fabricate its
provenance. If writing the combined manifest fails, leave that publication's prior
outputs and manifest untouched, log the Case and reason, and continue with other
Cases. Do not let one Case's provenance failure abort remaining Case work.

The publication step must write `provenance.json` once per Case per run, after
collecting all successful artifact records for that Case. Do not have workers
repeatedly read/replace the same Case file. Publish staged outputs and the staged
manifest as one best-effort Case-level operation; if a commit step fails, restore
the previous files where possible, log the Case failure, and do not continue
publishing that Case's staged outputs. Output writers must therefore support
staging rather than writing directly over published files. A process crash during
multi-file replacement cannot be made fully atomic without a transaction log; this
limitation must be explicit. Ordinary reported write/rename failures must not
leave a known partial update in place.

## 11. Parallelism and scale

- Capture the run context before workers fork; never invent a new operator/time
  independently in every worker.
- Preserve the current public lifter return values. Internal worker results can
  carry both the original result and per-artifact records.
- Ensure records for successful writes survive an error later in the same Case's
  work. Existing combined error boundaries need small targeted changes for this.
- Accumulate all records for a Case, then write that Case's `provenance.json`
  once per run. In `lift_all()`, nested lifters share the same per-Case collector;
  a standalone lifter flushes once for each Case it processes. No repeated or
  concurrent writes to the same Case document.
- Parallelize across Cases only. One worker owns a Case's output staging and
  manifest publication; same-Case lifters run serially within that owner. Separate
  processes/jobs writing the same Case require external serialization. An atomic
  rename alone does not solve concurrent read-modify-write conflicts.
- No extra nested pool is needed merely to stat or serialize small records.
- Stat only declared bindings, not every string in the Cohort or every directory.
- Batch filesystem observations where practical. Share observations for repeated
  references within a stable conversion batch; do not use an indefinite cache.
- Re-observe paths across separate lifter invocations. A path or symlink can change
  between them even when they share a run context.
- Bound caches and pending records by Case or processing batch. Do not require an
  all-Cohort provenance document in memory or on disk.
- Do not read source contents for provenance, compute full-file checksums, issue
  network calls, or spawn one `stat`/`ls` process per path.

The scaling target is work proportional to selected file bindings and written
outputs, with reuse for shared references. Any exact stat cache policy must honor
the stated assumption that input files remain stable during the conversion.

## 12. Frontend and inspection behavior

- Fetch the Case's provenance document once when needed, not one copy per component.
- Look up the component's actual output key. Do not assume a lifter produces only
  one output or that all files under a directory share the same sources.
- Show operator, run start, supplied source paths and available file observations.
- Use the same labels for production and backfill entries. `Run started` refers to
  the latest provenance-generation run, including when that run backfilled an
  existing output; it is not guaranteed to be the artifact's original publish time.
- Label dates explicitly: **Run started**, **Source modified**, and, when known,
  **Local file created**. Do not label file mtime as analysis or retrieval time.
- Source paths are labels. They are not automatically browser-downloadable links.
- If a manifest or entry is absent, show **Provenance unavailable**. Do not infer
  sources from current metadata or another component's entry.
- Serve under the same authorization boundary as the Case. Paths and operator
  names can expose sensitive information; never include credentials, signed URL
  tokens or full private configuration. Decide any path-display policy before
  deployment.
- Existing gOS data contracts remain unchanged. No full provenance object is
  embedded into metadata, Arrow payloads, event arrays or cohort catalogs.

A later CLI can inspect these files without loading the data. A new CLI framework
or frontend implementation is not required to deliver the package-side export.

## 13. Verification plan

### Pure unit tests

1. The same explicit inputs produce the same record/document without consulting
   clock, account, files or working directory.
2. Artifact merges preserve unrelated outputs and their original run contexts.
3. An input replacement affects only the corresponding output entry.
4. Invalid Case identities, artifact keys and unsupported schemas are rejected.
5. Empty input arrays and empty output objects retain their JSON types.
6. Unknown observations remain `null`; role/column bindings are not lost.

### Filesystem tests with temporary synthetic files

7. Regular file: correct supplied/absolute/resolved paths, size and mtime.
8. Symlink: preserve its supplied path and record the target observed at production.
   Retargeting it later does not mutate the stored record.
9. Relative paths and paths with spaces are handled without shell quoting tricks.
10. Missing, broken-link and inaccessible observations do not become fake zero-byte
    files or invented dates. Skip platform-specific permission tests when needed.
11. Unsupported birth time remains `null`; POSIX ctime is never used as creation.
12. Overwriting a file after production does not silently update saved provenance.
13. A failed temporary JSON write/rename does not corrupt the existing document.
14. Unavailable source stat details do not erase other successfully recorded facts.

### Lifter integration and compatibility tests

15. A successful simple output gets the selected source binding and the caller's
    run context; the data payload and public return value stay unchanged.
16. Metadata's annotation preference records the selected MAF rather than also
    claiming an unused fallback VCF. Selected graph and directly discovered
    resource paths are represented using the existing selection behavior.
17. Changing the Cohort after a conversion does not rewrite that conversion's
    captured bindings during export.
18. Valid empty events get a provenance entry.
19. A skipped existing artifact keeps its old operator, date and sources.
20. A failed replacement before touching old data preserves its old association.
21. A failed provenance stage/write prevents replacing published data for that
    Case, logs the Case and reason, and continues processing other Cases.
22. A failed commit restores prior Case outputs and provenance where possible;
    injected failures test every stage/rename and rollback boundary.
23. A successful events write followed by a highlights failure keeps an accurate
    events record and does not claim highlights succeeded.
24. Different outputs for one Case can retain different successful run contexts.
25. One-worker and two-worker runs give equivalent results across Cases; each Case
    has one manifest writer and no two workers update the same Case file.
26. No full provenance payload appears in metadata or the combined catalog.
27. Explicitly disabling provenance preserves legacy output content and returns.
28. Paths outside the Cohort supplied as direct resource arguments are represented;
    upstream manifests are never opened to expand the source list.
29. Instrument I/O in a many-Case synthetic test: no source-content reads or hashes
    for provenance, bounded observation cache, and no per-file subprocesses.
30. Backfill from a saved synthetic `cohort.rds` and existing output directory
    creates manifests without invoking conversion functions or changing payloads.
31. Backfill uses the ordinary schema and current provenance timestamp without a
    backfill marker; it does not invent a historical output timestamp.
32. Backfill considers documented Skilift output keys only; unknown files are not
    labeled as outputs.
33. If a Case's `provenance.json` already exists, default backfill leaves it
    unchanged and logs a skip; `force = TRUE` permits replacing it.
34. Backfill can read legacy Cohort RDS exports with no provenance fields.

Use only invented fixtures. Existing fixture guidance is at
`/Users/user1/projects/skilift/tests/fixtures/README.md`; the synthetic corpus is
at `/Users/user1/projects/skilift/inst/extdata/test_data/synthetic_cohorts`.
Do not copy downloaded clinical data into tests.

The existing environment has not established a runnable full Skilift/testthat
suite. Keep pure/helper tests runnable with base R and jsonlite where possible;
report missing package dependencies separately from test failures. Do not claim
full conversion coverage based only on schema/helper tests.

## 14. Implementation plan after review

1. **Freeze the narrow contract.** Complete.
2. **Add the shared source file.** Complete: shared assembly, observations,
   per-Case staging/publication, public backfill API and focused tests are in place.
3. **Prove publication.** Coverage, graph, hetsnps, filtered-events and metadata
   writers use default-on per-Case collection. Coverage, multi-output publication,
   malformed manifests, rollback, backfill and payload preservation have manual
   smoke checks; testthat remains unavailable in the current environment.
4. **Prove backfill without conversion.** Complete in code and manual smoke checks:
   backfill does not call lifters or alter data files, skips manifests by default,
   and requires `force = TRUE` to replace one.
5. **Prove the hard boundaries.** Metadata and filtered-events bindings include
   their selected files and known resources. Automated tests for fallback selection,
   legacy pairing and empty outputs remain to be run.
6. **Expand deliberately.** In progress: inventory and integrate remaining per-Case
   writers, including dynamic outputs and BAM/index handling. Do not claim all
   outputs are covered until each supported writer has an explicit binding and test.
7. **Connect orchestration.** `lift_all()` shares one run context and regenerates
   the combined metadata catalog after staged metadata is published. Documentation
   still needs to call out incomplete writer coverage and same-Case concurrency.
8. **Verify default-on publication.** Manual smoke checks and parse checks pass.
   Full testthat/package checks, broader compatibility checks and stat-count checks
   remain outstanding because required dependencies are unavailable.

No implementation step requires importing the other branch, creating upstream
lineage tables or changing the biological conversion rules.

## 15. Decisions confirmed for implementation

1. **Default:** production provenance is on by default for integrated lifters and
   `lift_all()`. An explicit disable option may preserve legacy behavior.
2. **Operator:** use the OS account. If unavailable, omit `operator`; do not warn
   or substitute another identity.
3. **Publication failure:** do not replace published data when its Case provenance
   cannot be written. Log the Case and reason, continue with other Cases, and leave
   that Case's previous data/provenance in place.
4. **Path visibility:** full paths and OS account may be shown to all users who can
   access the Case. Never include credentials or signed URL tokens.
5. **Parallelism:** one Case is owned and written once per run; no repeated or
   parallel writes to its `provenance.json`. Parallelize across Cases. Separate
   jobs targeting the same Case must be serialized.
6. **Creation time:** `created_at` may be `null`; do not add a dependency or
   substitute ctime/mtime.
7. **Backfill:** use the ordinary schema with no special label. `operator` and
   `started_at` identify the latest provenance-generation run, not a verified
   original publication event. The saved Cohort supplies the paths.
8. **Output discovery:** backfill only documented Skilift output keys that exist.
9. **Existing manifest:** default backfill preserves an existing `provenance.json`
   and logs a skip. `force = TRUE` explicitly permits replacing it.

The remaining API naming and internal representation choices belong to
implementation; they do not need another design approval.

## 16. Acceptance summary

The feature is ready when a reviewer can open one Case document, select a
published output, and see its operator, run start, immediate selected input paths
and saved filesystem observations—without reading the data or upstream provenance.

A user can also load a legacy `cohort.rds`, point to an existing output directory
and generate provenance manifests without rerunning conversions or changing data
outputs. The same schema and labels apply; `started_at` is the provenance-generation
run, not a verified historical output-publication time.

Production entries remain truthful after a partial rerun, a skipped output, an
input-path change and a failed replacement. Existing gOS payloads, scientific
behavior and public return conventions remain unchanged. Coverage limitations are
explicit, not represented as complete provenance.
