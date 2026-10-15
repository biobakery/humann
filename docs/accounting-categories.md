# Protected accounting: MVP

`humann/config.py` defines the protected identifiers and the two raw-count metadata fields. Core output and utilities inherit those definitions.

New gene-family, reaction, and pathway tables begin with count comments followed by the feature header:

```text
# READS_INITIAL	10000000
# READS_UNMAPPED	2000000
# Gene Family HUMAnN v4.0.0 Adjusted CPMs	sample
UniRef90_example	600000
UniRef90_other	400000
```

`READS_INITIAL` counts the input reads; `READS_UNMAPPED` counts reads remaining after the enabled alignment steps. Both are logged. Counts are measured from the existing read store, not inferred from abundance. Unknown counts are `NA`, including initial counts absent from historical gene tables and both counts for alignment-only BLAST input. Historical gene-family `READS_UNMAPPED` rows provide a fallback raw unmapped count; explicit count comments take precedence. Historical pathway `UNMAPPED` values are scaled abundances and cannot supply raw counts.

The default route is reads → Adjusted RPKs → Adjusted CPMs. The CPM denominator contains only mapped gene-family RPKs. Nonempty mapped community totals sum to one million; species contributions use the same factor. Counts are metadata outside that denominator. Existing RPK and Counts modes remain available. Precomputed gene tables retain their supplied abundance units.

Reactions and pathways are derived from the normalized gene scores. Newly ungroupable genes contribute to the existing `UNGROUPED` abundance category. Genes may contribute to multiple groups, and pathway abundance uses pathway structure, so downstream totals need not sum to one million. No second total-sum normalization is applied.

New pathway outputs no longer compute or emit `UNMAPPED` or `UNINTEGRATED`. Historical `UNMAPPED`, `UNINTEGRATED`, `UNGROUPED`, and `READS_UNMAPPED` rows remain protected, including named and stratified forms. Every profile ingestion containing either `UNMAPPED` or `UNINTEGRATED` emits one nonfatal `DEPRECATION` warning identifying the file and recommending CPM-scale abundances from the main workflow. Warning detection uses the canonical legacy names in `config.py`; zero-valued rows also trigger it. Historical values remain supported. They cannot become biological gene assignments, custom groups, taxonomy assignments, or candidates for statistical association or strain profiling.

Table transformations preserve raw-count comments and protected input rows by default. Regrouping sums ungrouped residuals even with `--function mean`. Renormalization preserves protected rows unchanged and excludes them from biological normalization. RNA/DNA normalization excludes accounting from smoothing and expression calculations while retaining source accounting values. Explicit filtering options such as `--protected N`, `--special n`, or taxonomy removal retain their existing purpose. Sample splitting selects the appropriate count column; joining fills unavailable counts with `NA`. Reducing a single sample preserves its counts; a summary across samples reports `NA` because it has no single underlying read library. BIOM stores the two counts as sample metadata and restores them to comments on TSV conversion.

This MVP introduces no expanded RPK/CPM summaries, grouped/integrated summary fields, header-summary mode, or additional normalization outputs. `humann_renorm_table` remains available for compatibility.
