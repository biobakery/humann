HUMAnN2 Tools Scripts

This scripts enable simple operations on HUMAnN2 input/output files.

Split a table:
 
* ./split_table.py -i <genetable.{tsv,biom}> -o <output_dir>
* Input: gene/pathway table (tsv or biom format)
* Output: gene/pathway tables (one per sample, in biom format if input is biom format)

Join tables:

* ./join_tables.py -i <input_dir> -o <genetable.{tsv,biom}>
* Input: a directory containing gene/pathway tables (tsv or biom format)
* Optional Input: "--file_name $STR" will only join gene tables with $STR in file name
* Output: a single gene table (biom format if input is biom format)

Rename table features:

* ./rename_table -i <old_table.tsv> -n <name_mapping.tsv> -o <new_table.tsv>
* Input: gene/pathway table (tsv format)
* Output: gene/pathway table with feature IDs mapped to human-readable names

Renormalize table:

* ./renorm_table -i <old_table.tsv> -o <new_table.tsv>
* Input: gene/pathway table (tsv format)
* Output: gene/pathway table normalized to new units (default=copies per million, cpm)

Regroup table features:

* ./regroup_table -i <old_table.tsv> -g <groupings.tsv> -o <new_table.tsv>
* Input: gene/pathway table (tsv format)
* Output: gene/pathway table with regrouped rows
