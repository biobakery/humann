HUMAnN2 Utility Scripts


Split a gene table
 
* ./split_gene_table.py -i <genetable.{tsv,biom}> -o <output_dir>
* Input: gene table (tsv or biom format)
* Output: gene tables (one per sample, in biom format if input is biom format)


Join gene tables

* ./join_gene_tables.py -i <input_dir> -o <genetable.{tsv,biom}>
* Input: a directory containing gene tables (tsv or biom format)
* Optional Input: "--file_name $STR" will only join gene tables with $STR in file name
* Output: a single gene table (biom format if input is biom format)
 
