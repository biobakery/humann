# command 1
humann_merge_abundance_tables --input-genes demo_genefamilies.tsv --input-pathways demo_pathabundance.tsv --output merge_output.tsv

# command 2
humann_merge_abundance_tables --input-genes demo_genefamilies.tsv --input-pathways demo_pathabundance.tsv --output merge_output_remove_taxonomy.tsv --remove-taxonomy
