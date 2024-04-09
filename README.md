# Coast_Seq_18S
Following the basic processing done for each seasonal dataset, the following scripts were run:

1. make_metadata.R (output: metadata_both.txt)
2. clean_up_ASV_wise.R (output: cleaned_otu_table_ASV_wise.txt)
3. no_sing_ASV_wise.R (output: no_sing_cleaned_otu_table_ASV_wise.txt)


workflow.py (input: DADA2_nochim_clean.otus, output: classified.txt)
new_clean_tax_ASV_ncbi_15_11_22.R (input: classified.txt, output: cleaned_tax_ncbi_12_12_22.txt)
new_final_trim_rarefying_ncbi_12_12_22.R (input: cleaned_tax_ncbi_12_12_22.txt, output: j_readsASVcssh_log.pdf ...)
find_depth_cutoff_ncbi.R (input: cleaned_tax_ncbi.txt, output: summary_depth_cutoff_ncbi.txt)

dist_by_sea_18S.R (input: f_otu_ncbi_5_01_22.txt, dist_by_sea.txt, output: distance_text_ssh_bray.pdf)
draft_keep_media.R (output: keep_below_ord_bc_root_ncbi.pdf)
summary_combined_tax.R (input: f_otu_ncbi.txt ..., output: barplot_stacked_supergroup_comparison2.pdf ...)