# Coast_Seq_18S
Following the basic processing done for each seasonal dataset, the following scripts were run:

1. make_metadata.R 
2. clean_up_ASV_wise.R
3. no_sing_ASV_wise.R (output: no_sing_cleaned_otu_table_ASV_wise.txt)
4. ...
5. workflow.py (used the script taxonomy_v0_1_1.R) (input: DADA2_nochim_clean.otus)
6. new_clean_tax_ASV_ncbi_15_11_22.R
7. new_final_trim_rarefying_ncbi_12_12_22.R