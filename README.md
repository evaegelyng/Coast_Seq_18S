Title of the study: Will be added upon acceptance of the manuscript, as the repository is public

Short summary of the study: Will be added upon acceptance of the manuscript, as the repository is public

Package versions: cutadapt 4.1, dada2 1.22.0, sickle-trim 1.33, phyloseq 1.38.0, ggplot2 3.3.6, vegan 2.6_2, reshape2 1.4.4, plyr 1.8.7, dplyr 1.0.10, scales 1.2.1, stringr 1.4.1, blast 2.12.0, tidyverse 2.0.0, biostrings 2.66.0, seqinr 4.2_30, taxizedb 0.3.0, RColorBrewer 1.1_3

Overview of folders/files and their contents: The folder "Basic processing" contains the gwf workflow file (workflow.py) and a folder with scripts for the basic processing of raw amplicon sequencing data, from demultiplexing of fastq files to the generation of an ASV table with a corresponding file of ASV sequences. The folder "Final_scripts" contains scripts for further processing, including taxonomic assignment, metadata formatting, ASV filtering based on negative controls and prevalence, and normalization by rarefying.

Instructions: The scripts for basic processing are all run from the workflow.py file. For the further processing, the scripts were run in the following order:

1. workflow.py
2. scripts/make_metadata.R 
3. scripts/clean_up_ASV_wise.R
4. scripts/no_sing_ASV_wise.R
5. scripts/Count_reads.r
6. ncbi_nt_tax/scripts/add_count_fa.r
7. ncbi_nt_tax/workflow.py (calling the scripts odin_240903.r taxonomy_v0_1_1_EES.R)
8. scripts/loki.r
9. scripts/normalize_250210.r