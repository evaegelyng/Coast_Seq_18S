#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 96G
#SBATCH -c 1
#SBATCH -t 12:00:00

# Remember to load conda environment "hmsc"
Rscript ../GitHub_repo/Both_seasons/scripts/new_clean_tax_pident90.R