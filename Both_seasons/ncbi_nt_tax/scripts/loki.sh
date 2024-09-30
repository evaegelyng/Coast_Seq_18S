#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1
#SBATCH -t 4:00:00

Rscript ~/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/scripts/loki.r