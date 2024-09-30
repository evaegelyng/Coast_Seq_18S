#In terminal run:
#conda activate hmsc

library(tidyverse)
library(seqinr)
library("Biostrings")

#load OTU table and fasta file
fastafile <- readDNAStringSet("results/DADA2_nochim_clean.otus")
OTUtable <- read.delim("results/no_sing_cleaned_otu_table_ASV_wise.txt")

#match seq name in OTU file with seq name in original fasta file, get sequence number and rowsums
seq_df <- OTUtable[match(row.names(OTUtable),names(fastafile)),] %>% 
  mutate(names = row.names(OTUtable), size = rowSums(OTUtable)) %>% 
  select(c(names, size))

#Add "size" and read counts to names
seq_name <- paste(seq_df$names, ";size=", seq_df$size, ";", sep = "")

#create a vector of sequences
sequence = as.list(paste(fastafile))

#write the FASTA file
write.fasta(sequences = sequence, names = seq_name, 
            as.string=FALSE, file.out="results/COSQ.fa")