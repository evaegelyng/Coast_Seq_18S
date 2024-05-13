# Counting reads at each step of 18S pipeline
# Run from Autumn/Gua18S/both_seasons

# Count reads before chimera removal
raw<-read.table("results/DADA2_raw.table", sep="\t",row.names = 1,header=TRUE)
sum(colSums(raw))
# 2023131440 = 2,023 M reads

# Count reads after chimera removal
nochim<-read.table("results/DADA2_nochim.table", sep="\t",row.names = 1,header=TRUE)
sum(colSums(nochim))
# 1945062804 = 1,945 M reads

# Check that the file DADA2_nochim_fixed.table has the same no. of reads as DADA2_nochim.table
fixed<-read.table("results/DADA2_nochim.table_fixed", sep="\t",row.names = 1,header=TRUE)
sum(colSums(fixed))
# 1945062804 Same as DADA2_nochim.table

# Count reads after cleaning based on controls
clean<-read.table("results/cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(clean))
# 1934618784 = 1,935 M reads

# Count reads after removal per sample of sequences found in a single PCR replicate
no.sing<-read.table("results/no_sing_cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(no.sing))
# 1812625200 = 1,813 M reads

# Count reads after taxonomic assignment
tax<-read.table("results/cleaned_tax_ncbi_12_12_22.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(tax))
# 762137264 = 762 M reads
