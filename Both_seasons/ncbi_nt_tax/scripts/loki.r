
#Run in terminal
#conda activate mjolnir
#cd both_seasons/ncbi_nt_tax/results
#sed "s/;.*;//" COSQ_SWARM_seeds.fasta > COSQ_SWARM_seeds_nosize.fasta

# LOKI: LULU Overseeing with Kinship Identification 

# This function is a convenient wrapper of LULU for the MJOLNIR metabarcoding pipeline.
# It starts from the combined dataset of abundances and taxonomy from the previous step (FRIGGA): LIBR.All_MOTUs.tsv
# Then a match list of representative MOTU sequences is created using VSEARCH and saved as a txt file.
# Then MOTUs which are potential pseudogenes are labelled and removed using LULU.
# The output includes the LIBR.match_list.txt file and 3 CSV files:
# - A final curated metabarcoding dataset: LIBR.Curated_LULU.tsv
# - A dataset of discarded MOTUs: LIBR.Discarded_LULU.tsv
# - A file with informaton on the fate of discarded MOTUs (with IDs of putative mother sequences): LIBR.Deleted_LULU_fate.tsv  


# lib is the name of the library to be processed. Usually a four-character uppercase name.
# Input file name must be in the format: LIBR.All_MOTUs.tsv. Then lib must be = "LIBR"
# min_id is the minimum identity between two sequences to be kept in the match_list output. Default: 0.84 
 
lib <- "COSQ"
#min_id <- 0.84
#min_id <- 0.90
#min_id <- 0.95
min_id <- 0.97

message("LOKI will produce a pairwise match list for LULU.")
#system(paste0("vsearch --usearch_global ",lib,"_SWARM_seeds_nosize.fasta --db ",lib,"_SWARM_seeds_nosize.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
#system(paste0("vsearch --usearch_global ",lib,"_SWARM_seeds_nosize.fasta --db ",lib,"_SWARM_seeds_nosize.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_match_list_90.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
#system(paste0("vsearch --usearch_global ",lib,"_SWARM_seeds_nosize.fasta --db ",lib,"_SWARM_seeds_nosize.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_match_list_95.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
system(paste0("vsearch --usearch_global ",lib,"_SWARM_seeds_nosize.fasta --db ",lib,"_SWARM_seeds_nosize.fasta --self --id ",min_id," --iddef 1 --userout ",lib,"_match_list_97.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"),intern=T,wait=T)
message("LOKI will now remove the pseudogenes with LULU.")
  
suppressPackageStartupMessages(library(lulu))
  
#Load the dataset
db <- read.table("~/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/cleaned_otu_pident90.txt",sep="\t",head=T,stringsAsFactors = F,check.names=F)
  
#Load the matchlist
#matchlist_name <- read.csv(paste0(lib,"_match_list.txt"),sep="\t",head=F,stringsAsFactors = F)
#matchlist_name <- read.csv(paste0(lib,"_match_list_90.txt"),sep="\t",head=F,stringsAsFactors = F)
#matchlist_name <- read.csv(paste0(lib,"_match_list_95.txt"),sep="\t",head=F,stringsAsFactors = F)
matchlist_name <- read.csv(paste0(lib,"_match_list_97.txt"),sep="\t",head=F,stringsAsFactors = F)
  
#Run LULU
curated_result <- lulu(db, matchlist_name)
  
#Get discarded table:
discarded_db <- db[rownames(db) %in% curated_result$discarded_otus,]
#write.table(discarded_db,paste0(lib,"_Discarded_LULU.tsv"),row.names = T,sep="\t",quote = F)
#write.table(discarded_db,paste0(lib,"_Discarded_LULU_90.tsv"),row.names = T,sep="\t",quote = F)
#write.table(discarded_db,paste0(lib,"_Discarded_LULU_95.tsv"),row.names = T,sep="\t",quote = F)
write.table(discarded_db,paste0(lib,"_Discarded_LULU_97.tsv"),row.names = T,sep="\t",quote = F)

#Get curated table:
curated_db <- db[rownames(db) %in% curated_result$curated_otus,]
curated_db <- curated_db[order(rownames(curated_db)),]
curated_db[,] <- curated_result$curated_table
curated_db$total_reads <- rowSums(curated_result$curated_table)
#write.table(curated_db,paste0(lib,"_Curated_LULU.tsv"),row.names = T,sep="\t",quote = F)
#write.table(curated_db,paste0(lib,"_Curated_LULU_90.tsv"),row.names = T,sep="\t",quote = F)
#write.table(curated_db,paste0(lib,"_Curated_LULU_95.tsv"),row.names = T,sep="\t",quote = F)
write.table(curated_db,paste0(lib,"_Curated_LULU_97.tsv"),row.names = T,sep="\t",quote = F)
  
#Get fate of deleted taxa 
deleted_otu_fate <- (curated_result$otu_map[curated_result$otu_map$curated=="merged",])
#write.table(deleted_otu_fate,paste0(lib,"_Deleted_LULU_fate.tsv"),row.names = T,sep="\t",quote = F)
#write.table(deleted_otu_fate,paste0(lib,"_Deleted_LULU_fate_90.tsv"),row.names = T,sep="\t",quote = F)
#write.table(deleted_otu_fate,paste0(lib,"_Deleted_LULU_fate_95.tsv"),row.names = T,sep="\t",quote = F)
write.table(deleted_otu_fate,paste0(lib,"_Deleted_LULU_fate_97.tsv"),row.names = T,sep="\t",quote = F)

#Save all results in RDS file
#saveRDS(curated_result,file=paste0(lib,"_LULU.rds"))
#saveRDS(curated_result,file=paste0(lib,"_LULU_90.rds"))
#saveRDS(curated_result,file=paste0(lib,"_LULU_95.rds"))
saveRDS(curated_result,file=paste0(lib,"_LULU_97.rds"))

message("LOKI is done. He kept ",nrow(curated_db)," MOTUs in the curated database, which He stored in file: ",paste0(lib,"_Curated_LULU.tsv"))
message("LOKI discarded ",nrow(discarded_db)," MOTUs, which He stored in file: ",paste0(lib,"_Discarded_LULU.tsv"))
message("LOKI stored the fate of the discarded MOTUs in file: ",paste0(lib,"_Deleted_LULU_fate.tsv"))

# Check for overmerging of MOTUs
##Load the taxonomy table
tax <- read.table("~/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/cleaned_tax_pident90.txt",sep="\t",head=T,stringsAsFactors = F,check.names=F)
## Make copy of deleted_otu_fate
fate <- deleted_otu_fate
## Add columns with the taxonomy of parent and daughter OTU
fate$class <- tax$class[match(rownames(fate),rownames(tax))]
fate$class_parent <- tax$class[match(fate$parent_id,rownames(tax))]
fate$same <- ifelse(fate$class==fate$class_parent,"yes","no")
## Count number of pairs with same vs different taxonomic assignments
message("Number of parent/daughter pairs with same class assignment")
nrow(fate[fate$same=="yes",]) # 1937 with min_id=0.84
message("Number of parent/daughter pairs with different class assignment")
nrow(fate[fate$same=="no",]) # 189 (9%) with min_id=0.84