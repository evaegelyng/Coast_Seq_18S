library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

#correction
#import data with check.names=F
otu_mat<-as.matrix(read.table("DADA2_nochim.table_fixed", sep="\t", header=T, row.names=1,check.names=F))

###Make phyloseq object from raw data
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
p_DADAwang = phyloseq(OTU)

#From sample_names generate metadata
sample_ID<-sample_names(p_DADAwang)

###Get sample source
source<-NA
metadata<-data.frame(cbind(sample_ID, source))
rownames(metadata)<-metadata$sample_ID
metadata$p_source<- sapply(strsplit(as.character(metadata$sample_ID), "_"), head, 1)
levels(factor(metadata$p_source))

metadata$source<-ifelse(metadata$p_source=="CNE"|metadata$p_source=="2CCNE", 
as.character("CNE"), ifelse(metadata$p_source=="control", as.character("control"), ifelse(grepl("2SN", metadata$p_source, fixed=T)|grepl("2WN", metadata$p_source, fixed=T), as.character("NTC"), as.character("Field_sample"))))

levels(factor(metadata$source))

#Get root
head(metadata)

metadata$root<-gsub("_+[^_]+$", "",metadata$sample_ID)

levels(factor(metadata$root))

###Get sample PCR replicate and seq run
metadata$po<- sapply(strsplit(as.character(metadata$sample_ID), "_"), tail, 1)
levels(factor(metadata$po))

metadata$PCR_replicate<-as.integer(gsub("\\D", "", metadata$po))

metadata$seq_run<-ifelse(grepl("^[0-9]+$", metadata$po, perl = T), as.integer(1), as.integer(2))

head(metadata)

###Get sample cluster
metadata$pn<-gsub('\\D','_', metadata$sample_ID)
metadata$pn2<-gsub(".*_(.+)__.*", "\\1", metadata$pn)
metadata[,c("root","pn2")]
levels(factor(metadata$pn2))


metadata$pn3<-gsub(".*[C]([^.]+)[_].*", "\\1", metadata$sample_ID)
metadata[,c("root","pn3")]

metadata$cluster<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn2), ifelse(metadata$source=="control", as.integer(gsub('\\D','', metadata$pn3)), NA))

metadata[,c("root","cluster")]
levels(factor(metadata$cluster))
str(metadata)

###Get field replicate
metadata$pn6<- sapply(strsplit(as.character(metadata$pn), "__"), tail, 1)
metadata$pn7<-sapply(strsplit(as.character(metadata$pn6), "_"), head, 1)

metadata$field_replicate<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn7), NA)

metadata[,c("root","field_replicate")]
levels(factor(metadata$field_replicate))
str(metadata)

###Get habitat
metadata$habitat<-ifelse(metadata$source=="Field_sample", as.character(gsub('\\d','', metadata$pn3)), NA)

metadata[,c("root","habitat")]
levels(factor(metadata$habitat))
str(metadata)

###Get extraction refs
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
extraction_refs<-read.table("extraction_refs_both_seasons.txt", sep="\t", header=T, row.names=1)
extraction_refs$sample_root<-row.names(extraction_refs)
extraction_refs$extraction_refs<-sapply(strsplit(as.character(extraction_refs$extraction_number), "_"), head, 1)

metadata$extraction_refs<-extraction_refs$extraction_refs[match(metadata$root, extraction_refs$sample_root)]

metadata[,c("root","extraction_refs")]
levels(factor(metadata$extraction_refs))
head(metadata)

###Get PSU refs
PSU_refs<-read.table("PSU_refs_both.txt", sep="\t", header=T, row.names=1)
PSU_refs$sample_root<-row.names(PSU_refs)
head(PSU_refs)

metadata$PSU_refs<-PSU_refs$PSU[match(metadata$root, PSU_refs$sample_root)]

metadata[,c("root","PSU_refs")]
levels(factor(metadata$PSU_refs))
head(metadata)

###Get substrate_type
metadata$substrate_type<-PSU_refs$substrate_type[match(metadata$root, PSU_refs$sample_root)]

metadata[,c("root","substrate_type")]
levels(factor(metadata$substrate_type))
str(metadata)

###Get season
metadata$season<-PSU_refs$season[match(metadata$root, PSU_refs$sample_root)]

metadata[,c("root","season")]
levels(factor(metadata$season))
head(metadata)

###Discard tmp variables and inspect metadata sheet
colnames(metadata)
test_md<-metadata[,c("sample_ID", "root", "source", "season", "seq_run", "substrate_type", "cluster", "habitat", "field_replicate", "extraction_refs", "PSU_refs", "PCR_replicate")]
head(test_md)
colnames(test_md)[2]<-"sample_root"


write.table(test_md, "metadata_both.txt", sep="\t", row.names = T)


