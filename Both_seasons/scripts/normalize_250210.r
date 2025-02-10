# Script to normalize phyloseq objects from COSQ data
# Should be run from the both_seasons folder, using the metabar_2021 environment

## Load packages
library(phyloseq)
library(tibble)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
library(vegan)

###Load OTU table
otu_tab<-read.table("~/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/COSQ_Curated_LULU_97.tsv", sep="\t", header=T, row.names=1,check.names=F)
###Remove total reads column 
otu_tab<-within(otu_tab,rm(total_reads))

## Count no. of ASVs before removing non-marine classes
nrow(otu_tab) # 10884
## Count no. of reads before removing non-marine classes
sum(colSums(otu_tab)) # 1,492,688,352

#Load taxonomy table
taxonomy_nt<-read.table("ncbi_nt_tax/results/COSQ_classified.tsv", sep='\t', header=T, quote="", fill=T, stringsAsFactors = FALSE)
### Remove read count from qseqid
head(taxonomy_nt$qseqid)
taxonomy_nt$qseqid<-sub("\\;.*","",taxonomy_nt$qseqid)
head(taxonomy_nt$qseqid)
rownames(taxonomy_nt)<-taxonomy_nt$qseqid

## Load table of marine/non-marine classes
env<-read.table("ncbi_nt_tax/results/18S_and_COI_classes_environment.txt",sep="\t", header=T, check.names=F)
## Check if any classes have not been assigned to marine/non-marine
taxonomy_nt_na <- taxonomy_nt[!taxonomy_nt$class %in% env$class,]
unique(taxonomy_nt_na$class) # NA, so all classes have been assigned

## Remove non-marine classes from taxonomy table
tax_mar <- taxonomy_nt[taxonomy_nt$class %in% env$class[env$marine=="yes"],]
###Convert to matrix for phyloseq
tax_mat<-as.matrix(tax_mar)

## Remove non-marine classes from OTU table
otu_mar <- otu_tab[rownames(otu_tab) %in% rownames(tax_mar),]

## Count no. of ASVs after removing non-marine classes
nrow(otu_mar) # 10673
## Count no. of reads after removing non-marine classes
sum(colSums(otu_mar)) # 1,491,982,988

###Convert to matrix for phyloseq
otu_mat<-as.matrix(otu_mar)

###Summarize no. of reads per PCR replicate
reads<-colSums(otu_mat)
mean(reads) # 396909.5
sd(reads) # 214975

#Load metadata
metadata<-read.table("results/metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

## Make phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
p_ncbi = phyloseq(OTU, TAX)

##SITE INFO
c_s<-read.table("results/metadata/cluster_site.txt", sep="\t", header=T)
metadata$Location<-c_s$Site_name[match(metadata$cluster, c_s$cluster)]
metadata$cluster<-as.integer(metadata$cluster)
metadata$cl_se<-as.character(paste(metadata$cluster,metadata$season,sep="_"))
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

COSQ_final = merge_phyloseq(p_ncbi, sampledata)
COSQ_final

##As rarefying removes data, first try normalizing by simply transforming to relative abundances
### First, merge PCR replicates from the same field sample
merged_1 = merge_samples(COSQ_final, "sample_root")
### Transform to rel abund
otu_tab_stand <- decostand(otu_table(merged_1), "total")
### Hellinger transformation (square-root) to reduce influence of a few very abundant species 
otu_tab_trans <- decostand(otu_tab_stand, "hellinger")
### Save table
write.table(otu_tab_trans,"otu_rel_hel_97.txt", sep="\t", row.names=TRUE)

## Rarefy PCR replicates to median depth, keeping replicates with lower depth
### Remove PCR replicates with zero reads
COSQ_final<-prune_samples(sample_sums(COSQ_final)>0, COSQ_final)

### Make a table with a column indicating which PCR replicates have a read depth above the median
readsi<-sample_sums(COSQ_final)
combinedi<-cbind(readsi, sample_data(COSQ_final))
combinedi<-data.frame(combinedi)
threshold<-round(median(combinedi$readsi))
combinedi$q<-combinedi$readsi>threshold

### Make histogram of raw read counts
mui <- ddply(combinedi, .(season, substrate_type), summarise, grp.mean=mean(readsi))
ggplot(combinedi, aes(x=readsi)) +
geom_histogram(aes(fill=q), position="identity", alpha=0.6, binwidth=2500) + geom_density(alpha=0.6) + geom_vline(data=mui, aes(xintercept=grp.mean), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(substrate_type ~ season, ncol=2, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Reads histogram plot", x ="Reads", y = "Count", fill = paste("Total reads > ",threshold,sep="")) + scale_x_continuous(limits=c(0,1000000)) + scale_y_continuous(limits=c(0,50))
ggsave("reads_hist_raw.pdf")

### Transfer the column generated above to the phyloseq object
sample_data(COSQ_final)$over_median<-combinedi$q[match(sample_data(COSQ_final)$sample_ID, combinedi$sample_ID)]

### Extract and then rarefy the PCR replicates with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(COSQ_final, over_median==TRUE), sample.size=as.numeric(threshold), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

# 1129OTUs were removed because they are no longer 
# present in any sample after random subsampling

### Extract the PCR replicates with a read depth at or below the median
below_t<-subset_samples(COSQ_final, over_median==FALSE)

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_merge <-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples.
COSQ_rare = filter_taxa(COSQ_merge, function(x) sum(x) > 0, TRUE)
COSQ_rare

## Rarefy samples to median read depth
### First, merge PCR replicates from the same field sample
merged = merge_samples(COSQ_rare, "sample_root")

## Rebuild sample data, as the merge_samples function only handles merging of the OTU table
d<-data.frame(sample_data(merged)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])

d$po<- sapply(strsplit(as.character(rownames(d)), "2C"), tail, 1)
d$pn<-gsub('\\d','', d$po)
d$pn1<-gsub(".*C(.+).*", "\\1", d$pn)
d$habitat<-ifelse(d$pn1=="EW"|d$pn1=="EB", "eelgrass", ifelse(d$pn1=="RW"|d$pn1=="RB", "rocks", "sand"))
d$substrate_type<-ifelse(grepl("B", d$pn1, fixed=T), "sediment", "water")
d$season<-ifelse(grepl("2C", as.character(rownames(d)), fixed=T), "autumn", "spring")
d$sample_root<-rownames(d)
d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)
d$cluster<-as.integer(d$pn2)

sample_data(merged)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

### Make a table with a column indicating which samples have a read depth above the median
reads<-sample_sums(merged)
combined<-cbind(reads, sample_data(merged))
combined<-data.frame(combined)
thres<-round(median(combined$reads))
combined$q<-combined$reads>thres

### Transfer the column generated above to the phyloseq object
sample_data(merged)$over_median<-combined$q[match(sample_data(merged)$sample_root, combined$sample_root)]

### Extract and then rarefy the samples with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(merged, over_median==TRUE), sample.size=as.numeric(thres), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

#1396OTUs were removed

### Extract the samples with a read depth at or below the median
below_t<-subset_samples(merged, over_median==FALSE)


### Merge the rarefied samples with the low-depth samples
COSQ_merge2<-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples
COSQ_rare2 = filter_taxa(COSQ_merge2, function(x) sum(x) > 0, TRUE)

### Check how many OTUs are left
COSQ_rare2

## Save final files
tax_m<-data.frame(tax_table(COSQ_rare2))
otu_m<-data.frame(otu_table(COSQ_rare2),check.names=F)

write.table(data.frame(sample_data(COSQ_rare2), check.names=F), "metadata/metadata_rarefy_97.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(otu_m, "otu_rarefy_97.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "tax_rarefy_97.txt", sep="\t", quote=FALSE, row.names=TRUE)