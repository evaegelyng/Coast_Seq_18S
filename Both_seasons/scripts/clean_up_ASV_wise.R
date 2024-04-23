
library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

#Load tables

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")
###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("DADA2_nochim.table_fixed", sep="\t", header=T, row.names=1,check.names=F))

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

p_DADAwang = phyloseq(OTU)


#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")

p_metadata<-read.table("metadata_both.txt", sep="\t", header=T)
p_metadata2<-p_metadata[,-5]
sampledata = sample_data(data.frame(p_metadata2, row.names=p_metadata2$sample_ID, stringsAsFactors=FALSE))
DADAwang1 = merge_phyloseq(p_DADAwang, sampledata)
DADAwang1 = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)

##Merging samples from different runs
sample_data(DADAwang1)$sr_PCR_r<-paste(sample_data(DADAwang1)$
sample_root, sample_data(DADAwang1)$PCR_replicate, sep="_")
final_table<-merge_samples(DADAwang1, "sr_PCR_r", fun=sum)

############Fix metadata --- (re-build it from 0)

sample_ID<-sample_names(final_table)

###Get sample source
source<-NA
metadata<-data.frame(cbind(sample_ID, source))
rownames(metadata)<-metadata$sample_ID
metadata$p_source<- sapply(strsplit(as.character(metadata$sample_ID), "_"), head, 1)

metadata$source<-ifelse(metadata$p_source=="CNE"|metadata$p_source=="2CCNE", 
as.character("CNE"), ifelse(metadata$p_source=="control", as.character("control"), ifelse(grepl("2SN", metadata$p_source, fixed=T)|grepl("2WN", metadata$p_source, fixed=T), as.character("NTC"), as.character("Field_sample"))))

metadata$root<-gsub("_+[^_]+$", "",metadata$sample_ID)
metadata$po<- sapply(strsplit(as.character(metadata$sample_ID), "_"), tail, 1)
metadata$PCR_replicate<-as.integer(gsub("\\D", "", metadata$po))
metadata$pn<-gsub('\\D','_', metadata$sample_ID)
metadata$pn2<-gsub(".*_(.+)__.*", "\\1", metadata$pn)
metadata$pn3<-gsub(".*[C]([^.]+)[_].*", "\\1", metadata$sample_ID)

metadata$cluster<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn2), ifelse(metadata$source=="control", as.integer(gsub('\\D','', metadata$pn3)), NA))

metadata$pn6<- sapply(strsplit(as.character(metadata$pn), "__"), tail, 1)
metadata$pn7<-sapply(strsplit(as.character(metadata$pn6), "_"), head, 1)
metadata$field_replicate<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn7), NA)
metadata$habitat<-ifelse(metadata$source=="Field_sample", as.character(gsub('\\d','', metadata$pn3)), NA)

extraction_refs<-read.table("extraction_refs_both_seasons.txt", sep="\t", header=T, row.names=1)
extraction_refs$sample_root<-row.names(extraction_refs)
extraction_refs$extraction_refs<-sapply(strsplit(as.character(extraction_refs$extraction_number), "_"), head, 1)

metadata$extraction_refs<-extraction_refs$extraction_refs[match(metadata$root, extraction_refs$sample_root)]

PSU_refs<-read.table("PSU_refs_both.txt", sep="\t", header=T, row.names=1)
PSU_refs$sample_root<-row.names(PSU_refs)
metadata$PSU_refs<-PSU_refs$PSU[match(metadata$root, PSU_refs$sample_root)]
metadata$substrate_type<-PSU_refs$substrate_type[match(metadata$root, PSU_refs$sample_root)]
metadata$season<-PSU_refs$season[match(metadata$root, PSU_refs$sample_root)]

new_metadata<-metadata[,c("sample_ID", "root", "source", "season", "substrate_type", "cluster", "habitat", "field_replicate", "extraction_refs", "PSU_refs", "PCR_replicate")]
colnames(new_metadata)[2]<-"sample_root"

###Setup new phyloseq object
new_sampledata = sample_data(data.frame(new_metadata, row.names=new_metadata$sample_ID, stringsAsFactors=FALSE))
sample_data(final_table)<-new_sampledata
otu_table(final_table)<-t(otu_table(final_table))
sum(sample_sums(DADAwang1))==sum(sample_sums(final_table))

##backup main table
str(otu_mat)

new_otu_mat<-as.matrix(data.frame(otu_table(final_table), check.names=F))
str(new_otu_mat)

utus<-new_otu_mat
sum(colSums(utus))
final_table

###SPRING

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/cleanup_ASV_wise/spring")

final_table_s<-subset_samples(final_table, season=="spring")

#Sed decon by PSU set
cat("spring_sed_NTC")
#Subset substrate type (sediment)
data_sedi<-subset_samples(final_table_s, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate  __!!!__ CNEs were used as NTCs for spring (NTCs missing)
  contam<-subset_samples(data_sedi, source=="CNE"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)  
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="CNE"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs 
  write.table(s_ct[,-1], file=paste("cont_list_sed_ntc_spring",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(new_otu_mat))-sum(colSums(utus))
after_ntc_decon_sed_spring<-utus

#Water decon by PSU set
cat("spring_wat_NTC")
#Subset substrate type (water)
data_wate<-subset_samples(final_table_s, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate  
  contam<-subset_samples(data_wate, source=="CNE"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="CNE"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_ntc_spring",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_ntc_decon_sed_spring))-sum(colSums(utus))
after_ntc_decon_wat_spring<-utus

#Sed decon by extraction set
cat("spring_sed_CNE")
#Subset substrate type (sediment)
####exclude missing CNE refs
otu_table(final_table_s) = otu_table(after_ntc_decon_sed_spring, taxa_are_rows = TRUE)
pdata_sedi<-subset_samples(final_table_s, !extraction_refs=="S2")
data_sedi<-subset_samples(pdata_sedi, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_sedi, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_sed_extr_spring",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ntc_decon_wat_spring))-sum(colSums(utus))
after_ext_decon_sed_spring<-utus

#Wat decon by extraction set
cat("spring_wat_CNE")
#Subset substrate type (water)
otu_table(final_table_s) = otu_table(after_ntc_decon_wat_spring, taxa_are_rows = TRUE)
data_wate<-subset_samples(final_table_s, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_wate, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_extr_spring",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_sed_spring))-sum(colSums(utus))
after_ext_decon_wate_spring<-utus

#Wat decon by field control
cat("spring_wat_field_control")
#Subset substrate type (water)
otu_table(final_table_s) = otu_table(after_ext_decon_wate_spring, taxa_are_rows = TRUE)
pdata_wate<-subset_samples(final_table_s, !cluster==3)
data_wate<-subset_samples(pdata_wate, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
clst<-levels(as.factor(tm[,"cluster"]))

for (ps in 1:length(clst)){
#summarize contaminants per cluster  
  contam<-subset_samples(data_wate, source=="control"&cluster==clst[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(clst[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, source=="Field_sample"&cluster==clst[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(clst[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_clst_spring",clst[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_wate_spring))-sum(colSums(utus))
after_clst_decon_wate_spring<-utus


#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#Season2 - keep fixing same utus object#

###AUTUMN

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/cleanup_ASV_wise/autumn")

final_table_a<-subset_samples(final_table, season=="autumn")

#Sed decon by PSU set
cat("autumn_sed_NTC")
#Subset substrate type (sediment)
data_sedi<-subset_samples(final_table_a, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate
  contam<-subset_samples(data_sedi, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs 
  write.table(s_ct[,-1], file=paste("cont_list_sed_ntc_autumn",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_ntc_decon_sed_spring))-sum(colSums(utus))
after_ntc_decon_sed_autumn<-utus

#Water decon by PSU set
cat("autumn_wat_NTC")
#Subset substrate type (water)
data_wate<-subset_samples(final_table_a, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate  
  contam<-subset_samples(data_wate, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")],   ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_ntc_autumn",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_ntc_decon_sed_autumn))-sum(colSums(utus))
after_ntc_decon_wat_autumn<-utus

#Sed decon by extraction set
cat("autumn_sed_CNE")
#Subset substrate type (sediment)
####exclude missing CNE refs
otu_table(final_table_a) = otu_table(after_ntc_decon_sed_autumn, taxa_are_rows = TRUE)
data_sedi<-subset_samples(final_table_a, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_sedi, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_sed_extr_autumn",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ntc_decon_wat_autumn))-sum(colSums(utus))
after_ext_decon_sed_autumn<-utus

#Wat decon by extraction set
cat("autumn_wat_CNE")
#Subset substrate type (water)
otu_table(final_table_a) = otu_table(after_ntc_decon_wat_autumn, taxa_are_rows = TRUE)
data_wate<-subset_samples(final_table_a, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_wate, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_extr_autumn",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_sed_autumn))-sum(colSums(utus))
after_ext_decon_wate_autumn<-utus

#Wat decon by field control
cat("autumn_wat_field_control")
#Subset substrate type (water)
otu_table(final_table_a) = otu_table(after_ext_decon_wate_autumn, taxa_are_rows = TRUE)
data_wate<-subset_samples(final_table_a, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
clst<-levels(as.factor(tm[,"cluster"]))

for (ps in 1:length(clst)){
#summarize contaminants per cluster  
  contam<-subset_samples(data_wate, source=="control"&cluster==clst[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(clst[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, sum)
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, source=="Field_sample"&cluster==clst[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, sum)
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(clst[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cont_list_wat_clst_autumn",clst[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_wate_autumn))-sum(colSums(utus))
after_clst_decon_wate_autumn<-utus

###Summarize and write new ASV table

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/cleanup_ASV_wise/summary")

new_otu_s<-subset_samples(final_table, season=="spring")
new_otu_a<-subset_samples(final_table, season=="autumn")

new_m_otu<-as.matrix(data.frame(otu_table(final_table), check.names=F))
new_m_otu_s<-as.matrix(data.frame(otu_table(new_otu_s), check.names=F))
new_m_otu_a<-as.matrix(data.frame(otu_table(new_otu_a), check.names=F))

andss<-after_ntc_decon_sed_spring[,colnames(new_m_otu_s)]
andws<-after_ntc_decon_wat_spring[,colnames(new_m_otu_s)]
aedss<-after_ext_decon_sed_spring[,colnames(new_m_otu_s)]
aedws<-after_ext_decon_wate_spring[,colnames(new_m_otu_s)]
acdws<-after_clst_decon_wate_spring[,colnames(new_m_otu_s)]

andsa<-after_ntc_decon_sed_autumn[,colnames(new_m_otu_a)]
andwa<-after_ntc_decon_wat_autumn[,colnames(new_m_otu_a)]
aedsa<-after_ext_decon_sed_autumn[,colnames(new_m_otu_a)]
aedwa<-after_ext_decon_wate_autumn[,colnames(new_m_otu_a)]
acdwa<-after_clst_decon_wate_autumn[,colnames(new_m_otu_a)]

summary_clean_up_spring<-cbind(colSums(new_m_otu_s),colSums(andss),colSums(andws),colSums(aedss),colSums(aedws),colSums(acdws))
colnames(summary_clean_up_spring)<-c("raw","ntc_sed","ntc_wat","cne_sed","cne_wat","cntrl_wat")
write.table(summary_clean_up_spring,"summary_clean_up_spring.txt")

summary_clean_up_autumn<-cbind(colSums(new_m_otu_a),colSums(andsa),colSums(andwa),colSums(aedsa),colSums(aedwa),colSums(acdwa))
colnames(summary_clean_up_autumn)<-c("raw","ntc_sed","ntc_wat","cne_sed","cne_wat","cntrl_wat")
write.table(summary_clean_up_autumn,"summary_clean_up_autumn.txt")

summary_clean_up_both_seasons<-cbind(colSums(new_m_otu),colSums(after_ntc_decon_sed_spring),colSums(after_ntc_decon_wat_spring),colSums(after_ext_decon_sed_spring),colSums(after_ext_decon_wate_spring),colSums(after_clst_decon_wate_spring),colSums(after_ntc_decon_sed_autumn),colSums(after_ntc_decon_wat_autumn),colSums(after_ext_decon_sed_autumn),colSums(after_ext_decon_wate_autumn),colSums(after_clst_decon_wate_autumn))
colnames(summary_clean_up_both_seasons)<-c("raw","ntc_sed_s","ntc_wat_s","cne_sed_s","cne_wat_s","cntrl_wat_s","ntc_sed_a","ntc_wat_a","cne_sed_a","cne_wat_a","cntrl_wat_a")
write.table(summary_clean_up_both_seasons,"summary_clean_up_both_seasons.txt")

###Write cleaned otu table

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")
write.table(utus, "cleaned_otu_table_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)

###Write cleaned metadata

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
write.table(data.frame(sample_data(final_table), check.names=F), "merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)

###Test cleaned table

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")
otu_matt<-as.matrix(read.table("cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

OTUt = otu_table(otu_matt, taxa_are_rows = TRUE)
p_DADAwangt = phyloseq(OTUt)

p_DADAwang
cat("total_reads_before_cleanup")
sum(sample_sums(p_DADAwang))
cat("\n")
p_DADAwangt
cat("total_reads_after_cleanup")
sum(sample_sums(p_DADAwangt))
