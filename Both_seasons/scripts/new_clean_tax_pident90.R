# First, run the following in terminal (interactive job)
#cd Gua18s/both_seasons/
#conda activate hmsc
#R

library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("dplyr")
library("scales")

#Load taxonomy tables
taxonomy_nt<-read.table("ncbi_nt_tax/results/COSQ_classified.tsv", sep='\t', header=T, quote="", fill=T, stringsAsFactors = FALSE)

### Remove read count from qseqid
head(taxonomy_nt$qseqid)
taxonomy_nt$qseqid<-sub("\\;.*","",taxonomy_nt$qseqid)
head(taxonomy_nt$qseqid)

nrow(taxonomy_nt) # 16025
levels(factor(taxonomy_nt[,8]))
head(taxonomy_nt)

#Removing taxa that contain only NA from kingdom to class
pwer<-subset(taxonomy_nt, !(is.na(taxonomy_nt$kingdom)&is.na(taxonomy_nt$phylum)&is.na(taxonomy_nt$class)))
nrow(pwer)
# 14461, so xx ASVs were removed

#Subset based on pident
wer2p<-subset(pwer, pident>=90)
nrow(wer2p) # 14461, so no ASVs removed. 
sum(is.na(wer2p$kingdom)) # 10626
sum(is.na(wer2p$phylum)) # 3271
sum(is.na(wer2p$class)) # 2866

length(levels(factor(wer2p$kingdom))) # 3
length(levels(factor(wer2p$phylum))) # 52
length(levels(factor(wer2p$class))) # 179

#Subset based on scores
wer3p<-subset(wer2p, kingdom_score>=90&phylum_score>=90&class_score>=90)
nrow(wer3p)
# 13526, so xx ASVs removed
sum(is.na(wer3p$kingdom)) # 9921
sum(is.na(wer3p$phylum)) # 3093
sum(is.na(wer3p$class)) # 2638

length(levels(factor(wer3p$kingdom))) # 3
length(levels(factor(wer3p$phylum))) # 52
length(levels(factor(wer3p$class))) # 175

#Removing "alternatives" column
wer3p[1,22]
wer3p<-wer3p[,-22]

#Removing taxa that contains NA at class level (thus, impossible to import kingdom and phylum)
wer4p<-subset(wer3p, !(is.na(wer3p$class)))
nrow(wer4p)
# 10888, so xx ASVs removed
sum(is.na(wer4p$kingdom)) # 7477
sum(is.na(wer4p$phylum)) # 3090
sum(is.na(wer4p$class)) # 0

length(levels(factor(wer4p$kingdom))) # 3
length(levels(factor(wer4p$phylum))) # 36
length(levels(factor(wer4p$class))) # 175

#Fill unassigned labels using taxallnomy
#Phylum based on class, trying with taxIDs
wcm<-subset(wer4p, is.na(wer4p$phylum) )
cm<-unique(wcm$staxid)
length(cm) # 569
t_t_f= noquote(paste(cm, collapse = ','))

#save list
write.table(t_t_f, file="ncbi_nt_tax/results/taxallnomy/list_phylum_from_class_pident90.txt",row.names=FALSE,quote=FALSE)

#Run through Taxallnomy web tool and load resulting table
#save list
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_phylum_from_class_taxall.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$phylum))

#filter based on unwanted terms
phyof<-subset(t_now, grepl("Phy_of_", t_now[,"phylum"], fixed=F)==FALSE)
list_truephylum<-unique(phyof$phylum)
list_truephylum
phyof$f_phylum<-sapply(strsplit(as.character(phyof$phylum), "_"), tail, 1)

#Fill unassigned phylum labels using taxallnomy table
##Only if class names match
for(e in 1:nrow(phyof))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"phylum"])&&wer4p[i,"class"]==phyof[e,"class"])
{
wer4p[i,"phylum"]<-phyof[e,"f_phylum"]
}
}
}

sum(is.na(wer4p$phylum)) # 3028
length(levels(factor(wer4p$phylum))) # 38

#Now, repeat it for kingdom based on class
wcm<-subset(wer4p, is.na(wer4p$kingdom))
cm<-unique(wcm$staxid)
length(cm) # 1727. Too big for Taxallnomy, so dividing into two datasets
cm_1 = cm[1:900]
cm_2 = cm[900:1727]
t_t_f= noquote(paste(cm_1, collapse = ','))
t_t_f_2= noquote(paste(cm_2, collapse = ','))

#save lists
write.table(t_t_f, file="ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_pident90_1.txt",row.names=FALSE,quote=FALSE)
write.table(t_t_f_2, file="ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_pident90_2.txt",row.names=FALSE,quote=FALSE)

#First chunk
#Run through Taxallnomy web tool and load resulting table
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_taxall_1.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$kingdom))

#filter based on unwanted terms
kinof<-subset(t_now, grepl("Kin_of_", t_now[,"kingdom"], fixed=F)==FALSE)
list_truekin<-unique(kinof$kingdom)
list_truekin
kinof$f_kingdom<-sapply(strsplit(as.character(kinof$kingdom), "_"), tail, 1)

#Fill unassigned kingdom labels using taxallnomy table
##Only if class names match
for(e in 1:nrow(kinof))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"kingdom"])&&wer4p[i,"class"]==kinof[e,"class"])
{
wer4p[i,"kingdom"]<-kinof[e,"f_kingdom"]
}
}
}

sum(is.na(wer4p$kingdom)) # 829
length(levels(factor(wer4p$kingdom))) # 7

#Second chunk
#Run through Taxallnomy web tool and load resulting table
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_taxall_2.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$kingdom))

#filter based on unwanted terms
kinof<-subset(t_now, grepl("Kin_of_", t_now[,"kingdom"], fixed=F)==FALSE)
list_truekin<-unique(kinof$kingdom)
list_truekin
kinof$f_kingdom<-sapply(strsplit(as.character(kinof$kingdom), "_"), tail, 1)

#Fill unassigned kingdom labels using taxallnomy table
##Only if class names match
for(e in 1:nrow(kinof))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"kingdom"])&&wer4p[i,"class"]==kinof[e,"class"])
{
wer4p[i,"kingdom"]<-kinof[e,"f_kingdom"]
}
}
}

sum(is.na(wer4p$kingdom)) # 746
length(levels(factor(wer4p$kingdom))) # 7

#Check taxIDs still missing kingdom or phylum label manually!!!!

wwo<-subset(wer4p, is.na(wer4p$kingdom) | is.na(wer4p$phylum))
sortwwo <- with(wwo, wwo[order(kingdom, phylum, class),])
sortwwo$kpc<-paste(sortwwo$kingdom,sortwwo$phylum,sortwwo$class)
unique(sortwwo$kpc)
write.table(unique(sortwwo$kpc), file="ncbi_nt_tax/results/taxallnomy/list_kin_phy_to_manual_pident90.txt",row.names=FALSE,quote=FALSE)

#Work in an excel sheet, import missing labels, but don't change lower/higher levels (to avoid creating synonyms and unnecessary future adjustments)

#load resulting table
t_nowf<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kin_phy_manual_edited.txt", sep="\t", header=T, stringsAsFactors = FALSE))

levels(factor(t_nowf$kingdom))
levels(factor(wer4p$kingdom))

levels(factor(t_nowf$phylum))
levels(factor(wer4p$phylum))

#Fill unassigned kingdom labels using manually edited table
##Only if class labels match
for(e in 1:nrow(t_nowf))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"kingdom"])&&wer4p[i,"class"]==t_nowf[e,"class"])
{
wer4p[i,"kingdom"]<-t_nowf[e,"kingdom"]
}
}
}

sum(is.na(wer4p$kingdom)) # 17
length(levels(factor(wer4p$kingdom))) # 10

#Inspect the remaining NAs
wwo<-subset(wer4p, is.na(wer4p$kingdom))
wwo$alt<-wer2p$alternatives[match(rownames(wwo), wer2p$qseqid)]
wwo

#Fill unassigned phylum labels using manually edited table
##Only if class labels match
for(e in 1:nrow(t_nowf))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"phylum"])&&wer4p[i,"class"]==t_nowf[e,"class"])
{
wer4p[i,"phylum"]<-t_nowf[e,"phylum"]
}
}
}

sum(is.na(wer4p$phylum)) # 735
length(levels(factor(wer4p$phylum))) # 46

## Simplify kingdom names to Plantae, Chromista, Metazoa, Protista, Fungi
unique(wer4p$kingdom)
wer4p$kingdom_5<-"NA"

wer4p<-wer4p %>%
  mutate(kingdom_5 = case_when(kingdom %in% c("Plantae","Viridiplantae") ~ "Plantae",
                             kingdom %in% c("Chromista","Haptista") ~ 'Chromista',
                             kingdom %in% c("Protozoa","Discoba","Sar","Amoebozoa") ~ 'Protozoa',
                             kingdom ==  "Metazoa" ~ 'Metazoa',
                             kingdom ==  "Fungi" ~ 'Fungi',TRUE ~ NA_character_),
         kingdom_5 = factor(kingdom_5, levels = c('Plantae',"Chromista","Protozoa","Metazoa","Fungi")))

unique(wer4p$kingdom_5)

wer5p<-wer4p
wer5p$kingdom<-wer5p$kingdom_5
wer5p<-within(wer5p,rm(kingdom_5))

##Inspect table row-by-row 
sortwer5p <- with(wer5p, wer5p[order(kingdom, phylum, class),])
f_wer5p<-unique(sortwer5p[,c("kingdom","phylum","class")])

write.table(f_wer5p, file="ncbi_nt_tax/results/taxallnomy/flist_to_screen_pident90.txt",row.names=FALSE,quote=FALSE, sep="\t")


#########END OF tHE SCRIPT###########

tax_mat_n <- wer5p
rownames(tax_mat_n)<-tax_mat_n$qseqid
tax_mat_n<-as.matrix(tax_mat_n[,c("kingdom","phylum","class")])

nrow(tax_mat_n) # 10888
levels(factor(tax_mat_n[,1]))
head(tax_mat_n)
summary(tax_mat_n)

#Make phyloseq object
otu_mat<-read.table("ncbi_nt_tax/results/COSQ_SWARM_output_counts.tsv", sep="\t", header=T, row.names=1,check.names=F)
otu_mat<-within(otu_mat,rm(cluster_weight,total_reads))
otu_mat<-as.matrix(otu_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat_n)

p_DADAwang = phyloseq(OTU, TAX)
p_DADAwang

#Load metadata
metadata<-read.table("results/metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_DADAwang, sampledata)
DADAwang1

# Remove any empty taxa or samples
t100 = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(t100)>0,t100)

tudao
sum(sample_sums(tudao))

tax_mat_b4<-data.frame(tax_table(tudao))
new_otu_mat<-data.frame(otu_table(tudao),check.names=F)

write.table(data.frame(sample_data(tudao), check.names=F), "results/metadata/cleaned_metadata_pident90.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(new_otu_mat, "results/cleaned_otu_pident90.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_mat_b4, "results/cleaned_tax_pident90.txt", sep="\t", quote=FALSE, row.names=TRUE)