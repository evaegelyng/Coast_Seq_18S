library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")

# Run from Gua18S/both_seasons

#Load tables

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("results/no_sing_cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_nt<-read.table("ncbi_nt_tax/results/classified.txt", sep='\t', header=T, quote="", fill=T, stringsAsFactors = FALSE)

cat("\n")
cat("raw_nt_tax_summary")
nrow(taxonomy_nt)
cat("\n")
levels(factor(taxonomy_nt[,8]))
cat("\n")
head(taxonomy_nt)
cat("\n")

#Select sequences according to cleaned ASV table
n_tax<- taxonomy_nt[taxonomy_nt$qseqid %in% rownames(otu_mat),]
rownames(n_tax)<-n_tax$qseqid
nrow(n_tax)
# 39610

#Removing taxa that contain only NA from kingdom to class
pwer<-subset(n_tax, !(is.na(n_tax$kingdom)&is.na(n_tax$phylum)&is.na(n_tax$class)))
nrow(pwer)
# 35686, so 3924 ASVs were removed

#Subset based on pident
wer2p<-subset(pwer, pident>=90)
nrow(wer2p) 
# 35686, so no ASVs removed. 
#sum(is.na(wer2p$kingdom))
#sum(is.na(wer2p$phylum))
#sum(is.na(wer2p$class))

#length(levels(factor(wer2p$kingdom)))
#length(levels(factor(wer2p$phylum)))
#length(levels(factor(wer2p$class)))


#Subset based on scores
wer3p<-subset(wer2p, kingdom_score>=90&phylum_score>=90&class_score>=90)
nrow(wer3p)
# 33194, so 2492 ASVs removed
sum(is.na(wer3p$kingdom))
sum(is.na(wer3p$phylum))
sum(is.na(wer3p$class))

length(levels(factor(wer3p$kingdom)))
length(levels(factor(wer3p$phylum)))
length(levels(factor(wer3p$class)))

#Removing "alternatives" column
wer3p<-wer3p[,-22]

#Removing taxa that contains NA at class level (thus, impossible to import kingdom and phylum)
wer4p<-subset(wer3p, !(is.na(wer3p$class)))
nrow(wer4p)
# 27355, so 5839 ASVs removed
sum(is.na(wer4p$kingdom))
sum(is.na(wer4p$phylum))
sum(is.na(wer4p$class))

length(levels(factor(wer4p$kingdom)))
length(levels(factor(wer4p$phylum)))
length(levels(factor(wer4p$class)))


#Fill unassigned labels using taxallnomy
#Phylum based on class, trying with taxIDs
wcm<-subset(wer4p, is.na(wer4p$phylum) )
cm<-unique(wcm$staxid)
length(cm) # 716
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
##Only if both families match
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

sum(is.na(wer4p$phylum))
length(levels(factor(wer4p$phylum)))

#Now, repeat it for kingdom based on class
wcm<-subset(wer4p, is.na(wer4p$kingdom))
cm<-unique(wcm$staxid)
length(cm) # 2221. Too big for Taxallnomy, so dividing into three datasets
cm_1 = cm[1:750]
cm_2 = cm[751:1500]
cm_3 = cm[1501:2221]
t_t_f= noquote(paste(cm_1, collapse = ','))
t_t_f_2= noquote(paste(cm_2, collapse = ','))
t_t_f_3= noquote(paste(cm_3, collapse = ','))

#save lists
write.table(t_t_f, file="ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_pident90_1.txt",row.names=FALSE,quote=FALSE)
write.table(t_t_f_2, file="ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_pident90_2.txt",row.names=FALSE,quote=FALSE)
write.table(t_t_f_3, file="ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_pident90_3.txt",row.names=FALSE,quote=FALSE)

#First chunk
#Run through Taxallnomy web tool and load resulting table
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_taxall_1.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$kingdom))

#filter based on unwanted terms
kinof<-subset(t_now, grepl("Kin_of_", t_now[,"kingdom"], fixed=F)==FALSE)
list_truekin<-unique(kinof$kingdom)
list_truekin
kinof$f_kingdom<-sapply(strsplit(as.character(kinof$kingdom), "_"), tail, 1)

#Fill unassigned order labels using taxallnomy table
##Only if both families match
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

sum(is.na(wer4p$kingdom))
length(levels(factor(wer4p$kingdom)))

#Second chunk
#Run through Taxallnomy web tool and load resulting table
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_taxall_2.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$kingdom))

#filter based on unwanted terms
kinof<-subset(t_now, grepl("Kin_of_", t_now[,"kingdom"], fixed=F)==FALSE)
list_truekin<-unique(kinof$kingdom)
list_truekin
kinof$f_kingdom<-sapply(strsplit(as.character(kinof$kingdom), "_"), tail, 1)

#Fill unassigned order labels using taxallnomy table
##Only if both families match
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

sum(is.na(wer4p$kingdom))
length(levels(factor(wer4p$kingdom)))

# Chunk 3
#Run through Taxallnomy web tool and load resulting table
t_now<-data.frame(read.table("ncbi_nt_tax/results/taxallnomy/list_kingdom_from_class_taxall_3.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$kingdom))

#filter based on unwanted terms
kinof<-subset(t_now, grepl("Kin_of_", t_now[,"kingdom"], fixed=F)==FALSE)
list_truekin<-unique(kinof$kingdom)
list_truekin
kinof$f_kingdom<-sapply(strsplit(as.character(kinof$kingdom), "_"), tail, 1)

#Fill unassigned order labels using taxallnomy table
##Only if both families match
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

sum(is.na(wer4p$kingdom))
length(levels(factor(wer4p$kingdom)))

#Check taxIDs still missing order label manually!!!!

wwo<-subset(wer4p, is.na(wer4p$kingdom) | is.na(wer4p$phylum))
sortwwo <- with(wwo, wwo[order(kingdom, phylum, class, order, family, genus),])
sortwwo$kpc<-paste(sortwwo$kingdom,sortwwo$phylum,sortwwo$class)
unique(sortwwo$kpc)
write.table(unique(sortwwo$kpc), file="ncbi_nt_tax/results/taxallnomy/list_kin_phy_to_manual_pident90.txt",row.names=FALSE,quote=FALSE)
