library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")

#Load tables
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("no_sing_cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results")

taxonomy_nt<-read.table("classified.txt", sep='\t', header=T, quote="", fill=T, stringsAsFactors = FALSE)

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

#Removing taxa that contains only NA from kin to genus
pwer<-subset(n_tax, !(is.na(n_tax$kingdom)&is.na(n_tax$phylum)&is.na(n_tax$class)&is.na(n_tax$order)&is.na(n_tax$family)&is.na(n_tax$genus)))
nrow(pwer)

#Subset based on pident
wer2p<-subset(pwer, pident==100)
nrow(wer2p)
sum(is.na(wer2p$kingdom))
sum(is.na(wer2p$phylum))
sum(is.na(wer2p$class))
sum(is.na(wer2p$order))
sum(is.na(wer2p$family))
sum(is.na(wer2p$genus))
length(levels(factor(wer2p$kingdom)))
length(levels(factor(wer2p$phylum)))
length(levels(factor(wer2p$class)))
length(levels(factor(wer2p$order)))
length(levels(factor(wer2p$family)))
length(levels(factor(wer2p$genus)))

#Subset based on scores
wer3p<-subset(wer2p, kingdom_score>=99&phylum_score>=99&class_score>=99&order_score>=99)
nrow(wer3p)
sum(is.na(wer3p$kingdom))
sum(is.na(wer3p$phylum))
sum(is.na(wer3p$class))
sum(is.na(wer3p$order))
sum(is.na(wer3p$family))
sum(is.na(wer3p$genus))
length(levels(factor(wer3p$kingdom)))
length(levels(factor(wer3p$phylum)))
length(levels(factor(wer3p$class)))
length(levels(factor(wer3p$order)))
length(levels(factor(wer3p$family)))
length(levels(factor(wer3p$genus)))

#Removing "alternatives" column
wer3p<-wer3p[,-22]

#Removing taxa that contains only NA from order to genus (thus, impossible to import order)
wer4p<-subset(wer3p, !(is.na(wer3p$order)&is.na(wer3p$family)&is.na(wer3p$genus)))
nrow(wer4p)
sum(is.na(wer4p$kingdom))
sum(is.na(wer4p$phylum))
sum(is.na(wer4p$class))
sum(is.na(wer4p$order))
sum(is.na(wer4p$family))
sum(is.na(wer4p$genus))
length(levels(factor(wer4p$kingdom)))
length(levels(factor(wer4p$phylum)))
length(levels(factor(wer4p$class)))
length(levels(factor(wer4p$order)))
length(levels(factor(wer4p$family)))
length(levels(factor(wer4p$genus)))


#Fill unassigned labels using taxallnomy
#First, order based on family, trying with taxIDs
wcm<-subset(wer4p, is.na(wer4p$order)&!is.na(wer4p$family)&wer4p$family_score>=99)
cm<-unique(wcm$staxid)
t_t_f= noquote(paste(cm, collapse = ','))

#save list
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
write.table(t_t_f, file="list_order_from_family.txt",row.names=FALSE,quote=FALSE)

#load resulting table
#save list
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
t_now<-data.frame(read.table("list_order_from_family_taxall.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$order))

#filter based on unwanted terms
ordof.a<-subset(t_now, grepl("Ord_of_", t_now[,"order"], fixed=F)==FALSE)
ordof.b<-subset(ordof.a, grepl("incertae sedis", ordof.a[,"order"], fixed=F)==FALSE)
ordof<-subset(ordof.b, grepl("unclassified", ordof.b[,"order"], fixed=F)==FALSE)
list_trueorder<-unique(ordof$order)
list_trueorder
ordof$f_order<-sapply(strsplit(as.character(ordof$order), "_"), tail, 1)

#Fill unassigned order labels using taxallnomy table
##Only if both families match
for(e in 1:nrow(ordof))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"order"])&&!is.na(wer4p[i,"family"])&&wer4p[i,"family_score"]>=99&&wer4p[i,"family"]==ordof[e,"family"])
{
wer4p[i,"order"]<-ordof[e,"f_order"]
}
}
}

sum(is.na(wer4p$order))
length(levels(factor(wer4p$order)))


#Now, repeat it for order based on genus
wcm<-subset(wer4p, is.na(wer4p$order)&!is.na(wer4p$genus)&wer4p$genus_score>=99)
cm<-unique(wcm$staxid)
t_t_f= noquote(paste(cm, collapse = ','))

#save list
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
write.table(t_t_f, file="list_order_from_genus.txt",row.names=FALSE,quote=FALSE)

#load resulting table
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
t_now<-data.frame(read.table("list_order_from_genus_taxall.txt", sep="\t", header=T, stringsAsFactors = FALSE))
levels(factor(t_now$order))

#filter based on unwanted terms
ordof.a<-subset(t_now, grepl("Ord_of_", t_now[,"order"], fixed=F)==FALSE)
ordof.b<-subset(ordof.a, grepl("incertae sedis", ordof.a[,"order"], fixed=F)==FALSE)
ordof<-subset(ordof.b, grepl("unclassified", ordof.b[,"order"], fixed=F)==FALSE)
list_trueorder<-unique(ordof$order)
list_trueorder
ordof$f_order<-sapply(strsplit(as.character(ordof$order), "_"), tail, 1)

#Fill unassigned order labels using taxallnomy table
##Only if both families match
for(e in 1:nrow(ordof))
{
for(i in 1:nrow(wer4p))
{
if(is.na(wer4p[i,"order"])&&!is.na(wer4p[i,"genus"])&&wer4p[i,"genus_score"]>=99&&wer4p[i,"genus"]==ordof[e,"genus"])
{
wer4p[i,"order"]<-ordof[e,"f_order"]
}
}
}

sum(is.na(wer4p$order))
length(levels(factor(wer4p$order)))

#Check taxIDs still missing order label manually!!!!

wwo<-subset(wer4p, is.na(wer4p$order))
sortwwo <- with(wwo, wwo[order(kingdom, phylum, class, order, family, genus),])
sortwwo$kpcofg<-paste(sortwwo$kingdom,sortwwo$phylum,sortwwo$class,sortwwo$order,sortwwo$family,sortwwo$genus)
unique(sortwwo$kpcofg)
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
write.table(unique(sortwwo$kpcofg), file="list_order_to_manual.txt",row.names=FALSE,quote=FALSE)

#Work in an excel sheet, import missing Orders, ?but don't? change lower/higher levels (to avoid creating synonyms and unnecessary future adjustments)

#load resulting table
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
t_nowf<-data.frame(read.table("list_order_manual_edited.txt", sep="\t", header=T, stringsAsFactors = FALSE))

levels(factor(t_nowf$order))
levels(factor(wer4p$order))

t_nowf$sts<-ifelse(is.na(t_nowf$family)&is.na(t_nowf$genus),"n",ifelse(!is.na(t_nowf$family)&is.na(t_nowf$genus),"f",ifelse(is.na(t_nowf$family)&!is.na(t_nowf$genus),"g","b")))

#Fill unassigned order labels using taxallnomy table
##Only if genus/family match
for(e in 1:nrow(t_nowf))
{
for(i in 1:nrow(wer4p))
{
wer4p[i,"order"]<-ifelse(is.na(wer4p[i,"order"])&&
!is.na(wer4p[i,"genus"])&&
t_nowf[e,"sts"]=="g"&&
wer4p[i,"genus_score"]>=99&&
wer4p[i,"genus"]==t_nowf[e,"genus"],t_nowf[e,"order"],
ifelse(is.na(wer4p[i,"order"])&&
!is.na(wer4p[i,"family"])&&
(t_nowf[e,"sts"]=="b"||
t_nowf[e,"sts"]=="f")&&
wer4p[i,"family_score"]>=99&&
wer4p[i,"family"]==t_nowf[e,"family"],
t_nowf[e,"order"],wer4p[i,"order"]))
}
}

sum(is.na(wer4p$order))
length(levels(factor(wer4p$order)))

#Inspect the remainning NA order
wwo<-subset(wer4p, is.na(wer4p$order))
wwo$alt<-wer2p$alternatives[match(rownames(wwo), wer2p$qseqid)]
wwo

#If it's due to low SCOREs (fam & gen), without an obvious fix, then just eliminate it
wer5p<-subset(wer4p, !is.na(wer4p$order))
nrow(wer5p)
sum(is.na(wer5p$order))
length(levels(factor(wer5p$order)))

#Proceed with adding missing ranks..
sortwwo <- with(wer5p, wer5p[order(kingdom, phylum, class, order, family, genus),])
sortwwo$kpcofg<-paste(sortwwo$kingdom,sortwwo$phylum,sortwwo$class,sortwwo$order,sortwwo$family,sortwwo$genus)
unique(sortwwo$kpcofg)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
write.table(unique(sortwwo$kpcofg), file="list_just_to_screen.txt",row.names=FALSE,quote=FALSE)

#Now finalize imputation/ corrections using final edited table
#load resulting table
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
f_tax_t<-data.frame(read.table("list_just_to_screen_edited_081222.txt", sep="\t", header=T, stringsAsFactors = FALSE))

#Create one new column for supergroup and another for marking changes, than incorporate based on unique orders

wer5p$supergroup<-NA
wer5p$altered<-NA

wer5p<-wer5p[,c("margin","qseqid","sseqid","staxid","pident","score","qcovs",          "supergroup","kingdom","phylum","class","order","family","genus","species","altered",       "species_score","genus_score","family_score","order_score","class_score",     "phylum_score","kingdom_score","pident.max.best","possible.misid","ssciname","score.id")]

#Overwrite all labels higher than order using final edited table
##Only if both order matches

f_tax_tuo<-unique(f_tax_t[,c("supergroup","kingdom","phylum","class","order")])
nrow(f_tax_tuo)

for(e in 1:nrow(f_tax_tuo))
{
for(i in 1:nrow(wer5p))
{
if(!is.na(wer5p[i,"order"])&&wer5p[i,"order"]==f_tax_tuo[e,"order"])
{
wer5p[i,8:11]<-f_tax_tuo[e,1:4]
wer5p[i,16]<-T
}
}
}

subset(wer5p, is.na(wer5p$altered))

##For orders not matching, use family/genus
f_tax_tf<-na.omit(f_tax_t[,c("supergroup","kingdom","phylum","class","order","family")])
f_tax_tuof<-unique(f_tax_tf[,c("supergroup","kingdom","phylum","class","order","family")])
nrow(f_tax_tuof)

for(e in 1:nrow(f_tax_tuof))
{
for(i in 1:nrow(wer5p))
{
if(is.na(wer5p[i,"altered"])&&
!is.na(wer5p[i,"family"])&&
wer5p[i,"family"]==f_tax_tuof[e,"family"])
{
wer5p[i,8:12]<-f_tax_tuof[e,1:5]
wer5p[i,16]<-T
}
}
}

subset(wer5p, is.na(wer5p$altered))

f_tax_to<-na.omit(f_tax_t[,c("supergroup","kingdom","phylum","class","order","genus")])
f_tax_tuofg<-unique(f_tax_to[,c("supergroup","kingdom","phylum","class","order","genus")])
nrow(f_tax_tuofg)

for(e in 1:nrow(f_tax_tuofg))
{
for(i in 1:nrow(wer5p))
{
if(is.na(wer5p[i,"altered"])&&
!is.na(wer5p[i,"genus"])&&
wer5p[i,"genus"]==f_tax_tuofg[e,"genus"])
{
wer5p[i,8:12]<-f_tax_tuofg[e,1:5]
wer5p[i,16]<-T
}
}
}

subset(wer5p, is.na(wer5p$altered))

##Inspect table row-by-row 

sortwer5p <- with(wer5p, wer5p[order(supergroup,kingdom, phylum, class, order, family, genus),])
f_wer5p<-unique(sortwer5p[,c("supergroup","kingdom","phylum","class","order","family","genus")])

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/ncbi_nt_tax/results/taxallnomy")
write.table(f_wer5p, file="flist_to_screen.txt",row.names=FALSE,quote=FALSE, sep="\t")

####Final Adjustments:

#Within Apicomplexa -> fix Margolisiella -> Marosporida, Marosporida_IS (c,o) taxID=1118043

wer5p$order[wer5p$staxid==1118043]<-"Marosporida_IS"


#########END OF tHE sCRIPT###########


wer5ps <- wer5p
wer5ps$seq_order<-sapply(strsplit(rownames(wer5ps), "seq"), tail, 1)
wer5ps$seq_order<-as.numeric(wer5ps$seq_order)
wer5ps <- with(wer5ps, wer5ps[order(seq_order),])
tax_mat_n<-as.matrix(wer5ps[,c("supergroup","kingdom","phylum","class","order","family","genus")])
colnames(tax_mat_n)<-c("Supergroup","Division","Phylum","Class","Order","Family","Genus")

cat("\n")
cat("nt_tax_summary")
nrow(tax_mat_n)
cat("\n")
levels(factor(tax_mat_n[,1]))
cat("\n")
head(tax_mat_n)
cat("\n")
summary(tax_mat_n)
cat("\n")

#Filter otu_table according to taxonomy
n_otu_mat<-otu_mat[rownames(otu_mat) %in% rownames(tax_mat_n),]

cat("\n")
cat("matched_otu_table")
nrow(n_otu_mat)

OTU = otu_table(n_otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat_n)

p_DADAwang = phyloseq(OTU, TAX)
p_DADAwang

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_DADAwang, sampledata)
DADAwang1

t100 = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(t100)>0,t100)

tudao
sum(sample_sums(tudao))

tax_mat_b4<-data.frame(tax_table(tudao))
new_otu_mat<-data.frame(otu_table(tudao),check.names=F)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")

write.table(data.frame(sample_data(tudao), check.names=F), "cleaned_ncbi_metadata_12_12_22.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

write.table(new_otu_mat, "cleaned_tax_ncbi_12_12_22.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_mat_b4, "cleaned_ncbi_12_12_22.txt", sep="\t", quote=FALSE, row.names=TRUE)

