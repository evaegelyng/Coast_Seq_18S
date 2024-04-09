library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
sessionInfo()

#Load tables
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

###Make phyloseq object from raw data
otu_mat_ncbi<-as.matrix(read.table("f_otu_ncbi.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("f_tax_ncbi.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat_ncbi, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_ncbi_metadata.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")

metadata_ncbi<-metadata
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

NCBIZAO = merge_phyloseq(p_SILVA, sampledata)
NCBIZAO

#Load tables
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

###Make phyloseq object from raw data
otu_mat_silva<-as.matrix(read.table("f_otu_silva.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_silva<-read.table("f_tax_silva.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_silva)

OTU = otu_table(otu_mat_silva, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_silva_metadata.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")

metadata_silva<-metadata
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

SILVAZAO = merge_phyloseq(p_SILVA, sampledata)
SILVAZAO

#Load tables
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

###Make phyloseq object from raw data
otu_mat_pr2<-as.matrix(read.table("f_otu_pr2.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_pr2<-read.table("f_tax_pr2.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_pr2)

OTU = otu_table(otu_mat_pr2, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_pr2_metadata.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")

metadata_pr2<-metadata
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

PR2ZAO = merge_phyloseq(p_SILVA, sampledata)
PR2ZAO

#####Plotting supergroup community composition
###Individual plots
#NCBI
#Rel abund
sample_data(NCBIZAO)$all<-"all"
NCBImgd<-merge_samples(NCBIZAO, "all")
datag = tax_glom(NCBImgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabNCBI<-tab[,-2:-7]
colnames(tabNCBI)[2]<-"ncbi"
#Rich
taxg<-data.frame(tax_table(NCBImgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(NCBImgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richNCBI<-z[,c("Supergroup","RR")]
colnames(richNCBI)[2]<-"ncbi"
total_reads_NCBI<-sum(sample_sums(NCBImgd))
total_taxa_NCBI<-ntaxa(NCBImgd)

#SILVA
sample_data(SILVAZAO)$all<-"all"
SILVAmgd<-merge_samples(SILVAZAO, "all")
datag = tax_glom(SILVAmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabSILVA<-tab[,-2:-7]
colnames(tabSILVA)[2]<-"silva"
#Rich
taxg<-data.frame(tax_table(SILVAmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(SILVAmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richSILVA<-z[,c("Supergroup","RR")]
colnames(richSILVA)[2]<-"silva"
total_reads_SILVA<-sum(sample_sums(SILVAmgd))
total_taxa_SILVA<-ntaxa(SILVAmgd)

#PR2
sample_data(PR2ZAO)$all<-"all"
pr2mgd<-merge_samples(PR2ZAO, "all")
datag = tax_glom(pr2mgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabpr2<-tab[,-2:-6]
colnames(tabpr2)[2]<-"pr2"
#Rich
taxg<-data.frame(tax_table(pr2mgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(pr2mgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-6]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richpr2<-z[,c("Supergroup","RR")]
colnames(richpr2)[2]<-"pr2"
total_reads_pr2<-sum(sample_sums(pr2mgd))
total_taxa_pr2<-ntaxa(pr2mgd)

#merging abund
erty<-merge(tabNCBI, tabpr2, by="Supergroup", all=T)
poiu<-merge(erty, tabSILVA, by="Supergroup", all=T)
poiu
for(i in 1:nrow(poiu)){
for(e in 1:ncol(poiu)){
if(is.na(poiu[i,e])){poiu[i,e]<-0}
}}

#merging rich
ertyr<-merge(richNCBI, richpr2, by="Supergroup", all=T)
poiur<-merge(ertyr, richSILVA, by="Supergroup", all=T)
poiur
for(i in 1:nrow(poiur)){
for(e in 1:ncol(poiur)){
if(is.na(poiur[i,e])){poiur[i,e]<-0}
}}

###Exclusive taxa plots
otu_n<-t(data.frame(otu_mat_ncbi,check.names=F))
otu_s<-t(data.frame(otu_mat_silva,check.names=F))
otu_p<-t(data.frame(otu_mat_pr2,check.names=F))
#NCBI
n_w_s<-subset(otu_n, !(rownames(otu_n) %in% rownames(otu_s)))
n_w_s_p<-subset(n_w_s, !(rownames(n_w_s) %in% rownames(otu_p)))
nrow(n_w_s_p)
n_tax_nwsp<- taxonomy_ncbi[rownames(taxonomy_ncbi) %in% rownames(n_w_s_p),]
tax_nwsp<-as.matrix(n_tax_nwsp)
otu_nwsp<-as.matrix(n_w_s_p)
OTU = otu_table(otu_nwsp, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_nwsp)
p_nwsp = phyloseq(OTU, TAX_b)
sampledata = sample_data(data.frame(metadata_ncbi, row.names=metadata_ncbi$sample_root, stringsAsFactors=FALSE))
p_nwsp2 = merge_phyloseq(p_nwsp, sampledata)
sample_data(p_nwsp2)$all<-"all"
nwspmgd<-merge_samples(p_nwsp2, "all")
datag = tax_glom(nwspmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabnwsp<-tab[,-2:-7]
colnames(tabnwsp)[2]<-"nwsp"
#Rich
taxg<-data.frame(tax_table(nwspmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(nwspmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richnwsp<-z[,c("Supergroup","RR")]
colnames(richnwsp)[2]<-"nwsp"
total_reads_nwsp<-sum(sample_sums(nwspmgd))
total_taxa_nwsp<-ntaxa(nwspmgd)

#SILVA
s_w_n<-subset(otu_s, !(rownames(otu_s) %in% rownames(otu_n)))
s_w_n_p<-subset(s_w_n, !(rownames(s_w_n) %in% rownames(otu_p)))
nrow(s_w_n_p)
n_tax_swnp<- taxonomy_silva[rownames(taxonomy_silva) %in% rownames(s_w_n_p),]
tax_swnp<-as.matrix(n_tax_swnp)
otu_swnp<-as.matrix(s_w_n_p)
OTU = otu_table(otu_swnp, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_swnp)
p_swnp = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_silva, row.names=metadata_silva$sample_root, stringsAsFactors=FALSE))
p_swnp2 = merge_phyloseq(p_swnp, sampledata)
sample_data(p_swnp2)$all<-"all"
swnpmgd<-merge_samples(p_swnp2, "all")
datag = tax_glom(swnpmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabswnp<-tab[,-2:-7]
colnames(tabswnp)[2]<-"swnp"
#Rich
taxg<-data.frame(tax_table(swnpmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(swnpmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richswnp<-z[,c("Supergroup","RR")]
colnames(richswnp)[2]<-"swnp"
total_reads_swnp<-sum(sample_sums(swnpmgd))
total_taxa_swnp<-ntaxa(swnpmgd)

#PR2
p_w_s<-subset(otu_p, !(rownames(otu_p) %in% rownames(otu_s)))
p_w_s_p<-subset(p_w_s, !(rownames(p_w_s) %in% rownames(otu_n)))
nrow(p_w_s_p)
n_tax_pwsn<- taxonomy_pr2[rownames(taxonomy_pr2) %in% rownames(p_w_s_p),]
tax_pwsn<-as.matrix(n_tax_pwsn)
otu_pwsn<-as.matrix(p_w_s_p)
OTU = otu_table(otu_pwsn, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_pwsn)
p_pwsn = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_pr2, row.names=metadata_pr2$sample_root, stringsAsFactors=FALSE))
p_pwsn2 = merge_phyloseq(p_pwsn, sampledata)
sample_data(p_pwsn2)$all<-"all"
pwsnmgd<-merge_samples(p_pwsn2, "all")
datag = tax_glom(pwsnmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabpwsn<-tab[,-2:-6]
colnames(tabpwsn)[2]<-"pwsn"
#Rich
taxg<-data.frame(tax_table(pwsnmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(pwsnmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-6]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richpwsn<-z[,c("Supergroup","RR")]
colnames(richpwsn)[2]<-"pwsn"
total_reads_pwsn<-sum(sample_sums(pwsnmgd))
total_taxa_pwsn<-ntaxa(pwsnmgd)

#merging abund
cvbn<-merge(tabnwsp, tabswnp, by="Supergroup", all=T)
jhgf<-merge(cvbn, tabpwsn, by="Supergroup", all=T)
jhgf
for(i in 1:nrow(jhgf)){
for(e in 1:ncol(jhgf)){
if(is.na(jhgf[i,e])){jhgf[i,e]<-0}
}}
#merging rich
cvbnr<-merge(richnwsp, richswnp, by="Supergroup", all=T)
jhgfr<-merge(cvbnr, richpwsn, by="Supergroup", all=T)
jhgfr
for(i in 1:nrow(jhgfr)){
for(e in 1:ncol(jhgfr)){
if(is.na(jhgfr[i,e])){jhgfr[i,e]<-0}
}}


###Overlap all plot
ns_tax<- taxonomy_ncbi[rownames(taxonomy_ncbi) %in% rownames(taxonomy_silva),]
nsp_tax<- ns_tax[rownames(ns_tax) %in% rownames(taxonomy_pr2),]
nrow(nsp_tax)
n_otu_nsp<- otu_n[rownames(otu_n) %in% rownames(nsp_tax),]
otu_nsp<-as.matrix(n_otu_nsp)
tax_nsp<-as.matrix(nsp_tax)
OTU = otu_table(otu_nsp, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_nsp)
p_nsp = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_ncbi, row.names=metadata_ncbi$sample_root, stringsAsFactors=FALSE))
p_nsp2 = merge_phyloseq(p_nsp, sampledata)
sample_data(p_nsp2)$all<-"all"
nspmgd<-merge_samples(p_nsp2, "all")
datag = tax_glom(nspmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabnsp<-tab[,-2:-7]
colnames(tabnsp)[2]<-"nsp"
#Rich
taxg<-data.frame(tax_table(nspmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(nspmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richnsp<-z[,c("Supergroup","RR")]
colnames(richnsp)[2]<-"nsp"
total_reads_nsp<-sum(sample_sums(nspmgd))
total_taxa_nsp<-ntaxa(nspmgd)


###Overlap pair plots
#NCBI_SILVA
ns_tax<- taxonomy_ncbi[rownames(taxonomy_ncbi) %in% rownames(taxonomy_silva),]
nrow(ns_tax)
n_otu_ns<- otu_n[rownames(otu_n) %in% rownames(ns_tax),]
otu_ns<-as.matrix(n_otu_ns)
tax_ns<-as.matrix(ns_tax)
OTU = otu_table(otu_ns, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_ns)
p_ns = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_ncbi, row.names=metadata_ncbi$sample_root, stringsAsFactors=FALSE))
p_ns2 = merge_phyloseq(p_ns, sampledata)
sample_data(p_ns2)$all<-"all"
nsmgd<-merge_samples(p_ns2, "all")
datag = tax_glom(nsmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabns<-tab[,-2:-7]
colnames(tabns)[2]<-"ns"
#Rich
taxg<-data.frame(tax_table(nsmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(nsmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richns<-z[,c("Supergroup","RR")]
colnames(richns)[2]<-"ns"
total_reads_ns<-sum(sample_sums(nsmgd))
total_taxa_ns<-ntaxa(nsmgd)


#NCBI_PR2
np_tax<-taxonomy_ncbi[rownames(taxonomy_ncbi) %in% rownames(taxonomy_pr2),]
nrow(np_tax)
n_otu_np<-otu_n[rownames(otu_n) %in% rownames(np_tax),]
otu_np<-as.matrix(n_otu_np)
tax_np<-as.matrix(np_tax)
OTU = otu_table(otu_np, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_np)
p_np = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_ncbi, row.names=metadata_ncbi$sample_root, stringsAsFactors=FALSE))
p_np2 = merge_phyloseq(p_np, sampledata)
sample_data(p_np2)$all<-"all"
npmgd<-merge_samples(p_np2, "all")
datag = tax_glom(npmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabnp<-tab[,-2:-7]
colnames(tabnp)[2]<-"np"
#Rich
taxg<-data.frame(tax_table(npmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(npmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richnp<-z[,c("Supergroup","RR")]
colnames(richnp)[2]<-"np"
total_reads_np<-sum(sample_sums(npmgd))
total_taxa_np<-ntaxa(npmgd)


#SILVA_PR2
sp_tax<- taxonomy_silva[rownames(taxonomy_silva) %in% rownames(taxonomy_pr2),]
nrow(sp_tax)
n_otu_sp<- otu_n[rownames(otu_n) %in% rownames(sp_tax),]
otu_sp<-as.matrix(n_otu_sp)
tax_sp<-as.matrix(sp_tax)
OTU = otu_table(otu_sp, taxa_are_rows = TRUE)
TAX_s = tax_table(tax_sp)
p_sp = phyloseq(OTU, TAX_s)
sampledata = sample_data(data.frame(metadata_silva, row.names=metadata_silva$sample_root, stringsAsFactors=FALSE))
p_sp2 = merge_phyloseq(p_sp, sampledata)
sample_data(p_sp2)$all<-"all"
spmgd<-merge_samples(p_sp2, "all")
datag = tax_glom(spmgd, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
tabsp<-tab[,-2:-7]
colnames(tabsp)[2]<-"sp"
#Rich
taxg<-data.frame(tax_table(spmgd), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Supergroup))
otug<-otu_table(spmgd)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)
tabr<-tabr[,-2:-7]
z<-expand.grid(Supergroup=clades, Richness=NA, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
z[z$Supergroup==clades[i],"Richness"]<-sum(gtr[,"all"] != 0)
}
z$Total<-sum(z$Richness)
z$RR<-z$Richness/z$Total
richsp<-z[,c("Supergroup","RR")]
colnames(richsp)[2]<-"sp"
total_reads_sp<-sum(sample_sums(spmgd))
total_taxa_sp<-ntaxa(spmgd)


#merging abund
cde<-merge(tabns, tabnp, by="Supergroup", all=T)
ihf<-merge(cde, tabsp, by="Supergroup", all=T)
ihf
for(i in 1:nrow(ihf)){
for(e in 1:ncol(ihf)){
if(is.na(ihf[i,e])){ihf[i,e]<-0}
}}

#merging rich
cder<-merge(richns, richnp, by="Supergroup", all=T)
ihfr<-merge(cder, richsp, by="Supergroup", all=T)
ihfr
for(i in 1:nrow(ihfr)){
for(e in 1:ncol(ihfr)){
if(is.na(ihfr[i,e])){ihfr[i,e]<-0}
}}

#######Merging all abund
ppqq<-merge(poiu, jhgf, all=T)
kkff<-merge(ihf, tabnsp, all=T)
final_t<-merge(ppqq, kkff, all=T)
for(i in 1:nrow(final_t)){
for(e in 1:ncol(final_t)){
if(is.na(final_t[i,e])){final_t[i,e]<-0}
}}
final_t

#######Merging all rich
ppqqr<-merge(poiur, jhgfr, all=T)
kkffr<-merge(ihfr, richnsp, all=T)
final_tr<-merge(ppqqr, kkffr, all=T)
for(i in 1:nrow(final_tr)){
for(e in 1:ncol(final_tr)){
if(is.na(final_tr[i,e])){final_tr[i,e]<-0}
}}
final_tr
colnames(final_tr)<-c("Supergroup","ncbi_r","silva_r","pr2_r","nwsp_r","swnp_r","pwsn_r","ns_r","np_r","sp_r","nsp_r")
final<-merge(final_t, final_tr, by="Supergroup")

##Plotting
rownames(final)<-final$Supergroup
final<-final[,-1]
combined<-t(final)
combined2<-data.frame(combined)
combined2$type<-c("single","single","single","exclusive","exclusive","exclusive","shared","shared","shared", "shared","single","single","single","exclusive","exclusive","exclusive","shared","shared","shared", "shared")
combined2$db<-c("NCBI","SILVA","PR2","NCBI","SILVA","PR2","NCBI + SILVA", "NCBI + PR2", "SILVA + PR2", "NCBI + SILVA + PR2","NCBI","SILVA","PR2","NCBI","SILVA","PR2","NCBI + SILVA", "NCBI + PR2", "SILVA + PR2", "NCBI + SILVA + PR2")
combined2$measure<-c("abundance","abundance","abundance","abundance","abundance","abundance","abundance","abundance",
"abundance","abundance","richness","richness","richness","richness","richness",
"richness","richness","richness","richness","richness")
combined2$total<-c(total_reads_NCBI,total_reads_SILVA,total_reads_pr2,total_reads_nwsp,
total_reads_swnp,total_reads_pwsn,total_reads_ns,total_reads_np,total_reads_sp,total_reads_nsp,
total_taxa_NCBI,total_taxa_SILVA,total_taxa_pr2,total_taxa_nwsp,total_taxa_swnp,total_taxa_pwsn,
total_taxa_ns,total_taxa_np,total_taxa_sp,total_taxa_nsp)
b<-melt(combined2, id=c("type","db","measure","total"), measure=colnames(combined))
b$label<-paste(b$db,"(",b$total,")")
colourCount = length(unique(b$variable))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/comparison")

ggplot(data=b, aes(x=as.factor(label), y=value, fill=variable)) + 
  geom_bar(stat="identity", size=0.2, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_wrap(~measure+type, scales="free")+ labs(title="Supergroup Relative Abundance and Richness", x ="Database comparison", y = "Relative abundance/richness", fill = "Supergroup") + theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 90, hjust = 1, size=6, vjust=0.5), strip.text = element_text(size=7), legend.title=element_text(size=8), legend.text=element_text(size=7), axis.ticks.length=unit(.03, "cm"), legend.key.size = unit(0.5, "cm"))
ggsave("barplot_stacked_supergroup_comparison2.pdf")





