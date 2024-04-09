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
otu_mat<-as.matrix(read.table("f_otu_ncbi_5_01_22.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("f_tax_ncbi_5_01_22.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_ncbi_metadata_5_01_22.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")
metadata$ssh<-paste(metadata$substrate_type, metadata$season, metadata$habitat, sep="_")

sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1<-subset_samples(DADAwang1, !cluster==2)

#Load distance matrix
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
distsea<-read.table("dist_by_sea.txt", sep="\t", header=T)

mds<-acast(distsea, SiteA ~ SiteB)

#prepare data
tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$ssh<-paste(d$substrate_type, d$season, d$habitat, sep="_")
d$sshc<-rownames(d)
sample_data(tudao)<-d[,c("ssh","cluster","season","habitat","substrate_type","sshc")]

####Create plots ~season+habitat -wise
dfmd<-data.frame(sample_data(tudao), stringsAsFactors=FALSE)
fal<-unique(dfmd$ssh)

#Changig order manually
fal<-c("sediment_spring_eelgrass","sediment_autumn_eelgrass","water_spring_eelgrass",
"water_autumn_eelgrass","sediment_spring_sand","sediment_autumn_sand",
"water_spring_sand","water_autumn_sand",
"sediment_spring_rocks","sediment_autumn_rocks",
"water_autumn_rocks","water_spring_rocks")
toplotd<-as.list(fal)



# Jaccard
for (i in 1:length(fal))
{

tcssh2<-subset_samples(tudao, ssh==fal[i])
tcssh0 = filter_taxa(tcssh2, function(x) sum(x) > 0, TRUE)
otug<-data.frame(otu_table(tcssh0), check.names=F)
rownames(otug)<-dfmd$cluster[match(rownames(otug), dfmd$sshc)]
dist_nj<-vegdist(otug, method="jaccard", upper=T)
sites<-as.integer(rownames(otug))

#matching dist matrix
distsea2<-distsea[distsea$SiteA %in% sites, ]
distsea3<-distsea2[distsea2$SiteB %in% sites, ]
comp<-rbind(c(min(distsea3$SiteA),min(distsea3$SiteA),NA), c(c(max(distsea3$SiteB),max(distsea3$SiteB),NA)))
colnames(comp)<-colnames(distsea3)
distsea4<-rbind(distsea3, comp)
mdsx<-acast(distsea4, SiteA ~ SiteB)

mmdsx<-as.matrix(mdsx)
mmdsx <- mmdsx[, order(as.integer(colnames(mmdsx)))]
mmdsx <- mmdsx[order(as.integer(rownames(mmdsx))), ]
mdist_nj<-as.matrix(dist_nj)
mdist_nj <- mdist_nj[, order(as.integer(colnames(mdist_nj)))]
mdist_nj <- mdist_nj[order(as.integer(rownames(mdist_nj))), ]

#storing data
toplotd[[i]]<- data.frame(km=mmdsx[upper.tri(mmdsx, diag = FALSE)], jacc=mdist_nj[upper.tri(mdist_nj, diag = FALSE)])
}

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/dist_by_sea")
pdf("distance_text_ssh.pdf")
par(mfrow=c(3,4))
for (t in 1:length(toplotd))
{
plot(toplotd[[t]]$jacc, toplotd[[t]]$km, 
   xlab="Jaccard", ylab="kilometers", main=fal[t], pch=1, cex =0.8, cex.main=0.9, cex.lab=0.7, cex.axis=0.7)
}
dev.off()




#Bray

for (i in 1:length(fal))
{

tcssh2<-subset_samples(tudao, ssh==fal[i])
tcssh0 = filter_taxa(tcssh2, function(x) sum(x) > 0, TRUE)
otug<-data.frame(otu_table(tcssh0), check.names=F)
rownames(otug)<-dfmd$cluster[match(rownames(otug), dfmd$sshc)]
dist_nj<-vegdist(otug, method="bray", upper=T)
sites<-as.integer(rownames(otug))

#matching dist matrix
distsea2<-distsea[distsea$SiteA %in% sites, ]
distsea3<-distsea2[distsea2$SiteB %in% sites, ]
comp<-rbind(c(min(distsea3$SiteA),min(distsea3$SiteA),NA), c(c(max(distsea3$SiteB),max(distsea3$SiteB),NA)))
colnames(comp)<-colnames(distsea3)
distsea4<-rbind(distsea3, comp)
mdsx<-acast(distsea4, SiteA ~ SiteB)

mmdsx<-as.matrix(mdsx)
mmdsx <- mmdsx[, order(as.integer(colnames(mmdsx)))]
mmdsx <- mmdsx[order(as.integer(rownames(mmdsx))), ]
mdist_nj<-as.matrix(dist_nj)
mdist_nj <- mdist_nj[, order(as.integer(colnames(mdist_nj)))]
mdist_nj <- mdist_nj[order(as.integer(rownames(mdist_nj))), ]

#storing data
toplotd[[i]]<- data.frame(km=mmdsx[upper.tri(mmdsx, diag = FALSE)], jacc=mdist_nj[upper.tri(mdist_nj, diag = FALSE)])
}

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/dist_by_sea")
pdf("distance_text_ssh_bray.pdf")
par(mfrow=c(3,4))
for (t in 1:length(toplotd))
{
plot(toplotd[[t]]$jacc, toplotd[[t]]$km, 
   xlab="Bray", ylab="kilometers", main=fal[t], pch=1, cex =0.8, cex.main=0.9, cex.lab=0.7, cex.axis=0.7)
}
dev.off()


