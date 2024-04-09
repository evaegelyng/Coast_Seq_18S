library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
library("corrplot")
library("gllvm")
library("gclus")
library(ALDEx2)
library(propr)
library(zCompositions)

sessionInfo()

#Load tables
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("f_otu_ncbi.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("f_tax_ncbi.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_ncbi_metadata.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")
metadata$snch<-paste(metadata$season, metadata$cluster, metadata$habitat,sep="_")

sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1
do<-data.frame(sample_data(DADAwang1))


#######################################################

#Phylum plots
####
#Richness
otuo<-data.frame(otu_table(DADAwang1))
taxo<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(DADAwang1))

clades<-levels(factor(taxo$Phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(DADAwang1, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

#Filtering by rich, abund and sites occur
pabund<-as.data.frame(cbind(log10(colMeans(ttab)), colnames(ttab)))
colnames(pabund)[2]<-"clades_p"
prich<-as.data.frame(cbind(log10(colMeans(rich_asv)), colnames(rich_asv)))
colnames(prich)[2]<-"clades_p"
trich<-t(rich_asv)
prich$sqrt_sites_occur <-sqrt(rowSums(trich != 0))
abund_rich_summary <- as.data.frame(merge(prich, pabund, by="clades_p"))
colnames(abund_rich_summary)<-c("clades_p","log10_mean_rich","sqrt_sites_occur","log10_mean_rel_abund")
abund_rich_summary$log10_mean_rich<-as.numeric(abund_rich_summary$log10_mean_rich)
abund_rich_summary$log10_mean_rel_abund<-as.numeric(abund_rich_summary$log10_mean_rel_abund)
summary(abund_rich_summary)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")

abund_rich_summary$clades_p = with(abund_rich_summary, reorder(clades_p, sqrt_sites_occur))
ggplot(abund_rich_summary, aes(clades_p, sqrt_sites_occur)) + 
geom_point(size=1, stroke=0.1, alpha=0.6) + theme_bw() + 
theme(axis.text.x = element_text(size=7, angle = 90, hjust = 1))
ggsave("summplot_sites_ncbi.pdf")

abund_rich_summary$clades_p = with(abund_rich_summary, reorder(clades_p, log10_mean_rich))
ggplot(abund_rich_summary, aes(clades_p, log10_mean_rich)) + 
geom_point(size=1, stroke=0.1, alpha=0.6) + theme_bw() + 
theme(axis.text.x = element_text(size=7, angle = 90, hjust = 1))
ggsave("summplot_rich_ncbi.pdf")

abund_rich_summary$clades_p = with(abund_rich_summary, reorder(clades_p, log10_mean_rel_abund))
ggplot(abund_rich_summary, aes(clades_p, log10_mean_rel_abund)) + 
geom_point(size=1, stroke=0.1, alpha=0.6) + theme_bw() + 
theme(axis.text.x = element_text(size=7, angle = 90, hjust = 1))
ggsave("summplot_abund_ncbi.pdf")

#all

summplot<-melt(abund_rich_summary, id="clades_p", measure=c("log10_mean_rel_abund","log10_mean_rich","sqrt_sites_occur"))
summplot$value<-as.numeric(summplot$value)

ggplot(summplot, aes(clades_p, value)) + 
geom_point(size=1, stroke=0.1, alpha=0.6) + geom_line() +
facet_wrap(~variable, ncol=1, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1))
ggsave("summplot_abund_rich_sites_ncbi.pdf")

abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

list_top15<-c(abundsorted$clades_p[1:15],richsorted$clades_p[1:15],sitessorted$clades_p[1:15])
phy_to_remove<-as.character(unique(list_top15))

##Making plot of all phyla
#Merge rich an freq

combinedr<-cbind(rich_asv, do)
design<-colnames(do)

br<-melt(combinedr, id=design, measure=clades)

a<-rownames(ttab)
b<-rownames(do)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, do)}

taxa<-colnames(ttab)
design<-colnames(do)

b<-melt(combined, id=design, measure=taxa)

b$merger<-paste(b$sshc, b$variable, sep="_")
br$merger<-paste(br$sshc, br$variable, sep="_")
colnames(br)[12]<-"rich"
b2<-b[,c("merger","value")]
br2<-br[,c("season","habitat","substrate_type","variable","merger","rich")]
dg<-merge(b2, br2, by="merger")

colnames(dg)[2]<-c("rel_abund")
dg$shsv<-paste(dg$season,dg$substrate_type,dg$habitat,dg$variable,sep="_")


setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")

ggplot(dg, aes(rich, log10(rel_abund), color=habitat, shape=substrate_type)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.6) +
facet_wrap(~variable, ncol=8, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

###
#Making plot of filtered phyla
`%ni%` <- Negate(`%in%`)
combined2<- subset(combined, select = names(combined) %ni% phy_to_remove)
combinedr2<- subset(combinedr, select = names(combinedr) %ni% phy_to_remove)
clades2<- as.character(colnames(combined2)[1:47])

br2<-melt(combinedr2, id=design, measure=clades2)
b2<-melt(combined2, id=design, measure=clades2)

b2$merger<-paste(b2$sshc, b2$variable, sep="_")
br2$merger<-paste(br2$sshc, br2$variable, sep="_")
colnames(br2)[12]<-"rich"
b2<-b2[,c("merger","value")]
br2<-br2[,c("season","habitat","substrate_type","variable","merger","rich")]
dg2<-merge(b2, br2, by="merger")

colnames(dg2)[2]<-c("rel_abund")
dg2$shsv<-paste(dg2$season,dg2$substrate_type,dg2$habitat,dg2$variable,sep="_")
dg2$Division<-tax$Division[match(dg2$variable, tax$Phylum)]

#Separate Metazoa/by groups of 4
mdg2<-subset(dg2,Division=="Metazoa")
clades_m<-unique(mdg2$variable)

mdg2.1<-mdg2[mdg2$variable %in% clades_m[1:4],]
mdg2.2<-mdg2[mdg2$variable %in% clades_m[5:8],]
mdg2.3<-mdg2[mdg2$variable %in% clades_m[9:12],]
mdg2.4<-mdg2[mdg2$variable %in% clades_m[13:16],]

nmdg2<-subset(dg2,!Division=="Metazoa")
nmdg2$variable = with(nmdg2, reorder(variable, rel_abund))
clades_nm<-unique(nmdg2$variable)

nmdg2.1<-nmdg2[nmdg2$variable %in% clades_nm[1:6],]
nmdg2.2<-nmdg2[nmdg2$variable %in% clades_nm[7:12],]
nmdg2.3<-nmdg2[nmdg2$variable %in% clades_nm[13:18],]
nmdg2.4<-nmdg2[nmdg2$variable %in% clades_nm[19:24],]
nmdg2.5<-nmdg2[nmdg2$variable %in% clades_nm[25:31],]


##Plotting meta and non-meta
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")

ggplot(mdg2, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.1, stroke=0.01, alpha=0.6, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.6, se=TRUE, fullrange=FALSE, level=0.95) +
facet_grid(variable~substrate_type, scale="free", space="free") + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Metaz_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(nmdg2, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.1, stroke=0.01, alpha=0.6, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.6, se=TRUE, fullrange=FALSE, level=0.95) +
facet_grid(variable~substrate_type, scale="free", space="free") + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")


#Plotting groups 4/5 separately

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/meta")

ggplot(mdg2.1, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Metaz1_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(mdg2.2, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Metaz2_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(mdg2.3, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Metaz3_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(mdg2.4, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("Metaz4_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

#Continue
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/nmeta")

ggplot(nmdg2.1, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz1_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(nmdg2.2, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz2_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")


ggplot(nmdg2.3, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz3_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")


ggplot(nmdg2.4, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz4_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")


ggplot(nmdg2.5, aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=0.5, stroke=0.1, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(variable~substrate_type, scale="free", ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=4, angle = 90, hjust = 1), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4)) +
theme(axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title="Phylum - Abundance vs. Richness", x ="ASV Richness", y = "log10 Relative Abundance")
ggsave("NMetaz5_Phylum_Rel_abund_rich_Clusters_ncbi.pdf")


##Plotting individually
#Meta
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/meta")

for (i in 1:length(clades_m))
{
ggplot(subset(mdg2, variable==clades_m[i]), aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=1, stroke=0.5, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(~substrate_type+season, ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1), axis.text.y = element_text(size=8), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),
legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(1, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title=paste(clades_m[i],"Abundance vs. Richness"), x ="ASV Richness", y = "log10 Relative Abundance")
ggsave(paste(clades_m[i],"Rel_abund_rich_Clusters_ncbi.pdf",sep="_"))
}

#Non-meta
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/nmeta")

for (i in 1:length(clades_nm))
{
ggplot(subset(nmdg2, variable==clades_nm[i]), aes(rich, log10(rel_abund), color=habitat)) + 
geom_point(size=1, stroke=0.5, alpha=0.5, shape=20) +
geom_smooth(aes(fill = habitat), method=loess, size=0.2, alpha=0.4, se=TRUE, fullrange=FALSE, level=0.95) +
facet_wrap(~substrate_type+season, ncol=2) + theme_bw() + 
theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1), axis.text.y = element_text(size=8), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),
legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(1, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank()) +
labs(title=paste(clades_nm[i],"Abundance vs. Richness"), x ="ASV Richness", y = "log10 Relative Abundance")
ggsave(paste(clades_nm[i],"Rel_abund_rich_Clusters_ncbi.pdf",sep="_"))
}

#Making overall - rich, abund and sites - relationship plots
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("ssh","variable","season","substrate_type", "habitat"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value))

richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("ssh","variable","season","substrate_type", "habitat"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$ssh, sep="_")
cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$ssh, sep="_")

f_pl_ra<-merge(cdata_rich, cdata_abund, by="refc")
pl_ra<-f_pl_ra[,c("refc","variable.x","season.x","substrate_type.x","habitat.x","mean_freq","sd_freq","mean_rich","sd_rich","sites")]
pl_ra$Division<-tax$Division[match(pl_ra$variable.x, tax$Phylum)]
pl_ra$met<-ifelse(pl_ra$Division=="Metazoa", "Metazoan", "Non-Metazoan")

#Overall plot
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")

ggplot(pl_ra, aes(mean_rich, mean_freq, color=met, group=met)) +
geom_point(size=1, shape=1, color="black", stroke=0.5) +
geom_errorbarh(aes(xmax = mean_rich + sd_rich, xmin = mean_rich - sd_rich), height = 0, size=0.5, alpha=0.8) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq), width = 0, size=0.5, alpha=0.8) + 
geom_smooth(aes(color=met),method = loess, se = FALSE, size=0.5) +
facet_wrap(season.x+substrate_type.x~habitat.x, ncol=3, scales="free") + theme_bw() +
labs(title="Phylum - Abundance vs. Richness", x ="Richness", y = "Relative Abundance") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_overall_ncbi.pdf")

#Plotting met and nmet phyla by season
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")
pl_ra_met<-subset(pl_ra, met=="Metazoan")
pl_ra_met$label_phy<-NA
pl_ra_met$label_phy<-ifelse(pl_ra_met$season.x=="spring",as.character(pl_ra_met$variable.x),NA)

ggplot(pl_ra_met, aes(mean_rich, mean_freq, color=season.x, group=season.x)) +
geom_point(aes(size=sites), shape=20, stroke=0.5, alpha=0.5) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.5,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.2, color="black") +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(habitat.x~substrate_type.x, ncol=2, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5), legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_met_ncbi.pdf")

ggplot(pl_ra_met, aes(mean_rich, mean_freq, color=season.x, group=season.x)) +
geom_point(size=1.2, shape=21, stroke=0.5, alpha=0.8) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.5,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.2, color="black") +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(habitat.x~substrate_type.x, ncol=2, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5), legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_met_ncbi2.pdf")

pl_ra_nmet<-subset(pl_ra, met=="Non-Metazoan")
pl_ra_nmet$label_phy<-NA
pl_ra_nmet$label_phy<-ifelse(pl_ra_nmet$season.x=="spring",as.character(pl_ra_nmet$variable.x),NA)

abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

list_top25<-c(abundsorted$clades_p[1:25],richsorted$clades_p[1:25],sitessorted$clades_p[1:25])

phy_to_remove<-as.character(unique(list_top25))
`%ni%` <- Negate(`%in%`)
pl_ra_nmet2<- pl_ra_nmet[!(pl_ra_nmet$variable.x %in% phy_to_remove),]

ggplot(pl_ra_nmet2, aes(mean_rich, mean_freq, color=season.x, group=season.x)) +
geom_point(aes(size=sites), shape=20, stroke=0.5, alpha=0.5) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.5,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.2, color="black") +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(habitat.x~substrate_type.x, ncol=2, scales="free") + theme_bw() +
labs(title="Non-Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5), legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_nonmet_ncbi.pdf")


ggplot(pl_ra_nmet2, aes(mean_rich, mean_freq, color=season.x, group=season.x)) +
geom_point(size=1.2, shape=21, stroke=0.5, alpha=0.8) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.5,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.2, color="black") +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(habitat.x~substrate_type.x, ncol=2, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5), legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_nonmet_ncbi2.pdf")

#Plotting met and nmet phyla by habitat
do$sh<-paste(do$substrate_type,do$habitat,sep="_")
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("sh","substrate_type", "habitat"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("sh","variable","substrate_type", "habitat"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value))

richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("sh","substrate_type", "habitat"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("sh","variable","substrate_type", "habitat"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$sh, sep="_")
cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$sh, sep="_")

f_pl_ra<-merge(cdata_rich, cdata_abund, by="refc")
pl_ra<-f_pl_ra[,c("refc","variable.x","substrate_type.x","habitat.x","mean_freq","sd_freq","mean_rich","sd_rich","sites")]
pl_ra$Division<-tax$Division[match(pl_ra$variable.x, tax$Phylum)]
pl_ra$met<-ifelse(pl_ra$Division=="Metazoa", "Metazoan", "Non-Metazoan")

#Plotting by habitat
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")
pl_ra_met<-subset(pl_ra, met=="Metazoan")
pl_ra_met$label_phy<-NA
pl_ra_met$label_phy<-ifelse(pl_ra_met$habitat.x=="rocks",as.character(pl_ra_met$variable.x),NA)

ggplot(pl_ra_met, aes(mean_rich, mean_freq, color=habitat.x, group=habitat.x)) +
geom_point(shape=1, size=0.2, stroke=0.1, alpha=0.5, color="black") +
geom_point(aes(size=sites), shape=16, stroke=0.5, alpha=0.5) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.9,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.1, color="black", alpha=0.6) +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(~substrate_type.x, ncol=1, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=8),
legend.title=element_text(size=7), legend.text=element_text(size=7), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
ggsave("Phylum_Rel_abund_rich_met_by_hab_ncbi.pdf")

ggplot(pl_ra_met, aes(mean_rich, mean_freq, color=habitat.x, group=habitat.x)) +
geom_point(shape=1, size=0.2, stroke=0.1, alpha=0.5, color="black") +
geom_point(size=1.4, shape=21, stroke=0.5, alpha=0.8) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 2.1,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.1, color="black", alpha=0.6) +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(~substrate_type.x, ncol=1, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=8),
legend.title=element_text(size=7), legend.text=element_text(size=7), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
ggsave("Phylum_Rel_abund_rich_met_by_hab_ncbi2.pdf")


pl_ra_nmet<-subset(pl_ra, met=="Non-Metazoan")
pl_ra_nmet$label_phy<-NA
pl_ra_nmet$label_phy<-ifelse(pl_ra_nmet$habitat.x=="rocks",as.character(pl_ra_nmet$variable.x),NA)

abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

list_top25<-c(abundsorted$clades_p[1:25],richsorted$clades_p[1:25],sitessorted$clades_p[1:25])

phy_to_remove<-as.character(unique(list_top25))
`%ni%` <- Negate(`%in%`)
pl_ra_nmet2<- pl_ra_nmet[!(pl_ra_nmet$variable.x %in% phy_to_remove),]

ggplot(pl_ra_nmet2, aes(mean_rich, mean_freq, color=habitat.x, group=habitat.x)) +
geom_point(shape=1, size=0.2, stroke=0.1, alpha=0.5, color="black") +
geom_point(aes(size=sites), shape=16, stroke=0.5, alpha=0.5) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 1.9,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.1, color="black", alpha=0.6) +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(~substrate_type.x, ncol=1, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=8),
legend.title=element_text(size=7), legend.text=element_text(size=7), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
ggsave("Phylum_Rel_abund_rich_nonmet_by_hab_ncbi.pdf")

ggplot(pl_ra_nmet2, aes(mean_rich, mean_freq, color=habitat.x, group=habitat.x)) +
geom_point(shape=1, size=0.2, stroke=0.1, alpha=0.5, color="black") +
geom_point(size=1.4, shape=21, stroke=0.5, alpha=0.8) +
geom_smooth(method = lm, se = FALSE, size=0.3) +
geom_text(aes(label=label_phy),size = 2.1,hjust = 1, vjust=0.5, nudge_x = -0.05, color="black", alpha=0.6) +
geom_line(aes(group=variable.x), size=0.1, color="black", alpha=0.6) +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
scale_radius() +
facet_wrap(~substrate_type.x, ncol=1, scales="free") + theme_bw() +
labs(title="Metazoan phyla - Abundance vs. Richness", x ="log10 Richness", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=8),
legend.title=element_text(size=7), legend.text=element_text(size=7), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
ggsave("Phylum_Rel_abund_rich_nonmet_by_hab_ncbi2.pdf")














#################################################################

######Overall salinity plot
#Making overall - rich, abund and sites - relationship plots

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Salinity),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")

#############################################

######### Without habitat

#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","Salinity"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("variable","season","substrate_type","Salinity"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season,cdata_abund$substrate_type,cdata_abund$Salinity, sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$Salinity, sep="_")

#FOr both abund and rich
f_pl_ra<-merge(cdata_rich, cdata_abund, by="refc")
f_pl_ra2<-subset(f_pl_ra, sites>=1)
pl_ra<-f_pl_ra2[,c("refc","variable.x","season.x","substrate_type.x","Salinity.x","mean_freq","sd_freq","mean_rich","sd_rich","sites")]
pl_ra$Division<-tax$Division[match(pl_ra$variable.x, tax$Phylum)]
pl_ra$met<-ifelse(pl_ra$Division=="Metazoa", "Metazoan", "Non-Metazoan")

pl_ra2<-melt(pl_ra, id=c("refc","variable.x","season.x","substrate_type.x","Salinity.x","met"), measure=c("mean_rich","sd_rich","mean_freq","sd_freq","sites"))

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-cdata_rich
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable,pl_r$season,pl_r$Salinity, sep="_")
pl_r$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r$ic2, bricho4$ic2)]

#Richness with log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=season, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_salinity_substrate_ratio.pdf",sep="_"))
}

#Richness classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity")
pg<-as.character(unique(pl_r$variable))

for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_rich, fill=season, group=season, shape=fjord)) +
geom_point(aes(shape=fjord, fill=season), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=season),method = loess, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(~substrate_type, ncol=1) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_salinity_season.pdf",sep="_"))
}

###############################sk
#Phyla individual plots

#Abundance
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/abund_salinity")
pg<-as.character(unique(pl_ra$variable.x))
for (i in 1:length(pg))
{
ggplot(subset(pl_ra, variable.x==pg[i]), aes(Salinity.x, mean_freq, color=season.x, group=season.x)) +
geom_point(size=1, shape=1, color="black", stroke=0.5) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq), width = 0, size=0.5, alpha=0.8) + 
geom_smooth(aes(color=season.x),method = loess, se = FALSE, size=0.5) +
facet_wrap(~substrate_type.x, ncol=1, scales="free_y") + theme_bw() +
labs(title="Phylum - Abundance vs. Salinity", x ="Salinity", y = "Relative Abundance") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_Rel_abund_salinity_ncbi.pdf",sep="_"))
}

#Richness
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity")
pg<-as.character(unique(pl_ra$variable.x))
for (i in 1:length(pg))
{
ggplot(subset(pl_ra, variable.x==pg[i]), aes(Salinity.x, mean_rich, color=season.x, group=season.x)) +
geom_point(size=1, shape=1, color="black", stroke=0.5) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.5, alpha=0.8) + 
geom_smooth(aes(color=season.x),method = loess, se = FALSE, size=0.5) +
facet_wrap(~substrate_type.x, ncol=1, scales="free_y") + theme_bw() +
labs(title="Phylum - Richness vs. Salinity", x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_Rel_rich_salinity_ncbi.pdf",sep="_"))
}


#Phyla individual plots

#Abundance
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/log10_abund_salinity")
pg<-as.character(unique(pl_ra$variable.x))
for (i in 1:length(pg))
{
ggplot(subset(pl_ra, variable.x==pg[i]), aes(Salinity.x, mean_freq, color=season.x, group=season.x)) +
geom_point(size=1, shape=1, color="black", stroke=0.5) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq), width = 0, size=0.5, alpha=0.8) + 
geom_smooth(aes(color=season.x),method = loess, se = FALSE, size=0.5) +
facet_wrap(~substrate_type.x, ncol=1, scales="free_y") +
scale_y_continuous(trans="log10") + theme_bw() +
labs(title="Phylum - Abundance vs. Salinity", x ="Salinity", y = "Relative Abundance") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_Rel_log10_abund_salinity_ncbi.pdf",sep="_"))
}

#Richness
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/log10_rich_salinity")
pg<-as.character(unique(pl_ra$variable.x))
for (i in 1:length(pg))
{
ggplot(subset(pl_ra, variable.x==pg[i]), aes(Salinity.x, mean_rich, color=season.x, group=season.x)) +
geom_point(size=1, shape=1, color="black", stroke=0.5) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.5, alpha=0.8) + 
geom_smooth(aes(color=season.x),method = loess, se = FALSE, size=0.5) +
facet_wrap(~substrate_type.x, ncol=1, scales="free_y") +
scale_y_continuous(trans="log10") + theme_bw() +
labs(title="Phylum - Richness vs. Salinity", x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_Rel_log10_rich_salinity_ncbi.pdf",sep="_"))
}


###################################################bk

######### With habitat

#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","Salinity"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("variable","season","substrate_type","Salinity"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season,cdata_abund$substrate_type,cdata_abund$Salinity, sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Salinity, sep="_")

#FOr both abund and rich
f_pl_ra<-merge(cdata_rich, cdata_abund, by="refc")
f_pl_ra2<-subset(f_pl_ra, sites>=1)
pl_ra<-f_pl_ra2[,c("refc","variable.x","season.x","substrate_type.x","Salinity.x","mean_freq","sd_freq","mean_rich","sd_rich","sites")]
pl_ra$Division<-tax$Division[match(pl_ra$variable.x, tax$Phylum)]
pl_ra$met<-ifelse(pl_ra$Division=="Metazoa", "Metazoan", "Non-Metazoan")

pl_ra2<-melt(pl_ra, id=c("refc","variable.x","season.x","substrate_type.x","Salinity.x","met"), measure=c("mean_rich","sd_rich","mean_freq","sd_freq","sites"))

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(cdata_rich, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable,pl_r$season,pl_r$habitat,pl_r$Salinity, sep="_")
pl_r$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Richness classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity")

pg<-as.character(unique(pl_r$variable))

for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_rich, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = loess, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_salinity_habitat_season.pdf",sep="_"))
}



######## Phylum ~ Class ~ Salinity
#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Salinity),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Class
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Class))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-3]
tabr<-tabr[,-2:-4]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Class")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Salinity, sep="_")

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Class")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Class)]
pl_r$Phylum<-tax$Phylum[match(pl_r$variable, tax$Class)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

##Inspect mean richness per Class so you decide on a cutoff
avg_rich<-ddply(pl_r, c("variable","season","habitat","Phylum"), dplyr::summarize, mean_class_rich = mean(mean_rich), length_class_h3= length(variable[mean_rich >=3]))

avg_rich2<-ddply(avg_rich, c("variable","Phylum"), dplyr::summarize, mean_class_rich = mean(mean_class_rich), mean_length_class_h3=mean(length_class_h3))

#Here arbitrary cutoff implemented
p_good_class<-subset(avg_rich2, mean_class_rich>=2&mean_length_class_h3>=2)

avg_rich_for_plot<-ddply(p_good_class, c("Phylum"), dplyr::summarize, num_classs=length(unique(as.character(variable))))

avg_rich_for_plot2<-subset(avg_rich_for_plot, num_classs>=2)

#Subset
p_pl_r_top<-subset(pl_r, (variable %in% p_good_class$variable))
pl_r_top<-subset(p_pl_r_top, (Phylum %in% avg_rich_for_plot2$Phylum))

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(pl_r_top, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r_top$ic2<-paste(pl_r_top$variable,pl_r_top$season,pl_r_top$habitat,pl_r_top$Salinity, sep="_")
pl_r_top$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r_top$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity/Class")
avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_classs>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_classs>=10,3.5,ifelse(avg_rich_for_plot2$num_classs>=6,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_classs<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_classs>=7,0.1,ifelse(avg_rich_for_plot2$num_classs>=6,1,2))

pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=Phylum, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.2, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Class_rich_salinity_phylum_substrate_ratio.pdf",sep="_"))
}


avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_classs>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_classs>=10,3.5,ifelse(avg_rich_for_plot2$num_classs>=6,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_classs<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_classs>=7,0.1,ifelse(avg_rich_for_plot2$num_classs>=6,1,2))

pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=habitat, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(aes(color=habitat), method = loess, se = FALSE, size=0.2, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Class_rich_salinity_phylum_habitat.pdf",sep="_"))
}



######## Phylum ~ Order ~ Salinity

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Order
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Order))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-4]
tabr<-tabr[,-2:-3]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Order==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Order")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Order
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Salinity, sep="_")

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Order")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Order)]
pl_r$Phylum<-tax$Phylum[match(pl_r$variable, tax$Order)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

##Inspect mean richness per Order so you decide on a cutoff
avg_rich<-ddply(pl_r, c("variable","season","habitat","Phylum"), dplyr::summarize, mean_order_rich = mean(mean_rich), length_order_h3= length(variable[mean_rich >=3]))

avg_rich2<-ddply(avg_rich, c("variable","Phylum"), dplyr::summarize, mean_order_rich = mean(mean_order_rich), mean_length_order_h3=mean(length_order_h3))

#Here arbitrary cutoff implemented
p_good_order<-subset(avg_rich2, mean_order_rich>=3&mean_length_order_h3>=5)

avg_rich_for_plot<-ddply(p_good_order, c("Phylum"), dplyr::summarize, num_orders=length(unique(as.character(variable))))

avg_rich_for_plot2<-subset(avg_rich_for_plot, num_orders>=2)

#Subset
p_pl_r_top<-subset(pl_r, (variable %in% p_good_order$variable))
pl_r_top<-subset(p_pl_r_top, (Phylum %in% avg_rich_for_plot2$Phylum))

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(pl_r_top, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r_top$ic2<-paste(pl_r_top$variable,pl_r_top$season,pl_r_top$habitat,pl_r_top$Salinity, sep="_")
pl_r_top$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r_top$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity/Order")
avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_orders>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_orders>=15,3,ifelse(avg_rich_for_plot2$num_orders>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_orders<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_orders>=7,0.1,2)
pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=Phylum, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.2, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Order_rich_salinity_phylum_substrate_ratio.pdf",sep="_"))
}


avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_orders>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_orders>=15,3,ifelse(avg_rich_for_plot2$num_orders>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_orders<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_orders>=7,0.1,2)
pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=habitat, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(aes(color=habitat), method = loess, se = FALSE, size=0.2, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Order_rich_salinity_phylum_habitat.pdf",sep="_"))
}




######## Phylum ~ Family ~ Salinity

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Family
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Family))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-5]
tabr<-tabr[,-2]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Family==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Family")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Family
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Salinity, sep="_")

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Family")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Family)]
pl_r$Phylum<-tax$Phylum[match(pl_r$variable, tax$Family)]
pl_r$Order<-tax$Order[match(pl_r$variable, tax$Family)]
pl_r$PO<-paste(pl_r$Phylum, pl_r$Order, sep="_")
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

##Inspect mean richness per Family so you decide on a cutoff
avg_rich<-ddply(pl_r, c("variable","season","habitat","PO"), dplyr::summarize, mean_family_rich = mean(mean_rich), length_family_h3= length(variable[mean_rich >=1]))

avg_rich2<-ddply(avg_rich, c("variable","PO"), dplyr::summarize, mean_family_rich = mean(mean_family_rich), mean_length_family_h3=mean(length_family_h3))

#Here arbitrary cutoff implemented
p_good_family<-subset(avg_rich2, mean_family_rich>=2&mean_length_family_h3>=5)

avg_rich_for_plot<-ddply(p_good_family, c("PO"), dplyr::summarize, num_familys=length(unique(as.character(variable))))

avg_rich_for_plot2<-subset(avg_rich_for_plot, num_familys>=2)

#Subset
p_pl_r_top<-subset(pl_r, (variable %in% p_good_family$variable))
pl_r_top<-subset(p_pl_r_top, (PO %in% avg_rich_for_plot2$PO))

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(pl_r_top, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r_top$ic2<-paste(pl_r_top$variable,pl_r_top$season,pl_r_top$habitat,pl_r_top$Salinity, sep="_")
pl_r_top$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r_top$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity/Family")

avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_familys>=6,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_familys>=15,3,ifelse(avg_rich_for_plot2$num_familys>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_familys<=6,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_familys>=6,1,2)
pg<-as.character(unique(pl_r_top$PO))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)

for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, PO==avg_rich_for_plot2[i,"PO"]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=PO, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.2, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"PO"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"PO"],"Family_rich_salinity_phylum_substrate_ratio.pdf",sep="_"))
}


avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_familys>=6,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_familys>=15,3,ifelse(avg_rich_for_plot2$num_familys>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_familys<=6,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_familys>=6,1,2)
pg<-as.character(unique(pl_r_top$PO))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, PO==avg_rich_for_plot2[i,"PO"]), aes(Salinity, mean_rich, fill=habitat, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(aes(color=habitat), method = loess, se = FALSE, size=0.2, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"PO"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"PO"],"Family_rich_salinity_phylum_habitat.pdf",sep="_"))
}



##################################################
##Chlorophyll
######### With habitat
##################################################
##Chlorophyll
######### With habitat
#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Chlorophyll<-merged_wat$Chlorophyll[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Chlorophyll),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Chlorophyll<-merged_wat$Chlorophyll[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")


#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","Chlorophyll"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("variable","season","substrate_type","Chlorophyll"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season,cdata_abund$substrate_type,cdata_abund$Chlorophyll, sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Chlorophyll"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Chlorophyll"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Chlorophyll, sep="_")

#FOr both abund and rich
f_pl_ra<-merge(cdata_rich, cdata_abund, by="refc")
f_pl_ra2<-subset(f_pl_ra, sites>=1)
pl_ra<-f_pl_ra2[,c("refc","variable.x","season.x","substrate_type.x","Chlorophyll.x","mean_freq","sd_freq","mean_rich","sd_rich","sites")]
pl_ra$Division<-tax$Division[match(pl_ra$variable.x, tax$Phylum)]
pl_ra$met<-ifelse(pl_ra$Division=="Metazoa", "Metazoan", "Non-Metazoan")

pl_ra2<-melt(pl_ra, id=c("refc","variable.x","season.x","substrate_type.x","Chlorophyll.x","met"), measure=c("mean_rich","sd_rich","mean_freq","sd_freq","sites"))

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Chlorophyll","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(cdata_rich, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Chlorophyll.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Chlorophyll.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Chlorophyll.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable,pl_r$season,pl_r$habitat,pl_r$Chlorophyll, sep="_")
pl_r$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_chlorophyll")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(log10(Chlorophyll), mean_rich, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. log10(Chlorophyll)"), x ="log10(Chlorophyll)", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_log10chlorophyll_habitat_substrate_ratio.pdf",sep="_"))
}


#Richness classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_chlorophyll")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(log10(Chlorophyll), mean_rich, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = loess, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - Richness vs. log10(Chlorophyll)"), x ="log10(Chlorophyll)", y = "Richness") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_rich_log10chlorophyll_habitat_season.pdf",sep="_"))
}



##################################################
##ABund Salinity
######### With habitat
#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Salinity),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

######
#Rel_abund
datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))
do<-data.frame(sample_data(fdt))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")


#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(ttab))
cdata_abund <- ddply(babundo, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value), sites=sum(value !=0))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season, cdata_abund$substrate_type, cdata_abund$habitat, cdata_abund$Salinity, sep="_")

#FOr abund only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_abund, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(cdata_abund, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_freq.x","mean_freq.y")]
bricho4$water_sediment_ratio<-bricho4$mean_freq.x/bricho4$mean_freq.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_freq.x/bricho4$mean_freq.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable,pl_r$season,pl_r$habitat,pl_r$Salinity, sep="_")
pl_r$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r$ic2, bricho4$ic2)]

#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/abund_salinity")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_freq, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log10 Relative Abundance vs. Salinity"), x ="Salinity", y = "log10 Relative Abundance") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_log10abund_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/abund_salinity")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, mean_freq, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = loess, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log 10 Relative Abundance vs. Salinity"), x ="Salinity", y = "log 10 Relative Abundance") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_log10abund_salinity_habitat_season.pdf",sep="_"))
}







##################################################
##reads_asv Salinity
######### With habitat
#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Salinity),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]


######
#Reads
datag = tax_glom(fdt, "Phylum")

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))
do<-data.frame(sample_data(fdt))

##################3
#OR RELABUND
#Rel_abund
datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))
do<-data.frame(sample_data(fdt))



###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")
do$sshf<-paste(do$substrate_type,do$season,do$habitat,do$field_replicate,sep="_")


#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$cluster, cdata_rich$Salinity, sep="_")

#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(ttab))

cdata_abund <- ddply(babundo, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value), sites=sum(value !=0))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season, cdata_abund$substrate_type, cdata_abund$habitat, cdata_abund$cluster, cdata_abund$Salinity, sep="_")

#FOr both abund and rich - SPECIFY INDEX
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("sshf","ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))

abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("sshf","ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(ttab))


babundo$refc<-paste(babundo$sshf, babundo$cluster, babundo$variable, sep="_")

bricho$refc<-paste(bricho$sshf, bricho$cluster, bricho$variable, sep="_")

f_pl_ra<-merge(babundo, bricho, by="refc")

#INDEX
f_pl_ra$avg_reads_asv<-f_pl_ra$value.y/sqrt(f_pl_ra$value.x)

f_pl_rabeta<-f_pl_ra[,c("season.x","substrate_type.x","cluster.x","habitat.x",
"Salinity.x","variable.x","avg_reads_asv","value.y")]

f_pl_rabeta2<-subset(f_pl_rabeta, value.y>=1)

cdata_ra <- ddply(f_pl_rabeta2, c("variable.x","season.x", "habitat.x","substrate_type.x","cluster.x","Salinity.x"), dplyr::summarize, mean_freq = mean(avg_reads_asv), sd_freq   = sd(avg_reads_asv), sites=sum(value.y !=0))

f_pl_ra2<-subset(cdata_ra, sites>=1)

pl_r<-f_pl_ra2[,c("variable.x","season.x","substrate_type.x","habitat.x","cluster.x","Salinity.x","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
#add beta diversity legendre measures as indicators of community change!!!

f_pl_rabeta2$ic<-paste(f_pl_rabeta2$cluster.x,f_pl_rabeta2$season.x,f_pl_rabeta2$variable.x,f_pl_rabeta2$habitat.x,sep="_")

f_pl_rabeta2w<-subset(f_pl_rabeta2, substrate_type.x=="water")
f_pl_rabeta2s<-subset(f_pl_rabeta2, substrate_type.x=="sediment")


f_pl_rabeta3<-merge(f_pl_rabeta2w, f_pl_rabeta2s, by="ic")

f_pl_rabeta4<-f_pl_rabeta3[,c("ic","season.x.x","cluster.x.x","variable.x.x","habitat.x.x","Salinity.x.x",
"value.y.x","value.y.y")]
f_pl_rabeta4$water_sediment_ratio<-f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y
f_pl_rabeta4$log_water_sediment_ratio<-log10(f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y)
f_pl_rabeta4$ic2<-paste(f_pl_rabeta4$variable.x.x,f_pl_rabeta4$season.x.x,f_pl_rabeta4$habitat.x.x,f_pl_rabeta4$Salinity.x.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable.x,pl_r$season.x,pl_r$habitat.x,pl_r$Salinity.x, sep="_")
pl_r$log_water_sediment_ratio<-f_pl_rabeta4$log_water_sediment_ratio[match(pl_r$ic2, f_pl_rabeta4$ic2)]
pl_r$rich<-f_pl_rabeta4$value.y.y[match(pl_r$ic2, f_pl_rabeta4$ic2)]

colnames(pl_r)<-c("variable","season","substrate_type"
,"habitat","cluster","Salinity","avg_reads_asv",
"sd_avg_reads_asv","sites","Division","met","fjord",
"ic2","log_water_sediment_ratio", "rich")


#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margalef_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


for (i in 1:length(pg))
{
vg<-subset(pl_r, variable==pg[i])
min_l<-min(vg$rich)
max_l<-max(vg$rich)
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=rich, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=rich), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(low = "yellow", high = "red",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margalef_perM_reads_salinity_habitat_rich.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margalef_salinity_habitat_season.pdf",sep="_"))
}

#Abund with habitat and log_water_sediment_ratio linear
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm,se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


for (i in 1:length(pg))
{
vg<-subset(pl_r, variable==pg[i])
min_l<-min(vg$rich)
max_l<-max(vg$rich)
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=rich, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=rich), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm,se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(low = "yellow", high = "red",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_perM_reads_salinity_habitat_rich.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm,se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_salinity_habitat_season.pdf",sep="_"))
}


#Abund with habitat and log_water_sediment_ratio linear
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log10 avg_reads_asv vs. Salinity"), x ="Salinity", y = "log10 avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10avg_reads_asv_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log 10 avg_reads_asv vs. Salinity"), x ="Salinity", y = "log 10 avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10avg_reads_asv_salinity_habitat_season.pdf",sep="_"))
}


######Same but with SQRT reads_asvs
#FOr both abund and rich
babundo$refc<-paste(babundo$ssh, babundo$cluster, babundo$variable, sep="_")

bricho$refc<-paste(bricho$ssh, bricho$cluster, bricho$variable, sep="_")

f_pl_ra<-merge(babundo, bricho, by="refc")
f_pl_ra$avg_reads_asv<-sqrt(f_pl_ra$value.x/f_pl_ra$value.y)
f_pl_rabeta<-f_pl_ra[,c("season.x","substrate_type.x","cluster.x","habitat.x",
"Salinity.x","variable.x","avg_reads_asv","value.y")]

f_pl_rabeta2<-subset(f_pl_rabeta, value.y>=1)

cdata_ra <- ddply(f_pl_rabeta2, c("variable.x","season.x", "habitat.x","substrate_type.x","cluster.x","Salinity.x"), dplyr::summarize, mean_freq = mean(avg_reads_asv), sd_freq   = sd(avg_reads_asv), sites=sum(value.y !=0))

f_pl_ra2<-subset(cdata_ra, sites>=1)

pl_r<-f_pl_ra2[,c("variable.x","season.x","substrate_type.x","habitat.x","cluster.x","Salinity.x","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)

f_pl_rabeta2$ic<-paste(f_pl_rabeta2$cluster.x,f_pl_rabeta2$season.x,f_pl_rabeta2$variable.x,f_pl_rabeta2$habitat.x,sep="_")

f_pl_rabeta2w<-subset(f_pl_rabeta2, substrate_type.x=="water")
f_pl_rabeta2s<-subset(f_pl_rabeta2, substrate_type.x=="sediment")


f_pl_rabeta3<-merge(f_pl_rabeta2w, f_pl_rabeta2s, by="ic")

f_pl_rabeta4<-f_pl_rabeta3[,c("ic","season.x.x","cluster.x.x","variable.x.x","habitat.x.x","Salinity.x.x",
"value.y.x","value.y.y")]
f_pl_rabeta4$water_sediment_ratio<-f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y
f_pl_rabeta4$log_water_sediment_ratio<-log10(f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y)
f_pl_rabeta4$ic2<-paste(f_pl_rabeta4$variable.x.x,f_pl_rabeta4$season.x.x,f_pl_rabeta4$habitat.x.x,f_pl_rabeta4$Salinity.x.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable.x,pl_r$season.x,pl_r$habitat.x,pl_r$Salinity.x, sep="_")
pl_r$log_water_sediment_ratio<-f_pl_rabeta4$log_water_sediment_ratio[match(pl_r$ic2, f_pl_rabeta4$ic2)]

colnames(pl_r)<-c("variable","season","substrate_type"
,"habitat","cluster","Salinity","avg_reads_asv",
"sd_avg_reads_asv","sites","Division","met","fjord",
"ic2","log_water_sediment_ratio")


#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log10 sqrt_avg_reads_asv vs. Salinity"), x ="Salinity", y = "log10 sqrt_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_logsqrt10avg_reads_asv_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log 10 sqrt_avg_reads_asv vs. Salinity"), x ="Salinity", y = "log 10 sqrt_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10sqrtavg_reads_asv_salinity_habitat_season.pdf",sep="_"))
}


#Abund with habitat and log_water_sediment_ratio, sqrt - no log10
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() +
labs(title=paste(pg[i]," - sqrt_avg_reads_asv vs. Salinity"), x ="Salinity", y = "sqrt_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_sqrt_avg_reads_asv_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - sqrt_avg_reads_asv vs. Salinity"), x ="Salinity", y = "sqrt_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_sqrt_avg_reads_asv_salinity_habitat_season.pdf",sep="_"))
}


######Same but with log10_relabund_t_sqrt_rich   !!!!!!!go up and set RELABUND
#FOr both abund and rich
babundo$refc<-paste(babundo$ssh, babundo$cluster, babundo$variable, sep="_")

bricho$refc<-paste(bricho$ssh, bricho$cluster, bricho$variable, sep="_")

f_pl_ra<-merge(babundo, bricho, by="refc")
f_pl_ra$avg_reads_asv<-sqrt(log10(f_pl_ra$value.x)*(-1*(f_pl_ra$value.y^2)))
f_pl_rabeta<-f_pl_ra[,c("season.x","substrate_type.x","cluster.x","habitat.x",
"Salinity.x","variable.x","avg_reads_asv","value.y")]

f_pl_rabeta2<-subset(f_pl_rabeta, value.y>=1)

cdata_ra <- ddply(f_pl_rabeta2, c("variable.x","season.x", "habitat.x","substrate_type.x","cluster.x","Salinity.x"), dplyr::summarize, mean_freq = mean(avg_reads_asv), sd_freq   = sd(avg_reads_asv), sites=sum(value.y !=0))

f_pl_ra2<-subset(cdata_ra, sites>=1)

pl_r<-f_pl_ra2[,c("variable.x","season.x","substrate_type.x","habitat.x","cluster.x","Salinity.x","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)

f_pl_rabeta2$ic<-paste(f_pl_rabeta2$cluster.x,f_pl_rabeta2$season.x,f_pl_rabeta2$variable.x,f_pl_rabeta2$habitat.x,sep="_")

f_pl_rabeta2w<-subset(f_pl_rabeta2, substrate_type.x=="water")
f_pl_rabeta2s<-subset(f_pl_rabeta2, substrate_type.x=="sediment")


f_pl_rabeta3<-merge(f_pl_rabeta2w, f_pl_rabeta2s, by="ic")

f_pl_rabeta4<-f_pl_rabeta3[,c("ic","season.x.x","cluster.x.x","variable.x.x","habitat.x.x","Salinity.x.x",
"value.y.x","value.y.y")]
f_pl_rabeta4$water_sediment_ratio<-f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y
f_pl_rabeta4$log_water_sediment_ratio<-log10(f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y)
f_pl_rabeta4$ic2<-paste(f_pl_rabeta4$variable.x.x,f_pl_rabeta4$season.x.x,f_pl_rabeta4$habitat.x.x,f_pl_rabeta4$Salinity.x.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable.x,pl_r$season.x,pl_r$habitat.x,pl_r$Salinity.x, sep="_")
pl_r$log_water_sediment_ratio<-f_pl_rabeta4$log_water_sediment_ratio[match(pl_r$ic2, f_pl_rabeta4$ic2)]

colnames(pl_r)<-c("variable","season","substrate_type"
,"habitat","cluster","Salinity","avg_reads_asv",
"sd_avg_reads_asv","sites","Division","met","fjord",
"ic2","log_water_sediment_ratio")


#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() +
labs(title=paste(pg[i]," - log10_relabund_t_sqrt_rich vs. Salinity"), x ="Salinity", y = "log10_relabund_t_sqrt_rich") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10_relabund_t_sqrt_rich_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - log10_relabund_t_sqrt_rich vs. Salinity"), x ="Salinity", y = "log10_relabund_t_sqrt_rich") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10_relabund_t_sqrt_rich_salinity_habitat_season.pdf",sep="_"))
}



######Same but with SQRT ONLY READS reads_asvs
#FOr both abund and rich
babundo$refc<-paste(babundo$ssh, babundo$cluster, babundo$variable, sep="_")

bricho$refc<-paste(bricho$ssh, bricho$cluster, bricho$variable, sep="_")

f_pl_ra<-merge(babundo, bricho, by="refc")
f_pl_ra$avg_reads_asv<-sqrt(f_pl_ra$value.x)/f_pl_ra$value.y
f_pl_rabeta<-f_pl_ra[,c("season.x","substrate_type.x","cluster.x","habitat.x",
"Salinity.x","variable.x","avg_reads_asv","value.y")]

f_pl_rabeta2<-subset(f_pl_rabeta, value.y>=1)

cdata_ra <- ddply(f_pl_rabeta2, c("variable.x","season.x", "habitat.x","substrate_type.x","cluster.x","Salinity.x"), dplyr::summarize, mean_freq = mean(avg_reads_asv), sd_freq   = sd(avg_reads_asv), sites=sum(value.y !=0))

f_pl_ra2<-subset(cdata_ra, sites>=1)

pl_r<-f_pl_ra2[,c("variable.x","season.x","substrate_type.x","habitat.x","cluster.x","Salinity.x","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)

f_pl_rabeta2$ic<-paste(f_pl_rabeta2$cluster.x,f_pl_rabeta2$season.x,f_pl_rabeta2$variable.x,f_pl_rabeta2$habitat.x,sep="_")

f_pl_rabeta2w<-subset(f_pl_rabeta2, substrate_type.x=="water")
f_pl_rabeta2s<-subset(f_pl_rabeta2, substrate_type.x=="sediment")


f_pl_rabeta3<-merge(f_pl_rabeta2w, f_pl_rabeta2s, by="ic")

f_pl_rabeta4<-f_pl_rabeta3[,c("ic","season.x.x","cluster.x.x","variable.x.x","habitat.x.x","Salinity.x.x",
"value.y.x","value.y.y")]
f_pl_rabeta4$water_sediment_ratio<-f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y
f_pl_rabeta4$log_water_sediment_ratio<-log10(f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y)
f_pl_rabeta4$ic2<-paste(f_pl_rabeta4$variable.x.x,f_pl_rabeta4$season.x.x,f_pl_rabeta4$habitat.x.x,f_pl_rabeta4$Salinity.x.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable.x,pl_r$season.x,pl_r$habitat.x,pl_r$Salinity.x, sep="_")
pl_r$log_water_sediment_ratio<-f_pl_rabeta4$log_water_sediment_ratio[match(pl_r$ic2, f_pl_rabeta4$ic2)]

colnames(pl_r)<-c("variable","season","substrate_type"
,"habitat","cluster","Salinity","avg_reads_asv",
"sd_avg_reads_asv","sites","Division","met","fjord",
"ic2","log_water_sediment_ratio")


#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log10 sqrt_or_avg_reads_asv vs. Salinity"), x ="Salinity", y = "log10 sqrt_or_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10sqrt_or_avg_reads_asv_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() + scale_y_continuous(trans="log10") +
labs(title=paste(pg[i]," - log 10 sqrt_or_avg_reads_asv vs. Salinity"), x ="Salinity", y = "log 10 sqrt_or_avg_reads_asv") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_ln_log10sqrt_or_avg_reads_asv_salinity_habitat_season.pdf",sep="_"))
}


######## Phylum ~ Order ~ Salinity SQRT reads_asv

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Order
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Order))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-4]
tabr<-tabr[,-2:-3]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Order==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]

######
#Rel_abund
datag = tax_glom(fdt, "Order")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Order
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")

#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$Salinity, sep="_")

#FOr rich only - adding fjord, water/sediment ratio, habitat
f_pl_r2<-subset(cdata_rich, sites>=1)
pl_r<-f_pl_r2[,c("refc","variable","season","substrate_type","habitat","cluster","Salinity","mean_rich","sd_rich","sites")]

datag = tax_glom(fdt, "Order")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Order)]
pl_r$Phylum<-tax$Phylum[match(pl_r$variable, tax$Order)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

##Inspect mean richness per Order so you decide on a cutoff
avg_rich<-ddply(pl_r, c("variable","season","habitat","Phylum"), dplyr::summarize, mean_order_rich = mean(mean_rich), length_order_h3= length(variable[mean_rich >=3]))

avg_rich2<-ddply(avg_rich, c("variable","Phylum"), dplyr::summarize, mean_order_rich = mean(mean_order_rich), mean_length_order_h3=mean(length_order_h3))

#Here arbitrary cutoff implemented
p_good_order<-subset(avg_rich2, mean_order_rich>=3&mean_length_order_h3>=5)

avg_rich_for_plot<-ddply(p_good_order, c("Phylum"), dplyr::summarize, num_orders=length(unique(as.character(variable))))

avg_rich_for_plot2<-subset(avg_rich_for_plot, num_orders>=2)

#Subset
p_pl_r_top<-subset(pl_r, (variable %in% p_good_order$variable))
pl_r_top<-subset(p_pl_r_top, (Phylum %in% avg_rich_for_plot2$Phylum))

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
bricho2<-subset(pl_r_top, sites>=1)
bricho2$ic<-paste(bricho2$cluster,bricho2$season,bricho2$variable,bricho2$habitat,sep="_")
bricho2w<-subset(bricho2, substrate_type=="water")
bricho2s<-subset(bricho2, substrate_type=="sediment")

bricho3<-merge(bricho2w, bricho2s, by="ic")

bricho4<-bricho3[,c("ic","season.x","cluster.x","variable.x","habitat.x","Salinity.x",
"mean_rich.x","mean_rich.y")]
bricho4$water_sediment_ratio<-bricho4$mean_rich.x/bricho4$mean_rich.y
bricho4$log_water_sediment_ratio<-log10(bricho4$mean_rich.x/bricho4$mean_rich.y)
bricho4$ic2<-paste(bricho4$variable.x,bricho4$season.x,bricho4$habitat.x,bricho4$Salinity.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r_top$ic2<-paste(pl_r_top$variable,pl_r_top$season,pl_r_top$habitat,pl_r_top$Salinity, sep="_")
pl_r_top$log_water_sediment_ratio<-bricho4$log_water_sediment_ratio[match(pl_r_top$ic2, bricho4$ic2)]

#Richness with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/rich_salinity/Order")
avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_orders>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_orders>=15,3,ifelse(avg_rich_for_plot2$num_orders>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_orders<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_orders>=7,0.1,2)
pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=log_water_sediment_ratio, group=Phylum, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(method = loess, se = FALSE, size=0.2, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Order_rich_salinity_phylum_substrate_ratio.pdf",sep="_"))
}


avg_rich_for_plot2$szp<-ifelse(avg_rich_for_plot2$num_orders>=7,0.8,1)
avg_rich_for_plot2$sz<-ifelse(avg_rich_for_plot2$num_orders>=15,3,ifelse(avg_rich_for_plot2$num_orders>=7,4,5))
avg_rich_for_plot2$nc<-ifelse(avg_rich_for_plot2$num_orders<=7,4,8)
avg_rich_for_plot2$fsx<-0.1
avg_rich_for_plot2$fsy<-ifelse(avg_rich_for_plot2$num_orders>=7,0.1,2)
pg<-as.character(unique(pl_r_top$Phylum))

min_l<-min(pl_r_top$log_water_sediment_ratio)
max_l<-max(pl_r_top$log_water_sediment_ratio)
for (i in 1:nrow(avg_rich_for_plot2))
{
ggplot(subset(pl_r_top, Phylum==avg_rich_for_plot2[i,"Phylum"]), aes(Salinity, mean_rich, fill=habitat, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=avg_rich_for_plot2[i,"szp"], stroke=0.1) +
geom_errorbar(aes(ymin = mean_rich-sd_rich, ymax = mean_rich+sd_rich), width = 0, size=0.1, alpha=0.8) + 
geom_smooth(aes(color=habitat), method = loess, se = FALSE, size=0.2, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_grid(variable~substrate_type+season, scales="free_y") + theme_bw() +
labs(title=paste(avg_rich_for_plot2[i,"Phylum"]," - Richness vs. Salinity"), x ="Salinity", y = "Richness") +
theme(strip.text = element_text(size=avg_rich_for_plot2[i,"sz"]), strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")),
legend.title=element_text(size=avg_rich_for_plot2[i,"sz"]),legend.text=element_text(size=avg_rich_for_plot2[i,"sz"]), axis.title=element_text(size=avg_rich_for_plot2[i,"sz"]), plot.title = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.x = element_text(size=avg_rich_for_plot2[i,"sz"]), axis.text.y = element_text(size=avg_rich_for_plot2[i,"sz"]), panel.spacing.x = unit(avg_rich_for_plot2[i,"fsx"], "lines"), panel.spacing.y = unit(avg_rich_for_plot2[i,"fsy"], "lines"))
ggsave(paste(avg_rich_for_plot2[i,"Phylum"],"Order_salinity_phylum_habitat.pdf",sep="_"))
}





##################################################
##shannon salinity
######### With habitat
#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]

##Update OTU table
with_NA1<-rownames(do[is.na(do$Salinity),])

newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

####
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(fdt))

clades<-levels(factor(taxo$Phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)

tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
ch<-do$sample_root
z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]



######
#Reads
datag = tax_glom(fdt, "Phylum")

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))
do<-data.frame(sample_data(fdt))

##################3
#OR RELABUND
#Rel_abund
datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))
do<-data.frame(sample_data(fdt))



###
#again

#Import water data
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)

merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")

do$Salinity<-merged_wat$Salinity[match(do$snch, merged_wat$snch)]
do$ssh<-paste(do$substrate_type,do$season,do$habitat,sep="_")
do$sshf<-paste(do$substrate_type,do$season,do$habitat,do$field_replicate,sep="_")


#rich (some changes made - cluster, habitat added, temperature removed)
richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))
cdata_rich <- ddply(bricho, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_rich = mean(value), sd_rich = sd(value), sites=sum(value !=0))

cdata_rich$refc<-paste(cdata_rich$variable, cdata_rich$season,cdata_rich$substrate_type, cdata_rich$habitat, cdata_rich$cluster, cdata_rich$Salinity, sep="_")

#abund
abundo<-cbind(ttab, do)
babundo<-melt(abundo, id=c("ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(ttab))

cdata_abund <- ddply(babundo, c("variable","season", "habitat","substrate_type","cluster","Salinity"), dplyr::summarize, mean_freq = mean(value), sd_freq   = sd(value), sites=sum(value !=0))

cdata_abund$refc<-paste(cdata_abund$variable, cdata_abund$season, cdata_abund$substrate_type, cdata_abund$habitat, cdata_abund$cluster, cdata_abund$Salinity, sep="_")


#SHANNON

z<-expand.grid("sample_root"=ch, stringsAsFactors=FALSE)
for(i in 1:length(clades))
{
sgsg<-subset_taxa(fdt, Phylum==clades[i])
sgsg = filter_taxa(sgsg, function(x) sum(x) > 0, TRUE)
sgsg = prune_samples(sample_sums(sgsg)>0,sgsg)
otuosg<-as.matrix(data.frame(otu_table(sgsg)))
braw<-diversity(otuosg, index = "shannon", MARGIN = 1)
z$nada<-NA
t<-1+i
colnames(z)[t]<-clades[i]
z[,t]<-braw[match(z$sample_root, names(braw))]
}

z<-z[,-1]

#FOr both abund and rich - SPECIFY INDEX

shano<-cbind(z, do)
bshano<-melt(shano, id=c("sshf","ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(z))

richo<-cbind(rich_asv, do)
bricho<-melt(richo, id=c("sshf","ssh","season","substrate_type", "habitat","cluster","Salinity"), measure=colnames(rich_asv))

bshano$refc<-paste(bshano$sshf, bshano$cluster, bshano$variable, sep="_")

bricho$refc<-paste(bricho$sshf, bricho$cluster, bricho$variable, sep="_")

f_pl_ra<-merge(bshano, bricho, by="refc")

#INDEX

f_pl_rabeta<-f_pl_ra[,c("season.x","substrate_type.x","cluster.x","habitat.x",
"Salinity.x","variable.x","value.x","value.y")]

f_pl_rabeta2<-subset(f_pl_rabeta, value.y>=1)

cdata_ra <- ddply(f_pl_rabeta2, c("variable.x","season.x", "habitat.x","substrate_type.x","cluster.x","Salinity.x"), dplyr::summarize, mean_freq = mean(value.x), sd_freq   = sd(value.x), sites=sum(value.y !=0))

f_pl_ra2<-subset(cdata_ra, sites>=1)

pl_r<-f_pl_ra2[,c("variable.x","season.x","substrate_type.x","habitat.x","cluster.x","Salinity.x","mean_freq","sd_freq","sites")]

datag = tax_glom(fdt, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)

pl_r$Division<-tax$Division[match(pl_r$variable, tax$Phylum)]
pl_r$met<-ifelse(pl_r$Division=="Metazoa", "Metazoan", "Non-Metazoan")
pl_r$fjord<-ifelse(pl_r$cluster==2|pl_r$cluster==12|pl_r$cluster==16|pl_r$cluster==29, "fjord","open")

#water/sediment ratio (keep it for later when using higher taxonomic resolution)
#add beta diversity legendre measures as indicators of community change!!!

f_pl_rabeta2$ic<-paste(f_pl_rabeta2$cluster.x,f_pl_rabeta2$season.x,f_pl_rabeta2$variable.x,f_pl_rabeta2$habitat.x,sep="_")

f_pl_rabeta2w<-subset(f_pl_rabeta2, substrate_type.x=="water")
f_pl_rabeta2s<-subset(f_pl_rabeta2, substrate_type.x=="sediment")


f_pl_rabeta3<-merge(f_pl_rabeta2w, f_pl_rabeta2s, by="ic")

f_pl_rabeta4<-f_pl_rabeta3[,c("ic","season.x.x","cluster.x.x","variable.x.x","habitat.x.x","Salinity.x.x",
"value.y.x","value.y.y")]
f_pl_rabeta4$water_sediment_ratio<-f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y
f_pl_rabeta4$log_water_sediment_ratio<-log10(f_pl_rabeta4$value.y.x/f_pl_rabeta4$value.y.y)
f_pl_rabeta4$ic2<-paste(f_pl_rabeta4$variable.x.x,f_pl_rabeta4$season.x.x,f_pl_rabeta4$habitat.x.x,f_pl_rabeta4$Salinity.x.x, sep="_")

#tresg<-subset(bricho4, !water_sediment_ratio==Inf|!water_sediment_ratio==NaN)
#avg_ratio<-ddply(tresg, c("variable.x","season.x","habitat.x","cluster.x","Salinity.x"), dplyr::summarize, mean_ratio = mean(water_sediment_ratio), sd_ratio = sd(water_sediment_ratio))

pl_r$ic2<-paste(pl_r$variable.x,pl_r$season.x,pl_r$habitat.x,pl_r$Salinity.x, sep="_")
pl_r$log_water_sediment_ratio<-f_pl_rabeta4$log_water_sediment_ratio[match(pl_r$ic2, f_pl_rabeta4$ic2)]
pl_r$rich<-f_pl_rabeta4$value.y.y[match(pl_r$ic2, f_pl_rabeta4$ic2)]

colnames(pl_r)<-c("variable","season","substrate_type"
,"habitat","cluster","Salinity","avg_shannon",
"sd_shannon","sites","Division","met","fjord",
"ic2","log_water_sediment_ratio", "rich")


#Abund with habitat and log_water_sediment_ratio
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/shannon")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_shannon, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_shannon-sd_shannon, ymax = avg_shannon+sd_shannon), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - shannon vs. Salinity"), x ="Salinity", y = "shannon") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_shannon_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


for (i in 1:length(pg))
{
vg<-subset(pl_r, variable==pg[i])
min_l<-min(vg$rich)
max_l<-max(vg$rich)
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_shannon, fill=rich, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=rich), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_shannon-sd_shannon, ymax = avg_shannon+sd_shannon), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(low = "yellow", high = "red",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - shannon vs. Salinity"), x ="Salinity", y = "shannon") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_shannon_reads_salinity_habitat_rich.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/shannon")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_shannon, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_shannon-sd_shannon, ymax = avg_shannon+sd_shannon), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - shannon vs. Salinity"), x ="Salinity", y = "shannon") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_shannon_salinity_habitat_season.pdf",sep="_"))
}

#Abund with habitat and log_water_sediment_ratio linear
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")
pg<-as.character(unique(pl_r$variable))

min_l<-min(pl_r$log_water_sediment_ratio)
max_l<-max(pl_r$log_water_sediment_ratio)
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=log_water_sediment_ratio, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=log_water_sediment_ratio), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm,se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_salinity_habitat_substrate_ratio.pdf",sep="_"))
}


for (i in 1:length(pg))
{
vg<-subset(pl_r, variable==pg[i])
min_l<-min(vg$rich)
max_l<-max(vg$rich)
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=rich, group=habitat, shape=fjord)) +
geom_point(aes(shape=fjord, fill=rich), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(method = lm,se = FALSE, size=0.3, color="black", alpha=0.7) +
scale_fill_gradient2(low = "yellow", high = "red",na.value = "grey50", limits=c(min_l,max_l)) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(substrate_type~season+habitat, ncol=3) + theme_bw() + 
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_perM_reads_salinity_habitat_rich.pdf",sep="_"))
}


#Abund classic
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich/reads_asvs")

pg<-as.character(unique(pl_r$variable))
pl_r$hss<-paste(pl_r$habitat,pl_r$substrate_type, pl_r$season, sep="_")
for (i in 1:length(pg))
{
ggplot(subset(pl_r, variable==pg[i]), aes(Salinity, avg_reads_asv, fill=habitat, group=hss, shape=fjord)) +
geom_point(aes(shape=fjord, fill=habitat), color="black", size=2, stroke=0.4) +
geom_errorbar(aes(ymin = avg_reads_asv-sd_avg_reads_asv, ymax = avg_reads_asv+sd_avg_reads_asv), width = 0, size=0.4, alpha=0.8) + 
geom_smooth(aes(color=habitat),method = lm,se = FALSE, size=0.3, alpha=0.7) +
scale_shape_manual(values=c(21,22)) +
facet_wrap(season~substrate_type, ncol=2) + theme_bw() +
labs(title=paste(pg[i]," - margalef vs. Salinity"), x ="Salinity", y = "margalef") +
theme(strip.text = element_text(size=5),
legend.title=element_text(size=5),legend.text=element_text(size=5), axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6))
ggsave(paste(pg[i],"Phylum_margaleflm_salinity_habitat_season.pdf",sep="_"))
}


######
#Corr rich rel
######
#Phylum
#Prepare merged object
tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sshc<-rownames(d)
head(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc")]

#Rich
taxg<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Phylum))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-rowSums(combined[,1:length(clades)],)
design<-colnames(d)

br<-melt(combined, id=c(design, "Total"), measure=clades)
br$rel_rich<-br$value/br$Total

#Freq
####get taxa abund
datag = tax_glom(tudao, "Phylum")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)


#Merge rich an freq
b$merger<-paste(b$sshc, b$variable, sep="_")
br$merger<-paste(br$sshc, br$variable, sep="_")
b2<-b[,c("merger","value")]
br2<-br[,c("season","habitat","substrate_type","variable","merger","rel_rich")]
dg<-merge(b2, br2, by="merger")

colnames(dg)[2]<-c("rel_freq")
dg$shsv<-paste(dg$season,dg$substrate_type,dg$habitat,dg$variable,sep="_")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_rich")

ggplot(dg, aes(rich, log10(rel_abund), color=habitat, shape=substrate_type)) + 
geom_point(size=1, stroke=0.2, alpha=0.6) +
facet_wrap(~variable, ncol=8, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Phylum - Abundance vs. Richness", x ="log10 Relative Richness", y = "log10 Relative Abundance")
ggsave("Phylum_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=season)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Phylum - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Phylum_Rel_abund_rich_season_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Phylum - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Phylum_Rel_abund_rich_habitat_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=season), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Phylum - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Phylum_Rel_abund_rich_habitat_season_ncbi.pdf")


#Cross plots
cdata_rich_freq1 <- ddply(dg, c("season","substrate_type", "variable"), summarise,
               N_freq    = length(rel_freq),
               mean_freq = mean(rel_freq),
               sd_freq   = sd(rel_freq),
               se_freq   = sd_freq / sqrt(N_freq),
               N_rich    = length(rel_rich),
               mean_rich = mean(rel_rich),
               sd_rich   = sd(rel_rich),
               se_rich   = sd_rich / sqrt(N_rich)
)

ggplot(cdata_rich_freq1, aes(mean_rich, mean_freq, color=season)) +
geom_point(size=0.2) +
geom_errorbarh(aes(xmax = mean_rich + sd_rich, xmin = mean_rich - sd_rich, color= season, linetype=substrate_type),height = 0, size=0.5, alpha=0.8) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = season, linetype=substrate_type),width = 0, size=0.5, alpha=0.8) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() +
labs(title="Phylum - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_st_season_mean_cross_ncbi.pdf")


cdata_rich_freq2 <- ddply(dg, c("substrate_type", "habitat", "variable"), summarise,
               N_freq    = length(rel_freq),
               mean_freq = mean(rel_freq),
               sd_freq   = sd(rel_freq),
               se_freq   = sd_freq / sqrt(N_freq),
               N_rich    = length(rel_rich),
               mean_rich = mean(rel_rich),
               sd_rich   = sd(rel_rich),
               se_rich   = sd_rich / sqrt(N_rich)
)

ggplot(cdata_rich_freq2, aes(mean_rich, mean_freq, color=habitat)) +
geom_point(size=0.2) +
geom_errorbarh(aes(xmax = mean_rich + sd_rich, xmin = mean_rich - sd_rich, color= habitat, linetype=substrate_type),height = 0, size=0.5, alpha=0.8) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = habitat, linetype=substrate_type),width = 0, size=0.5, alpha=0.8) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() +
labs(title="Phylum - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Phylum_Rel_abund_rich_habitat_st_mean_cross_ncbi.pdf")

