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

sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1

#####
#####
#Species accumulation plots
#Cluster ~ substrate_type + season + habitat

ju<-expand.grid("Sample_root"=c("NA"), "Richness"=c("NA"), "SD"=c("NA"), "Sites"=c("NA"), "ssc"=c("NA"), "sshc"=c("NA"), stringsAsFactors=FALSE)
metadata<-data.frame(sample_data(DADAwang1))

grupao<-levels(as.factor(sample_data(DADAwang1)$ssc))
for (i in 1:length(grupao))
{
wat10<-subset_samples(DADAwang1, ssc==grupao[i])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat10 = prune_samples(sample_sums(wat10)>0,wat10)

gru<-levels(as.factor(sample_data(wat10)$sshc))
for (e in 1:length(gru))
{
wat100<-subset_samples(wat10, sshc==gru[e])
wat100 = filter_taxa(wat100, function(x) sum(x) > 0, TRUE)
wat100 = prune_samples(sample_sums(wat100)>0,wat100)

rotu<-data.frame(otu_table(wat100),check.names=F)
if(nrow(rotu)<=1) next
sdfr<-specaccum(rotu,method="random")
er<-data.frame(rownames(rotu), sdfr$richness, sdfr$sd, sdfr$sites, grupao[i], gru[e])
colnames(er)<-c("Sample_root","Richness","SD", "Sites", "ssc","sshc")
ju<-rbind(ju, er)
}
}

ju<-ju[-1,]
ju$Substrate_type<-metadata$substrate_type[match(ju$Sample_root, metadata$sample_root)]
ju$Habitat<-metadata$habitat[match(ju$Sample_root, metadata$sample_root)]
ju$Season<-metadata$season[match(ju$Sample_root, metadata$sample_root)]
ju$Cluster<-metadata$cluster[match(ju$Sample_root, metadata$sample_root)]

ju3<-ju
ju3$Richness<-as.numeric(ju3$Richness)
ju3$SD<-as.numeric(ju3$SD)
ju3$Sites<-as.numeric(ju3$Sites)
ju3$Cluster<-as.numeric(as.character(ju3$Cluster))

ju3$Cluster <- factor(ju3$Cluster, levels = unique(ju3$Cluster[order(ju3$Cluster)]))
ju3$Season <- factor(ju3$Season, levels = c("spring", "autumn"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/accu/merged")

dodge_t <- position_dodge(width=0.3)

ggplot(subset(ju3, Substrate_type=="sediment"), aes(x = Sites, y = Richness, color = as.factor(Season), linetype=as.factor(Habitat), group=sshc)) +
geom_point(size=0.9, shape=20, alpha=0.6, stroke=0.2, position=dodge_t) +
geom_line(size=0.4, alpha=0.6, position=dodge_t) +
scale_linetype_manual(values=c("longdash", "dotted", "dotdash")) + scale_color_manual(values=c("black", "red")) +
geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), position=position_dodge(0.3), width=0.2, size=0.2, linetype="solid") +
facet_wrap(~Cluster, ncol=5) +
ggtitle("ASV accumulation - Sediment") +
theme_bw() +
theme(plot.title = element_text(size=8), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), axis.text.x = element_text(size=4.5), axis.text.y = element_text(size=4.5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=4), legend.key.size = unit(0.8, "cm")) +
scale_x_continuous(breaks=c(1,2,3))
ggsave("merged_ncbi_cluster_sediment.pdf")

ggplot(subset(ju3, Substrate_type=="water"), aes(x = Sites, y = Richness, color = as.factor(Season), linetype=as.factor(Habitat), group=sshc)) +
geom_point(size=0.9, shape=20, alpha=0.6, stroke=0.2, position=dodge_t) +
geom_line(size=0.4, alpha=0.6, position=dodge_t) +
scale_linetype_manual(values=c("longdash", "dotted", "dotdash")) + scale_color_manual(values=c("black", "red")) +
geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), position=position_dodge(0.3), width=0.2, size=0.2, linetype="solid") +
facet_wrap(~Cluster, ncol=5) +
ggtitle("ASV accumulation - Water") +
theme_bw() +
theme(plot.title = element_text(size=8), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), axis.text.x = element_text(size=4.5), axis.text.y = element_text(size=4.5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=4), legend.key.size = unit(0.8, "cm")) +
scale_x_continuous(breaks=c(1,2,3))
ggsave("merged_ncbi_cluster_water.pdf")

#Cluster ~ substrate_type + season

ju<-expand.grid("Sample_root"=c("NA"), "Richness"=c("NA"), "SD"=c("NA"), "Sites"=c("NA"), "ssc"=c("NA"), stringsAsFactors=FALSE)
metadata<-data.frame(sample_data(DADAwang1))

grupao<-levels(as.factor(sample_data(DADAwang1)$ssc))
for (i in 1:length(grupao))
{
wat10<-subset_samples(DADAwang1, ssc==grupao[i])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat10 = prune_samples(sample_sums(wat10)>0,wat10)

rotu<-data.frame(otu_table(wat10),check.names=F)
if(nrow(rotu)<=1) next
sdfr<-specaccum(rotu,method="random")
er<-data.frame(rownames(rotu), sdfr$richness, sdfr$sd, sdfr$sites, grupao[i])
colnames(er)<-c("Sample_root","Richness","SD", "Sites", "ssc")
ju<-rbind(ju, er)
}

ju<-ju[-1,]
ju$Substrate_type<-metadata$substrate_type[match(ju$Sample_root, metadata$sample_root)]
ju$Season<-metadata$season[match(ju$Sample_root, metadata$sample_root)]
ju$Cluster<-metadata$cluster[match(ju$Sample_root, metadata$sample_root)]

ju3<-ju
ju3$Richness<-as.numeric(ju3$Richness)
ju3$SD<-as.numeric(ju3$SD)
ju3$Sites<-as.numeric(ju3$Sites)
ju3$Cluster<-as.numeric(as.character(ju3$Cluster))

ju3$Cluster <- factor(ju3$Cluster, levels = unique(ju3$Cluster[order(ju3$Cluster)]))
ju3$Season <- factor(ju3$Season, levels = c("spring", "autumn"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/accu/merged")

dodge_t <- position_dodge(width=0.3)

ggplot(ju3, aes(x = Sites, y = Richness, color = as.factor(Season), linetype=as.factor(Substrate_type), group=ssc)) +
geom_point(size=0.9, shape=20, alpha=0.6, stroke=0.2, position=dodge_t) +
geom_line(size=0.4, alpha=0.6, position=dodge_t) +
scale_linetype_manual(values=c("longdash", "dotted")) + scale_color_manual(values=c("black", "red")) +
geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), position=position_dodge(0.3), width=0.2, size=0.2, linetype="solid") +
facet_wrap(~Cluster, ncol=4) +
ggtitle("ASV accumulation - all habitats") +
theme_bw() +
theme(plot.title = element_text(size=8), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), axis.text.x = element_text(size=4.5), axis.text.y = element_text(size=4.5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=4), legend.key.size = unit(0.8, "cm")) + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9))
ggsave("merged_ncbi_cluster_all_hab.pdf")

#Cluster ~ season

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_ncbi_metadata.txt", sep="\t", header=T)
metadata$sc<-paste(metadata$season, metadata$cluster, sep="_")
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))
DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1

ju<-expand.grid("Sample_root"=c("NA"), "Richness"=c("NA"), "SD"=c("NA"), "Sites"=c("NA"), "sc"=c("NA"), stringsAsFactors=FALSE)
metadata<-data.frame(sample_data(DADAwang1))

grupao<-levels(as.factor(sample_data(DADAwang1)$sc))
for (i in 1:length(grupao))
{
wat10<-subset_samples(DADAwang1, sc==grupao[i])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat10 = prune_samples(sample_sums(wat10)>0,wat10)

rotu<-data.frame(otu_table(wat10),check.names=F)
if(nrow(rotu)<=1) next
sdfr<-specaccum(rotu,method="random")
er<-data.frame(rownames(rotu), sdfr$richness, sdfr$sd, sdfr$sites, grupao[i])
colnames(er)<-c("Sample_root","Richness","SD", "Sites", "sc")
ju<-rbind(ju, er)
}

ju<-ju[-1,]
ju$Season<-metadata$season[match(ju$Sample_root, metadata$sample_root)]
ju$Cluster<-metadata$cluster[match(ju$Sample_root, metadata$sample_root)]

ju3<-ju
ju3$Richness<-as.numeric(ju3$Richness)
ju3$SD<-as.numeric(ju3$SD)
ju3$Sites<-as.numeric(ju3$Sites)
ju3$Cluster<-as.numeric(as.character(ju3$Cluster))

ju3$Cluster <- factor(ju3$Cluster, levels = unique(ju3$Cluster[order(ju3$Cluster)]))
ju3$Season <- factor(ju3$Season, levels = c("spring", "autumn"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/accu/merged")

dodge_t <- position_dodge(width=0.3)

ggplot(ju3, aes(x = Sites, y = Richness, linetype=as.factor(Season), group=sc)) +
geom_point(size=0.9, shape=20, alpha=0.6, stroke=0.2, position=dodge_t) +
geom_line(size=0.4, alpha=0.6, position=dodge_t) +
scale_linetype_manual(values=c("longdash", "dotted")) +
geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), position=position_dodge(0.3), width=0.2, size=0.2, linetype="solid") +
facet_wrap(~Cluster, ncol=4) +
ggtitle("ASV accumulation - all habitats and substrates") +
theme_bw() +
theme(plot.title = element_text(size=8), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), axis.text.x = element_text(size=4.5), axis.text.y = element_text(size=4.5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=4), legend.key.size = unit(0.8, "cm")) + scale_x_continuous(breaks=c(1,3,6,9,12,15,18))
ggsave("merged_ncbi_cluster_all_hab_sub.pdf")

# ~ Cluster

ju<-expand.grid("Sample_root"=c("NA"), "Richness"=c("NA"), "SD"=c("NA"), "Sites"=c("NA"), "Cluster"=c("NA"), stringsAsFactors=FALSE)
metadata<-data.frame(sample_data(DADAwang1))

grupao<-levels(as.factor(sample_data(DADAwang1)$cluster))
for (i in 1:length(grupao))
{
wat10<-subset_samples(DADAwang1, cluster==grupao[i])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat10 = prune_samples(sample_sums(wat10)>0,wat10)

rotu<-data.frame(otu_table(wat10),check.names=F)
if(nrow(rotu)<=1) next
sdfr<-specaccum(rotu,method="random")
er<-data.frame(rownames(rotu), sdfr$richness, sdfr$sd, sdfr$sites, grupao[i])
colnames(er)<-c("Sample_root","Richness","SD", "Sites", "Cluster")
ju<-rbind(ju, er)
}

ju<-ju[-1,]
ju$Cluster<-metadata$cluster[match(ju$Sample_root, metadata$sample_root)]

ju3<-ju
ju3$Richness<-as.numeric(ju3$Richness)
ju3$SD<-as.numeric(ju3$SD)
ju3$Sites<-as.numeric(ju3$Sites)
ju3$Cluster<-as.numeric(as.character(ju3$Cluster))

ju3$Cluster <- factor(ju3$Cluster, levels = unique(ju3$Cluster[order(ju3$Cluster)]))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/accu/merged")

dodge_t <- position_dodge(width=0.3)

ggplot(ju3, aes(x = Sites, y = Richness, group=Cluster)) +
geom_point(size=0.9, shape=20, alpha=0.6, stroke=0.2, position=dodge_t) +
geom_line(size=0.4, alpha=0.6, position=dodge_t) +
geom_errorbar(aes(ymin = Richness - SD, ymax = Richness + SD), position=position_dodge(0.3), width=0.2, size=0.2, linetype="solid") +
facet_wrap(~Cluster, ncol=4) +
ggtitle("ASV accumulation - all habitats, substrates and seasons") +
theme_bw() +
theme(plot.title = element_text(size=8), axis.title.x = element_text(size=6), axis.title.y = element_text(size=6), axis.text.x = element_text(size=4.5), axis.text.y = element_text(size=4.5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5)) + scale_x_continuous(breaks=c(1,6,12,18,24,30,36))
ggsave("merged_ncbi_cluster_all_hab_sub_sea.pdf")

#####
#####
#Rarefaction plots

#WHOLE DATASET

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = 1000)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = TRUE)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/rare/merged")

grupao<-levels(as.factor(sample_data(DADAwang1)$cluster))
for (e in 1:length(grupao))
{
wat10<-subset_samples(DADAwang1, cluster==grupao[e])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat101 = prune_samples(sample_sums(wat10)>0,wat10)
x<-data.frame(otu_table(wat101),check.names=F)
tot <- rowSums(x)
S <- rowSums(x > 0)
nr <- nrow(x)

out <- lapply(seq_len(nr), rarefun)
df <- do.call(rbind, out)

sdf<-data.frame(sample_data(wat101))
sdf$Sample <- rownames(sdf)
data <- merge(df, sdf, by = "Sample")
labels <- data.frame(x = tot, y = S, Sample = rownames(x))
labels <- merge(labels, sdf, by = "Sample")

colnames(data)[2]<-"ASVs"

ggplot(data = data, aes_string(x = "Size", y = "ASVs", group = "Sample")) + 
labs(x = "Reads", y = "ASVs") +
geom_text(data = labels, aes_string(x = "x", y = "y", label = "sample_root"), size = 4, hjust = 0) +
geom_line() + facet_wrap(substrate_type~habitat+season, ncol=2) + xlim(0, max(labels$x + 200000)) + theme_bw() + theme(axis.text.x = element_text(size=5), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5)) + theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste(grupao[e],"rare_ncbi.pdf",sep="_"))
}

#####
#####
#Taxonomy plots


#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("f_ncbi_metadata.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")

sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1

#For stacked plots, use "tudao" as starting object:
#Stacked

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


######
#Relative abundance
######
#Supergroup

d<-data.frame(sample_data(DADAwang1))
datarg = transform_sample_counts(DADAwang1, function(x) x/sum(x))
datag = tax_glom(datarg, "Supergroup")

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Supergroup
tab<-tab[,-1:-7]
ttab<-t(data.frame(round(tab, digits=3),check.names=F))

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Supergroup Relative Abundance", x ="Supergroup", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Supergorup_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Supergroup Relative Abundance", x ="Supergroup", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Supergorup_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Supergroup Relative Abundance", x ="Supergroup", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Supergorup_ncbi.pdf")

#Stacked

datag = tax_glom(tudao, "Supergroup")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Supergroup
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

cdata2$variable<-factor(cdata2$variable, levels=c("Alveolata","Amoebozoa","Archaeplastida","Excavata",
"Hacrobia","Opisthokonta","Rhizaria","SA13C06","Stramenopiles"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#D9D9D9")) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Supergroup Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Supergroup") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Supergroup_ncbi.pdf")


#Division

d<-data.frame(sample_data(DADAwang1))
datarg = transform_sample_counts(DADAwang1, function(x) x/sum(x))
datag = tax_glom(datarg, "Division")

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Division
tab<-tab[,-1:-6]
ttab<-t(data.frame(round(tab, digits=3),check.names=F))

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Division Relative Abundance", x ="Division", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Supergorup_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Division Relative Abundance", x ="Division", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Supergorup_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Division Relative Abundance", x ="Division", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Supergorup_ncbi.pdf")

#Stacked

datag = tax_glom(tudao, "Division")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Division
tab<-tab[,-1:-6]
ttab<-t(data.frame(tab, check.names=F))

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(16, "Paired"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Division Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Division") + theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Division_ncbi.pdf")

#Phylum - Non Metazoa
no<-subset_taxa(DADAwang1, !Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))
datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:15])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:15])
}


a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Phylum_nop_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Phylum_nop_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Phylum_nop_ncbi.pdf")

#Stacked
no<-subset_taxa(tudao, !Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:15])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:15])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(10, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Non-Metazoa Phylum Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Non-Metazoa Phylum") + theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Phylum_nop_ncbi_5_1_22.pdf")

#Phylum - Metazoa
no<-subset_taxa(DADAwang1, Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))
datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}


a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Phylum_m_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Phylum_m_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Phylum_m_ncbi.pdf")

#Stacked
no<-subset_taxa(tudao, Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(10, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Metazoa Phylum Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Metazoa Phylum") + theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Phylum_m_ncbi_5_1_22.pdf")


#Division - Opisthokonta
no<-subset_taxa(DADAwang1, Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Phylum")

d<-data.frame(sample_data(datag))

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Phylum_op_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Phylum_op_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Phylum Relative Abundance", x ="Phylum", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Phylum_op_ncbi.pdf")

#Stacked

datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, Division=="Metazoa")
datag = tax_glom(no, "Phylum")

d<-data.frame(sample_data(datag))

tax<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
otu<-otu_table(datag)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-6]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(10, "Set3"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Metazoa Phylum Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Metazoa Phylum") + theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Phylum_op_ncbi.pdf")


#Class - Non Metazoa

datarg = transform_sample_counts(DADAwang1, function(x) x/sum(x))
no<-subset_taxa(datarg, !Division=="Metazoa")
datag = tax_glom(no, "Class")

TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}


a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Class_nom_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Class_nom_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Class_nom_ncbi.pdf")

#Stacked

datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, !Division=="Metazoa")
datag = tax_glom(no, "Class")

TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(10, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Non-Metazoa Class Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Non-Metazoa Class") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Class_nom_ncbi.pdf")

#Class - Metazoa

datarg = transform_sample_counts(DADAwang1, function(x) x/sum(x))
no<-subset_taxa(datarg, Division=="Metazoa")
datag = tax_glom(no, "Class")

TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

#Side-by-side
cdata <- ddply(b, c("season", "substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_season_Class_m_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_st_Class_m_ncbi.pdf")

ggplot(cdata, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Class Relative Abundance", x ="Class", y = "Relative abundance", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_sbs_habitat_Class_m_ncbi.pdf")

#Stacked

datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, Division=="Metazoa")
datag = tax_glom(no, "Class")

TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:9])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]
pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA

for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:9])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

colourCount = length(unique(cdata2$variable))
getPalette = colorRampPalette(brewer.pal(10, "Set3"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund/stacked")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Metazoans Relative Abundance", x ="Cluster number", y = "Relative abundance", fill = "Metazoans") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_stacked_Class_m_ncbi.pdf")


#####Richness
#Supergroup

#Side-by-side
taxg<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Supergroup))
d<-data.frame(sample_data(DADAwang1))
ch<-d$sshc
otug<-otu_table(DADAwang1)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-2:-7]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)

fgr<-melt(combined, id=c("season","substrate_type","habitat"), measure=clades)
head(fgr)

cdata_rich <- ddply(fgr, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Supergroup Richness", x ="Supergroup", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Supergroup_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Supergroup Richness", x ="Supergroup", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Supergroup_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Supergroup Richness", x ="Supergroup", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Supergroup_ncbi.pdf")

#Stacked

taxg<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Supergroup))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-2:-7]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Supergroup==clades[i])
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
br$rel<-br$value/br$Total

br$variable<-factor(br$variable, levels=c("Alveolata","Amoebozoa","Archaeplastida","Excavata",
"Hacrobia","Opisthokonta","Rhizaria","SA13C06","Stramenopiles"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#D9D9D9")) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Supergroup Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Supergroup") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Supergroup_ncbi.pdf")

#Division

#Side-by-side
taxg<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(DADAwang1))
ch<-d$sshc
otug<-otu_table(DADAwang1)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-2:-7]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)

fgr<-melt(combined, id=c("season","substrate_type","habitat"), measure=clades)
head(fgr)

cdata_rich <- ddply(fgr, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Division_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Division_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Division_ncbi.pdf")

#Stacked

taxg<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-2:-7]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
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
br$rel<-br$value/br$Total

colourCount = length(clades)
getPalette = colorRampPalette(brewer.pal(length(clades), "Paired"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Division Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Division") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Division_ncbi.pdf")


#Phylum - Non metazoa
#skip for getting relative within non-opis
#Calculate totals for Division relative richness
#Rich
taxgt<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
cladest<-levels(factor(taxgt$Division))
d<-data.frame(sample_data(DADAwang1))
ch<-d$sshc
otugt<-otu_table(DADAwang1)
otugt<-t(data.frame(otugt, check.names=F))
tabrt<-cbind(taxgt, otugt)

#Tax-level-wise
tabrt<-tabrt[,-1]
tabrt<-tabrt[,-2:-6]
zt<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladest))
{
gtrt<-subset(tabrt, Division==cladest[i])
richt<-colSums(gtrt[,-1] != 0)
zt<-cbind(zt, richt)
tt<-1+i
colnames(zt)[tt]<-cladest[i]
}
zt<-zt[,-1]
combinedt<-cbind(zt, d)
combinedt$Total<-rowSums(combinedt[,1:length(cladest)],)


#Now, subset and go!
no<-subset_taxa(DADAwang1, !Division=="Metazoa")

#Side-by-side
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Phylum))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
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
z$Total<-rowSums(z)
combined<-cbind(z, d)
#Old way
#combined$Total<-combinedt$Total[match(combined$sshc, combinedt$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:15])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))

cdata_rich <- ddply(br, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Phylum_nom_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Phylum_nom_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Phylum_nom_ncbi.pdf")

#Stacked
#skip
#Calculate totals for Division relative richness
#Rich
taxgts<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
cladests<-levels(factor(taxgts$Division))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otugts<-otu_table(tudao)
otugts<-t(data.frame(otugts, check.names=F))
tabrts<-cbind(taxgts, otugts)

#Tax-level-wise
tabrts<-tabrts[,-1]
tabrts<-tabrts[,-2:-6]
zts<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladests))
{
gtrts<-subset(tabrts, Division==cladests[i])
richts<-colSums(gtrts[,-1] != 0)
zts<-cbind(zts, richts)
tts<-1+i
colnames(zts)[tts]<-cladests[i]
}
zts<-zts[,-1]
combinedts<-cbind(zts, d)
combinedts$Total<-rowSums(combinedts[,1:length(cladests)],)

#Now, subset and go!
no<-subset_taxa(tudao, !Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Phylum))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
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
z$Total<-rowSums(z)
combined<-cbind(z, d)
#Old way
#combined$Total<-combinedt$Total[match(combined$sshc, combinedt$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:15])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))
br$rel<-br$value/br$Total

colourCount = length(top9clades)+1
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Non-Metazoa Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Non-Metazoa") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Phylum_nom_ncbi.pdf")


#Phylum - Metazoa
#skip for getting relative within non-opis
#Calculate totals for Division relative richness
#Rich
taxgt<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
cladest<-levels(factor(taxgt$Division))
d<-data.frame(sample_data(DADAwang1))
ch<-d$sshc
otugt<-otu_table(DADAwang1)
otugt<-t(data.frame(otugt, check.names=F))
tabrt<-cbind(taxgt, otugt)

#Tax-level-wise
tabrt<-tabrt[,-1]
tabrt<-tabrt[,-2:-6]
zt<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladest))
{
gtrt<-subset(tabrt, Division==cladest[i])
richt<-colSums(gtrt[,-1] != 0)
zt<-cbind(zt, richt)
tt<-1+i
colnames(zt)[tt]<-cladest[i]
}
zt<-zt[,-1]
combinedt<-cbind(zt, d)
combinedt$Total<-rowSums(combinedt[,1:length(cladest)],)


#Now, subset and go!
no<-subset_taxa(DADAwang1, Division=="Metazoa")

#Side-by-side
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Phylum))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
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
z$Total<-rowSums(z)
combined<-cbind(z, d)
#Old way
#combined$Total<-combinedt$Total[match(combined$sshc, combinedt$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))

cdata_rich <- ddply(br, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Phylum_m_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Phylum_m_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Phylum_m_ncbi.pdf")

#Stacked
#skip
#Calculate totals for Division relative richness
#Rich
taxgts<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
cladests<-levels(factor(taxgts$Division))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otugts<-otu_table(tudao)
otugts<-t(data.frame(otugts, check.names=F))
tabrts<-cbind(taxgts, otugts)

#Tax-level-wise
tabrts<-tabrts[,-1]
tabrts<-tabrts[,-2:-6]
zts<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladests))
{
gtrts<-subset(tabrts, Division==cladests[i])
richts<-colSums(gtrts[,-1] != 0)
zts<-cbind(zts, richts)
tts<-1+i
colnames(zts)[tts]<-cladests[i]
}
zts<-zts[,-1]
combinedts<-cbind(zts, d)
combinedts$Total<-rowSums(combinedts[,1:length(cladests)],)

#Now, subset and go!
no<-subset_taxa(tudao, Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Phylum))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
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
z$Total<-rowSums(z)
combined<-cbind(z, d)
#Old way
#combined$Total<-combinedt$Total[match(combined$sshc, combinedt$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))
br$rel<-br$value/br$Total

colourCount = length(top9clades)+1
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Metazoa Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Metazoa") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Phylum_m_ncbi.pdf")



#Division - Opisthokonta

no<-subset_taxa(DADAwang1, Division=="Opisthokonta")

#Side-by-side
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1]
tabr<-tabr[,-2:-6]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)

fgr<-melt(combined, id=c("season","substrate_type","habitat"), measure=clades)
head(fgr)

cdata_rich <- ddply(fgr, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Division_op_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Division_op_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Division Richness", x ="Division", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Division_op_ncbi.pdf")

#Stacked
no<-subset_taxa(tudao, Division=="Opisthokonta")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1]
tabr<-tabr[,-2:-6]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedts$Total[match(combined$sshc, combinedts$sshc)]
design<-colnames(d)

br<-melt(combined, id=c(design, "Total"), measure=clades)
br$rel<-br$value/br$Total

colourCount = length(clades)
getPalette = colorRampPalette(brewer.pal(length(clades), "Set3"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Opisthokonta Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Opisthokonta") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Division_op_ncbi.pdf")

#Class - Non Metazoa
#Calculate Class totals 
#Rich
taxgj<-data.frame(tax_table(DADAwang1), stringsAsFactors=FALSE)
cladesj<-levels(factor(taxgj$Class))
d<-data.frame(sample_data(DADAwang1))
ch<-d$sshc
otugj<-otu_table(DADAwang1)
otugj<-t(data.frame(otugj, check.names=F))
tabrj<-cbind(taxgj, otugj)

#Tax-level-wise
tabrj<-tabrj[,-1:-2]
tabrj<-tabrj[,-2:-5]
zj<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladesj))
{
gtrj<-subset(tabrj, Class==cladesj[i])
richj<-colSums(gtrj[,-1] != 0)
zj<-cbind(zj, richj)
tj<-1+i
colnames(zj)[tj]<-cladesj[i]
}
zj<-zj[,-1]
combinedj<-cbind(zj, d)
combinedj$Total<-rowSums(combinedj[,1:length(cladesj)],)

#Now, subset and go
no<-subset_taxa(DADAwang1, !Division=="Metazoa")

#Side-by-side
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedj$Total[match(combined$sshc, combinedj$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=top9clades)

cdata_rich <- ddply(br, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Class_nom_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Class_nom_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Class_nom_ncbi.pdf")

#Stacked
#Calculate Class totals 
#Rich
taxgjw<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
cladesjw<-levels(factor(taxgjw$Class))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otugjw<-otu_table(tudao)
otugjw<-t(data.frame(otugjw, check.names=F))
tabrjw<-cbind(taxgjw, otugjw)

#Tax-level-wise
tabrjw<-tabrjw[,-1:-2]
tabrjw<-tabrjw[,-2:-5]
zjw<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladesjw))
{
gtrjw<-subset(tabrjw, Class==cladesjw[i])
richjw<-colSums(gtrjw[,-1] != 0)
zjw<-cbind(zjw, richjw)
tjw<-1+i
colnames(zjw)[tjw]<-cladesj[i]
}
zjw<-zjw[,-1]
combinedjw<-cbind(zjw, d)
combinedjw$Total<-rowSums(combinedjw[,1:length(cladesjw)],)

#Now, subset and go!
no<-subset_taxa(tudao, !Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedjw$Total[match(combined$sshc, combinedjw$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))
br$rel<-br$value/br$Total

colourCount = length(top9clades)+1
getPalette = colorRampPalette(brewer.pal(length(top9clades)+1, "Set1"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Non-Metazoa Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Non-Metazoa") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Class_nom_ncbi.pdf")

#Class - Metazoa

no<-subset_taxa(DADAwang1, Division=="Metazoa")

#Side-by-side
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedj$Total[match(combined$sshc, combinedj$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=top9clades)

cdata_rich <- ddply(br, c("season","substrate_type","habitat", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)

head(cdata_rich)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/sbs")
dodge_t <- position_dodge(width=0.8)
ggplot(cdata_rich, aes(variable, mean, fill=season)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(substrate_type~habitat) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "season") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_season_Class_m_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=substrate_type)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~habitat) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "substrate") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_st_Class_m_ncbi.pdf")

ggplot(cdata_rich, aes(variable, mean, fill=habitat)) + geom_errorbar(aes(ymax=mean+sd, ymin=mean), size=0.4, position = position_dodge(0.8), width=0.6) + geom_bar(size=0.4, width=0.6, color="black", stat="identity", position=dodge_t) + facet_wrap(season~substrate_type) + labs(title="Class Richness", x ="Class", y = "Richness", fill = "habitat") + theme_bw() + theme(axis.text.y = element_text(size=7), axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7)) + scale_fill_grey()
ggsave("barplot_richness_sbs_habitat_Class_m_ncbi.pdf")

#Stacked
no<-subset_taxa(tudao, Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
head(taxg)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedjw$Total[match(combined$sshc, combinedjw$sshc)]
top9clades<-names(sort(colSums(combined[,1:length(clades)],), TRUE) [1:9])
combined$Total_top<-rowSums(combined[,top9clades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=c(top9clades, "Other"))
br$rel<-br$value/br$Total

colourCount = length(top9clades)+1
getPalette = colorRampPalette(brewer.pal(length(top9clades)+1, "Set3"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/rich/stacked")

ggplot(data=br, aes(x=as.factor(cluster), y=rel, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values = getPalette(colourCount)) + facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="Metazoa Relative Richness", x ="Cluster number", y = "Relative Richness", fill = "Metazoa") + theme_bw() + scale_y_continuous(limits=c(0, 1.01), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("barplot_richness_stacked_Class_m_ncbi.pdf")

######
#Corr rich rel
######
#Division
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
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-2:-7]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
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
datag = tax_glom(tudao, "Division")
datarg = transform_sample_counts(datag, function(x) x/sum(x))

tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Division
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

ggplot(dg, aes(log10(rel_rich), log10(rel_freq), color=habitat, shape=substrate_type)) + 
geom_point(size=1, stroke=0.2, alpha=0.6) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division - Abundance vs. Richness", x ="log10 Relative Richness", y = "log10 Relative Abundance")
ggsave("Division_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=season)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_Rel_abund_rich_season_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_Rel_abund_rich_habitat_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=season), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_Rel_abund_rich_habitat_season_ncbi.pdf")


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
labs(title="Division - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Division_Rel_abund_rich_st_season_mean_cross_ncbi.pdf")


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
labs(title="Division - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Division_Rel_abund_rich_habitat_st_mean_cross_ncbi.pdf")


#Division non Opisthokonta
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

#Define groups - top rich != top abund
topclades<-c("Ochrophyta","Ciliophora","Cercozoa","Dinoflagellata","Pseudofungi","Chlorophyta","Cryptophyta","Haptophyta")

#Calculate totals for Division relative richness
#Rich
taxgt<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
cladest<-levels(factor(taxgt$Division))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otugt<-otu_table(tudao)
otugt<-t(data.frame(otugt, check.names=F))
tabrt<-cbind(taxgt, otugt)

#Tax-level-wise
tabrt<-tabrt[,-1]
tabrt<-tabrt[,-2:-6]
zt<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladest))
{
gtrt<-subset(tabrt, Division==cladest[i])
richt<-colSums(gtrt[,-1] != 0)
zt<-cbind(zt, richt)
tt<-1+i
colnames(zt)[tt]<-cladest[i]
}
zt<-zt[,-1]
combinedt<-cbind(zt, d)
combinedt$Total<-rowSums(combinedt[,1:length(cladest)],)

#Now, with subset
#Rich
no<-subset_taxa(tudao, !Division=="Opisthokonta")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1]
tabr<-tabr[,-2:-6]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Division==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
#
combined$Total<-combinedt$Total[match(combined$sshc, combinedt$sshc)]
combined$Total_top<-rowSums(combined[,topclades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=topclades)
br$rel_rich<-br$value/br$Total

#Freq
####get taxa abund
datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, !Division=="Opisthokonta")
datag = tax_glom(no, "Division")
p_topcladesseq<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
topcladesseq<-rownames(p_topcladesseq[p_topcladesseq$Division %in% topclades,])
data9 = prune_taxa(topcladesseq, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Division
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)

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

ggplot(dg, aes(log10(rel_rich), log10(rel_freq), color=habitat, shape=substrate_type)) + 
geom_point(size=1, stroke=0.2, alpha=0.6) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="log10 Relative Richness", y = "log10 Relative Abundance")
ggsave("Division_nop_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=season)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_nop_Rel_abund_rich_season_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_nop_Rel_abund_rich_habitat_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=season), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Division_nop_Rel_abund_rich_habitat_season_ncbi.pdf")


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
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Division_nop_Rel_abund_rich_st_season_mean_cross_ncbi.pdf")


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
labs(title="Division non Opisthokonta - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Division_nop_Rel_abund_rich_habitat_st_mean_cross_ncbi.pdf")


#Class Metazoa
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

#Define groups - top rich != top abund
topclades<-c("Mollusca","Arthropoda","Nematoda","Platyhelminthes","Cnidaria","Rotifera","Bryozoa","Craniata")

#Calculate Class totals 
#Rich
taxgj<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
cladesj<-levels(factor(taxgj$Class))
d<-data.frame(sample_data(tudao))
ch<-d$sshc
otugj<-otu_table(tudao)
otugj<-t(data.frame(otugj, check.names=F))
tabrj<-cbind(taxgj, otugj)

#Tax-level-wise
tabrj<-tabrj[,-1:-2]
tabrj<-tabrj[,-2:-5]
zj<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(cladesj))
{
gtrj<-subset(tabrj, Class==cladesj[i])
richj<-colSums(gtrj[,-1] != 0)
zj<-cbind(zj, richj)
tj<-1+i
colnames(zj)[tj]<-cladesj[i]
}
zj<-zj[,-1]
combinedj<-cbind(zj, d)
combinedj$Total<-rowSums(combinedj[,1:length(cladesj)],)

#Rich
no<-subset_taxa(tudao, Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedj$Total[match(combined$sshc, combinedj$sshc)]
combined$Total_top<-rowSums(combined[,topclades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=topclades)
br$rel_rich<-br$value/br$Total

#Freq
####get taxa abund
datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, Division=="Metazoa")
datag = tax_glom(no, "Class")
p_topcladesseq<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
topcladesseq<-rownames(p_topcladesseq[p_topcladesseq$Class %in% topclades,])
data9 = prune_taxa(topcladesseq, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)

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

ggplot(dg, aes(log10(rel_rich), log10(rel_freq), color=habitat, shape=substrate_type)) + 
geom_point(size=1, stroke=0.2, alpha=0.6) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class Metazoa - Abundance vs. Richness", x ="log10 Relative Richness", y = "log10 Relative Abundance")
ggsave("Class_m_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=season)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_m_Rel_abund_rich_season_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_m_Rel_abund_rich_habitat_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=season), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_m_Rel_abund_rich_habitat_season_ncbi.pdf")


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
labs(title="Class Metazoa - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Class_m_Rel_abund_rich_st_season_mean_cross_ncbi.pdf")


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
labs(title="Class Metazoa - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Class_m_Rel_abund_rich_habitat_st_mean_cross_ncbi.pdf")


#Class non-Metazoa
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

#Define groups - top rich != top abund
topclades<-c("Bacillariophyta","Dinophyceae","Filosa-Thecofilosea","Spirotrichea","Oligohymenophorea","Labyrinthulomycetes","Oomycota","Litostomatea","Syndiniales","Phaeophyceae","Mamiellophyceae","Cryptophyceae")

#Rich
no<-subset_taxa(tudao, !Division=="Metazoa")
taxg<-data.frame(tax_table(no), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Class))
d<-data.frame(sample_data(no))
ch<-d$sshc
otug<-otu_table(no)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#Tax-level-wise
tabr<-tabr[,-1:-2]
tabr<-tabr[,-2:-5]
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, Class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}
z<-z[,-1]
combined<-cbind(z, d)
combined$Total<-combinedj$Total[match(combined$sshc, combinedj$sshc)]
combined$Total_top<-rowSums(combined[,topclades],)
combined$Other<-combined$Total-combined$Total_top
design<-colnames(d)
br<-melt(combined, id=c(design, "Total"), measure=topclades)
br$rel_rich<-br$value/br$Total

#Freq
####get taxa abund
datarg = transform_sample_counts(tudao, function(x) x/sum(x))
no<-subset_taxa(datarg, !Division=="Metazoa")
datag = tax_glom(no, "Class")
p_topcladesseq<-data.frame(tax_table(datag), stringsAsFactors=FALSE)
topcladesseq<-rownames(p_topcladesseq[p_topcladesseq$Class %in% topclades,])
data9 = prune_taxa(topcladesseq, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Class
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)

#Careful with hyphen-separated taxa names
colnames(ttab)[8]<-"Filosa-Thecofilosea"

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

ggplot(dg, aes(log10(rel_rich), log10(rel_freq), color=habitat, shape=substrate_type)) + 
geom_point(size=1, stroke=0.2, alpha=0.6) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class Metazoa - Abundance vs. Richness", x ="log10 Relative Richness", y = "log10 Relative Abundance")
ggsave("Class_nom_Rel_abund_rich_Clusters_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=season)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class non-Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_nom_Rel_abund_rich_season_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=substrate_type), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class non-Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_nom_Rel_abund_rich_habitat_st_ncbi.pdf")

ggplot(dg, aes(rel_rich, rel_freq, color=habitat)) + 
stat_ellipse(aes(lty=season), type = "norm", size=0.5, level=.99) +
facet_wrap(~variable, ncol=3, scales="free") + theme_bw() + 
theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
theme(legend.key.size = unit(0.3, "cm")) +
labs(title="Class non-Metazoa - Abundance vs. Richness", x ="Relative Richness", y = "Relative Abundance")
ggsave("Class_nom_Rel_abund_rich_habitat_season_ncbi.pdf")


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
labs(title="Class non-Metazoa - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Class_nom_Rel_abund_rich_st_season_mean_cross_ncbi.pdf")


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
labs(title="Class non-Metazoa - Abundance vs. Richness", x ="Mean Relative Richness per Cluster", y = "Mean Relative Abundance per Cluster") +
theme(axis.title=element_text(size=7), plot.title = element_text(size=8), axis.text.x = element_text(size=6, angle=90, hjust=1), axis.text.y = element_text(size=6))
ggsave("Class_nom_Rel_abund_rich_habitat_st_mean_cross_ncbi.pdf")


######
#Prevalence
######
#Division
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

##Calculate prev.
taxg<-data.frame(tax_table(tudao), stringsAsFactors=FALSE)
clades<-levels(factor(taxg$Division))
d<-data.frame(sample_data(tudao))
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
tabr<-cbind(taxg, otug)

#prev.
recseq<-expand.grid("ASV"=rownames(otug), stringsAsFactors=FALSE)
d$ssth<-paste(d$season, d$substrate_type, d$habitat, sep="_")
d$sth<-paste(d$substrate_type, d$habitat, sep="_")
hbtao<-levels(factor(d$ssth))

for (i in 1:length(hbtao))
{
hbx<-unique(subset(d, ssth==hbtao[i])$sshc)
otughbr<-rowSums(otug[,hbx] != 0)/length(hbx)
otughbr2<-rowSums(otug[,hbx])/sum(rowSums(otug[,hbx]))
recseq<-cbind(recseq, otughbr, otughbr2)
t<-2*i
colnames(recseq)[t:(t+1)]<-c(hbtao[i], paste(hbtao[i],"rel",sep="_"))
}

recseq<-recseq[,-1]
combinedz<-cbind(recseq, taxg)
combinedz$ASV<-rownames(combinedz)
fseq1<-melt(combinedz, id=c("Division","ASV"), measure=hbtao)
colnames(fseq1)<-c("Division", "ASV", "prev_id", "prev")
fseq2<-melt(combinedz, id=c("Division", "ASV"), measure=paste(hbtao,"rel",sep="_"))
colnames(fseq2)<-c("Division", "ASV", "rel_id", "rel_freq")
fseq1$track<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$rel_id<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")
fseq2$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")

ffseq<-merge(fseq1, fseq2, by="tutu")
identical(ffseq$Division.x,ffseq$Division.y)
prev_data<-ffseq[,c("tutu", "Division.x", "ASV.x", "prev_id", "prev", "rel_freq")]
colnames(prev_data)<-c("ASV_group", "Division", "ASV", "group", "prev", "rel_freq")

#optional?
prev_data2<-subset(prev_data, prev>0)
#
prev_data3<-subset(prev_data, rel_freq>0)
prev_data4<-subset(prev_data3, !Division=="Excavata")
prev_data5<-subset(prev_data4, !Division=="Amoebozoa")

prev_data5$season<-d$season[match(prev_data5$group, d$ssth)]
prev_data5$sth<-d$sth[match(prev_data5$group, d$ssth)]

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_prev")

prev_summaryf <- ddply(prev_data5, c("Division","season","sth"), transform,
               l10_freq = round(log10(rel_freq), digits=2)
)
prev_summaryf$uyt<-paste(prev_summaryf$Division, prev_summaryf$group, sep="_")

ggplot(prev_summaryf, aes(prev, l10_freq, color=season, group=uyt)) + 
geom_point(size=0.4, stroke=0.01, alpha=0.5, shape=1, position=position_dodge(width = 0.05)) +
geom_smooth(method=loess, se=FALSE, size=0.3, alpha=0.8) +
facet_wrap(Division~sth, ncol=6) + 
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"), legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0.05)) +
labs(title="Division - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season")
ggsave("Division_Rel_abund_Prev_ASVs_ncbi.pdf")


#mean_SD
prev_summary <- ddply(prev_data5, c("Division","season","sth", "prev"), summarize,
               mean_freq = mean(log10(rel_freq)),
               sd_freq   = sd(log10(rel_freq))
)

ggplot(prev_summary, aes(prev, mean_freq)) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = season), position = position_dodge(0.2), width = 0.05, size=0.2, alpha=0.8) +
geom_point(aes(color = season), size=0.4, stroke=0.1, alpha=0.8, shape=21, position = position_dodge(0.2)) +
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
facet_wrap(Division~sth, ncol=6) + 
labs(title="Division - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season") +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), axis.ticks.length=unit(.04, "cm"), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"),legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0))
ggsave("Division_Rel_abund_Prev_mean_sd_ncbi.pdf")

#Division - non Opisthokonta

#Prepare merged object
tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sshc<-rownames(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc")]

##Calculate prev.
#save raw table for totals
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))

#get top6 taxa
no<-subset_taxa(tudao, !Division=="Opisthokonta")
TopNOTUs <- sort(tapply(taxa_sums(no), tax_table(no)[, "Division"], sum), TRUE)[1:6]
data6 <- subset_taxa(no, Division %in% names(TopNOTUs))

taxg6<-data.frame(tax_table(data6), stringsAsFactors=FALSE)
d<-data.frame(sample_data(tudao))
otug6<-otu_table(data6)
otug6<-t(data.frame(otug6, check.names=F))
tabr6<-cbind(taxg6, otug6)

#prev.
recseq6<-expand.grid("ASV"=rownames(otug6), stringsAsFactors=FALSE)
d$ssth<-paste(d$season, d$substrate_type, d$habitat, sep="_")
d$sth<-paste(d$substrate_type, d$habitat, sep="_")
hbtao<-levels(factor(d$ssth))

for (i in 1:length(hbtao))
{
hbx<-unique(subset(d, ssth==hbtao[i])$sshc)
otughbr6<-rowSums(otug6[,hbx] != 0)/length(hbx)
otughbr62<-rowSums(otug6[,hbx])/sum(rowSums(otug[,hbx]))
recseq6<-cbind(recseq6, otughbr6, otughbr62)
t<-2*i
colnames(recseq6)[t:(t+1)]<-c(hbtao[i], paste(hbtao[i],"rel",sep="_"))
}

recseq6<-recseq6[,-1]
combinedz<-cbind(recseq6, taxg6)
combinedz$ASV<-rownames(combinedz)
fseq1<-melt(combinedz, id=c("Division","ASV"), measure=hbtao)
colnames(fseq1)<-c("Division", "ASV", "prev_id", "prev")
fseq2<-melt(combinedz, id=c("Division", "ASV"), measure=paste(hbtao,"rel",sep="_"))
colnames(fseq2)<-c("Division", "ASV", "rel_id", "rel_freq")
fseq1$track<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$rel_id<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")
fseq2$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")

ffseq<-merge(fseq1, fseq2, by="tutu")
identical(ffseq$Division.x,ffseq$Division.y)
prev_data<-ffseq[,c("tutu", "Division.x", "ASV.x", "prev_id", "prev", "rel_freq")]
colnames(prev_data)<-c("ASV_group", "Division", "ASV", "group", "prev", "rel_freq")

#optional?
prev_data2<-subset(prev_data, prev>0)
#
prev_data3<-subset(prev_data, rel_freq>0)

prev_data3$season<-d$season[match(prev_data3$group, d$ssth)]
prev_data3$sth<-d$sth[match(prev_data3$group, d$ssth)]

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_prev")

prev_summaryf <- ddply(prev_data3, c("Division","season","sth"), transform,
               l10_freq = round(log10(rel_freq), digits=2)
)
prev_summaryf$uyt<-paste(prev_summaryf$Division, prev_summaryf$group, sep="_")

ggplot(prev_summaryf, aes(prev, l10_freq, color=season, group=uyt)) + 
geom_point(size=0.4, stroke=0.01, alpha=0.5, shape=1, position=position_dodge(width = 0.05)) +
geom_smooth(method=loess, se=FALSE, size=0.3, alpha=0.8) +
facet_wrap(Division~sth, ncol=6) + 
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"), legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0.05)) +
labs(title="Division (non Opisthokonta) - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season")
ggsave("Division_nop_Rel_abund_Prev_ASVs_ncbi.pdf")


#mean_SD
prev_summary <- ddply(prev_data3, c("Division","season","sth", "prev"), summarize,
               mean_freq = mean(log10(rel_freq)),
               sd_freq   = sd(log10(rel_freq))
)

ggplot(prev_summary, aes(prev, mean_freq)) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = season), position = position_dodge(0.2), width = 0.05, size=0.2, alpha=0.8) +
geom_point(aes(color = season), size=0.4, stroke=0.1, alpha=0.8, shape=21, position = position_dodge(0.2)) +
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
facet_wrap(Division~sth, ncol=6) + 
labs(title="Division (non Opisthokonta) - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season") +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), axis.ticks.length=unit(.04, "cm"), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"),legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0))
ggsave("Division_nop_Rel_abund_Prev_mean_sd_ncbi.pdf")

#Metazoan

#Prepare merged object
tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sshc<-rownames(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc")]

##Calculate prev.
#save raw table for totals
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))

#get top6 taxa
no<-subset_taxa(tudao, Division=="Metazoa")
TopNOTUs <- sort(tapply(taxa_sums(no), tax_table(no)[, "Class"], sum), TRUE)[1:6]
data6 <- subset_taxa(no, Class %in% names(TopNOTUs))

taxg6<-data.frame(tax_table(data6), stringsAsFactors=FALSE)
d<-data.frame(sample_data(tudao))
otug6<-otu_table(data6)
otug6<-t(data.frame(otug6, check.names=F))
tabr6<-cbind(taxg6, otug6)

#prev.
recseq6<-expand.grid("ASV"=rownames(otug6), stringsAsFactors=FALSE)
d$ssth<-paste(d$season, d$substrate_type, d$habitat, sep="_")
d$sth<-paste(d$substrate_type, d$habitat, sep="_")
hbtao<-levels(factor(d$ssth))

for (i in 1:length(hbtao))
{
hbx<-unique(subset(d, ssth==hbtao[i])$sshc)
otughbr6<-rowSums(otug6[,hbx] != 0)/length(hbx)
otughbr62<-rowSums(otug6[,hbx])/sum(rowSums(otug[,hbx]))
recseq6<-cbind(recseq6, otughbr6, otughbr62)
t<-2*i
colnames(recseq6)[t:(t+1)]<-c(hbtao[i], paste(hbtao[i],"rel",sep="_"))
}

recseq6<-recseq6[,-1]
combinedz<-cbind(recseq6, taxg6)
combinedz$ASV<-rownames(combinedz)
fseq1<-melt(combinedz, id=c("Class","ASV"), measure=hbtao)
colnames(fseq1)<-c("Class", "ASV", "prev_id", "prev")
fseq2<-melt(combinedz, id=c("Class", "ASV"), measure=paste(hbtao,"rel",sep="_"))
colnames(fseq2)<-c("Class", "ASV", "rel_id", "rel_freq")
fseq1$track<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$rel_id<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")
fseq2$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")

ffseq<-merge(fseq1, fseq2, by="tutu")
identical(ffseq$Class.x,ffseq$Class.y)
prev_data<-ffseq[,c("tutu", "Class.x", "ASV.x", "prev_id", "prev", "rel_freq")]
colnames(prev_data)<-c("ASV_group", "Class", "ASV", "group", "prev", "rel_freq")

#optional?
prev_data2<-subset(prev_data, prev>0)
#
prev_data3<-subset(prev_data, rel_freq>0)

prev_data3$season<-d$season[match(prev_data3$group, d$ssth)]
prev_data3$sth<-d$sth[match(prev_data3$group, d$ssth)]

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_prev")

prev_summaryf <- ddply(prev_data3, c("Class","season","sth"), transform,
               l10_freq = round(log10(rel_freq), digits=2)
)
prev_summaryf$uyt<-paste(prev_summaryf$Class, prev_summaryf$group, sep="_")

ggplot(prev_summaryf, aes(prev, l10_freq, color=season, group=uyt)) + 
geom_point(size=0.4, stroke=0.01, alpha=0.5, shape=1, position=position_dodge(width = 0.05)) +
geom_smooth(method=loess, se=FALSE, size=0.3, alpha=0.8) +
facet_wrap(Class~sth, ncol=6) + 
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"), legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0.05)) +
labs(title="Class Metazoa - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season")
ggsave("Class_m_Rel_abund_Prev_ASVs_ncbi.pdf")


#mean_SD
prev_summary <- ddply(prev_data3, c("Class","season","sth", "prev"), summarize,
               mean_freq = mean(log10(rel_freq)),
               sd_freq   = sd(log10(rel_freq))
)

ggplot(prev_summary, aes(prev, mean_freq)) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = season), position = position_dodge(0.2), width = 0.05, size=0.2, alpha=0.8) +
geom_point(aes(color = season), size=0.4, stroke=0.1, alpha=0.8, shape=21, position = position_dodge(0.2)) +
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
facet_wrap(Class~sth, ncol=6) + 
labs(title="Class Metazoa - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season") +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), axis.ticks.length=unit(.04, "cm"), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"),legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0))
ggsave("Class_m_Rel_abund_Prev_mean_sd_ncbi.pdf")



#non-Metazoan

#Prepare merged object
tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sshc<-rownames(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc")]

##Calculate prev.
#save raw table for totals
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))

#get top6 taxa
no<-subset_taxa(tudao, !Division=="Metazoa")
TopNOTUs <- sort(tapply(taxa_sums(no), tax_table(no)[, "Class"], sum), TRUE)[1:6]
data6 <- subset_taxa(no, Class %in% names(TopNOTUs))

taxg6<-data.frame(tax_table(data6), stringsAsFactors=FALSE)
d<-data.frame(sample_data(tudao))
otug6<-otu_table(data6)
otug6<-t(data.frame(otug6, check.names=F))
tabr6<-cbind(taxg6, otug6)

#prev.
recseq6<-expand.grid("ASV"=rownames(otug6), stringsAsFactors=FALSE)
d$ssth<-paste(d$season, d$substrate_type, d$habitat, sep="_")
d$sth<-paste(d$substrate_type, d$habitat, sep="_")
hbtao<-levels(factor(d$ssth))

for (i in 1:length(hbtao))
{
hbx<-unique(subset(d, ssth==hbtao[i])$sshc)
otughbr6<-rowSums(otug6[,hbx] != 0)/length(hbx)
otughbr62<-rowSums(otug6[,hbx])/sum(rowSums(otug[,hbx]))
recseq6<-cbind(recseq6, otughbr6, otughbr62)
t<-2*i
colnames(recseq6)[t:(t+1)]<-c(hbtao[i], paste(hbtao[i],"rel",sep="_"))
}

recseq6<-recseq6[,-1]
combinedz<-cbind(recseq6, taxg6)
combinedz$ASV<-rownames(combinedz)
fseq1<-melt(combinedz, id=c("Class","ASV"), measure=hbtao)
colnames(fseq1)<-c("Class", "ASV", "prev_id", "prev")
fseq2<-melt(combinedz, id=c("Class", "ASV"), measure=paste(hbtao,"rel",sep="_"))
colnames(fseq2)<-c("Class", "ASV", "rel_id", "rel_freq")
fseq1$track<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$rel_id<-paste(fseq1$prev_id,"rel",sep="_")
fseq1$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")
fseq2$tutu<-paste(fseq1$ASV,fseq1$rel_id,sep="_")

ffseq<-merge(fseq1, fseq2, by="tutu")
identical(ffseq$Class.x,ffseq$Class.y)
prev_data<-ffseq[,c("tutu", "Class.x", "ASV.x", "prev_id", "prev", "rel_freq")]
colnames(prev_data)<-c("ASV_group", "Class", "ASV", "group", "prev", "rel_freq")

#optional?
prev_data2<-subset(prev_data, prev>0)
#
prev_data3<-subset(prev_data, rel_freq>0)

prev_data3$season<-d$season[match(prev_data3$group, d$ssth)]
prev_data3$sth<-d$sth[match(prev_data3$group, d$ssth)]

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/tax/ncbi/abund_prev")

prev_summaryf <- ddply(prev_data3, c("Class","season","sth"), transform,
               l10_freq = round(log10(rel_freq), digits=2)
)
prev_summaryf$uyt<-paste(prev_summaryf$Class, prev_summaryf$group, sep="_")

ggplot(prev_summaryf, aes(prev, l10_freq, color=season, group=uyt)) + 
geom_point(size=0.4, stroke=0.01, alpha=0.5, shape=1, position=position_dodge(width = 0.05)) +
geom_smooth(method=loess, se=FALSE, size=0.3, alpha=0.8) +
facet_wrap(Class~sth, ncol=6) + 
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"), legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0.05)) +
labs(title="Class non-Metazoa - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season")
ggsave("Class_nom_Rel_abund_Prev_ASVs_ncbi.pdf")


#mean_SD
prev_summary <- ddply(prev_data3, c("Class","season","sth", "prev"), summarize,
               mean_freq = mean(log10(rel_freq)),
               sd_freq   = sd(log10(rel_freq))
)

ggplot(prev_summary, aes(prev, mean_freq)) +
geom_errorbar(aes(ymin = mean_freq-sd_freq, ymax = mean_freq+sd_freq, color = season), position = position_dodge(0.2), width = 0.05, size=0.2, alpha=0.8) +
geom_point(aes(color = season), size=0.4, stroke=0.1, alpha=0.8, shape=21, position = position_dodge(0.2)) +
geom_hline(yintercept=-3, linetype="dashed", color = "black", size=0.2) +
facet_wrap(Class~sth, ncol=6) + 
labs(title="Class non-Metazoa - ASVs Relative abundance vs. Prevalence", x ="Prevalence", y = "log10 ASV Relative Abundance", fill = "Season") +
theme_bw() + 
theme(axis.title=element_text(size=6), plot.title = element_text(size=7), axis.text.x = element_text(size=5, angle=90, hjust=1), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.03,0,0.03,0, "cm")), strip.text = element_text(size=5), legend.title=element_text(size=5), axis.ticks.length=unit(.04, "cm"), legend.text=element_text(size=5)) +
theme(legend.key.size = unit(0.3, "cm"),legend.position="bottom") +
scale_x_continuous(limits=c(0, 1), expand = c(0, 0))
ggsave("Class_nom_Rel_abund_Prev_mean_sd_ncbi.pdf")



