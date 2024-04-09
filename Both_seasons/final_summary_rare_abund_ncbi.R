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

