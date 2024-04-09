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
otu_mat<-as.matrix(read.table("cleaned_tax_ncbi_12_12_22.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("cleaned_ncbi_12_12_22.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_mat_b)
p_ncbi = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("cleaned_ncbi_metadata_12_12_22.txt", sep="\t", header=T)

##SITE INFO
c_s<-read.table("cluster_site.txt", sep="\t", header=T)
metadata$Location<-c_s$Site_name[match(metadata$cluster, c_s$cluster)]
metadata$cluster<-as.integer(metadata$cluster)
metadata$cl_se<-as.character(paste(metadata$cluster,metadata$season,sep="_"))
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_ncbi, sampledata)
DADAwang1

##Make histogram for raw read counts
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022")

readsi<-sample_sums(DADAwang1)
combinedi<-cbind(readsi, sample_data(DADAwang1))
combinedi<-data.frame(combinedi)

threshold<-round(quantile(combinedi$readsi, c(.05)), digits=-3)
combinedi$q025<-combinedi$readsi>threshold

mui <- ddply(combinedi, .(season, substrate_type), summarise, grp.mean=mean(readsi))

ggplot(combinedi, aes(x=readsi)) +
geom_histogram(aes(fill=q025), position="identity", alpha=0.6, binwidth=8000) + geom_density(alpha=0.6) + geom_vline(data=mui, aes(xintercept=grp.mean), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(substrate_type ~ season, ncol=2, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Reads histogram plot", x ="Reads", y = "Count", fill = paste("Total reads > ",threshold,sep="")) + scale_x_continuous(limits=c(0,1000000)) + scale_y_continuous(limits=c(0,50))
ggsave("reads_hist_raw_ncbi.pdf")

##Make histogram for merged read counts (~FIELD_replicate)
#prepare data
tudao<-merge_samples(DADAwang1, "sample_root")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

head(d)
sample_data(tudao)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

#make histogram

reads<-sample_sums(tudao)
combined<-cbind(reads, sample_data(tudao))
head(combined)
combined<-data.frame(combined)

threshold<-round(quantile(combined$reads, c(.05)), digits=-3)
combined$q025<-combined$reads>threshold

mu <- ddply(combined, .(season, substrate_type), summarise, grp.mean=mean(reads))

ggplot(combined, aes(x=reads)) +
geom_histogram(aes(fill=q025), position="identity", alpha=0.6, binwidth=50000)+
geom_density(alpha=0.6)+
geom_vline(data=mu, aes(xintercept=grp.mean),
           linetype="dashed")+
labs(title="Reads histogram plot", x ="Reads", y = "Count", fill = paste("Total reads > ",threshold,sep=""))+
theme_classic() + facet_wrap(substrate_type ~ season, ncol=2, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + scale_x_continuous(limits=c(0,3000000),labels = comma) + scale_y_continuous(limits=c(0,25),labels = comma)
ggsave("reads_hist_merged_ncbi.pdf")

###Make sediment read summary per cluster and per field_replicate
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment")

###Although water and sediment plots are built separately, 
#this first part addresses the whole dataset,
#aiming for the comparison between different substrate types

###Sediment PCR rarefaction - diagnosis

#Use threshold of whole dataset
threshold<-round(quantile(combinedi$readsi, c(.05)), digits=-3)
combinedi$q025<-combinedi$readsi>threshold

xreads_sed <- ddply(subset(combinedi, substrate_type=="sediment"), .(sample_root, season, cluster, habitat), summarize, Total_replicates = length(unique(readsi)), Total_overq025 = length(unique(readsi[q025==T])), Total_reads = sum(readsi), Min_reads = min(readsi), Mean_reads = mean(readsi), SD_reads = sd(readsi))

xreads_sed$hb_season<-paste(xreads_sed$season, xreads_sed$habitat, sep="_")
xreads_sed$hb_season <- factor(xreads_sed$hb_season, levels = c("spring_EB", "autumn_EB", "spring_RB", "autumn_RB","spring_SB","autumn_SB"))

ggplot(data = data.frame(xreads_sed), aes(hb_season, Total_reads, fill = as.factor(Total_overq025))) +  geom_bar(stat="identity", position = position_dodge2(), color="black", size=0.02) + geom_text(aes(label = round(Min_reads, digits=0), y = Total_reads + 500000), position = position_dodge2(1), angle=90, color="black", size=0.9) + facet_wrap(~cluster, ncol=4) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
    theme(legend.key.size = unit(0.3, "cm")) + scale_y_continuous(limits=c(0,3300000),labels = comma) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Sediment - reads per field replicate",
        x ="Field Replicates", y = "Reads", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("reads_sed_cluster_ncbi.pdf")

zcom<-ddply(combinedi, .(sample_root), transform, Total_overq025=length(unique(sample_ID[q025==T])))
dastr<-subset(zcom, !Total_overq025==4)
ruins<-levels(factor(dastr$sample_root))
xpmo<-combinedi[combinedi$sample_root %in% ruins, ]

xreadsi <- subset(xpmo, substrate_type=="sediment")
xreadsi$hb_season<-paste(xreadsi$season, xreadsi$habitat, sep="_")
xreadsi$hb_season <- factor(xreadsi$hb_season, levels = c("spring_EB", "autumn_EB", "spring_RB", "autumn_RB","spring_SB","autumn_SB"))

xreadsi2  <- ddply(xreadsi, .(sample_root, cl_se), transform, Total_overq025=length(unique(sample_ID[q025==T])))

xreadsi3<-xreadsi2[ order(xreadsi2$cluster, xreadsi2$sample_root), ]
xreadsi3$sample_root <- factor(xreadsi3$sample_root, levels = unique(xreadsi3$sample_root[order(xreadsi3$cluster)]))

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red", size=0.2) + facet_wrap(~sample_root, ncol=8, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Sediment - field replicates below threshold", x ="PCR replicate", y = "number of reads (millions)", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("sed_summary_q025_ncbi.pdf")

##alternative threshold (frtotal>=50000)
threshold<-50000
combinedi$q025<-combinedi$readsi>threshold

zcom<-ddply(combinedi, .(sample_root), transform, Total_overq025=length(unique(sample_ID[q025==T])))
dastr<-subset(zcom, !Total_overq025==4)
ruins<-levels(factor(dastr$sample_root))
xpmo<-combinedi[combinedi$sample_root %in% ruins, ]

xreadsi <- subset(xpmo, substrate_type=="sediment")
xreadsi$hb_season<-paste(xreadsi$season, xreadsi$habitat, sep="_")
xreadsi$hb_season <- factor(xreadsi$hb_season, levels = c("spring_EB", "autumn_EB", "spring_RB", "autumn_RB","spring_SB","autumn_SB"))

xreadsi2  <- ddply(xreadsi, .(sample_root, cl_se), transform, Total_overq025=length(unique(sample_ID[q025==T])))

xreadsi3<-xreadsi2[ order(xreadsi2$cluster, xreadsi2$sample_root), ]
xreadsi3$sample_root <- factor(xreadsi3$sample_root, levels = unique(xreadsi3$sample_root[order(xreadsi3$cluster)]))

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red", size=0.2) + facet_wrap(~sample_root, ncol=10, scales = "free_x") + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
labels = trans_format("log10", math_format(10^.x))) + theme_bw()  + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Sediment - field replicates below threshold", x ="PCR replicate", y = "number of reads", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("alt_sed_summary_50000_ncbi.pdf")


setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water")

###Water PCR rarefaction - diagnosis

#Use threshold of whole dataset
threshold<-round(quantile(combinedi$readsi, c(.05)), digits=-3)
combinedi$q025<-combinedi$readsi>threshold

xreads_sed <- ddply(subset(combinedi, substrate_type=="water"), .(sample_root, season, cluster, habitat), summarize, Total_replicates = length(unique(readsi)), Total_overq025 = length(unique(readsi[q025==T])), Total_reads = sum(readsi), Min_reads = min(readsi), Mean_reads = mean(readsi), SD_reads = sd(readsi))

xreads_sed$hb_season<-paste(xreads_sed$season, xreads_sed$habitat, sep="_")
xreads_sed$hb_season <- factor(xreads_sed$hb_season, levels = c("spring_EW", "autumn_EW", "spring_RW", "autumn_RW","spring_SW","autumn_SW"))

ggplot(data = data.frame(xreads_sed), aes(hb_season, Total_reads, fill = as.factor(Total_overq025))) +  geom_bar(stat="identity", position = position_dodge2(), color="black", size=0.02) + geom_text(aes(label = round(Min_reads, digits=0), y = Total_reads + 650000), position = position_dodge2(1), angle=90, color="black", size=0.9) + facet_wrap(~cluster, ncol=4) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +
    theme(legend.key.size = unit(0.3, "cm")) + scale_y_continuous(limits=c(0,5000000),labels = comma) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Water - reads per field replicate",
        x ="Field Replicates", y = "Reads", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("reads_water_cluster_ncbi.pdf")

zcom<-ddply(combinedi, .(sample_root), transform, Total_overq025=length(unique(sample_ID[q025==T])))
dastr<-subset(zcom, !Total_overq025==4)
ruins<-levels(factor(dastr$sample_root))
xpmo<-combinedi[combinedi$sample_root %in% ruins, ]

xreadsi <- subset(xpmo, substrate_type=="water")
xreadsi$hb_season<-paste(xreadsi$season, xreadsi$habitat, sep="_")
xreadsi$hb_season <- factor(xreadsi$hb_season, levels = c("spring_EW", "autumn_EW", "spring_RW", "autumn_RW","spring_SW","autumn_SW"))

xreadsi2  <- ddply(xreadsi, .(sample_root, cl_se), transform, Total_overq025=length(unique(sample_ID[q025==T])))

xreadsi3<-xreadsi2[ order(xreadsi2$cluster, xreadsi2$sample_root), ]
xreadsi3$sample_root <- factor(xreadsi3$sample_root, levels = unique(xreadsi3$sample_root[order(xreadsi3$cluster)]))

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red", size=0.2) + facet_wrap(~sample_root, ncol=7, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Water - field replicates below threshold", x ="PCR replicate", y = "number of reads (millions)", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("water_summary_q025_ncbi.pdf")

#alternatively

readsi<-sample_sums(DADAwang1)
combinedi<-cbind(readsi, sample_data(DADAwang1))
combinedi<-data.frame(combinedi)

threshold<-50000
combinedi$q025<-combinedi$readsi>threshold

zcom<-ddply(combinedi, .(sample_root), transform, Total_overq025=length(unique(sample_ID[q025==T])))
dastr<-subset(zcom, !Total_overq025==4)
ruins<-levels(factor(dastr$sample_root))
xpmo<-combinedi[combinedi$sample_root %in% ruins, ]

xreadsi <- subset(xpmo, substrate_type=="water")
xreadsi$hb_season<-paste(xreadsi$season, xreadsi$habitat, sep="_")
xreadsi$hb_season <- factor(xreadsi$hb_season, levels = c("spring_EW", "autumn_EW", "spring_RW", "autumn_RW","spring_SW","autumn_SW"))

xreadsi2  <- ddply(xreadsi, .(sample_root, cl_se), transform, Total_overq025=length(unique(sample_ID[q025==T])))

xreadsi3<-xreadsi2[ order(xreadsi2$cluster, xreadsi2$sample_root), ]
xreadsi3$sample_root <- factor(xreadsi3$sample_root, levels = unique(xreadsi3$sample_root[order(xreadsi3$cluster)]))

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red", size=0.2) + facet_wrap(~sample_root, ncol=10, scales = "free_x") + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
labels = trans_format("log10", math_format(10^.x))) + theme_bw()  + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Water - field replicates below threshold", x ="PCR replicate", y = "number of reads", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("alt_wat_summary_50000_ncbi.pdf")

###### Make rarefaction plots to support the thresholds

#WHOLE DATASET
frare<-DADAwang1
sample_data(frare)$fp<-paste(as.character(sample_data(frare)$field_replicate),
as.character(sample_data(frare)$PCR_replicate),sep=".")

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = 500)
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

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/rare")

grupao<-levels(as.factor(sample_data(frare)$cluster))
for (e in 1:length(grupao))
{
wat10<-subset_samples(frare, cluster==grupao[e])
wat10 = filter_taxa(wat10, function(x) sum(x) > 0, TRUE)
wat101 = prune_samples(sample_sums(wat10)>0,wat10)
x<-t(data.frame(otu_table(wat101),check.names=F))
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
geom_text(data = labels, aes_string(x = "x", y = "y", label = "fp"), size = 2, hjust = 0) +
geom_line() + facet_wrap(season~habitat+substrate_type, ncol=2) + scale_x_continuous(labels = function(x) format(x, big.mark = ",",scientific = FALSE),breaks=seq(0,max(labels$x),10000)) + theme_bw() + theme(axis.text.x = element_text(angle=90, size=4), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5)) + theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste(grupao[e],"rare_ncbi.pdf",sep="_"))
}

##Make the same for merged read counts (~FIELD_replicate)
#prepare data
tudao<-merge_samples(DADAwang1, "sample_root")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

head(d)
sample_data(tudao)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = 5000)
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

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/rare/merged")

grupao<-levels(as.factor(sample_data(tudao)$cluster))
for (e in 1:length(grupao))
{
wat10<-subset_samples(tudao, cluster==grupao[e])
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
geom_text(data = labels, aes_string(x = "x", y = "y", label = "sample_root"), size = 2, hjust = 0) +
geom_line(size=0.2) +
geom_vline(xintercept = c(50000,100000,150000), linetype = "dashed", size=0.5, color="red") +
geom_errorbar(aes(ymin=ASVs-.se, ymax=ASVs+.se), width=.1) +
 facet_wrap(season~habitat+substrate_type, ncol=2) + scale_x_continuous(labels = function(x) format(x, big.mark = ",",scientific = FALSE),breaks=seq(0,max(labels$x),25000)) + theme_bw() + theme(axis.text.x = element_text(angle=90, size=4), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5)) + theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste(grupao[e],"rare_merged_ncbi.pdf",sep="_"))
}





######
###
#######
####Before calculating stats, get rid of obvious outliers
# (here, defined as pcr replicates having a single replicate, ASVs < 10 &| reads < 1000)

sd1<-data.frame(sample_data(DADAwang1))
otug<-t(data.frame(otu_table(DADAwang1), check.names=F))
sd1$reads<-sample_sums(DADAwang1)
rich<-rowSums(otug != 0)
sd1$rich<-rich[match(sd1$sample_ID, names(rich))]
sd1$qt<-sd1$rich<10|sd1$reads<1000
sd0<-ddply(sd1, .(sample_root), transform, q_pcreps=length(unique(sample_ID[qt==T])))

badsam<-unique(as.character(sd0[which(sd0$q_pcreps>2),"sample_root"]))
badsam
DADAwang1

#remove frs with a single pcr replicate - which for 18S must match C33RW3,C32RW3
DADAwang2<-subset_samples(DADAwang1, !(sample_root %in% badsam))

badsam2<-as.character(sd0[which(sd0$q_pcreps<=2&sd0$qt==T),"sample_ID"])
badsam2
DADAwang2

#then, remove other pcreps below threshold
DADAwang3<-subset_samples(DADAwang2, !(sample_ID %in% badsam2))
DADAwang3

DADAwang3 = filter_taxa(DADAwang3, function(x) sum(x) > 0, TRUE)
DADAwang3

#Just check if all samples have at least 3 pcreps
sd3<-data.frame(sample_data(DADAwang3))
otug3<-t(data.frame(otu_table(DADAwang3), check.names=F))
sd3$reads3<-sample_sums(DADAwang3)
rich3<-rowSums(otug3 != 0)
sd3$rich3<-rich3[match(sd3$sample_ID, names(rich3))]
sd3$qt<-sd3$rich3<10|sd3$reads3<1000
sd03<-ddply(sd3, .(sample_root), transform, q_pcreps=length(unique(sample_ID[qt==T])))

subset(sd03, q_pcreps>0)
sd03$q_pcreps

###Define threshold based on quantile 0.5%

thrsh<-.005

#merge and estimate total per field replicates
tudao<-merge_samples(DADAwang3, "sample_root")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

head(d)
sample_data(tudao)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

#

sd9<-data.frame(sample_data(tudao))
otu9<-data.frame(otu_table(tudao), check.names=F)
sd9$reads<-sample_sums(tudao)
rich<-rowSums(otu9 != 0)
sd9$rich<-rich[match(sd9$sample_root, names(rich))]

#Set threshold
qreads_S<-quantile(sd9[which(sd9$substrate_type=="sediment"),"reads"], thrsh)
qrich_S<-quantile(sd9[which(sd9$substrate_type=="sediment"),"rich"], thrsh)

qreads_W<-quantile(sd9[which(sd9$substrate_type=="water"),"reads"], thrsh)
qrich_W<-quantile(sd9[which(sd9$substrate_type=="water"),"rich"], thrsh)

#register thresholds

qreads_S
qrich_S

qreads_W
qrich_W

sd9$qt<-ifelse(sd9$substrate_type=="sediment", sd9$reads<=qreads_S|sd9$rich<=qrich_S, ifelse(sd9$substrate_type=="water", sd9$reads<=qreads_W|sd9$rich<=qrich_W, NA))

dsm<-unique(as.character(sd9[which(sd9$qt==T),"sample_root"]))

#register frs below threshold

dsm

######MANUAL EDITION - compare threshold qt==T with visual inspection of rarefaction plots

#prevent removal of "false outliers" - samples with observed true low rich/ worth keeping read totals
tlr<-c("2C20RW2","2C25EB3","2C25RB2","C29EW1","C29EW2","C30EW2")
#add missing relative outliers - samples with "false high rich"
fhr<-c("C20RW1","2C1SB3","2C32RB3","C29EB2")

fo<-dsm[! dsm %in% tlr]
fo2<-c(fo,fhr)

##Make sure all samples are outliers or has biased totals
DADAwang4<-subset_samples(DADAwang3, !(sample_root %in% fo2))
DADAwang4

DADAwang5 = filter_taxa(DADAwang4, function(x) sum(x) > 0, TRUE)
DADAwang5

sd1<-data.frame(sample_data(DADAwang5))
otug<-t(data.frame(otu_table(DADAwang5), check.names=F))
sd1$reads<-sample_sums(DADAwang5)
rich<-rowSums(otug != 0)
sd1$rich<-rich[match(sd1$sample_ID, names(rich))]
sd0<-ddply(sd1, .(sample_root), transform, fcr_tot=sum(reads), q_pcreps=length(unique(sample_ID)))

## define threshold (min fcr_tot ~ substrate_type)
ts2<-unique(min(sd0[which(sd0$substrate_type=="sediment"),"fcr_tot"]))
tw2<-unique(min(sd0[which(sd0$substrate_type=="water"),"fcr_tot"]))

sd0$qt<-ifelse(sd0$substrate_type=="sediment", sd0$reads<=ts2, ifelse(sd0$substrate_type=="water", sd0$reads<=qreads_W, NA))

## check thresholds, min reads, min rich, min pcreps, min fcr total
ts2
tw2
summary(sd0)

zcom<-sd0

#Set a parameter to define normalization
zcom$Norm_rule<-ifelse(zcom$qt==T, "r_1",
ifelse(zcom$qt==F&sd0$substrate_type=="sediment", "r_2",
ifelse(zcom$qt==F&sd0$substrate_type=="water", "r_3", NA)))

#update object's sample data
sample_data(DADAwang5)$final_rule<-zcom$Norm_rule[match(sample_data(DADAwang5)$sample_ID, zcom$sample_ID)]
steps<-levels(factor(zcom$Norm_rule))

#Check trash samples and save it
#ssss<-subset(zcom, is.na(zcom$Norm_rule))
#write.table(levels(factor(ssss$sample_ID)), "list_discarded_PCR_replicates_ncbi.txt", sep=" ")

#Normalize each category according to the associated sample sizes
steps<-c("r_2","r_3")

sizes<-c(ts2,tw2)
sizes2<-round(as.numeric(sizes),digits=0)
prtc2<-as.character(steps)
gaba<-data.frame(cbind(sizes2, prtc2), stringsAsFactors = FALSE)

for(e in 1:length(steps)) {
  assign(paste0("samples", gaba[e,2]), rarefy_even_depth(subset_samples(DADAwang5, final_rule==gaba[e,2]), sample.size=as.numeric(gaba[e,1]), replace=FALSE, rngseed= 13072021))
}

#merge objects
myobs<-paste0("samples",gaba$prtc2)
tudao_norm<-do.call(merge_phyloseq, mget(myobs), quote = FALSE)
tudao_unnorm<-subset_samples(DADAwang5, final_rule=="r_1")
tudao_p_merge<-merge_phyloseq(tudao_norm, tudao_unnorm)

tudao_p_merge
sample_sums(tudao_p_merge)

sample_data(tudao_p_merge)<-sample_data(tudao_p_merge)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]
#merge pcreps
tudao_pos_merge<-merge_samples(tudao_p_merge, "sample_root")

#Re-build sample_data
tudao0 = filter_taxa(tudao_pos_merge, function(x) sum(x) > 0, TRUE)
p_ncbi = prune_samples(sample_sums(tudao0)>0,tudao0)

d<-data.frame(sample_data(p_ncbi)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-ifelse(d$habitat==1|d$habitat==4, "eelgrass", ifelse(d$habitat==2|d$habitat==5,"rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

sample_data(p_ncbi)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

#Now normalize again according to substrate_type
sed_pm<-subset_samples(p_ncbi, substrate_type=="sediment")
wat_pm<-subset_samples(p_ncbi, substrate_type=="water")

stf<-min(sample_sums(sed_pm))
wtf<-min(sample_sums(wat_pm))

cas<-rarefy_even_depth(sed_pm, sample.size=stf, replace=FALSE, rngseed= 13072021)

caw<-rarefy_even_depth(wat_pm, sample.size=wtf, replace=FALSE, rngseed= 13072021)

ncbi<-merge_phyloseq(cas, caw)
ncbi2 = filter_taxa(ncbi, function(x) sum(x) > 0, TRUE)
p_ncbi2 = prune_samples(sample_sums(ncbi2)>0,ncbi2)

#Save final files

tax_m<-data.frame(tax_table(p_ncbi2))
otu_m<-data.frame(otu_table(p_ncbi2),check.names=F)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
write.table(data.frame(sample_data(p_ncbi2), check.names=F), "f_ncbi_metadata_5_01_22.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")
write.table(otu_m, "f_otu_ncbi_5_01_22.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "f_tax_ncbi_5_01_22.txt", sep="\t", quote=FALSE, row.names=TRUE)


###########
###
####
###
####

#####  ALL LINES BELOW RELATE TO SUPPLEMENTARY FIGURES NOT CURRENTLY USED FOR NORMALIZATION 
#####  IT INCLUDES CENTROID-BASED THRESHOLDS COMPARISON, PCR/FIELD REPLICATES ORDINATIONS, ETC

##Distance-based screening
###Currently, dists jaccard and bray are included separately 

##Extract, for each PCR rep (row):
#col1 -> distance to fr centroid (frC)
#col2 -> distance to cssh centroid (csshC)
#col3 -> average dist to other PCR reps
#col4 -> average dist to other fr PCR reps
#col5 -> average dist among other PCR reps
#col6 -> average dist among fr PCR reps
#col7 -> sd dist among other PCR reps
#col8 -> sd dist among fr PCR reps
#col9 -> read count
#col10 -> FR avg read count
#col11 -> cssh avg read count
#...also richness

####First, bray

####Do it separately for each cluster (cssh) (to increase resolution and reduce time ;)

#First, create the result collector df:
result_dist_b<-expand.grid(samples=sample_names(DADAwang3), d_to_frC=NA, d_to_csshC=NA, avg_d_to_pcreps=NA, avg_d_to_frpcreps=NA, avg_d_other_pcreps=NA, avg_d_other_frpcreps=NA, sd_d_other_pcreps=NA, sd_d_other_frpcreps=NA, reads=NA, fr_avg_reads=NA, cssh_avg_reads=NA, ASVs=NA, fr_avg_ASVs=NA, cssh_avg_ASVs=NA)
rownames(result_dist_b)<-result_dist_b$samples
result_dist_b<-result_dist_b[,-1]

sample_data(DADAwang3)$cssh<-paste(sample_data(DADAwang3)$cluster,sample_data(DADAwang3)$substrate_type,
sample_data(DADAwang3)$season,sample_data(DADAwang3)$habitat,sep="_")
allg<-unique(sample_data(DADAwang3)$cssh)

#Here, the distances between PCR replicates and read counts are summarized
#First, cluster-wise (~ [season + habitat + substrate_type])
#Check distance
for (csample in 1:length(allg))
{
#subset cluster
tcssh<-subset_samples(DADAwang3, cssh==allg[csample])
tcssh0 = filter_taxa(tcssh, function(x) sum(x) > 0, TRUE)
#mk data table
p_rt<-data.frame(sample_data(tcssh0))
p_rt$reads<-sample_sums(tcssh0)
#calculate cluster dist data
otug<-t(data.frame(otu_table(tcssh0), check.names=F))
dist_n<-vegdist(otug, method="bray", upper=T)
rich<-rowSums(otug != 0)
p_rt$rich<-rich[match(p_rt$sample_ID, names(rich))]
klhao<-cbind(as.data.frame(as.matrix(dist_n)),p_rt)
zin<-betadisper(dist_n, p_rt$field_replicate, type = "centroid")
zin2<-betadisper(dist_n, p_rt$cssh, type = "centroid")
onezin<-cbind(zin$distances, zin2$distances)
slist<-rownames(p_rt)
#### Loop for collecting stats (replicate-wise)
for (pcrsample in 1:length(slist))
{
result_dist_b[slist[pcrsample],1:2]<-as.numeric(onezin[slist[pcrsample],1:2])
ga_in<-subset(klhao, field_replicate==p_rt[slist[pcrsample],"field_replicate"])
ga_out<-subset(klhao, !field_replicate==p_rt[slist[pcrsample],"field_replicate"])
bool<-rownames(klhao)==slist[pcrsample]
result_dist_b[slist[pcrsample],3]<-as.numeric(mean(na.omit(ga_in[!bool,slist[pcrsample]])))
result_dist_b[slist[pcrsample],4]<-as.numeric(mean(na.omit(ga_out[!bool,slist[pcrsample]])))
vvk<-rownames(ga_in)
llk <- vvk[!vvk %in% slist[pcrsample]]
st<-as.matrix(na.omit(ga_in[llk,llk]))
st[lower.tri(st, diag = FALSE)]<-NA
stf<-st[st!=0]
result_dist_b[slist[pcrsample],5]<-as.numeric(mean(na.omit(stf)))
yyk<-rownames(ga_out)
wt<-as.matrix(na.omit(ga_out[,yyk]))
wt[lower.tri(wt, diag = FALSE)]<-NA
wtf<-wt[wt!=0]
result_dist_b[slist[pcrsample],6]<-as.numeric(mean(na.omit(wtf)))
result_dist_b[slist[pcrsample],7]<-as.numeric(sd(na.omit(stf)))
result_dist_b[slist[pcrsample],8]<-as.numeric(sd(na.omit(wtf)))
result_dist_b[slist[pcrsample],9]<-as.numeric(p_rt[slist[pcrsample],"reads"])
result_dist_b[slist[pcrsample],10]<-as.numeric(round(mean(ga_in[,"reads"])), digits=0)
result_dist_b[slist[pcrsample],11]<-as.numeric(round(mean(p_rt[,"reads"])), digits=0)
result_dist_b[slist[pcrsample],12]<-as.numeric(p_rt[slist[pcrsample],"rich"])
result_dist_b[slist[pcrsample],13]<-as.numeric(round(mean(ga_in[,"rich"])), digits=0)
result_dist_b[slist[pcrsample],14]<-as.numeric(round(mean(p_rt[,"rich"])), digits=0)
}
}


####Then, jaccard
#First, create the result collector df:
result_dist_j<-expand.grid(samples=sample_names(DADAwang3), d_to_frC=NA, d_to_csshC=NA, avg_d_to_pcreps=NA, avg_d_to_frpcreps=NA, avg_d_other_pcreps=NA, avg_d_other_frpcreps=NA, sd_d_other_pcreps=NA, sd_d_other_frpcreps=NA, reads=NA, fr_avg_reads=NA, cssh_avg_reads=NA,
ASVs=NA, fr_avg_ASVs=NA, cssh_avg_ASVs=NA)
rownames(result_dist_j)<-result_dist_j$samples
result_dist_j<-result_dist_j[,-1]

#Here, the distances between PCR replicates and read counts are summarized
#First, cluster-wise (~ [season + habitat + substrate_type])
#Check distance
for (csample in 1:length(allg))
{
#subset cluster
tcssh<-subset_samples(DADAwang3, cssh==allg[csample])
tcssh0 = filter_taxa(tcssh, function(x) sum(x) > 0, TRUE)
#mk data table
p_rt<-data.frame(sample_data(tcssh0))
p_rt$reads<-sample_sums(tcssh0)
#calculate cluster dist data
otug<-t(data.frame(otu_table(tcssh0), check.names=F))
dist_n<-vegdist(otug, method="jaccard", upper=T)
rich<-rowSums(otug != 0)
p_rt$rich<-rich[match(p_rt$sample_ID, names(rich))]
klhao<-cbind(as.data.frame(as.matrix(dist_n)),p_rt)
zin<-betadisper(dist_n, p_rt$field_replicate, type = "centroid")
zin2<-betadisper(dist_n, p_rt$cssh, type = "centroid")
onezin<-cbind(zin$distances, zin2$distances)
slist<-rownames(p_rt)
#### Loop for collecting stats (replicate-wise)
for (pcrsample in 1:length(slist))
{
result_dist_j[slist[pcrsample],1:2]<-as.numeric(onezin[slist[pcrsample],1:2])
ga_in<-subset(klhao, field_replicate==p_rt[slist[pcrsample],"field_replicate"])
ga_out<-subset(klhao, !field_replicate==p_rt[slist[pcrsample],"field_replicate"])
bool<-rownames(klhao)==slist[pcrsample]
result_dist_j[slist[pcrsample],3]<-as.numeric(mean(na.omit(ga_in[!bool,slist[pcrsample]])))
result_dist_j[slist[pcrsample],4]<-as.numeric(mean(na.omit(ga_out[!bool,slist[pcrsample]])))
vvk<-rownames(ga_in)
llk <- vvk[!vvk %in% slist[pcrsample]]
st<-as.matrix(na.omit(ga_in[llk,llk]))
st[lower.tri(st, diag = FALSE)]<-NA
stf<-st[st!=0]
result_dist_j[slist[pcrsample],5]<-as.numeric(mean(na.omit(stf)))
yyk<-rownames(ga_out)
wt<-as.matrix(na.omit(ga_out[,yyk]))
wt[lower.tri(wt, diag = FALSE)]<-NA
wtf<-wt[wt!=0]
result_dist_j[slist[pcrsample],6]<-as.numeric(mean(na.omit(wtf)))
result_dist_j[slist[pcrsample],7]<-as.numeric(sd(na.omit(stf)))
result_dist_j[slist[pcrsample],8]<-as.numeric(sd(na.omit(wtf)))
result_dist_j[slist[pcrsample],9]<-as.numeric(p_rt[slist[pcrsample],"reads"])
result_dist_j[slist[pcrsample],10]<-as.numeric(round(mean(ga_in[,"reads"])), digits=0)
result_dist_j[slist[pcrsample],11]<-as.numeric(round(mean(p_rt[,"reads"])), digits=0)
result_dist_j[slist[pcrsample],12]<-as.numeric(p_rt[slist[pcrsample],"rich"])
result_dist_j[slist[pcrsample],13]<-as.numeric(round(mean(ga_in[,"rich"])), digits=0)
result_dist_j[slist[pcrsample],14]<-as.numeric(round(mean(p_rt[,"rich"])), digits=0)
}
}

######Creating summary (within field rep, within cluster, and read-count)
####Combine it afterwards 

sddf<-data.frame(sample_data(DADAwang3))
result_dist_j$sample_ID<-rownames(result_dist_j)
result_dist_b$sample_ID<-rownames(result_dist_b)
result_dist_j$sample_root<-sddf$sample_root[match(rownames(result_dist_j), sddf$sample_ID)]
result_dist_b$sample_root<-sddf$sample_root[match(rownames(result_dist_b), sddf$sample_ID)]
result_dist_j$cssh<-sddf$cssh[match(rownames(result_dist_j), sddf$sample_ID)]
result_dist_b$cssh<-sddf$cssh[match(rownames(result_dist_b), sddf$sample_ID)]
result_dist_j$substrate_type<-sddf$substrate_type[match(rownames(result_dist_j), sddf$sample_ID)]
result_dist_b$substrate_type<-sddf$substrate_type[match(rownames(result_dist_b), sddf$sample_ID)]

#####Adding summaries for field replicate (distances & ASVs+reads)

aflb<-ddply(result_dist_b, .(sample_root), transform, frC=(d_to_frC-mean(d_to_frC))/sd(d_to_frC), 
dpcr=(avg_d_to_pcreps-mean(avg_d_to_pcreps))/sd(avg_d_to_pcreps), sdpcr=(sd_d_other_pcreps-mean(sd_d_other_pcreps))/sd(sd_d_other_pcreps), readsfr=round((reads-mean(reads))/sd(reads), digits=2), ASVsfr=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

aflj<-ddply(result_dist_j, .(sample_root), transform, frC=(d_to_frC-mean(d_to_frC))/sd(d_to_frC), 
dpcr=(avg_d_to_pcreps-mean(avg_d_to_pcreps))/sd(avg_d_to_pcreps), sdpcr=(sd_d_other_pcreps-mean(sd_d_other_pcreps))/sd(sd_d_other_pcreps), readsfr=round((reads-mean(reads))/sd(reads), digits=2), ASVsfr=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

#####Adding summaries for cssh (distances & ASVs+reads)

aflb2<-ddply(aflb, .(cssh), transform, frC2=(d_to_frC-mean(d_to_frC))/sd(d_to_frC), csshC=(d_to_csshC-mean(d_to_csshC))/sd(d_to_csshC), 
dfr=(avg_d_to_frpcreps-mean(avg_d_to_frpcreps))/sd(avg_d_to_frpcreps), sdpcr2=(sd_d_other_pcreps-mean(sd_d_other_pcreps))/sd(sd_d_other_pcreps), readscssh=round((reads-mean(reads))/sd(reads), digits=2), ASVscssh=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

aflj2<-ddply(aflj, .(cssh), transform, frC2=(d_to_frC-mean(d_to_frC))/sd(d_to_frC), csshC=(d_to_csshC-mean(d_to_csshC))/sd(d_to_csshC), 
dfr=(avg_d_to_frpcreps-mean(avg_d_to_frpcreps))/sd(avg_d_to_frpcreps), sdpcr2=(sd_d_other_pcreps-mean(sd_d_other_pcreps))/sd(sd_d_other_pcreps), readscssh=round((reads-mean(reads))/sd(reads), digits=2), ASVscssh=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

#####Adding rich/read summary

aflb3<-ddply(aflb2, .(substrate_type), transform, readsz=round((reads-mean(reads))/sd(reads), digits=2), ASVsz=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

aflj3<-ddply(aflj2, .(substrate_type), transform, readsz=round((reads-mean(reads))/sd(reads),digits=2), ASVsz=round((ASVs-mean(ASVs))/sd(ASVs), digits=2))

#####Adding THRESHOLDS
#using "jaccard" output first
#Merge jaccard and bray
#Changing colnames

afl<-aflj3
colnames(afl)
colnames(afl)[1:8]<-c("d_to_frC_j","d_to_csshC_j","avg_d_to_pcreps_j",
"avg_d_to_frpcreps_j","avg_d_other_pcreps_j","avg_d_other_frpcreps_j",
"sd_d_other_pcreps_j","sd_d_other_frpcreps_j")

colnames(afl)[19:21]<-c("frC_j","dpcr_j","sdpcr_j")
colnames(afl)[24:27]<-c("frC2_j","csshC_j","dfr_j","sdpcr2_j")

pafb<-aflb3[,c("sample_ID","d_to_frC","d_to_csshC","avg_d_to_pcreps",
"avg_d_to_frpcreps","avg_d_other_pcreps","avg_d_other_frpcreps",
"sd_d_other_pcreps","sd_d_other_frpcreps","frC","dpcr","sdpcr","frC2","csshC","dfr","sdpcr2")]

colnames(pafb)<-c("sample_ID","d_to_frC_b","d_to_csshC_b","avg_d_to_pcreps_b",
"avg_d_to_frpcreps_b","avg_d_other_pcreps_b","avg_d_other_frpcreps_b",
"sd_d_other_pcreps_b","sd_d_other_frpcreps_b","frC_b","dpcr_b","sdpcr_b","frC2_b","csshC_b","dfr_b","sdpcr2_b")

hias<-merge(afl, pafb, by="sample_ID", all=T)

hias$cluster<-sample_data(DADAwang3)$cluster[match(hias$sample_ID, sample_data(DADAwang3)$sample_ID)]
hias$season<-sample_data(DADAwang3)$season[match(hias$sample_ID, sample_data(DADAwang3)$sample_ID)]
hias$ss<-paste(hias$season,hias$substrate_type, sep="_")

###Thresh loop for reads and ASVs

thr1<-hias

#set quantiles
selector1<-c("q01","q025","q05","q1")
selector2<-c(.01,.025,.05,.1)
#set variables
selector3<-c("reads","ASVs")
#digits
selector5<-c(-3,0)

result_thresh<-expand.grid(q=selector1,measure=selector3,substrate_type=c("sediment","water"),value=NA)

for (i in 1:length(selector3))
{
selector4<-paste(selector3[i],selector1,sep="_")
for (e in 1:length(selector1))
{
threshold_s<-round(quantile(thr1[which(hias$substrate_type=="sediment"),selector3[i]], selector2[e]), digits=selector5[i])
threshold_w<-round(quantile(thr1[which(hias$substrate_type=="water"),selector3[i]], selector2[e]), digits=selector5[i])

thr1[,selector4[e]]<-ifelse(thr1$substrate_type=="sediment", thr1[,selector3[i]]>threshold_s, ifelse(thr1$substrate_type=="water", thr1[,selector3[i]]>threshold_w, NA))
unique(is.na(thr1[,selector4[e]]))

result_thresh[which(result_thresh$substrate_type=="sediment"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_s
result_thresh[which(result_thresh$substrate_type=="water"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_w
}
}

afl2<-ddply(thr1, .(sample_root), transform, Total_overq01=length(unique(sample_ID[reads_q01==T])), Total_overq025=length(unique(sample_ID[reads_q025==T])),
Total_overq05=length(unique(sample_ID[reads_q05==T])),
Total_overq1=length(unique(sample_ID[reads_q1==T])),
Total_ASVsoverq01=length(unique(sample_ID[ASVs_q01==T])), Total_ASVsoverq025=length(unique(sample_ID[ASVs_q025==T])),
Total_ASVsoverq05=length(unique(sample_ID[ASVs_q05==T])),
Total_ASVsoverq1=length(unique(sample_ID[ASVs_q1==T])), Total_pcreps=length(unique(sample_ID)))

afl2$reads_t<-ifelse(afl2$reads_q01==T&afl2$reads_q025==T&afl2$reads_q05==T&afl2$reads_q1==T,
"good",ifelse(afl2$reads_q01==T&afl2$reads_q025==T&afl2$reads_q05==T&afl2$reads_q1==F,
"10%",ifelse(afl2$reads_q01==T&afl2$reads_q025==T&afl2$reads_q05==F&afl2$reads_q1==F,
"5%",ifelse(afl2$reads_q01==T&afl2$reads_q025==F&afl2$reads_q05==F&afl2$reads_q1==F,
"2.5%",ifelse(afl2$reads_q01==F&afl2$reads_q025==F&afl2$reads_q05==F&afl2$reads_q1==F,
"1%",NA)))))

afl2$ASVs_t<-ifelse(afl2$ASVs_q01==T&afl2$ASVs_q025==T&afl2$ASVs_q05==T&afl2$ASVs_q1==T,
"good",ifelse(afl2$ASVs_q01==T&afl2$ASVs_q025==T&afl2$ASVs_q05==T&afl2$ASVs_q1==F,
"10%",ifelse(afl2$ASVs_q01==T&afl2$ASVs_q025==T&afl2$ASVs_q05==F&afl2$ASVs_q1==F,
"5%",ifelse(afl2$ASVs_q01==T&afl2$ASVs_q025==F&afl2$ASVs_q05==F&afl2$ASVs_q1==F,
"2.5%",ifelse(afl2$ASVs_q01==F&afl2$ASVs_q025==F&afl2$ASVs_q05==F&afl2$ASVs_q1==F,
"1%",NA)))))

afl2$reads_t2<-ifelse(afl2$reads>=50000,
">50000",ifelse(afl2$reads>=20000&afl2$reads<50000,
"<50000",ifelse(afl2$reads>=10000&afl2$reads<20000,
"<20000",ifelse(afl2$reads>=5000&afl2$reads<10000, 
"<10000",ifelse(afl2$reads<5000,"<5000", NA)))))

afl2$ASVs_t2<-ifelse(afl2$ASVs>=100,
"good",ifelse(afl2$ASVs>=50&afl2$ASVs<100,
"<100",ifelse(afl2$ASVs<50,
"<50",NA)))

###PLOTs
######Make classic reads vs. ASVs log plot
setwd(paste("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022"))

afl2$reads_t2 <- factor(afl2$reads_t2, levels = c(">50000","<50000","<20000","<10000",
"<5000"))
ggplot(data=afl2, aes(x=ASVs, y=reads)) +
geom_point(aes(color=as.factor(reads_t2)),size=1.4, alpha=0.8, shape=1, stroke=1) + 
geom_hline(yintercept = 50000, linetype = "dashed", size=0.05) +
geom_vline(xintercept = 100, linetype = "dashed", size=0.05) +
scale_colour_manual(values = c("snow3","palegreen4","orange","red","red4")) +
facet_wrap(.~ss, ncol=2) +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
labels = trans_format("log10", math_format(10^.x))) +
theme_classic()
ggsave("ASVs_reads_logclassic.pdf")


#####Histogram plots
setwd(paste("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022"))

hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","reads_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$reads_t <- factor(hiasplotd$reads_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=reads_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("reads_distj_ncbi.pdf")

hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","reads_t","substrate_type"), measure=c("frC_b","dpcr_b","frC2_b","csshC_b","dfr_b"))
hiasplotd$reads_t <- factor(hiasplotd$reads_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=reads_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("reads_distb_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","ASVs_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$ASVs_t <- factor(hiasplotd$ASVs_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=ASVs_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("ASVs_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","ASVs_t","substrate_type"), measure=c("frC_b","dpcr_b","frC2_b","csshC_b","dfr_b"))
hiasplotd$ASVs_t <- factor(hiasplotd$ASVs_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=ASVs_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("ASVs_distb_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","reads_t2","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$reads_t2 <- factor(hiasplotd$reads_t2, levels = c(">50000","<50000","<20000","<10000",
"<5000"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=reads_t2), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("reads2_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","ASVs_t2","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$ASVs_t2 <- factor(hiasplotd$ASVs_t2, levels = c("good","<100",
"<50"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=ASVs_t2), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("ASVs2_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","reads_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$reads_t <- factor(hiasplotd$reads_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=reads_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("reads_asvread_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","ASVs_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$ASVs_t <- factor(hiasplotd$ASVs_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=ASVs_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("ASVs_asvread_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","reads_t2","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$reads_t2 <- factor(hiasplotd$reads_t2, levels = c(">50000","<50000","<20000","<10000",
"<5000"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=reads_t2), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("reads2_asvread_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","ASVs_t2","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$ASVs_t2 <- factor(hiasplotd$ASVs_t2, levels = c("good","<100",
"<50"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=ASVs_t2), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("ASVs2_asvread_ncbi.pdf")


###Thresh loop for dists (fcr and dpcr)

thr2<-hias

#set quantiles
selector1<-c("q90","q95","q975","q99")
selector2<-c(.90,.95,.975,.99)
#set variables
selector3<-c("frC_j","dpcr_j")
#digits
selector5<-c(2,2)

result_thresh<-expand.grid(q=selector1,measure=selector3,substrate_type=c("sediment","water"),value=NA)


for (i in 1:length(selector3))
{
selector4<-paste(selector3[i],selector1,sep="_")
for (e in 1:length(selector1))
{
threshold_s<-round(quantile(thr2[which(hias$substrate_type=="sediment"),selector3[i]], selector2[e],na.rm=T), digits=selector5[i])
threshold_w<-round(quantile(thr2[which(hias$substrate_type=="water"),selector3[i]], selector2[e],na.rm=T), digits=selector5[i])

thr2[,selector4[e]]<-ifelse(thr2$substrate_type=="sediment", thr2[,selector3[i]]>threshold_s, ifelse(thr2$substrate_type=="water", thr2[,selector3[i]]>threshold_w, NA))
unique(is.na(thr2[,selector4[e]]))

result_thresh[which(result_thresh$substrate_type=="sediment"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_s
result_thresh[which(result_thresh$substrate_type=="water"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_w
}
}

afl2<-ddply(thr2, .(sample_root), transform, Total_frC_joverq90=length(unique(sample_ID[frC_j_q90==T])), Total_frC_joverq95=length(unique(sample_ID[frC_j_q95==T])),
Total_frC_joverq975=length(unique(sample_ID[frC_j_q975==T])),
Total_frC_joverq99=length(unique(sample_ID[frC_j_q99==T])),
Total_dpcr_joverq90=length(unique(sample_ID[dpcr_j_q90==T])), Total_dpcr_joverq95=length(unique(sample_ID[dpcr_j_q95==T])),
Total_dpcr_joverq975=length(unique(sample_ID[dpcr_j_q975==T])),
Total_dpcr_joverq99=length(unique(sample_ID[dpcr_j_q99==T])), Total_pcreps=length(unique(sample_ID)))

afl2$frC_j_t<-ifelse(afl2$frC_j_q90==F&afl2$frC_j_q95==F&afl2$frC_j_q975==F&afl2$frC_j_q99==F,
"good",ifelse(afl2$frC_j_q90==T&afl2$frC_j_q95==F&afl2$frC_j_q975==F&afl2$frC_j_q99==F,
"10%",ifelse(afl2$frC_j_q90==T&afl2$frC_j_q95==T&afl2$frC_j_q975==F&afl2$frC_j_q99==F,
"5%",ifelse(afl2$frC_j_q90==T&afl2$frC_j_q95==T&afl2$frC_j_q975==T&afl2$frC_j_q99==F,
"2.5%",ifelse(afl2$frC_j_q90==T&afl2$frC_j_q95==T&afl2$frC_j_q975==T&afl2$frC_j_q99==T,
"1%",NA)))))

afl2$dpcr_j_t<-ifelse(afl2$dpcr_j_q90==F&afl2$dpcr_j_q95==F&afl2$dpcr_j_q975==F&afl2$dpcr_j_q99==F,
"good",ifelse(afl2$dpcr_j_q90==T&afl2$dpcr_j_q95==F&afl2$dpcr_j_q975==F&afl2$dpcr_j_q99==F,
"10%",ifelse(afl2$dpcr_j_q90==T&afl2$dpcr_j_q95==T&afl2$dpcr_j_q975==F&afl2$dpcr_j_q99==F,
"5%",ifelse(afl2$dpcr_j_q90==T&afl2$dpcr_j_q95==T&afl2$dpcr_j_q975==T&afl2$dpcr_j_q99==F,
"2.5%",ifelse(afl2$dpcr_j_q90==T&afl2$dpcr_j_q95==T&afl2$dpcr_j_q975==T&afl2$dpcr_j_q99==T,
"1%",NA)))))


#####Histogram plots
setwd(paste("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022"))

hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","frC_j_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$frC_j_t <- factor(hiasplotd$frC_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=frC_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("frcjT_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","dpcr_j_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$dpcr_j_t <- factor(hiasplotd$dpcr_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=dpcr_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("dpcrjT_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","frC_j_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$frC_j_t <- factor(hiasplotd$frC_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=frC_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("frcjT_readsasvs_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","dpcr_j_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$dpcr_j_t <- factor(hiasplotd$dpcr_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=dpcr_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("dpcrjT_readsasvs_ncbi.pdf")


###Thresh loop for sdpcr_j and sdpcr2_j

thr3<-subset(hias, !is.na(hias$sdpcr_j)&!is.na(hias$sdpcr2_j))

#set quantiles
selector1<-c("q01","q025","q05","q1")
selector2<-c(.01,.025,.05,.1)
#set variables
selector3<-c("sdpcr_j","sdpcr2_j")
#digits
selector5<-c(2,2)

result_thresh<-expand.grid(q=selector1,measure=selector3,substrate_type=c("sediment","water"),value=NA)


for (i in 1:length(selector3))
{
selector4<-paste(selector3[i],selector1,sep="_")
for (e in 1:length(selector1))
{
threshold_s<-round(quantile(thr3[which(thr3$substrate_type=="sediment"),selector3[i]], selector2[e]), digits=selector5[i])
threshold_w<-round(quantile(thr3[which(thr3$substrate_type=="water"),selector3[i]], selector2[e]), digits=selector5[i])

thr3[,selector4[e]]<-ifelse(thr3$substrate_type=="sediment", thr3[,selector3[i]]>threshold_s, ifelse(thr3$substrate_type=="water", thr3[,selector3[i]]>threshold_w, NA))
unique(is.na(thr3[,selector4[e]]))

result_thresh[which(result_thresh$substrate_type=="sediment"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_s
result_thresh[which(result_thresh$substrate_type=="water"&
result_thresh$measure==selector3[i]&
result_thresh$q==selector1[e]),"value"]<-threshold_w
}
}

afl2<-ddply(thr3, .(sample_root), transform, Total_sdpcr_joverq01=length(unique(sample_ID[sdpcr_j_q01==T])), Total_sdpcr_joverq025=length(unique(sample_ID[sdpcr_j_q025==T])),
Total_sdpcr_joverq05=length(unique(sample_ID[sdpcr_j_q05==T])),
Total_sdpcr_joverq1=length(unique(sample_ID[sdpcr_j_q1==T])),
Total_sdpcr2_joverq01=length(unique(sample_ID[sdpcr2_j_q01==T])), Total_sdpcr2_joverq025=length(unique(sample_ID[sdpcr2_j_q025==T])),
Total_sdpcr2_joverq05=length(unique(sample_ID[sdpcr2_j_q05==T])),
Total_sdpcr2_joverq1=length(unique(sample_ID[sdpcr2_j_q1==T])), Total_pcreps=length(unique(sample_ID)))

afl2$sdpcr_j_t<-ifelse(afl2$sdpcr_j_q01==T&afl2$sdpcr_j_q025==T&afl2$sdpcr_j_q05==T&afl2$sdpcr_j_q1==T,
"good",ifelse(afl2$sdpcr_j_q01==T&afl2$sdpcr_j_q025==T&afl2$sdpcr_j_q05==T&afl2$sdpcr_j_q1==F,
"10%",ifelse(afl2$sdpcr_j_q01==T&afl2$sdpcr_j_q025==T&afl2$sdpcr_j_q05==F&afl2$sdpcr_j_q1==F,
"5%",ifelse(afl2$sdpcr_j_q01==T&afl2$sdpcr_j_q025==F&afl2$sdpcr_j_q05==F&afl2$sdpcr_j_q1==F,
"2.5%",ifelse(afl2$sdpcr_j_q01==F&afl2$sdpcr_j_q025==F&afl2$sdpcr_j_q05==F&afl2$sdpcr_j_q1==F,
"1%",NA)))))

afl2$sdpcr2_j_t<-ifelse(afl2$sdpcr2_j_q01==T&afl2$sdpcr2_j_q025==T&afl2$sdpcr2_j_q05==T&afl2$sdpcr2_j_q1==T,
"good",ifelse(afl2$sdpcr2_j_q01==T&afl2$sdpcr2_j_q025==T&afl2$sdpcr2_j_q05==T&afl2$sdpcr2_j_q1==F,
"10%",ifelse(afl2$sdpcr2_j_q01==T&afl2$sdpcr2_j_q025==T&afl2$sdpcr2_j_q05==F&afl2$sdpcr2_j_q1==F,
"5%",ifelse(afl2$sdpcr2_j_q01==T&afl2$sdpcr2_j_q025==F&afl2$sdpcr2_j_q05==F&afl2$sdpcr2_j_q1==F,
"2.5%",ifelse(afl2$sdpcr2_j_q01==F&afl2$sdpcr2_j_q025==F&afl2$sdpcr2_j_q05==F&afl2$sdpcr2_j_q1==F,
"1%",NA)))))


#####Histogram plots
setwd(paste("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022"))

hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","sdpcr_j_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$sdpcr_j_t <- factor(hiasplotd$sdpcr_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=sdpcr_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("sdpcr_j_t_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","sdpcr_j_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$sdpcr_j_t <- factor(hiasplotd$sdpcr_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=sdpcr_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("sdpcr_j_t_asvread_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","sdpcr2_j_t","substrate_type"), measure=c("frC_j","dpcr_j","frC2_j","csshC_j","dfr_j"))
hiasplotd$sdpcr2_j_t <- factor(hiasplotd$sdpcr2_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=sdpcr2_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("sdpcr2_j_t_distj_ncbi.pdf")



hiasplotd<-melt(afl2, id=c("ss", "reads","ASVs","sample_ID","sdpcr2_j_t","substrate_type"), measure=c("readscssh","ASVscssh","readsfr","ASVsfr",
"readsz","ASVsz"))
hiasplotd$sdpcr2_j_t <- factor(hiasplotd$sdpcr2_j_t, levels = c("good","10%",
"5%","2.5%","1%"))

ggplot(hiasplotd, aes(x=value)) +
geom_histogram(aes(fill=sdpcr2_j_t), position="identity", alpha=0.6, bins=60, color="black", size=0.05)+
geom_vline(xintercept = 0, linetype = "dashed", size=0.05) +
theme_classic() + 
facet_wrap(variable~ss, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","gold","orange","red","red4"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave("sdpcr2_j_t_asvread_ncbi.pdf")



#########Make ordination plot for visual inspection
#using reads threshold

fodf<-afl2
fodf$reads_t3<-ifelse(fodf$reads>=50000,T,F)

head(sample_data(DADAwang3))
PD<-DADAwang3
sample_data(PD)$reads_q025<-fodf$reads_q025[match(sample_data(PD)$sample_ID, fodf$sample_ID)]
sample_data(PD)$reads_t3<-fodf$reads_t3[match(sample_data(PD)$sample_ID, fodf$sample_ID)]
head(sample_data(PD))

sample_data(PD)$css<-paste(sample_data(PD)$cluster,sample_data(PD)$season,sample_data(PD)$substrate_type,sep="_")
sample_data(PD)$fp<-paste(sample_data(PD)$field_replicate,sample_data(PD)$PCR_replicate,sep=".")
allg<-unique(sample_data(PD)$css)

#Set plot aspcts
bgvej<- c("transparent","red")
pcvec <- c(21,22,23)

for (csample in 1:length(allg))
{
#subset cluster
tcssh<-subset_samples(PD, css==allg[csample])
#create list for storing habitat-wise plot
toplotj<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotb<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotd<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotg<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotr<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotl<-as.list(levels(unique(sample_data(tcssh)$habitat)))
toplotgc<-as.list(levels(unique(sample_data(tcssh)$habitat)))
for (h in 1:length(toplotj))
{
tcssh2<-subset_samples(tcssh, habitat==toplotj[[h]])
tcssh0 = filter_taxa(tcssh2, function(x) sum(x) > 0, TRUE)
p_rt<-data.frame(sample_data(tcssh0))
p_rt$lab<-ifelse(p_rt$reads_t3==T,"",as.character(p_rt$fp))
toplotd[[h]]<-factor(p_rt$reads_t3, levels = c("TRUE","FALSE"))
toplotg[[h]]<-p_rt$field_replicate
toplotgc[[h]]<-p_rt$cssh
toplotr[[h]]<-as.character(p_rt$lab)
#calculate cluster dist data
otug<-t(data.frame(otu_table(tcssh0), check.names=F))
dist_nj<-vegdist(otug, method="jaccard", upper=T)
toplotj[[h]]<-cmdscale(dist_nj, k=2, eig=T)
dist_nb<-vegdist(otug, method="bray", upper=T)
toplotb[[h]]<-cmdscale(dist_nb, k=2, eig=T)
}
if (unique(p_rt[,"substrate_type"])=="sediment") {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment/ordination")
} else {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water/ordination")
}
pdf(paste(allg[csample],".pdf",sep=""))
par(mfrow=c(length(toplotj),2))
for (t in 1:length(toplotj))
{
pj<-ordiplot(toplotj[[t]], type = "n", display="sites", main=paste(toplotl[[t]],"jaccard",sep="_"))
points(pj, "sites", col= "black", pch = pcvec[as.factor(toplotg[[t]])], bg = bgvej[as.factor(toplotd[[t]])], cex=1.5)
text(pj, "sites", labels = toplotr[[t]], cex=0.4)
ordiellipse(pj, toplotg[[t]], display="sites", kind = "se", conf=0.95, draw = "polygon", col="grey", border = "lightgrey")
ordiellipse(pj, toplotgc[[t]], display="sites", kind = "se", conf=0.95, draw = "polygon", col="grey", border = "lightgrey")
pb<-ordiplot(toplotb[[t]], type = "n", display="sites", main=paste(toplotl[[t]],"bray",sep="_"))
points(pb, "sites", col= "black", pch = pcvec[as.factor(toplotg[[t]])], bg = bgvej[as.factor(toplotd[[t]])], cex=1.5)
text(pb, "sites", labels = toplotr[[t]], cex=0.4)
ordiellipse(pb, toplotg[[t]], display="sites", kind = "se", conf=0.95, draw = "polygon", col="grey", border = "lightgrey")
ordiellipse(pb, toplotgc[[t]], display="sites", kind = "se", conf=0.95, draw = "polygon", col="grey", border = "lightgrey")
}
dev.off()
}


####Make Cluster hists

fodf<-afl2
fodf$reads_t3<-ifelse(fodf$reads>=50000,T,F)
selector1<-c("frC_j","ASVsfr","dpcr_j")

allg<-unique(fodf$ss)

for (csample in 1:length(allg))
{
for (i in 1:length(selector1))
{

fodf2<-subset(fodf, ss==allg[csample])

if (unique(fodf2[,"substrate_type"])=="sediment") {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment/cluster_hists")
} else {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water/cluster_hists")
}
fodf2$reads_t3 <- factor(fodf2$reads_t3, levels = c("TRUE","FALSE"))
ggplot(fodf2, aes_string(x=selector1[i])) +
geom_histogram(aes(fill=reads_t3), position="identity", alpha=0.6, bins=15, color="black", size=0.05)+
theme_classic() + 
facet_wrap(.~cluster, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","red"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste(allg[csample],selector1[i],".pdf",sep="_"))
}
}

fodf<-afl2
fodf$reads_t3<-ifelse(fodf$reads>=50000,T,F)
selector1<-c("frC2_j","dfr_j","ASVscssh","sdpcr2_j")

allg<-unique(fodf$ss)

for (csample in 1:length(allg))
{
for (i in 1:length(selector1))
{

fodf2<-subset(fodf, ss==allg[csample])

if (unique(fodf2[,"substrate_type"])=="sediment") {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment/cluster_hists")
} else {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water/cluster_hists")
}
fodf2$reads_t3 <- factor(fodf2$reads_t3, levels = c("TRUE","FALSE"))
ggplot(fodf2, aes_string(x=selector1[i])) +
geom_histogram(aes(fill=reads_t3), position="identity", alpha=0.6, bins=15, color="black", size=0.05)+
theme_classic() + 
facet_wrap(.~cluster, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","red"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste("2",allg[csample],selector1[i],".pdf",sep="_"))
}
}

fodf<-afl2
fodf$ASVs_t3<-ifelse(fodf$ASVs>=100,T,F)
selector1<-c("frC_j","readsfr","dpcr_j")

allg<-unique(fodf$ss)

for (csample in 1:length(allg))
{
for (i in 1:length(selector1))
{

fodf2<-subset(fodf, ss==allg[csample])

if (unique(fodf2[,"substrate_type"])=="sediment") {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment/cluster_hists")
} else {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water/cluster_hists")
}
fodf2$ASVs_t3 <- factor(fodf2$ASVs_t3, levels = c("TRUE","FALSE"))
ggplot(fodf2, aes_string(x=selector1[i])) +
geom_histogram(aes(fill=ASVs_t3), position="identity", alpha=0.6, bins=15, color="black", size=0.05)+
theme_classic() + 
facet_wrap(.~cluster, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","red"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste("ASVs",allg[csample],selector1[i],".pdf",sep="_"))
}
}

fodf<-afl2
fodf$ASVs_t3<-ifelse(fodf$ASVs>=100,T,F)
selector1<-c("frC2_j","dfr_j","readscssh","sdpcr2_j")

allg<-unique(fodf$ss)

for (csample in 1:length(allg))
{
for (i in 1:length(selector1))
{

fodf2<-subset(fodf, ss==allg[csample])

if (unique(fodf2[,"substrate_type"])=="sediment") {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/sediment/cluster_hists")
} else {
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022/water/cluster_hists")
}
fodf2$ASVs_t3 <- factor(fodf2$ASVs_t3, levels = c("TRUE","FALSE"))
ggplot(fodf2, aes_string(x=selector1[i])) +
geom_histogram(aes(fill=ASVs_t3), position="identity", alpha=0.6, bins=15, color="black", size=0.05)+
theme_classic() + 
facet_wrap(.~cluster, ncol=4, scale="free_y") + 
scale_colour_manual(values = c("transparent","red"), aesthetics="fill") +
theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 6), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm"))
ggsave(paste("2ASVs",allg[csample],selector1[i],".pdf",sep="_"))
}
}

#Find a way to discard bad pcreps now before normalization
#Perhaps, use frC2+Rarefaction diagnose - bad frs
#Then, change the script below for an rarefaction threshold optimized version



##Create final categorization before normalization
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")

readsi<-sample_sums(DADAwang1)
combinedi<-cbind(readsi, sample_data(DADAwang1))
combinedi<-data.frame(combinedi)
threshold<-round(quantile(combinedi$readsi, .025), digits=-3)
combinedi$q<-combinedi$readsi>threshold
zcom<-ddply(combinedi, .(sample_root), transform, Total_overq=length(unique(sample_ID[q==T])))

#samples with all 4 PCR replicates above threshold
#r_et
zcom$Norm_rule<-ifelse(zcom$Total_overq==4, "r_et", NA)

#samples with 3 PCR replicates above threshold
#r_to_1
xcom<-ddply(zcom, .(sample_root), transform, 
Norm_rule=ifelse(Total_overq==4, "r_et",
ifelse(length(unique(sample_ID[q==T]))==3&&length(unique(sample_ID[readsi>=threshold+threshold/3]))==3, "r_to_1", NA)))

#r_to_2
ycom<-ddply(xcom, .(sample_root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&length(unique(sample_ID[readsi>=threshold+threshold/2]))==2,"r_to_2", NA)))

#r_to_3
vcom<-ddply(ycom, .(sample_root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&max(readsi)>=2*threshold,"r_to_3", NA)))

#r_to_4
ucom<-ddply(vcom, .(sample_root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.85*threshold&&length(unique(sample_ID[readsi>=1.05*threshold]))==3,"r_to_4", NA)))

#r_to_5
icom<-ddply(ucom, .(sample_root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.85*threshold&&length(unique(sample_ID[readsi>=1.075*threshold]))==2,"r_to_5", NA)))

#r_to_6
ecom<-ddply(icom, .(sample_root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.75*threshold&&max(readsi)>=1.25*threshold,"r_to_6", NA)))

#Calculate max depth per field replicate
p_norm_sample_data_h<-ddply(ecom, .(sample_root), transform, maxi=max(readsi))

#Define step-wise normalization
p_norm_sample_data<-ddply(p_norm_sample_data_h, .(sample_ID), transform, final_rule=
ifelse(is.na(Norm_rule), "trash", 
ifelse(Norm_rule=="r_et", "threshold", 
ifelse(Norm_rule=="r_to_1"&&q=="TRUE","threshold_plus_third", 
ifelse(Norm_rule=="r_to_1"&&q=="FALSE", "trash", ifelse(Norm_rule=="r_to_2"&&q=="TRUE"&&readsi>=threshold+threshold/2,"threshold_plus_half", ifelse(Norm_rule=="r_to_2"&&q=="TRUE"&&readsi<threshold+threshold/2, "threshold", ifelse(Norm_rule=="r_to_2"&&q=="FALSE", "trash", ifelse(Norm_rule=="r_to_3"&&q=="TRUE"&&readsi>=2*threshold,"threshold_twice", ifelse(Norm_rule=="r_to_3"&&q=="TRUE"&&readsi<2*threshold, "threshold", ifelse(Norm_rule=="r_to_3"&&q=="FALSE", "trash",
ifelse(Norm_rule=="r_to_4"&&q=="TRUE"&&readsi>=1.05*threshold,"threshold_1050",  ifelse(Norm_rule=="r_to_4"&&q=="FALSE", "keep_0850",
ifelse(Norm_rule=="r_to_5"&&q=="TRUE"&&readsi>=1.075*threshold,"threshold_1075", 
ifelse(Norm_rule=="r_to_5"&&q=="TRUE"&&readsi<1.075*threshold, "threshold", ifelse(Norm_rule=="r_to_5"&&q=="FALSE", "keep_0850",
ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi>=1.25*threshold&&readsi>=maxi,"threshold_1250",
ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi>=1.25*threshold&&readsi<maxi,"threshold", ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi<1.25*threshold, "threshold", ifelse(Norm_rule=="r_to_6"&&q=="FALSE", "keep_0750",
"issue"))))))))))))))))))))

#update object's sample data
sample_data(DADAwang1)$final_rule<-p_norm_sample_data$final_rule[match(sample_data(DADAwang1)$sample_ID, p_norm_sample_data$sample_ID)]
steps<-levels(factor(p_norm_sample_data$final_rule))

#Check trash samples and save it
ssss<-subset(p_norm_sample_data, final_rule=="trash")
write.table(levels(factor(ssss$sample_ID)), "list_discarded_PCR_replicates_ncbi.txt", sep=" ")

cccc<-subset(p_norm_sample_data, is.na(Norm_rule))
write.table(levels(factor(cccc$sample_root)), "list_discarded_field_replicates_ncbi.txt", sep=" ")

#Normalize each category according to the associated sample sizes
steps=steps[steps!="trash"]

sizes<-c(threshold,
0.85*threshold,
1.075*threshold,
1.050*threshold,
threshold+threshold/3,
threshold+threshold/2,
threshold*2,
0.75*threshold,
1.250*threshold
)

sizes2<-round(as.numeric(sizes),digits=0)

prtc<-c("threshold",
"keep_0850",
"threshold_1075",
"threshold_1050",
"threshold_plus_third",
"threshold_plus_half",
"threshold_twice",
"keep_0750",
"threshold_1250")

prtc2<-as.character(prtc)
gaba<-data.frame(cbind(sizes2, prtc2), stringsAsFactors = FALSE)
pepe<-gaba[match(steps, gaba$prtc2),]
stnd<-subset(pepe, !is.na(pepe$prtc2))

for(e in 1:nrow(stnd)) {
  assign(paste0("samples", stnd[e,2]), rarefy_even_depth(subset_samples(DADAwang1, final_rule==stnd[e,2]), sample.size=as.numeric(stnd[e,1]), replace=FALSE, rngseed= 13072021))
}

#merge objects
myobs<-paste0("samples",stnd$prtc2)
tudao_r_et<-do.call(merge_phyloseq, mget(myobs), quote = FALSE)

tudao_pos_merge<-merge_samples(tudao_r_et, "sample_root")

#Re-build sample_data
tudao0 = filter_taxa(tudao_pos_merge, function(x) sum(x) > 0, TRUE)
p_ncbi = prune_samples(sample_sums(tudao0)>0,tudao0)

d<-data.frame(sample_data(p_ncbi)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

sample_data(p_ncbi)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

#Save final files

tax_m<-data.frame(tax_table(p_ncbi))
otu_m<-data.frame(otu_table(p_ncbi),check.names=F)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
write.table(data.frame(sample_data(p_ncbi), check.names=F), "f_ncbi_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results")
write.table(otu_m, "f_otu_ncbi.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "f_tax_ncbi.txt", sep="\t", quote=FALSE, row.names=TRUE)





















###########trash


#######Determine cutoff here:
hias$lab<-ifelse(hias$readsfr<0.75&hias$readscssh<0.75&hias$ASVscssh<0.75, as.character(hias$sample_root), as.character(""))
hias$threshold<-ifelse(hias$readsfr<0.75&hias$readscssh<0.75&hias$ASVscssh<0.75, T, F)

cvflj1<-melt(hias, id=c("threshold","reads","ASVs","sample_ID","frc_j","frc2_j","q01","q025","q05","q1","lab","substrate_type"), measure=c("readsz","ASVsz"))

cvflb1<-melt(hias, id=c("reads","ASVs","sample_ID","frc_b","frc2_b","q01","q025","q05","q1","lab","substrate_type"), measure=c("readsz","ASVsz"))

cvflj2<-melt(hias, id=c("reads","ASVs","sample_ID","frc_j","frc2_j","q01","q025","q05","q1","lab","substrate_type"), measure=c("readsfr","readscssh","ASVsfr","ASVscssh"))

cvflb2<-melt(hias, id=c("reads","ASVs","sample_ID","frc_b","frc2_b","q01","q025","q05","q1","lab","substrate_type"), measure=c("readsfr","readscssh","ASVsfr","ASVscssh"))

cvfl3<-melt(hias, id=c("threshold","reads","ASVs","sample_ID","ASVsz","ASVsfr","ASVscssh","lab","substrate_type"), measure=c("reads","ASVs"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022")


ggplot(data=hias, aes(x=ASVs, y=reads)) +
geom_jitter(aes(color=log(readscssh)), size=0.8, alpha=0.4, shape=20) + 
geom_hline(aes(yintercept = reads), thresholds, linetype = "dashed", size=0.5) +
geom_vline(aes(xintercept = ASVs), thresholds, linetype = "dashed", size=0.5) +
facet_grid(.~substrate_type) +
theme_classic()
ggsave("reads_q.pdf")

ggplot(data=cvflj1, aes(x=frc2_j, y=reads)) +
geom_jitter(aes(color=threshold), size=0.8, alpha=0.4, shape=20) + 
geom_hline(aes(yintercept = reads), thresholds, linetype = "dashed", size=0.5) +
facet_grid(.~substrate_type) +
theme_classic()
ggsave("reads_frc2.pdf")



ggplot(data=cvflj1, aes(x=frc_j, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(-1,0,1), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("j_frc_reads_ASVs_z.pdf")

ggplot(data=cvflj2, aes(x=frc_j, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("j_frc_reads_ASVs_frc_cssh.pdf")

ggplot(data=cvflj1, aes(x=frc2_j, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(-1,0,1), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("j_frc2_reads_ASVs_z.pdf")

ggplot(data=cvflj2, aes(x=frc2_j, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("j_frc2_reads_ASVs_frc_cssh.pdf")

ggplot(data=cvflb1, aes(x=frc_b, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(-1,0,1), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("b_frc_reads_ASVs_z.pdf")

ggplot(data=cvflb2, aes(x=frc_b, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("b_frc_reads_ASVs_frc_cssh.pdf")

ggplot(data=cvflb1, aes(x=frc2_b, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(-1,0,1), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("b_frc2_reads_ASVs_z.pdf")

ggplot(data=cvflb2, aes(x=frc2_b, y=value)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("b_frc2_reads_ASVs_frc_cssh.pdf")


ggplot(data=cvfl3, aes(x=sqrt(ASVs), y=log(reads))) +
geom_jitter(aes(color=threshold), size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y")+
theme_classic()
ggsave("reads_q.pdf")

ggplot(data=cvfl3, aes(x=value, y=ASVsz)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(-1,0,1), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("reads_ASVsz.pdf")

ggplot(data=cvfl3, aes(x=value, y=ASVsfr)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("reads_ASVsfr.pdf")

ggplot(data=cvfl3, aes(x=value, y=ASVscssh)) +
geom_jitter(aes(color=q),size=0.8, alpha=0.4, shape=20) +
facet_grid(variable~substrate_type, scale="free_y") + geom_hline(yintercept = c(0.5,1,1.5), linetype = "dashed", size=0.5) +
theme_classic()
ggsave("reads_ASVscssh.pdf")


#Best plots

ggplot(data=subset(cvflj2, variable=="readsfr"&value>1), aes(x=frc2_j, y=value)) +
geom_jitter(aes(color=q, size=sqrt(ASVs)/10), alpha=0.5, shape=20) + geom_smooth(method="lm",se=F,color="black",size=0.5) +
facet_grid(q~substrate_type) +
theme_classic()
ggsave("j_frc2_readsfrc_pos.pdf")

ggplot(data=subset(cvflj2, variable=="readsfr"&value<1), aes(x=frc2_j, y=value)) +
geom_jitter(aes(color=q, size=sqrt(ASVs)/10), alpha=0.5, shape=20) + geom_smooth(method="lm",se=F,color="black",size=0.5) +
facet_grid(q~substrate_type, scale="free_y") +
theme_classic()
ggsave("j_frc2_readsfrc_neg.pdf")

ggplot(data=subset(cvflj2, variable=="readsfr"|variable=="readscssh"), aes(x=q, y=value)) +
geom_jitter(aes(color=q, size=sqrt(ASVs)/10), alpha=0.5, shape=20) + 
facet_grid(variable~substrate_type, scale="free_y") +
theme_classic()
ggsave("j_q_reads.pdf")

ggplot(data=subset(cvfl3, variable=="readscssh"), aes(x=log(ASVscssh
), y=value)) +
geom_jitter(aes(color=q, size=sqrt(reads)/10), alpha=0.5, shape=20) + geom_smooth(method="lm",se=F,color="black",size=0.5) +
facet_grid(q~substrate_type) +
theme_classic()
ggsave("j_readsASVcssh_log.pdf")



subset(afl2, q==F)[,c("sample_ID","Total_overq","j_dist_out_fr","j_dist_out_cssh","b_dist_out_fr","b_dist_out_cssh","reads")]
















ordiellipse(pj, toplotg[[t]], display="sites", kind = "se", conf=0.95, draw = "lines")




pdf(paste(allg[csample],".pdf",sep=""))
pj<-ordiplot(toplotj[[t]], type = "n", display="sites")
points(pj, "sites", col=toplotq[[t]], pch =toplotdj[[t]])
dev.off()

ordiplot(toplotj[[t]], type = "n", display="sites", pch = pcvej[toplotdj[[t]]], bg = bgvec[toplotg[[t]]])
unique(p_rt[,"substrate_type"])=="sediment"

toplotb[[h]]<-betadisper(dist_nb, p_rt$field_replicate, type = "centroid")
plot(toplotb[[t]],label=F, hull = FALSE, ellipse = TRUE, col= colvec[toplotq[[h]]], pch = pcveb[toplotdb[[h]]], bg = bgvec[toplotg[[h]]])

df.q <- ddply(df.example, .(model, type),
              summarize, q=quantile(value, c(.025, .975))) 



#Now, estimate mean and standard deviations per field replicate, make some test plots for overall dispersions

afl<-result_dist
sddf<-data.frame(sample_data(DADAwang1))
afl$sample_root<-sddf$sample_root[match(rownames(afl), sddf$sample_ID)]
sdafl<-ddply(afl, .(sample_root), summarize, sd_d_to_frC=sd(d_to_frC), sd_d_to_csshC=sd(d_to_csshC), sd_avg_d_to_pcreps=sd(avg_d_to_pcreps), sd_avg_d_to_frpcreps=sd(avg_d_to_frpcreps), sd_reads=sd(reads), avg_reads=mean(reads), cv_reads=sd(reads)/mean(reads)*100)

cvfl<-melt(sdafl, id=c("sample_root","cv_reads"), measure=c("sd_d_to_frC","sd_d_to_csshC"))

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary_2022")
ggplot(data=cvfl, aes(x=value, y=cv_reads)) +
geom_jitter(size=0.1) +
geom_smooth(method=lm, se=T, level=0.5) +
facet_wrap(.~variable, ncol=1)
ggsave("cv_reads_distances.pdf")

mod1 <- lm(sd_d_to_csshC~cv_reads,sdafl)
outl = order(-cooks.distance(mod1))[1:10]

pdf("test_out_cooks.pdf")
par(mfrow=c(1,2))
plot(mod1,which=1,labels.id=sdafl$sample_root)
plot(fitted(mod1),residuals(mod1))
panel.smooth(fitted(mod1),residuals(mod1))
text(fitted(mod1)[outl]-0.045,residuals(mod1)[outl],
sdafl$sample_root[outl],col="red",cex=0.8)
dev.off()


h <- hist(x, breaks=20, plot=F) # h$breaks and h$mids
cols <- c('grey', "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")  
k <- cols[findInterval(h$mids, quantile(x), rightmost.closed=T, all.inside=F) + 1]
# plot the histogram with the colours
plot(h, col=k)

dists<-colnames(sdafl)[2:5]
reads<-colnames(sdafl)[6:8]

for (i in 1:length(dists))
{
for (e in 1:length(reads))
{
algo<-sdafl[,c("sample_root",dists[i],reads[e])]
fit <- lm(as.formula(paste(reads[e], "~", dists[i])),data=algo)
dat <- predict(fit, interval="confidence")
algo$inside <- ifelse(algo[,reads[e]] < dat[,"upr"] & algo[,reads[e]] > dat[,"lwr"], "", as.character(algo$sample_root))
ggplot(algo, aes_string(y = reads[e],x = dists[i])) +
geom_point(size=0.1) +
geom_smooth(method=lm) + geom_text(aes(label=inside),size=2)
ggsave(paste(dists[i],reads[e],".pdf",sep=""))
}}




