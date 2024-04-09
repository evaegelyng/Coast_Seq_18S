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
otu_mat<-as.matrix(read.table("cleaned_tax_ncbi.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("cleaned_ncbi.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_mat_b)
p_ncbi = phyloseq(OTU, TAX_b)

#Load metadata
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/metadata")
metadata<-read.table("cleaned_ncbi_metadata.txt", sep="\t", header=T)

##SITE INFO
c_s<-read.table("cluster_site.txt", sep="\t", header=T)
metadata$Location<-c_s$Site_name[match(metadata$cluster, c_s$cluster)]
metadata$cluster<-as.integer(metadata$cluster)
metadata$cl_se<-as.character(paste(metadata$cluster,metadata$season,sep="_"))
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_ncbi, sampledata)
DADAwang1

##Make histogram for raw read counts
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary")

readsi<-sample_sums(DADAwang1)
combinedi<-cbind(readsi, sample_data(DADAwang1))
combinedi<-data.frame(combinedi)

threshold<-round(quantile(combinedi$readsi, c(.025)), digits=-3)
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

threshold<-round(quantile(combined$reads, c(.025)), digits=-3)
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
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary/sediment")

###Although water and sediment plots are built separately, 
#this first part addresses the whole dataset,
#aiming for the comparison between different substrate types

###Sediment PCR rarefaction - diagnosis

#Use threshold of whole dataset
threshold<-round(quantile(combinedi$readsi, c(.025)), digits=-3)
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

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red") + facet_wrap(~sample_root, ncol=6, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Sediment - field replicates below threshold", x ="PCR replicate", y = "number of reads (millions)", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("sed_summary_q025_ncbi.pdf")


setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/reads_summary/water")

###Water PCR rarefaction - diagnosis

#Use threshold of whole dataset
threshold<-round(quantile(combinedi$readsi, c(.025)), digits=-3)
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

ggplot(data = data.frame(xreadsi3), aes(sample_root, readsi, fill = as.factor(Total_overq025))) +  geom_bar(aes(sample_ID, readsi, fill = as.factor(Total_overq025)), stat="identity", position = position_dodge2(preserve='single'), color="black", size=0.02) + geom_hline(yintercept=threshold, linetype="dashed", color = "red") + facet_wrap(~sample_root, ncol=7, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=6), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + scale_fill_manual(values=c("red1", "darkorange1", "gold1", "forestgreen", "grey50")) + labs(title="Water - field replicates below threshold", x ="PCR replicate", y = "number of reads (millions)", fill = paste("Number of PCR replicates \n with total reads > ",threshold,sep=""))
ggsave("water_summary_q025_ncbi.pdf")

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

