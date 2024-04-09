library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/depth_cutoff")

#Inputs
pr2_t<-read.table("summary_depth_cutoff_pr2.txt", sep="\t", header=T)
silva_t<-read.table("summary_depth_cutoff_silva.txt", sep="\t", header=T)
ncbi_t<-read.table("summary_depth_cutoff_ncbi.txt", sep="\t", header=T)
pr2_t$db<-"pr2"
silva_t$db<-"silva"
ncbi_t$db<-"ncbi"

all_t<-rbind(pr2_t, silva_t, ncbi_t)
all_t$Quantile<-as.character(all_t$Quantile)
all_t[5,1]<-as.character(0.5)
all_t[11,1]<-as.character(0.5)
all_t[17,1]<-as.character(0.5)

dasf<-melt(all_t, id=c("Quantile","db"), measure=c("Threshold","Total_replicates_discarded",
"Total_reads","Total_ASVs","Field_replicates","Depth",
"Mean_richness","SD_richness","Mean_bray_centroid_field_rep", "SD_bray_centroid_field_rep"))
head(dasf)

ggplot(dasf, aes(as.factor(Quantile), value, fill=db)) + geom_bar(stat="identity", position="dodge") + facet_wrap(~variable, ncol=3, scales="free_y")+theme_bw() + theme(axis.text.x = element_text(size=7), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=7))
ggsave("depth_cutoff_comparison.pdf")
