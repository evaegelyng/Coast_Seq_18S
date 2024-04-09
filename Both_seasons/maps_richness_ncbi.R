library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
library("maps")
library("mapproj")
library("dplyr")
library("mapdata")
library("ggmap")
library('lwgeom')
library('viridis')
library("maptools")
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(gdistance)
library(ncdf4)
library(Matrix)
library(rnaturalearth)
library(rnaturalearthdata)

library(SDraw)

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

#Remove Eelgrass reads
#DADAwang2<-subset_taxa(DADAwang1, !Phylum=="Streptophyta")
#doi<-data.frame(sample_data(DADAwang2))

#For stacked plots, use "tudao", in which field replicates are merged as starting object:
#Alternatively, see below how to conserve field replicates allowing sd calculation
#Stacked

tudao<-merge_samples(DADAwang1, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
d$sshc<-rownames(d)
d$substrate_type<-metadata$substrate_type[match(d$sshc, metadata$sshc)]
d$habitat<-metadata$habitat[match(d$sshc, metadata$sshc)]
d$season<-metadata$season[match(d$sshc, metadata$sshc)]
d$snch<-paste(d$season,d$cluster,d$habitat,sep="_")
head(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc","snch")]

#Stacked

d<-data.frame(sample_data(tudao))
otug<-otu_table(tudao)
otug<-t(data.frame(otug, check.names=F))
taxom<-data.frame(tax_table(tudao))
ch<-d$sshc

#Overall richness

rich<-colSums(otug != 0)
combined<-cbind(rich, d)

#Creating separate tables for metazoans and non-metazoans
notug<-data.frame(otug)
notug$Division<-taxom$Division[match(rownames(notug), rownames(taxom))]
otugm<-subset(notug, Division=="Metazoa")
otugm<-subset(otugm, select=-c(Division))
otugnm<-subset(notug, !Division=="Metazoa")
otugnm<-subset(otugnm, select=-c(Division))

richm<-colSums(otugm != 0)
combinedm<-cbind(richm, d)

richnm<-colSums(otugnm != 0)
combinednm<-cbind(richnm, d)

###Retrieve coordinates/other variables

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

combined$Latitude<-merged_wat$Latitude[match(combined$snch, merged_wat$snch)]
combined$Longitude<-merged_wat$Longitude[match(combined$snch, merged_wat$snch)]

head(combined)

combined$habitat <- factor(combined$habitat, levels = c("sand", "rocks", "eelgrass"))
combined$season <- factor(combined$season, levels = c("spring", "autumn"))
combined$substrate_type <- factor(combined$substrate_type, levels = c("sediment", "water"))
combined$ss<-paste(combined$substrate_type, combined$season,  sep="_")
combined$ss <- factor(combined$ss, levels = c("water_autumn","sediment_autumn","sediment_spring","water_spring"))
combined$stsc<-paste(combined$substrate_type, combined$season, combined$cluster, sep="_")
combined$ssh<-paste(combined$substrate_type, combined$season, combined$habitat, sep="_")
combined$ssh <- factor(combined$ssh, levels = c("water_autumn_eelgrass","water_autumn_rocks","water_autumn_sand","sediment_autumn_eelgrass",
"sediment_autumn_rocks","sediment_autumn_sand","sediment_spring_eelgrass",
"sediment_spring_rocks","sediment_spring_sand","water_spring_eelgrass","water_spring_rocks",
"water_spring_sand"))

#Need to fix to keep right coordinate mapping

#combined<-ddply(combined, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combined$max_rich<-max(rich)

####Also for metazoans and non-metazoans

#Meta
combinedm$Latitude<-merged_wat$Latitude[match(combinedm$snch, merged_wat$snch)]
combinedm$Longitude<-merged_wat$Longitude[match(combinedm$snch, merged_wat$snch)]

head(combinedm)

combinedm$habitat <- factor(combinedm$habitat, levels = c("sand", "rocks", "eelgrass"))
combinedm$season <- factor(combinedm$season, levels = c("spring", "autumn"))
combinedm$substrate_type <- factor(combinedm$substrate_type, levels = c("sediment", "water"))
combinedm$ss<-paste(combinedm$substrate_type, combinedm$season,  sep="_")
combinedm$ss <- factor(combinedm$ss, levels = c("water_autumn","sediment_autumn","sediment_spring","water_spring"))
combinedm$stsc<-paste(combinedm$substrate_type, combinedm$season, combinedm$cluster, sep="_")
combinedm$ssh<-paste(combinedm$substrate_type, combinedm$season, combinedm$habitat, sep="_")
combinedm$ssh <- factor(combinedm$ssh, levels = c("water_autumn_eelgrass","water_autumn_rocks","water_autumn_sand","sediment_autumn_eelgrass",
"sediment_autumn_rocks","sediment_autumn_sand","sediment_spring_eelgrass",
"sediment_spring_rocks","sediment_spring_sand","water_spring_eelgrass","water_spring_rocks",
"water_spring_sand"))

combinedm$max_rich<-max(richm)

#Non-meta
combinednm$Latitude<-merged_wat$Latitude[match(combinednm$snch, merged_wat$snch)]
combinednm$Longitude<-merged_wat$Longitude[match(combinednm$snch, merged_wat$snch)]

head(combinednm)

combinednm$habitat <- factor(combinednm$habitat, levels = c("sand", "rocks", "eelgrass"))
combinednm$season <- factor(combinednm$season, levels = c("spring", "autumn"))
combinednm$substrate_type <- factor(combinednm$substrate_type, levels = c("sediment", "water"))
combinednm$ss<-paste(combinednm$substrate_type, combinednm$season,  sep="_")
combinednm$ss <- factor(combinednm$ss, levels = c("water_autumn","sediment_autumn","sediment_spring","water_spring"))
combinednm$stsc<-paste(combinednm$substrate_type, combinednm$season, combinednm$cluster, sep="_")
combinednm$ssh<-paste(combinednm$substrate_type, combinednm$season, combinednm$habitat, sep="_")
combinednm$ssh <- factor(combinednm$ssh, levels = c("water_autumn_eelgrass","water_autumn_rocks","water_autumn_sand","sediment_autumn_eelgrass",
"sediment_autumn_rocks","sediment_autumn_sand","sediment_spring_eelgrass",
"sediment_spring_rocks","sediment_spring_sand","water_spring_eelgrass","water_spring_rocks",
"water_spring_sand"))

combinednm$max_rich<-max(richnm)

##Working on Metazoan separated plots

tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")

myshapefile_df <- fortify(so_ne_10, region ="NAME_0")

colnames(combinedm)[1]<-"ASVs"
colnames(combinednm)[1]<-"ASVs"

colnames(combined)[1]<-"ASVs"

hbts<-as.character(unique(combinedm$habitat))
for (i in 1:length(hbts))
{
hobj<-subset(combinedm, habitat==hbts[i])
ggplot(data = myshapefile_df) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray90") + geom_point(data=hobj, aes(Longitude, Latitude, fill=ASVs), shape=21, size=4.8) + geom_text(data=hobj, aes(Longitude, Latitude, label=ASVs), size=2) + scale_fill_gradient(low = "yellow", high = "red", limits=c(1, unique(combinedm$max_rich))) + theme_bw() + facet_wrap(substrate_type~season, ncol=2)
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
ggsave(paste("Metazoan",hbts[i],"ASV_richness.pdf", sep="_"))
}

hbts<-as.character(unique(combinednm$habitat))
for (i in 1:length(hbts))
{
hobnj<-subset(combinednm, habitat==hbts[i])
ggplot(data = myshapefile_df) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray90") + geom_point(data=hobnj, aes(Longitude, Latitude, fill=ASVs), shape=21, size=5) + geom_text(data=hobnj, aes(Longitude, Latitude, label=ASVs), size=1.8) + scale_fill_gradient(low = "yellow", high = "red", limits=c(1, unique(combinednm$max_rich))) + theme_bw() + facet_wrap(substrate_type~season, ncol=2)
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
ggsave(paste("Non_Metazoan",hbts[i],"ASV_richness.pdf", sep="_"))
}

hbtsc<-as.character(unique(combined$habitat))
for (i in 1:length(hbtsc))
{
hobcj<-subset(combined, habitat==hbtsc[i])
ggplot(data = myshapefile_df) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray90") + geom_point(data=hobcj, aes(Longitude, Latitude, fill=ASVs), shape=21, size=5) + geom_text(data=hobcj, aes(Longitude, Latitude, label=ASVs), size=1.8) + scale_fill_gradient(low = "yellow", high = "red", limits=c(1, unique(combined$max_rich))) + theme_bw() + facet_wrap(substrate_type~season, ncol=2)
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
ggsave(paste(hbtsc[i],"ASV_richness.pdf", sep="_"))
}


#######
#Alternative for field replicates (ASV level here - see loop for alternative tax level below)

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

########

jt<-data.frame(sample_data(DADAwang1))
otugr<-otu_table(DADAwang1)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr != 0)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

#Create reference latitude for cluster

#Get map
tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))


ggplot(combinedr2spring, aes(x = ssh, y = rich, fill = ssh)) + geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.6) +
geom_col(position="identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = c("red2","indianred1","darkgreen","lightgreen","dodgerblue2","lightskyblue1")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing7_spring.pdf")

ggplot(combinedr2autumn, aes(x = ssh, y = rich, fill = ssh)) + geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.6) +
geom_col(position="identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = c("red2","indianred1","darkgreen","lightgreen","dodgerblue2","lightskyblue1")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing7_autumn.pdf")


color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))


df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)]) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.5, y = cl_Latitude - 0.5, xmax = cl_Longitude + 0.5, ymax = cl_Latitude + 0.5)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave("ASV_richness_spring.pdf")

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)]) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.5, y = cl_Latitude - 0.5, xmax = cl_Longitude + 0.5, ymax = cl_Latitude + 0.5)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave("ASV_richness_autumn.pdf")

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = color_table$Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("ASV_richness_legend.pdf")


##############
#Loop for geting richness plots at different taxonomic levels
##############
#scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)])

#Get map
tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)

tax_level<-c("Phylum","Class","Order","Family")
for (i in 1:length(tax_level))
{

tax_def<-tax_glom(DADAwang1, tax_level[i])
jt<-data.frame(sample_data(tax_def))
otugr<-otu_table(tax_def)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr != 0)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))

color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))

Color<-c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1")
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich/rich_tax_lev")

df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"lev_rich_spring.pdf",sep="_"))

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"lev_rich_autumn.pdf",sep="_"))

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave(paste(tax_level[i],"lev_rich_legend.pdf",sep="_"))
}


####
####
#Loop within Supergroups read totals
##############
#scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)])

#Get map

tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)
sg_def<-tax_glom(DADAwang1, "Supergroup")
tax_sg<-tax_table(sg_def)
tax_level<-as.vector(unique(tax_sg)[,"Supergroup"])

for (i in 1:length(tax_level))
{
tax_def<-subset_taxa(sg_def, Supergroup==tax_level[i])
jt<-data.frame(sample_data(tax_def))
otugr<-otu_table(tax_def)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))

color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))

Color<-c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1")
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich/reads_Supergroup")

size_ref<-ifelse(tax_level[i]=="Opisthokonta", 0.4, ifelse(tax_level[i]=="Stramenopiles"|tax_level[i]=="Alveolata"|tax_level[i]=="Archaeplastida", 0.6, ifelse(tax_level[i]=="Hacrobia", 0.8, 1)))
df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) +
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - size_ref, y = cl_Latitude - size_ref, xmax = cl_Longitude + size_ref, ymax = cl_Latitude + size_ref)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"spring_reads.pdf",sep="_"))

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - size_ref, y = cl_Latitude - size_ref, xmax = cl_Longitude + size_ref, ymax = cl_Latitude + size_ref)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"autumn_reads.pdf",sep="_"))

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave(paste(tax_level[i],"legend_reads.pdf",sep="_"))
}


####

####
####
#Loop within Supergroups ASV richness
##############
#scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)])

#Get map
tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)
sg_def<-DADAwang1
tax_sg<-tax_table(sg_def)
tax_level<-as.vector(unique(tax_sg[,"Supergroup"]))

for (i in 1:length(tax_level))
{
tax_def<-subset_taxa(sg_def, Supergroup==tax_level[i])
jt<-data.frame(sample_data(tax_def))
otugr<-otu_table(tax_def)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr != 0)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))

color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))

Color<-c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1")
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich/asv_rich_Supergroup")

df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"spring_ASV_rich.pdf",sep="_"))

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"autumn_ASV_rich.pdf",sep="_"))

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave(paste(tax_level[i],"legend_ASV_rich.pdf",sep="_"))
}


####




#Loop within Phylum read totals
##############
#scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)])

#Get map

tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)
po_data<-tax_glom(DADAwang1, "Phylum")

#Filter phyla with abundance lower than 10 in at least 30 % of samples
#Add this step to class modelling
po_data2 = filter_taxa(po_data, function(x) sum(x > 10) > (0.3*length(x)), TRUE)

#Filter the taxa using a variation coef cutoff of 3.0
po_data3F = filter_taxa(po_data2, function(x) sd(x)/mean(x) > 1.12, TRUE)

tudao1 = filter_taxa(po_data3F, function(x) sum(x) > 0, TRUE)
sg_def = prune_samples(sample_sums(tudao1)>0,tudao1)

tax_sg<-tax_table(sg_def)
tax_level<-as.vector(unique(tax_sg)[,"Phylum"])

for (i in 1:length(tax_level))
{
tax_def<-subset_taxa(sg_def, Phylum==tax_level[i])
jt<-data.frame(sample_data(tax_def))
otugr<-otu_table(tax_def)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))

color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))

Color<-c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1")
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/reads_Phylum")

#size_ref<-ifelse(tax_level[i]=="Opisthokonta", 0.4, ifelse(tax_level[i]=="Stramenopiles"|tax_level[i]=="Alveolata"|tax_level[i]=="Archaeplastida", 0.6, ifelse(tax_level[i]=="Hacrobia", 0.8, 1)))
df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) +
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.5, y = cl_Latitude - 0.5, xmax = cl_Longitude + 0.5, ymax = cl_Latitude + 0.5)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"spring_reads.pdf",sep="_"))

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.5, y = cl_Latitude - 0.5, xmax = cl_Longitude + 0.5, ymax = cl_Latitude + 0.5)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"autumn_reads.pdf",sep="_"))

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave(paste(tax_level[i],"legend_reads.pdf",sep="_"))
}


####

####
####
#Loop within Supergroups ASV richness
##############
#scale_fill_manual(values = color_table$Color[match(.$sth, color_table$sth)])

#Get map
tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 15.5, 54.5, 57.8)
so_ne_10 <- crop(ctry_shps, so_dn)
sg_def<-DADAwang1
tax_sg<-tax_table(sg_def)
tax_level<-as.vector(unique(tax_sg[,"Supergroup"]))

for (i in 1:length(tax_level))
{
tax_def<-subset_taxa(sg_def, Supergroup==tax_level[i])
jt<-data.frame(sample_data(tax_def))
otugr<-otu_table(tax_def)
otugr<-t(data.frame(otugr, check.names=F))
chr<-jt$sample_root

richr<-colSums(otugr != 0)
combinedr<-cbind(richr, jt)
combinedr$snch<-paste(combinedr$season,combinedr$cluster,combinedr$habitat,sep="_")

combinedr2<-ddply(combinedr, .(sshc, snch), summarise,  rich=mean(richr), rich_sd=sd(richr))

combinedr2$Latitude<-merged_wat$Latitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$Longitude<-merged_wat$Longitude[match(combinedr2$snch, merged_wat$snch)]
combinedr2$cluster<-combinedr$cluster[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$season<-combinedr$season[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$habitat<-combinedr$habitat[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$substrate_type<-combinedr$substrate_type[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssc<-combinedr$ssc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ss<-combinedr$ss[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$ssh<-combinedr$ssh[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$stsc<-combinedr$stsc[match(combinedr2$sshc, combinedr$sshc)]
combinedr2$snch<-combinedr$snch[match(combinedr2$sshc, combinedr$sshc)]

combinedr2$habitat <- factor(combinedr2$habitat, levels = c("eelgrass", "rocks", "sand"))
combinedr2$season <- factor(combinedr2$season, levels = c("spring", "autumn"))
combinedr2$ss<-paste(combinedr2$substrate_type, combinedr2$season,  sep="_")
combinedr2$stsc<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$cluster, sep="_")
combinedr2$ssh<-paste(combinedr2$substrate_type, combinedr2$season, combinedr2$habitat, sep="_")
combinedr2$ssh <- factor(combinedr2$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass","sediment_autumn_rocks",
"water_autumn_rocks","sediment_autumn_sand","water_autumn_sand",
"sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2<-ddply(combinedr2, .(cluster), transform, cl_Latitude=mean(Latitude), cl_Longitude=mean(Longitude))

combinedr2$max_rich<-max(combinedr2$rich)+max((combinedr2$rich_sd[-238]))

combinedr2spring<-subset(combinedr2, season=="spring")

combinedr2autumn<-subset(combinedr2, season=="autumn")

combinedr2spring$ssh <- factor(combinedr2spring$ssh, levels = c("sediment_spring_eelgrass","water_spring_eelgrass",
"sediment_spring_rocks","water_spring_rocks","sediment_spring_sand","water_spring_sand"))

combinedr2autumn$ssh <- factor(combinedr2autumn$ssh, levels = c("sediment_autumn_eelgrass","water_autumn_eelgrass",
"sediment_autumn_rocks","water_autumn_rocks","sediment_autumn_sand",
"water_autumn_sand"))

combinedr2spring$sth<-paste(combinedr2spring$substrate_type, combinedr2spring$habitat,  sep="_")
combinedr2autumn$sth<-paste(combinedr2autumn$substrate_type, combinedr2autumn$habitat,  sep="_")

combinedr2spring$sth <- factor(combinedr2spring$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand","water_sand"))

combinedr2autumn$sth <- factor(combinedr2autumn$sth, levels = c("sediment_eelgrass","water_eelgrass","sediment_rocks",
"water_rocks","sediment_sand",
"water_sand"))

color_table <- tibble(sth = c("sediment_eelgrass","water_eelgrass",
"sediment_rocks","water_rocks","sediment_sand","water_sand"),
Color = c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1"))

Color<-c("red2","indianred1","springgreen4","lightgreen","dodgerblue2","lightskyblue1")
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich/asv_rich_Supergroup")

df.grobs <- combinedr2spring %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"spring_ASV_rich.pdf",sep="_"))

df.grobs <- combinedr2autumn %>% group_by(cluster, cl_Longitude, cl_Latitude, sshc) %>% mutate(sth = factor(sth, levels = color_table$sth)) %>%
do(subplots = ggplot(.,
aes(x = sth, y = rich, fill = sth)) +
geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.5) + 
geom_col(position = "identity", width = 0.9, color = "black",
alpha = 1, size=0.3, show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = unique(color_table$Color[match(.$sth, color_table$sth)])) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.4, y = cl_Latitude - 0.4, xmax = cl_Longitude + 0.4, ymax = cl_Latitude + 0.4)))

d<-ggplot(so_ne_10) + geom_path(aes(long, lat, group=group), size=0.1) + geom_polygon(aes(long, lat, group=group), size=0.5, fill="gray85") + theme_bw() +
theme(axis.ticks = element_blank(), panel.background = element_rect(colour = "gray95", fill="gray95"), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"))
d+df.grobs$subgrobs
ggsave(paste(tax_level[i],"autumn_ASV_rich.pdf",sep="_"))

ggplot(combinedr2autumn, aes(x = sth, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = Color) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave(paste(tax_level[i],"legend_ASV_rich.pdf",sep="_"))
}


####



ggplot() +
geom_path(data=mymaps, aes(long, lat, group=group), size=0.2) +
geom_point(data=subset(combined, substrate_type=="sediment"), aes(Longitude, Latitude, color=habitat, size=rich), alpha=0.6, stroke=0.8) +
scale_shape_manual(values=c(21)) +
coord_map("polyconic") +
facet_wrap(habitat~season, ncol=2) +
labs(x=NULL, y=NULL) +
theme_bw() +
theme(axis.text=element_blank(),strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
theme(legend.position="bottom") +
theme(axis.ticks = element_blank()) +
scale_size_binned()

ggsave("map_richness_overall_sediment_ncbi.pdf")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")
ggplot() +
geom_path(data=mymaps, aes(long, lat, group=group), size=0.2) +
geom_point(data=subset(combined, substrate_type=="water"), aes(Longitude, Latitude, color=habitat, size=rich), alpha=0.6, stroke=0.8) +
scale_shape_manual(values=c(21)) +
coord_map("polyconic") +
facet_wrap(habitat~season, ncol=2) +
labs(x=NULL, y=NULL) +
theme_bw() +
theme(axis.text=element_blank(),strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
theme(legend.position="bottom")+
theme(axis.ticks = element_blank()) +
scale_size_binned()

ggsave("map_richness_overall_water_ncbi.pdf")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")
ggplot() +
geom_path(data=mymaps, aes(long, lat, group=group), size=0.2) +
geom_point(data=combined, aes(Longitude, Latitude, color=habitat, size=rich, shape=substrate_type), alpha=0.6, stroke=0.8) +
scale_shape_manual(values=c(0,1)) +
coord_map("polyconic") +
facet_wrap(habitat~season, ncol=2) +
labs(x=NULL, y=NULL) +
theme_bw() +
theme(axis.text=element_blank(),strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
theme(legend.position="bottom")+
theme(axis.ticks = element_blank()) +
scale_size_binned()

ggsave("map_richness_overall_ncbi.pdf")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")

df.grobs <- combined %>% group_by(cluster, cl_Longitude, cl_Latitude, max_rich) %>%
do(subplots = ggplot(.,
aes(x = ssh, y = rich, fill = habitat)) +
geom_col(position = "identity", width = 1, color = "black",
alpha = 0.7,
show.legend = FALSE) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0, unique(.$max_rich))) +
scale_fill_manual(values = c("royalblue4", "snow", "sandybrown")) +
coord_polar() + 
theme_void()) %>% mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots), x = cl_Longitude - 0.3, y = cl_Latitude - 0.3, xmax = cl_Longitude + 0.3, ymax = cl_Latitude + 0.3)))

d<- ggplot() + geom_path(data=mymaps, aes(long, lat, group=group), size=0.2)
d+df.grobs$subgrobs
ggsave("testing5.pdf")

########

dodge_t <- position_dodge(width=0.8)

ggplot(combinedr2, aes(x = ssh, y = rich, fill = ssh)) + geom_errorbar(aes(ymax=rich+rich_sd, ymin=rich), size=0.2, width=0.6) +
geom_col(position="identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = c("white","black","white","black","white","black","white","black","white","black","white","black")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing6.pdf")

ggplot(combined, aes(x = ss, y = rich, fill = habitat)) +
geom_col(position = "identity", width = 1, color = "black", alpha=0.5) +
scale_fill_manual(values = c("grey", "green", "red")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing.pdf")

ggplot(combined, aes(x = ss, y = rich, fill = habitat)) +
geom_col(position = "identity", width = 1, color = "black", alpha=0.8) +
scale_fill_manual(values = c("yellow", "green", "red")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing1.pdf")

ggplot(combined, aes(x = ssh, y = rich, fill = habitat)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = c("royalblue4", "snow", "sandybrown")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing3.pdf")

ggplot(combined, aes(x = ssh, y = rich, fill = ssh)) +
geom_col(position = "identity", width = 1, color = "black", alpha=1) +
scale_fill_manual(values = c("lightskyblue1","darkolivegreen1","coral","royalblue2",
"seagreen3","orangered4","royalblue2", "seagreen3","orangered4", "lightskyblue1","darkolivegreen1","coral")) +
coord_polar() + facet_wrap(~cluster, ncol=5) + theme_void()
ggsave("testing4.pdf")

########
#########

#LEFTOVERS


p_md_wat<- subset(na.omit(doi), substrate_type=="water")
md_wat<-unique(p_md_wat[c("season", "cluster", "habitat")])
md_wat$Sample_ID<-ifelse(md_wat$season=="spring",paste("C",md_wat$cluster,md_wat$habitat,sep=""),
paste("2C",md_wat$cluster,md_wat$habitat,sep=""))
doi$Salinity<-merged_wat$Salinity[match(doi$snch, merged_wat$snch)]
doi$PO4<-merged_wat$PO4[match(doi$snch, merged_wat$snch)]
for(i in 1:nrow(doi)){doi[i,"PO4"]<-ifelse(doi[i,"PO4"]<=0,0,doi[i,"PO4"])}

#Import sediment organic content
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Environmental_data/data")
TP_sed_autunm<-read.table("TP_org_inorg_dens_sed_autumn.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
TP_sed_spring<-read.table("TP_org_inorg_dens_sed_spring.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
TP_sed_both<-rbind(TP_sed_autunm, TP_sed_spring)
TP_sed_both$po<- sapply(strsplit(as.character(TP_sed_both$Sample_ID), "2C"), tail, 1)
TP_sed_both$po2<- sapply(strsplit(as.character(TP_sed_both$po), "C"), tail, 1)
TP_sed_both$pcl<-gsub('\\D','_', TP_sed_both$po2)
TP_sed_both$cl<- sapply(strsplit(as.character(TP_sed_both$pcl), "__"), head, 1)
TP_sed_both$pn<-gsub('\\d','', TP_sed_both$po2)
TP_sed_both$hb<-ifelse(TP_sed_both$pn=="EB", "eelgrass", ifelse(TP_sed_both$pn=="RB", "rocks", "sand"))
TP_sed_both$poi<- sapply(strsplit(as.character(TP_sed_both$Sample_ID), "2C"), head, 1)
TP_sed_both$sn<-ifelse(TP_sed_both$poi=="", "autumn", "spring")
TP_sed_both$snch<-paste(TP_sed_both$sn,TP_sed_both$cl,TP_sed_both$hb,sep="_")

OC_f_wat <- ddply(TP_sed_both, .(snch), summarise, grp.mean=mean(Organic_content))

doi$Organic_Content<-ifelse(doi$substrate_type=="sediment", TP_sed_both$Organic_content[match(doi$sample_root, TP_sed_both$Sample_ID)], OC_f_wat$grp.mean[match(doi$snch, OC_f_wat$snch)])
doi$Organic_Content<-round(doi$Organic_Content, digits=4)

##Update OTU table
with_NA<-rownames(doi[is.na(doi$Salinity),])
newddw<-subset_samples(DADAwang1, !(sample_root %in% with_NA))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

o_data<-tax_glom(fdt, "Class")
otuo<-data.frame(otu_table(o_data))
taxo<-data.frame(tax_table(o_data), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Class
ncol(otuo)
do<-data.frame(sample_data(o_data))
do$Salinity<-doi$Salinity[match(do$snch, doi$snch)]
do$PO4<-doi$PO4[match(do$snch, doi$snch)]
do$Organic_Content<-doi$Organic_Content[match(do$sample_root, doi$sample_root)]
do$habitat <- factor(do$habitat, levels = c("sand", "rocks", "eelgrass"))
do$season <- factor(do$season, levels = c("spring", "autumn"))
do$cluster<-as.character(do$cluster)
do[,11:13]<-scale(do[,11:13])

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

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/ncbi/rich")

ggplot(data = world) + geom_sf() +  coord_sf(xlim = c(54, 58), ylim = c(7, 16), expand = FALSE)
ggsave("testingmap.pdf")


pdf("testingmap.pdf")
plot(basemap)
dev.off()

basemap2<-getData('SRTM', lon=10, lat=55)

basemap <-getData('GADM', country='DNK', level=1)
basemap <- raster("eea_r_3035_100_m_clc12_V18_5_land_mask.tif")

lakes50 <- ne_download(scale = 50, type = 'lakes', category = 'physical')
rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
urbanareas50 <- ne_download(scale = 50, type = 'urban_areas', category = 'cultural')
populatedplaces50 <- ne_download(scale = 50, type = 'populated_places', category = 'cultural')
coastline50 <- ne_download(scale = 50, type = 'coastline', category = 'physical')

st_agr(lakes50) <- 'constant'


so_lakes50 <- crop(lakes50, so_dn)
so_rivers50 <- crop(rivers50, so_dn)
so_urbanareas50 <- crop(urbanareas50, so_dn)
so_populatedplaces50 <- crop(populatedplaces50, so_dn)
so_coastline50 <- crop(coastline50, so_dn)


tre = list()
tre$countries = c("DNK","DEU","SWE")
ctry_shps = do.call("bind", lapply(tre$countries, function(x) getData('GADM', country=x, level=0)))

so_dn <- extent(7, 16, 54, 58)
so_ne_10 <- crop(ctry_shps, so_dn)

pdf("testingmap.pdf")
plot(so_ne_10)
dev.off()
