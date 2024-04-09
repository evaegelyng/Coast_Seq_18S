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

#Remove Eelgrass reads

#Subset by substrate type & remove cluster 2
pDADAwang1s<-subset_samples(DADAwang1, substrate_type=="sediment")
DADAwang1s<-subset_samples(pDADAwang1s, !cluster==2)
doi<-data.frame(sample_data(DADAwang1s))

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

#what is this for?
p_md_wat<- subset(na.omit(doi), substrate_type=="water")
md_wat<-unique(p_md_wat[c("season", "cluster", "habitat")])
md_wat$Sample_ID<-ifelse(md_wat$season=="spring",paste("C",md_wat$cluster,md_wat$habitat,sep=""),
paste("2C",md_wat$cluster,md_wat$habitat,sep=""))
#

doi$Salinity<-merged_wat$Salinity[match(doi$snch, merged_wat$snch)]

#Import sediment data
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

doi$Organic_content<-ifelse(doi$substrate_type=="sediment", TP_sed_both$Organic_content[match(doi$sample_root, TP_sed_both$Sample_ID)], OC_f_wat$grp.mean[match(doi$snch, OC_f_wat$snch)])
doi$Organic_content<-round(doi$Organic_content, digits=5)

##Update OTU table
with_NA1<-rownames(doi[is.na(doi$Salinity),])
newddw<-subset_samples(DADAwang1s, !(sample_root %in% with_NA1))
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE)
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)
po_dataF<-tax_glom(fdt, "Phylum")

#Remove "rare" and low-frequency phyla, using mean total phyla richness, abundance and number of sites present
#Richness
otuo<-data.frame(otu_table(fdt))
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
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

#Rel_abund
datarg = transform_sample_counts(po_dataF, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]
ttab<-t(data.frame(tab, check.names=F))

#Filtering by rich, abund and sites occurring
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
abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

list_top15<-c(abundsorted$clades_p[1:15],richsorted$clades_p[1:15],sitessorted$clades_p[1:15])
phy_to_remove<-as.character(unique(list_top15))

po_data1F = subset_taxa(po_dataF, !(Phylum %in% phy_to_remove))

#Remove Picozoa
#po_data1 = subset_taxa(po_data1F, !Phylum=="Picozoa")

#Filter phyla with abundance lower than 10 in at least 10 % of samples
#Add this step to class modelling
po_data2 = filter_taxa(po_data1F, function(x) sum(x > 10) > (0.1*length(x)), TRUE)

#Filter the taxa using a variation coef cutoff of 3.0
#po_data3F = filter_taxa(po_data2, function(x) sd(x)/mean(x) > 1.2, TRUE)

#Keep only top 30 phyla (richness)
toprichsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"], decreasing=T ),]
list_top30<-c(toprichsorted$clades_p[1:30])
po_data2.1 = subset_taxa(po_data2, (Phylum %in% list_top30))

tudao1 = filter_taxa(po_data2.1, function(x) sum(x) > 0, TRUE)
o_dataF = prune_samples(sample_sums(tudao1)>0,tudao1)

taxoF<-data.frame(tax_table(o_dataF), stringsAsFactors=FALSE)
phy_to_keep<-taxoF$Phylum

po_data2 = subset_taxa(fdt, Phylum %in% phy_to_keep)

tudao1 = filter_taxa(po_data2, function(x) sum(x) > 0, TRUE)
o_data = prune_samples(sample_sums(tudao1)>0,tudao1)

otuo<-data.frame(otu_table(o_data))
taxo<-data.frame(tax_table(o_data), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$Phylum
ncol(otuo)
do<-data.frame(sample_data(o_data))
do$Salinity<-doi$Salinity[match(do$snch, doi$snch)]
do$Organic_content<-doi$Organic_content[match(do$sample_root, doi$sample_root)]

do$habitat <- factor(do$habitat, levels = c("sand", "rocks", "eelgrass"))
do$season <- factor(do$season, levels = c("spring", "autumn"))
do$cluster<-as.character(do$cluster)
do[,11:12]<-scale(do[,11:12])

#Mk dummy variables
head(do)
do$season<-as.character(do$season)

for (h in 1:nrow(do))
{
if(do[h,"season"]=="spring") {
do[h,"season"]<-as.numeric(0)
} else {
do[h,"season"]<-as.numeric(1)}

if(do[h,"habitat"]=="rocks") {
do[h,"habitat_rocks"]<-as.numeric(1)
} else {
do[h,"habitat_rocks"]<-as.numeric(0)}

if(do[h,"habitat"]=="eelgrass") {
do[h,"habitat_eel"]<-as.numeric(1)
} else {
do[h,"habitat_eel"]<-as.numeric(0)}

}

do$season<-as.numeric(do$season)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/regs/gllvm/sed")

#tranformation
#otuo<-log1p(otuo)
#otuo<-round(otuo,digits=4)

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

######Test to gaussian
#rich0<-cmultRepl(rich_asv, output = "p-counts")
#pr <- propr(counts = rich0, metric = "rho", ivar = "clr")
#trans_rich<-as.data.frame(pr@logratio)
######


fit_0<-gllvm(rich_asv, family = "negative.binomial", num.lv = 3, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")

fit_0.2<-gllvm(rich_asv, family = "negative.binomial", num.lv = 2, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")

fit_0.1<-gllvm(rich_asv, family = "negative.binomial", num.lv = 1, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")


#num.lv = 2, method = "EVA", sd.errors=F, control = list(max.iter = 500, maxit = 4000), #control.start = list(starting.val = "res", n.init = 3, jitter.var=0.4))
###!!###$$###!!###

#Best model
fit_e <- gllvm(rich_asv, do, num.lv = 3,
formula = ~ season + habitat_rocks + habitat_eel + Salinity + Organic_content + season * (habitat_rocks + habitat_eel + Salinity + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#Best model
fit_e.2 <- gllvm(rich_asv, do, num.lv = 2,
formula = ~ season + habitat_rocks + habitat_eel + Salinity + Organic_content + season * (habitat_rocks + habitat_eel + Salinity + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#Best model
fit_e.1 <- gllvm(rich_asv, do, num.lv = 1,
formula = ~ season + habitat_rocks + habitat_eel + Salinity + Organic_content + season * (habitat_rocks + habitat_eel + Salinity + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e2 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season + habitat_rocks + habitat_eel + Salinity + season * (habitat_rocks + habitat_eel + Salinity), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e3 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season + habitat_rocks + habitat_eel + Organic_content + season * (habitat_rocks + habitat_eel + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e4 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season + Salinity + Organic_content + season * (Salinity + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e5 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season + habitat_rocks + habitat_eel + Salinity + Organic_content, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e6 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ habitat_rocks + habitat_eel + Salinity + Organic_content, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e7 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season + habitat_rocks + habitat_eel + Salinity, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e8 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ habitat_rocks + habitat_eel + season, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e9 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ habitat_rocks + habitat_eel, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e10 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ season, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#
fit_e11 <- gllvm(rich_asv, do, num.lv = 3,formula = ~ Salinity, row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

fit_0
rescov0 <- getResidualCov(fit_0)
rescov0$trace
rescov0$var.q

fit_0.2
rescov0.2 <- getResidualCov(fit_0.2)
rescov0.2$trace
rescov0.2$var.q

fit_0.1
rescov0.1 <- getResidualCov(fit_0.1)
rescov0.1$trace
rescov0.1$var.q


fit_e
rescov <- getResidualCov(fit_e)
rescov$trace
rescov$var.q
1 - rescov$trace / rescov0$trace

fit_e.2
rescov.2 <- getResidualCov(fit_e.2)
rescov.2$trace
rescov.2$var.q
1 - rescov.2$trace / rescov0$trace

fit_e.1
rescov.1 <- getResidualCov(fit_e.1)
rescov.1$trace
rescov.1$var.q
1 - rescov.1$trace / rescov0$trace

fit_e2
rescov2 <- getResidualCov(fit_e2)
rescov2$trace
rescov2$var.q
1 - rescov2$trace / rescov0$trace

fit_e3
rescov3 <- getResidualCov(fit_e3)
rescov3$trace
rescov3$var.q
1 - rescov3$trace / rescov0$trace

fit_e4
rescov4 <- getResidualCov(fit_e4)
rescov4$trace
rescov4$var.q
1 - rescov4$trace / rescov0$trace

fit_e5
rescov5 <- getResidualCov(fit_e5)
rescov5$trace
rescov5$var.q
1 - rescov5$trace / rescov0$trace

fit_e6
rescov6 <- getResidualCov(fit_e6)
rescov6$trace
rescov6$var.q
1 - rescov6$trace / rescov0$trace

fit_e7
rescov7 <- getResidualCov(fit_e7)
rescov7$trace
rescov7$var.q
1 - rescov7$trace / rescov0$trace

fit_e8
rescov8 <- getResidualCov(fit_e8)
rescov8$trace
rescov8$var.q
1 - rescov8$trace / rescov0$trace

fit_e9
rescov9 <- getResidualCov(fit_e9)
rescov9$trace
rescov9$var.q
1 - rescov9$trace / rescov0$trace

fit_e10
rescov10 <- getResidualCov(fit_e10)
rescov10$trace
rescov10$var.q
1 - rescov10$trace / rescov0$trace

fit_e11
rescov11 <- getResidualCov(fit_e11)
rescov11$trace
rescov11$var.q
1 - rescov11$trace / rescov0$trace

cat("NOW BACK TO BEST MODEL")

summary(fit_e)
getLV(fit_e)
coef(fit_e)
fit_e$sd

summary(fit_e.2)
summary(fit_e.1)

adf<-fit_e$lvs

pdf("phylum_sed_diagn_nb_rich_null.pdf")
plot(fit_0, mfrow = c(3,2),var.colors = 1) 
dev.off()

pdf("phylum_sed_diagn_nb_rich_full.pdf")
plot(fit_e, mfrow = c(3,2),var.colors = 1) 
dev.off()

#New ugly (but correct) plots

rbPal <- c("#D2B48C", "#5F9EA0")
# "#8B4513", "#D2691E", "#006400", "#32CD32","#B266FF", "#FF66FF")
hbPal <- c("#6666FF","#FF8000","#00CC00")
sbPal2 <- c("#66FFFF","#66CC00","#CCCC00","#FF9933","#FF0000")
obPal <- c("#DEB887","#F4A460","#D2691E","#A0522D","#8B4513")

#Inspect colour breaks in continuous variables
jajas<-doi$Salinity[match(do$snch, doi$snch)]
binx<-as.numeric(cut(jajas, breaks = 5))
jajao<-log10(doi$Organic_content[match(do$snch, doi$snch)])
binx<-as.numeric(cut(jajao, breaks = 5))
sort(jajao)
sort(binx)

dp<-do
dp$legend<-NA
dp$legend<-ifelse(dp$season==0,as.numeric(0),as.numeric(1))
dp$pchh<-ifelse(dp$habitat=="sand",as.numeric(0),ifelse(dp$habitat=="rocks",as.numeric(1),as.numeric(2)))

pdf("ordi_phylum_nb_sed_rich_null3r_LVs12_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0, symbols = T, s.colors = Col, main = "season", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colo <- obPal[as.numeric(cut(jajao, breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Colo, main = "organic content", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("0-0.004","0.004-0.008","0.008-0.017","0.017-0.035","0.035-0.07"), col=obPal, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_sed_rich_null3r_LVs13_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0, symbols = T, s.colors = Col, main = "season", cex.spp=0.4, which.lvs = c(1,3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.4, which.lvs = c(1,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colo <- obPal[as.numeric(cut(jajao, breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Colo, main = "organic content", cex.spp=0.4, which.lvs = c(1,3))
legend("topleft", legend = c("0-0.004","0.004-0.008","0.008-0.017","0.017-0.035","0.035-0.07"), col=obPal, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.4, which.lvs = c(1,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 


pdf("ordi_phylum_nb_sed_rich_null3r_LVs23_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0, symbols = T, s.colors = Col, main = "season", cex.spp=0.4, which.lvs = c(2,3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.4, which.lvs = c(2,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colo <- obPal[as.numeric(cut(jajao, breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Colo, main = "organic content", cex.spp=0.4, which.lvs = c(2,3))
legend("topleft", legend = c("0-0.004","0.004-0.008","0.008-0.017","0.017-0.035","0.035-0.07"), col=obPal, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.4, which.lvs = c(2,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 


#ggplots
caca<-cbind(adf,do)
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segss <- merge(caca, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segss$sh<-paste(segss$season, segss$habitat)

ggplot(data=segss, aes(x=LV1, y=LV2)) +
geom_point(aes(x=LV1, y=LV2, color=habitat, shape=factor(season)), size=0.8, alpha=0.7) + 
geom_segment(aes(xend = cLV1, yend = cLV2, colour = habitat), alpha=0.9, size=0.2) + geom_point(aes(cLV1, cLV2, colour = habitat, shape=factor(season)), size = 2.5, alpha=0.7) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_colour_manual(values = c("blue1","darkgreen","red1")) +
theme_bw() +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca$LV1), max(caca$LV1)), ylim = c(min(caca$LV2),max(caca$LV2)))
ggsave("ordi_phylum_nb_sed_rich_full_old.pdf")

ggplot(data=segss, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = habitat), alpha=0.3, size=0.05) + geom_point(aes(cLV1, cLV2, colour = habitat, shape=factor(season)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_colour_manual(values = c("blue1","darkgreen","red1")) +
theme_bw() +
scale_shape_manual(values = c(1,2)) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca$LV1), max(caca$LV1)), ylim = c(min(caca$LV2),max(caca$LV2)))
ggsave("ordi_phylum_nb_sed_rich_full_old2.pdf")

ad0<-fit_0$lvs
caca0<-cbind(ad0,do)
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_point(aes(x=LV1, y=LV2, color=habitat, shape=factor(season)), size=0.8, alpha=0.7) + 
geom_segment(aes(xend = cLV1, yend = cLV2, colour = habitat), alpha=0.9, size=0.2) + geom_point(aes(cLV1, cLV2, colour = habitat, shape=factor(season)), size = 2.5, alpha=0.7) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_colour_manual(values = c("blue1","darkgreen","red1")) +
theme_bw() +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_sed_rich_null_old.pdf")

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = habitat), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV2, colour = habitat, shape=factor(season)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_colour_manual(values = c("blue1","darkgreen","red1")) +
theme_bw() +
scale_shape_manual(values = c(1,2)) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_sed_rich_null_old2.pdf")


ds<-names(sort(rescov0$var.q, decreasing=T))[1:2]
pad0<-fit_0$lvs
ad0<-pad0[,c(ds)]
caca0<-cbind(ad0,do)
colnames(caca0)[1:2]<-c("LV1","LV3")
cent2 <- aggregate(cbind(LV1,LV3)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV3")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV3")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

ppexpvar<-sort(rescov0$var.q, decreasing=T)[1:2]
pexpvar<-100*(ppexpvar/rescov0$trace)
expvar<-round(pexpvar, 2)

ggplot(data=segs, aes(x=LV1, y=LV3)) +
geom_segment(aes(xend = cLV1, yend = cLV3, colour = sh), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV3, colour = sh), shape=1, size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV3, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
theme_bw() +
scale_colour_manual(values = c("red3","darkgreen","blue3","indianred1","lightgreen","lightskyblue1")) +
xlab(paste("LV1","-",expvar[1],"%",sep=" ")) + ylab(paste("LV3","-",expvar[2],"%",sep=" ")) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_sed_rich_null.pdf")

#ggsave("sed_sal_ordi_phylum_neg_b_gllvm_subs_habs_random_full_0_va.pdf")

pdf("sed_sal_phylum_env_effect_nb_rich_va_top30.pdf")
coefplot(fit_e, cex.ylab = 0.5, mar = c(4, 5.5, 2, 1), mfrow=c(1,2))
dev.off()

cr <- getResidualCor(fit_e)

pdf("sed_sal_phylum_nb_rich_corr_va_top30.pdf")
corrplot(cr[order.single(cr), order.single(cr)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()

cr0 <- getResidualCor(fit_0)

pdf("sed_sal_phylum_nb_rich_corr_0_va_top30.pdf")
corrplot(cr0[order.single(cr0), order.single(cr0)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()
# Environmental correlation
modelmat <- model.matrix(fit_e$formula, data = do) %>% 
   as.matrix %>% 
   {.[,-1]} # Remove intercept
linpred <- tcrossprod(modelmat, fit_e$params$Xcoef)
envircor <- cov2cor(cov(linpred))

pdf("sed_sal_phylum_nb_rich_envicorr_va_top30.pdf")
corrplot(envircor[order.single(envircor), order.single(envircor)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()

#Extracting coefs gllvm
p_summary<-as.data.frame(summary(fit_e)$Coef.tableX)
p_pred<-strsplit(rownames(summary(fit_e)$Coef.tableX), ":")
interC<-summary(fit_e, spp.intercepts = T)$Coefficients[,1]
pred<-sapply(p_pred, "[[", 1)
pred_int<-sapply(p_pred, "[[", 2)

#Sorted predictors
spred<-c("season","habitat_rocks","habitat_eel","Salinity","Organic_content",
"season*habitat_rocks","season*habitat_eel","season*Salinity","season*Organic_content")

#Here relies on the response variable table (rich|abund)
alltaxa<-pred_int[1:ncol(rich_asv)]
tlength<-length(spred)*ncol(rich_asv)
p_summary$taxa <- rep_len(alltaxa, length.out=tlength)
p_summary$covariate<-NA
tlength

#The "number of non-interaction terms" * number of taxa 
p_summary$covariate[1:150]<-pred[1:150]
p_summary$covariate[151:tlength]<-paste(pred[151:tlength], pred_int[151:tlength], sep="*")
colnames(p_summary)[4]<-"P_value"
p_summary$intercept <- rep_len(interC, length.out=tlength)

summary_f<-p_summary
summary_f$est_if_sig<-ifelse(summary_f$P_value <.05, summary_f$Estimate, NA)

#Sort covariate and taxa
summary_f$covariate <- factor(summary_f$covariate, levels = spred)
summary_f$taxa <- factor(summary_f$taxa, levels = c(rownames(envircor[order.single(envircor), order.single(envircor)])))
min_l<-min(summary_f$Estimate)
max_l<-max(summary_f$Estimate)

summary_f$taxa_int<-paste(summary_f$taxa," (",round(summary_f$intercept,2),")",sep=" ")

#Plot
ggplot(summary_f, aes(covariate, taxa, fill=Estimate, label=round(est_if_sig,2))) +
geom_tile(color = "black") +
labs(x = NULL, y = NULL, fill = "Covariate effect") +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3", limits=c(min_l,max_l)) +
geom_text(size=2) +
theme_classic() +
scale_x_discrete(expand=c(0,0)) +
scale_y_discrete(expand=c(0,0)) +
theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text=element_text(size=8))
ggsave("sed_sal_summary_heatmap_env_eff_phylum_nb_rich_va_top30.pdf")

#Division plot
summary_f$Division<-taxo$Division[match(summary_f$taxa, taxo$Phylum)]

#Metazoan plot
summary_fm<-subset(summary_f, Division=="Metazoa")

ggplot(summary_fm, aes(covariate, taxa, fill=Estimate, label=round(est_if_sig,2))) +
geom_tile(color = "black") +
labs(x = NULL, y = NULL, fill = "Covariate effect") +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3", limits=c(min_l,max_l)) +
geom_text(size=2) +
theme_classic() +
scale_x_discrete(expand=c(0,0)) +
scale_y_discrete(expand=c(0,0)) +
theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text=element_text(size=8))
ggsave("sed_sal_Meta_summary_heatmap_env_eff_phylum_nb_rich_va_top30.pdf")

#Non_Metazoan plot
summary_fnm<-subset(summary_f, !Division=="Metazoa")

ggplot(summary_fnm, aes(covariate, taxa, fill=Estimate, label=round(est_if_sig,2))) +
geom_tile(color = "black") +
labs(x = NULL, y = NULL, fill = "Covariate effect") +
scale_fill_gradient2(mid="white",low="firebrick3",high="royalblue3", limits=c(min_l,max_l)) +
geom_text(size=2) +
theme_classic() +
scale_x_discrete(expand=c(0,0)) +
scale_y_discrete(expand=c(0,0)) +
theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text=element_text(size=8))
ggsave("sed_sal_Non_meta_summary_heatmap_env_eff_phylum_nb_rich_va_top30.pdf")


anova(fit_0, fit_e)


