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

pDADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1<-subset_samples(pDADAwang1, !cluster==2)
DADAwang1

#Remove Eelgrass reads
doi<-data.frame(sample_data(DADAwang1))

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
doi$Temperature<-merged_wat$Temperature[match(doi$snch, merged_wat$snch)]
doi$PO4<-merged_wat$PO4[match(doi$snch, merged_wat$snch)]
doi$Chlorophyll<-merged_wat$Chlorophyll[match(doi$snch, merged_wat$snch)]
for(i in 1:nrow(doi)){doi[i,"PO4"]<-ifelse(doi[i,"PO4"]<=0,0,doi[i,"PO4"])}

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

TP_f_wat <- ddply(TP_sed_both, .(snch), summarise, grp.mean=mean(TP))

doi$TP<-ifelse(doi$substrate_type=="sediment", TP_sed_both$TP[match(doi$sample_root, TP_sed_both$Sample_ID)], TP_f_wat$grp.mean[match(doi$snch, TP_f_wat$snch)])
doi$TP<-round(doi$TP, digits=4)

WC_f_wat <- ddply(TP_sed_both, .(snch), summarise, grp.mean=mean(Watercontent))

doi$Water_content<-ifelse(doi$substrate_type=="sediment", TP_sed_both$Watercontent[match(doi$sample_root, TP_sed_both$Sample_ID)], WC_f_wat$grp.mean[match(doi$snch, WC_f_wat$snch)])
doi$Water_content<-round(doi$Water_content, digits=4)

IC_f_wat <- ddply(TP_sed_both, .(snch), summarise, grp.mean=mean(Inorganic_content))

doi$Inorganic_content<-ifelse(doi$substrate_type=="sediment", TP_sed_both$Inorganic_content[match(doi$sample_root, TP_sed_both$Sample_ID)], IC_f_wat$grp.mean[match(doi$snch, IC_f_wat$snch)])
doi$Inorganic_content<-round(doi$Inorganic_content, digits=5)

##Update OTU table
with_NA1<-rownames(doi[is.na(doi$Salinity),])
with_NA2<-rownames(doi[is.na(doi$Chlorophyll),])
newddw<-subset_samples(DADAwang1, !(sample_root %in% c(with_NA1, with_NA2)))
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
rich_asvF<-z[,-1]

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
prich<-as.data.frame(cbind(log10(colMeans(rich_asvF)), colnames(rich_asvF)))
colnames(prich)[2]<-"clades_p"
trich<-t(rich_asvF)
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

#Remove some more hard to model phyla

plusP<-c("Ctenophora","Hyphochytriomycota","Picozoa")
po_data1 = subset_taxa(po_dataF, !(Phylum %in% c(phy_to_remove, plusP)))

#Filter phyla with abundance lower than 10 in at least 10 % of samples
#Add this step to class modelling
po_data2F = filter_taxa(po_data1, function(x) sum(x > 10) > (0.1*length(x)), TRUE)

#Filter the taxa using a variation coef cutoff of 3.0
#po_data3F = filter_taxa(po_data2F, function(x) sd(x)/mean(x) > 1.12, TRUE)
#po_data2

tudao1F = filter_taxa(po_data2F, function(x) sum(x) > 0, TRUE)
o_dataF = prune_samples(sample_sums(tudao1F)>0,tudao1F)

#o_data

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
do$PO4<-doi$PO4[match(do$snch, doi$snch)]
do$Temperature<-doi$Temperature[match(do$snch, doi$snch)]
do$Chlorophyll<-doi$Chlorophyll[match(do$snch, doi$snch)]

do$Organic_content<-doi$Organic_content[match(do$sample_root, doi$sample_root)]
do$TP<-doi$TP[match(do$sample_root, doi$sample_root)]
do$Water_content<-doi$Water_content[match(do$sample_root, doi$sample_root)]
do$Inorganic_content<-doi$Inorganic_content[match(do$sample_root, doi$sample_root)]

do$habitat <- factor(do$habitat, levels = c("sand", "rocks", "eelgrass"))
do$season <- factor(do$season, levels = c("spring", "autumn"))
do$cluster<-as.character(do$cluster)
do[,11:18]<-scale(do[,11:18])

#Mk dummy variables
head(do)
do$season<-as.character(do$season)
do$substrate_type<-as.character(do$substrate_type)

for (h in 1:nrow(do))
{
if(do[h,"season"]=="spring") {
do[h,"season"]<-as.numeric(0)
} else {
do[h,"season"]<-as.numeric(1)}

if(do[h,"substrate_type"]=="sediment") {
do[h,"substrate_type"]<-as.numeric(0)
} else {
do[h,"substrate_type"]<-as.numeric(1)}

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
do$substrate_type<-as.numeric(do$substrate_type)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/regs/gllvm/both")

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
#doa<-data.frame(sample_data(o_data))
#cond<-as.vector(doa$substrate_type)
#x <-aldex.clr(otuo2, cond, denom = "iqlr")
#x <- aldex.clr(synth2, blocks, denom="lvha")
#x <- aldex.clr(otuo2, cond, denom = "all")
#rho <- propr::aldex2propr(x, how = "rho")
#trans_otuo<-as.data.frame(rho@logratio)
#control = list(reltol = 1e-7, optimizer = "nlminb", max.iter = 2000, maxit=6000)
#control = list(max.iter = 2000, maxit = 4000, optim.method = "BFGS")

fit_0<-gllvm(rich_asv, family = "negative.binomial", num.lv = 2, method = "VA", control.start = list(starting.val ="zero"))

fit_0.r<-gllvm(rich_asv, family = "negative.binomial", num.lv = 2, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")

fit_0.3<-gllvm(rich_asv, family = "negative.binomial", num.lv = 3, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")

fit_0.4<-gllvm(rich_asv, family = "negative.binomial", num.lv = 4, method = "VA", control.start = list(starting.val ="zero"), row.eff="random")

#num.lv = 2, method = "EVA", sd.errors=F, control = list(max.iter = 500, maxit = 4000), #control.start = list(starting.val = "res", n.init = 3, jitter.var=0.4))

#Best current model
#fit_e <- gllvm(rich_asv, do, num.lv = 2, formula = ~ season + habitat_rocks + habitat_eel + substrate_type + Salinity + Chlorophyll + Temperature + Water_content + Organic_content + substrate_type * ( habitat_rocks + habitat_eel + season + Salinity + Chlorophyll + Temperature + Water_content + Organic_content), row.eff = ~(1 | cluster), family = "negative.binomial", method = "VA", control.start = list(starting.val ="zero"))
###!!###$$###!!###

#fit_e
#rescov <- getResidualCov(fit_e)
#rescov$trace
#rescov$var.q

fit_0
rescov0 <- getResidualCov(fit_0)
rescov0$trace
rescov0$var.q

fit_0.r
rescov0.r <- getResidualCov(fit_0.r)
rescov0.r$trace
rescov0.r$var.q

fit_0.3
rescov3 <- getResidualCov(fit_0.3)
rescov3$trace
rescov3$var.q

fit_0.4
rescov4 <- getResidualCov(fit_0.4)
rescov4$trace
rescov4$var.q

pdf("phylum_diagn_nb_rich_null2.pdf")
plot(fit_0, mfrow = c(3,2),var.colors = 1) 
dev.off()

pdf("phylum_diagn_nb_rich_null2r.pdf")
plot(fit_0.r, mfrow = c(3,2),var.colors = 1) 
dev.off()

pdf("phylum_diagn_nb_rich_null3r.pdf")
plot(fit_0.3, mfrow = c(3,2),var.colors = 1) 
dev.off()

pdf("phylum_diagn_nb_rich_null4r.pdf")
plot(fit_0.4, mfrow = c(3,2),var.colors = 1) 
dev.off()


#New ugly (but correct) plots
stbPal<- c("#808000","#00FA9A")
rbPal <- c("#D2B48C", "#5F9EA0")
# "#8B4513", "#D2691E", "#006400", "#32CD32","#B266FF", "#FF66FF")
hbPal <- c("#6666FF","#FF8000","#00CC00")
sbPal2 <- c("#66FFFF","#66CC00","#CCCC00","#FF9933","#FF0000")
cbPal <- c("#B2FF66","#80FF00","#00CC00","#009900","#336600")

#Inspect colour breaks in continuous variables
jajas<-doi$Salinity[match(do$snch, doi$snch)]
binx<-as.numeric(cut(jajas, breaks = 5))
jajac<-log10(doi$Chlorophyll[match(do$snch, doi$snch)])
binx<-as.numeric(cut(jajac, breaks = 5))

dp<-do
dp$legend<-NA
dp$legend<-ifelse(dp$season==0,as.numeric(0),as.numeric(1))
dp$legendst<-ifelse(dp$substrate_type==0,as.numeric(0),as.numeric(1))
dp$pchh<-ifelse(dp$habitat=="sand",as.numeric(0),ifelse(dp$habitat=="rocks",as.numeric(1),as.numeric(2)))

pdf("ordi_phylum_nb_rich_null_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(2, 1))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_nullr_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.r, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(1, 2))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.r, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(1, 2))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.r, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(1,2))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.r, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(1,2))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null3r_LVs12_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(2, 1))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(2, 1))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.3, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.3, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null3r_LVs13_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(1, 3))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(1, 3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.3, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(1,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.3, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(1,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null3r_LVs23_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(2, 3))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.3, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(2, 3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.3, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(2,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.3, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(2,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null4r_LVs12_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(2, 1))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(2, 1))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.4, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.4, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(2,1))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null4r_LVs13_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(1, 3))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(1, 3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.4, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(1,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.4, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(1,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 

pdf("ordi_phylum_nb_rich_null4r_LVs23_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Colst <- stbPal[as.numeric(cut(dp[,"legendst"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Colst, main = "substrate type", cex.spp=0.3, pch=1, which.lvs = c(2, 3))
legend("topleft", legend = c("Sediment", "Water"), col=stbPal, pch=1, bty = "n")
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_0.4, symbols = T, s.colors = Col, main = "season", cex.spp=0.3, which.lvs = c(2, 3))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_0.4, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.3, which.lvs = c(2,3))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")
Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_0.4, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.3, which.lvs = c(2,3))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
dev.off() 


#ggplots

ad0<-fit_0$lvs
caca0<-cbind(ad0,do)
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = sh), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV2, colour = sh, shape=factor(substrate_type)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
theme_bw() +
scale_colour_manual(values = c("red3","darkgreen","blue3","indianred1","lightgreen","lightskyblue1")) +
scale_shape_manual(values = c(1,2)) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_rich_null2.pdf")

ad0<-fit_0.r$lvs
caca0<-cbind(ad0,do)
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = sh), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV2, colour = sh, shape=factor(substrate_type)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
theme_bw() +
scale_colour_manual(values = c("red3","darkgreen","blue3","indianred1","lightgreen","lightskyblue1")) +
scale_shape_manual(values = c(1,2)) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_rich_null2r.pdf")


ds<-names(sort(rescov3$var.q, decreasing=T))[1:2]
pad0<-fit_0.3$lvs
ad0<-pad0[,c(ds)]
caca0<-cbind(ad0,do)
colnames(caca0)[1:2]<-c("LV1","LV2")
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

pexpvar<-sort(rescov3$var.q, decreasing=T)[1:2]
expvar<-round(100*(pexpvar/rescov3$trace), 2)

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = sh), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV2, colour = sh, shape=factor(substrate_type)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
theme_bw() +
scale_colour_manual(values = c("red3","darkgreen","blue3","indianred1","lightgreen","lightskyblue1")) +
scale_shape_manual(values = c(1,2)) +
xlab(paste("LV1","-",expvar[1],"%",sep=" ")) + ylab(paste("LV2","-",expvar[2],"%",sep=" ")) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_rich_null3r.pdf")


ds<-names(sort(rescov3$var.q, decreasing=T))[1:2]
pad0<-fit_0.4$lvs
ad0<-pad0[,c(ds)]
caca0<-cbind(ad0,do)
colnames(caca0)[1:2]<-c("LV1","LV2")
cent2 <- aggregate(cbind(LV1,LV2)~sshc,caca0,mean)
colnames(cent2)[2:3]<-c("cLV1","cLV2")
segs <- merge(caca0, setNames(cent2, c("sshc","cLV1","cLV2")), by = "sshc", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

pexpvar<-sort(rescov3$var.q, decreasing=T)[1:2]
expvar<-round(100*(pexpvar/rescov3$trace), 2)

ggplot(data=segs, aes(x=LV1, y=LV2)) +
geom_segment(aes(xend = cLV1, yend = cLV2, colour = sh), alpha=0.3, size=0.05) + 
geom_point(aes(cLV1, cLV2, colour = sh, shape=factor(substrate_type)), size = 2.5, alpha=0.8) +
geom_text(aes(cLV1, cLV2, label=cluster), size=1.7, color="gray20") +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
theme_bw() +
scale_colour_manual(values = c("red3","darkgreen","blue3","indianred1","lightgreen","lightskyblue1")) +
scale_shape_manual(values = c(1,2)) +
xlab(paste("LV1","-",expvar[1],"%",sep=" ")) + ylab(paste("LV2","-",expvar[2],"%",sep=" ")) +
theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
coord_cartesian(xlim = c(min(caca0$LV1), max(caca0$LV1)), ylim = c(min(caca0$LV2),max(caca0$LV2)))
ggsave("ordi_phylum_nb_rich_null4r.pdf")


#Phyla correlation structure

cr0 <- getResidualCor(fit_0)

pdf("phylum_nb_corr_rich_null2.pdf")
corrplot(cr0[order.single(cr0), order.single(cr0)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()

cr0r <- getResidualCor(fit_0.r)

pdf("phylum_nb_corr_rich_null2r.pdf")
corrplot(cr0r[order.single(cr0r), order.single(cr0r)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()

cr03 <- getResidualCor(fit_0.3)

pdf("phylum_nb_corr_rich_null3.pdf")
corrplot(cr03[order.single(cr03), order.single(cr03)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()

cr04 <- getResidualCor(fit_0.4)

pdf("phylum_nb_corr_rich_null4.pdf")
corrplot(cr04[order.single(cr04), order.single(cr04)], diag = FALSE, type = "lower", method = "square", tl.cex = 0.55, tl.srt = 45, tl.col = "black")
dev.off()


anova(fit_0, fit_0.r, fit_0.3, fit_0.4)


