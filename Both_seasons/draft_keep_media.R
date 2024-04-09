yul<-DADAwang1
sample_data(yul)$over_median<-p_norm_sample_data$q[match(sample_data(yul)$sample_ID, p_norm_sample_data$sample_ID)]

above_t<-rarefy_even_depth(subset_samples(yul, over_median==TRUE), sample.size=as.numeric(round(median(combinedi$readsi))), replace=FALSE, rngseed= 13072021)

below_t<-subset_samples(yul, over_median==FALSE)

#merge objects
tudao_r_new<-merge_phyloseq(above_t, below_t)
tudao_pos_merge_new<-merge_samples(tudao_r_new, "sample_root")

#Re-build sample_data
tudaok = filter_taxa(tudao_pos_merge_new, function(x) sum(x) > 0, TRUE)
p_ncbi = prune_samples(sample_sums(tudaok)>0,tudaok)

otug<-otu_table(p_ncbi)
otugi<-t(data.frame(otug, check.names=F))
rich<-colSums(otugi != 0)

#Calculate beta_diversity (as PERMANOVA results)
dist_n<-vegdist(data.frame(otug), method="bray")
d$s_st_h_c<-paste(d$cluster,d$season,d$substrate_type,d$habitat, sep="_")
zin<-betadisper(dist_n, as.factor(d$s_st_h_c))
tre<-with(zin, tapply(distances, as.factor(d$s_st_h_c), "mean"))

#Write info
beta_br[6,"Quantile"]<-"keep_all"
beta_br[6,"Total_reads"]<-sum(sample_sums(p_ncbi))
beta_br[6,"Total_ASVs"]<-ntaxa(p_ncbi)
beta_br[6,"Field_replicates"]<-nsamples(p_ncbi)
beta_br[6,"Depth"]<-ifelse(length(levels(as.factor(sample_sums(p_ncbi))))==1, levels(as.factor(sample_sums(p_ncbi))), mean(sample_sums(p_ncbi)))
beta_br[6,"Mean_richness"]<-round(mean(rich), digits=0)
beta_br[6,"SD_richness"]<-round(sd(rich), digits=0)
beta_br[6,"Total_replicates_discarded"]<-length(levels(factor(combinedi$sample_ID)))-nsamples(tudao_r_new)
beta_br[6,"Mean_bray_centroid_field_rep"]<-mean(tre)
beta_br[6,"SD_bray_centroid_field_rep"]<-sd(tre)

d<-data.frame(sample_data(p_ncbi)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])
d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$sample_root<-rownames(d)

d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)

d$cluster<-as.integer(d$pn2)

sample_data(p_ncbi)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

p_ncbi.2 <- transform_sample_counts(p_ncbi, function(x) sqrt(x))
ordb<-ordinate(p_ncbi.2, "PCoA", "bray")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/depth_cutoff")

plot<-plot_ordination(p_ncbi.2, ordb, type = "sample", color="habitat", shape="season") + geom_text(aes(label=cluster), size=1.5, color="black", alpha=0.9) + geom_point(size=1.6, alpha=0.2) + stat_ellipse(aes(lty=season), type = "norm") + facet_wrap(~substrate_type, ncol=1) + theme_bw()
ggsave("keep_below_ord_bc_root_ncbi.pdf")


