#### Calculating waterway distances ### 

# Based on this tutorial: #
# https://gis.stackexchange.com/questions/305070/finding-shortest-path-without-barriers-using-the-sf-package-in-r

# set wd
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/dist_sites")

# Add libraries
library(sp)
library(raster)
library(dplyr)
library(rgeos)
library(rgdal)
library(gdistance)
library(SDraw)
library(ncdf4)
library(Matrix)


###### Load site locations
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

coord<-subset(merged_wat, sn=="spring"&hb=="sand")[,c("Sample_ID","Longitude","Latitude","cl")]
coord<-coord[order(coord$cl),]
rownames(coord) <- NULL


#### Step 1:

# Choose some projections
# You will need to choose a metric one so you can calculate distances later...
# PROJECTION IN is the standard flat map projection
# PROJECTION_OUT is a metric one (gives you results in meters/kms)

PROJECTION_IN = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#PROJECTION_OUT = "+proj=lcc +lat_1=-54 +lat_2=-54.75 +lat_0=-55 +lon_0=-37 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
#PROJECTION_OUT = "+init=EPSG:3035"
PROJECTION_OUT = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

#Localities
#Agernæs Havn	9.986035242	55.19995896
#Allinge	14.80350748	55.27889500
#Anholt havn	11.50908377	56.71370855
#Assens Havn	9.884813328	55.27208598
#Bagenkop havn	10.66949082	54.75283305
#...

#### Step 2: Open landmask & create distance layer ######
# set wd
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/dist_sites")

# Landmask is from here: https://sdi.eea.europa.eu/catalogue/EEA_Reference_Catalogue/api/records/228e9bf2-0b81-4506-ab72-d42dd2ad19b6

p_landmask<- raster('eea_r_3035_100_m_clc12_V18_5_land_mask.tif')

#bornholm
#raster.ext<-extent(as.data.frame(point1_trans)$Lon - 350000, as.data.frame(point1_trans)$Lon + 32000, as.data.frame(point1_trans)$Lat - 95000, as.data.frame(point1_trans)$Lat + 225000)

#raster.ext<-extent(as.data.frame(point1_trans)$Lon - 100000, as.data.frame(point1_trans)$Lon + 100000, as.data.frame(point1_trans)$Lat - 100000, as.data.frame(point1_trans)$Lat + 100000)

##
# Define your points and project them
br_dist<-data.frame(coord)
br_dist$dist<-NA
for (i in 18:nrow(coord))
{

point1<-data.frame(coord[6,2:3]) ## Location 1 input GPS
colnames(point1)<-c("Lon", "Lat")
coordinates(point1)<-~Lon + Lat
proj4string(point1)<-PROJECTION_IN
point1_trans<-spTransform(point1, PROJECTION_OUT)

#Problem with cluster16
if(coord[i,"cl"]==16)
{point2<-data.frame(t(c(10.56285,55.52186)))} else
{point2<-data.frame(coord[i,2:3])}

fig_name<-coord[i,"cl"]
colnames(point2)<-c("Lon", "Lat")
coordinates(point2)<-~Lon + Lat
proj4string(point2)<-PROJECTION_IN
point2_trans<-spTransform(point2, PROJECTION_OUT)

#Making fixed limits for lat and long omitted paths
landmask<-p_landmask
maxlat<-ifelse(coord[i,"cl"]>11,3860000,ifelse(point1_trans$Lat>point2_trans$Lat,as.data.frame(point1_trans)$Lat + 50000,as.data.frame(point2_trans)$Lat + 50000))

#for cluster < 5
#maxlong<-ifelse(coord[i,"cl"]==15|coord[i,"cl"]==17|coord[i,"cl"]==18,4400000,ifelse(point1_trans$Lat>point2_trans$Lon,as.data.frame(point1_trans)$Lon + 50000,as.data.frame(point2_trans)$Lon + 50000))

#for cluster 5:7
maxlong<-ifelse(coord[i,"cl"]<20,4400000,ifelse(point1_trans$Lat>point2_trans$Lon,as.data.frame(point1_trans)$Lon + 50000,as.data.frame(point2_trans)$Lon + 50000))

minlong<-4200000


#defining limits
raster.ext<-extent(ifelse(point1_trans$Lon<point2_trans$Lon,as.data.frame(point1_trans)$Lon - 50000,as.data.frame(point2_trans)$Lon - 50000),maxlong,ifelse(point1_trans$Lat<point2_trans$Lat,as.data.frame(point1_trans)$Lat - 50000,as.data.frame(point2_trans)$Lat - 50000),maxlat)

land.crop<-crop(landmask, raster.ext) # Crop straight away as this file is massive
plot(land.crop)
land.crop[land.crop == 1] <- NA # Change land to NA so that it's impossible to cross
land.crop[land.crop == 0] <- 1 # Now sea is 1 and land is NA
land.new<-aggregate(land.crop, fact=5) # change resolution otherwise impossible to do anything

#Check whether points are within geo zone
land.new
coordinates(point1_trans)
coordinates(point2_trans)

remove(landmask)
r <- transition(land.new, mean, directions = 16) ## Has to be 4, 8 or 16 - the higher the better
tr <- geoCorrection(r, "c") # this corrects some kind of map distortion (see: https://www.rdocumentation.org/packages/gdistance/versions/1.3-6/topics/geoCorrection)
   
distance <- gdistance::shortestPath(tr, coordinates(point1_trans), coordinates(point2_trans), output ="SpatialLines")

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Gua18s/both_seasons/results/maps/dist_sites/paths_fig")
pdf(paste("dist_7","_",fig_name,".pdf",sep=""))
plot(land.new)
plot(distance, add=T)
plot(point1_trans, add=T)
plot(point2_trans, add=T)
dev.off()

# Measure length of path in km
dist_km<-as.data.frame(lineLength(distance, byid=FALSE)/1000) ### 
colnames(dist_km)<-c("dist_km")
dist_km
br_dist[i,"dist"]<-dist_km
}

br_dist
