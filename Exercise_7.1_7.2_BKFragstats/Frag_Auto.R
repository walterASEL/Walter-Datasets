library(SDMTools)
library(raster)
library(plyr)
library(maptools)
library(rgdal)

#Load vegetation raster layer clipped in ArcMap 
crops <-raster("crop2012utm12.tif")
plot(crops)
class(crops)
as.matrix(table(values(crops)))
proj4string(crops)
crops

# reclassify the values into 9 groups
# all values between 0 and 20 equal 1, etc.
m <- c(-Inf,0,NA,2, 7, 2, 20, 60, 3, 60, 70, 4, 110, 132, 5, 133, 150, 6, 151, 191, 7, 192,Inf,NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(crops, rclmat)
plot(rc)
rc
as.matrix(table(values(rc)))

#Now we get into Landscape Metrics with the SDTM Tool
#Calculate the Patch statistics
ps.data = PatchStat(rc)
ps.data

#Calculate the Class statistics
cl.data = ClassStat(rc)
cl.data

############################################################################
############################################################################
#Some research designs may need landscape metrics for several areas that may 
#be available as a shapefile or some other polygon layer.  

rm(list=ls())

library(SDMTools)
library(raster)
library(rgdal)
library(maptools)
library(plyr)

### load raster file into R
raster <- raster("county_hab")
raster

### load PA shapefile into R
HareCounties <- readOGR(dsn=".", layer="Hare_Counties")
HareCounties
plot(HareCounties)
proj4string(HareCounties)
proj4string(raster)
image(raster)
plot(HareCounties, add=T)

#Let's project Counties to habitat just to be safe
new.crs <-CRS("+proj=lcc +lat_1=41.95 +lat_2=40.88333333333333 +lat_0=40.16666666666666 +lon_0=-77.75 +x_0=600000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
county <- spTransform(HareCounties, CRS=new.crs)
HareCounties <- county
proj4string(HareCounties)
proj4string(raster)
#Matching projections successful!
row.names(HareCounties)<-as.character(HareCounties$COUNTY_NAM)
names.polygons<-sapply(HareCounties@polygons, function(x) slot(x,"ID")) 
text(coordinates(HareCounties), labels=sapply(slot(HareCounties, "polygons"), function(i) slot(i, "ID")), cex=0.8)

#Now we want to export by County name (i.e., COUNTY_NAM) as individual shapefiles
indata <- HareCounties
innames <- unique(HareCounties@data$COUNTY_NAM)
innames <- innames[1:2]
outnames <- innames

# set up output table
#output <- as.data.frame(matrix(0,nrow=length(innames),ncol=38))

# begin loop to create separate county shapefiles 
for (i in 1:length(innames)){
  data <- indata[which(indata$COUNTY_NAM==innames[i]),]
  if(dim(data)[1] != 0){
    writePolyShape(data,fn=paste(outnames[i],sep="/"),factor2char=TRUE)
    write.table(innames, "List.txt", eol=".shp\n", col.names=FALSE, quote=FALSE, row.names=FALSE)
}
}

#Read in a list of shapefiles files from above
Listshps<-read.table("List.txt",sep="\t",header=F)
#colnames(Listshps) <- c("id")
Listshps

shape <- function(Listshps) {
file <- as.character(Listshps[1,])
shp <- readShapeSpatial(file)
mask <- mask(raster,shp)
### Calculate the Class statistics in each county
cl.data <- ClassStat(mask)
}

results <- ddply(Listshps, 1, shape)
results
write.table(results, "FragCounty.txt")







