library(sp)
library(lattice)
library(rgdal)
library(rgeos)
library(raster)

muleys <-read.csv("muleysexample.csv", header=T)
summary(muleys$id)

muley8 <- subset(muleys, id=="D8")
str(muley8)
summary <- table(muley8$UTM_Zone,muley8$id)
summary(muley8$id)
muley8$id <- factor(muley8$id)

#Remove outlier locations if needed
summary(muley8$Long)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -111.8  -108.9  -108.9  -108.9  -108.9  -108.8 
#NOTE: Min. of -111.8 is an outlier so remove
summary(muley8$Lat)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  33.38   37.84   37.84   37.83   37.85   37.86 
#NOTE: Min. of 33.38 is an outlier so remove
newmuley8 <-subset(muley8, muley8$Long > -111.7 & muley8$Lat > 37.80)
str(newmuley8)
muley8 <- newmuley8

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muley8$Long, y = muley8$Lat)
crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
head(coords)

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muley8, proj4string = CRS(crs))
head(deer.spdf)
class(deer.spdf)
proj4string(deer.spdf)
plot(deer.spdf, axes=T)

#Again let's project the deer.spdf to Albers 
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
deer.albers <-spTransform(deer.spdf, CRS=Albers.crs)
class(deer.albers)
proj4string(deer.albers)
head(deer.spdf)
head(deer.albers)

bbox(deer.albers)
bb1 <- cbind(x=c(-1115562,-1115562,-1120488,-1120488, -1115562), y=c(1718097,1722611,1722611,1718097,1718097))
AlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(deer.albers)))
plot(AlbersSP)
points(deer.albers, col="red")

#Load vegetation raster layer textfile clipped in ArcMap 
veg <-raster("cropnlcd.tif")
plot(veg)
plot(AlbersSP,add=T)
points(deer.albers, col="red")

#Clip using the raster imported with "raster" package 
bbclip <- crop(veg, AlbersSP)
plot(bbclip)
plot(AlbersSP,add=T)
points(deer.albers, col="red")

#Let's expand the polygon created in R to encompass more of an area around locations
bb1 <- cbind(x=c(-1115000,-1115000,-1121000,-1121000, -1115000), y=c(1717000,1723000,1723000,1717000,1717000))
AlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(deer.albers)))

bbclip <- crop(veg, AlbersSP)
plot(bbclip)
plot(AlbersSP,lwd=2, add=T)
points(deer.albers, col="red")

extract(bbclip,deer.albers)
settbuff=gBuffer(deer.albers,width=100)
plot(bbclip)
plot(settbuff, add=T, lty=2)
table(extract(bbclip,settbuff))

#Cell size of raster layer
res(bbclip)

 30^2
#[1] 900
900*37
#[1] 33300
(900*37)/1000000
#[1] 0.0333 km2

#Let's create the buffers a different way so we don't
#have one continuous polygon
settbuff=gBuffer(deer.albers, width=100, byid=TRUE)
plot(bbclip)
points(deer.albers, col="blue")
plot(settbuff, add=TRUE,lty=8)

e= extract(bbclip,settbuff)
et=lapply(e,table)

et[[328]]

