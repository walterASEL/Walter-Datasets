library(sp)
library(lattice)
library(rgdal)
library(rgeos)
library(raster)

study.counties<-readOGR(dsn=".",layer="MDcounties")
str(study.counties)
class(study.counties)
proj4string(study.counties)
plot(study.counties)
study.counties@data$StateCO
#Labels each county with @plotOrder of each polygon (i.e., county)
text(coordinates(study.counties), labels=sapply(slot(study.counties, "polygons"), function(i) slot(i, "ID")), cex=0.8)

muleys <-read.csv("muleysexample.csv", header=T)
str(muleys)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.50 & muleys$Lat > 37.8 & muleys$Long < -107)
muleys <- newmuleys

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$Long, y = muleys$Lat)
crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
head(coords)
plot(coords)

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(crs))
head(deer.spdf)
class(deer.spdf)
proj4string(deer.spdf)
points(deer.spdf,col="red")

#Used to rename labels by county name otherwise plot order would be used
#because duplicate counties within each state (i.e., CO, UT)
row.names(study.counties)<-as.character(study.counties$StateCO)
str(study.counties@polygons[3], max.level=3)
names.polygons<-sapply(study.counties@polygons, function(x) slot(x,"ID")) 
#Now add labels of State and County to Map
plot(study.counties)
text(coordinates(study.counties), labels=sapply(slot(study.counties, "polygons"), function(i) slot(i, "ID")), cex=0.3)

#Now lets extract counties within the extent of the mule deer locations
int <- gIntersection(study.counties,deer.spdf)#requires rgeos library
clipped <- study.counties[int,]
MDclip <- as(clipped, "SpatialPolygons")
plot(MDclip,pch=16)
#Now add labels of State and County to Map
text(coordinates(MDclip), labels=sapply(slot(MDclip, "polygons"), function(i) slot(i, "ID")), cex=0.8)
bbox(MDclip)

HexPts <-spsample(MDclip,type="hexagonal", n=1000, offset=c(0,0))
HexPols <- HexPoints2SpatialPolygons(HexPts)
proj4string(HexPols) <- CRS(crs)
plot(HexPols, add=T)

study.zoom<-readOGR(dsn=".",layer="MDzoom")
plot(study.zoom,pch=16)
points(deer.spdf,col="red")

HexPts2 <-spsample(study.zoom,type="hexagonal", n=500, offset=c(0,0))
HexPols2 <- HexPoints2SpatialPolygons(HexPts2)
proj4string(HexPols2) <- CRS(crs)
plot(HexPols2, add=T)
#Now add labels to each hexagon for unique ID
text(coordinates(HexPols2), labels=sapply(slot(HexPols2, "polygons"), function(i) slot(i, "ID")), cex=0.3)

#Intersect the points with the polygon shapefile they occur in and add to deer locations
o = over(deer.spdf,study.counties)
new = cbind(deer.spdf@data, o)
head(o)
head(deer.spdf)
head(new)

#Used to rename labels by hexagonal grid ID only otherwise plot order with "IDxx" would be used
#and would throw an error (i.e., ID2, ID3)
row.names(HexPols2)<-as.character(HexPols2@plotOrder)
str(HexPols2@plotOrder[3], max.level=3)
names.hex<-sapply(HexPols2@polygons, function(x) slot(x,"ID"))

o2 = over(deer.spdf,HexPols2)
new2 = cbind(deer.spdf@data,o2)
head(new2)
deer.spdf@data[1:10,]
HexPols2

#Now plot with new grid IDs
plot(study.zoom,pch=16)
points(deer.spdf,col="red")
plot(HexPols2, add=T)
#Now add labels of State and County to Map
text(coordinates(HexPols2), labels=sapply(slot(HexPols2, "polygons"), function(i) slot(i, "ID")), cex=0.5)

#We can also create and clip with a rectangular grid and do the same thing
proj4string(deer.spdf)
bbox(deer.spdf@coords)
bb <- cbind(x=c(-108.83966,-108.83966,-108.9834,-108.9834, -108.83966), y=c(37.8142, 37.86562,37.86562,37.8142,37.8142))
SP <- SpatialPolygons(list(Polygons(list(Polygon(bb)),"1")), proj4string=CRS(proj4string(MDclip)))
plot(SP)
proj4string(SP)
points(deer.spdf,col="red")

#CODE below will create polygon and keep original metadata of MDclip polygon if needed
gI <- gIntersects(MDclip, SP, byid=TRUE)
out <- vector(mode="list", length=length(which(gI)))
ii <- 1
for (i in seq(along=gI)) if (gI[i]) {

out[[ii]] <- 
gIntersection(MDclip[i,], SP); row.names(out[[ii]]) <- 
row.names(MDclip)[i]; ii <- ii+1
}

out1 <- do.call("rbind", out)
plot(out1, col = "khaki", bg = "azure2", add = TRUE)

#but do remember to reset par() to defaults with:

oldpar <- par(plt = c(0.57, 0.87, 0.4, 0.7), new = TRUE)

#before the insert and:

par(oldpar)

#afterwards!

plot(out1)
points(deer.spdf,col="red")

#Load vegetation raster layer textfile clipped in ArcMap 
veg <-raster("extentnlcd2.txt")
plot(veg)
class(veg)

#Clip using the raster imported with "raster" package 
bbclip <- crop(veg, SP)
veg
#WON'T WORK because projections are not the same, WHY?
#Let's check projections of layers we are working with now.
proj4string(MDclip)
proj4string(deer.spdf)
proj4string(SP)
proj4string(veg)

#Will not work so let's project the deer.spdf to Albers 
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
deer.albers <-spTransform(deer.spdf, CRS=Albers.crs)
points(deer.albers, col="red")
class(deer.albers)
proj4string(deer.albers)
head(deer.spdf)
head(deer.albers)

bbox(deer.albers)
bb1 <- cbind(x=c(-1115562,-1115562,-1127964,-1127964,-1115562), y=c(1718097, 1724867,1724867,1718097,1718097))
AlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(deer.albers)))
plot(AlbersSP)

#Check to see all our layers are now in Albers projection
proj4string(veg)
proj4string(deer.albers)
proj4string(AlbersSP)

plot(veg)
points(deer.albers, col="red")

#Clip using the raster imported with "raster" package 
bbclip <- crop(veg, AlbersSP)
plot(bbclip)
points(deer.albers, col="red")
plot(AlbersSP, lwd=5, add=T)
text(coordinates(AlbersSP), labels="Colorado Mule Deer")

