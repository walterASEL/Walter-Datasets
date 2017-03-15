library(sp)
library(lattice)
library(rgdal)
library(rgeos)
library(raster)
library(spatstat)
library(polyCub)#replaces gpclib function
library(pgirmess)

muleys <-read.csv("muleysexample.csv", header=T)
summary(muleys$id)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.80)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$X, y = muleys$Y)
crs<-"+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(crs))
deer.spdf[1:5,]
class(deer.spdf)
proj4string(deer.spdf)

roads<-readOGR(dsn=".",layer="AlbersRoads")
rivers<-readOGR(dsn=".",layer="AlbersRivers")
plot(roads,pch=16)
points(deer.spdf, col="red")
plot(rivers,add=T, col="blue",pch=16)

#Will not work so let's project the deer.spdf to Albers 
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
deer.albers <-spTransform(deer.spdf, CRS=Albers.crs)
points(deer.albers, col="red")
class(deer.albers)
proj4string(deer.albers)
deer.spdf[1:5,]
deer.albers[1:5,]

bbox(deer.albers)
bb1 <- cbind(x=c(-1106865,-1106865,-1145027,-1145027, -1106865), y=c(1695607, 1729463,1729463,1695607,1695607))
AlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(deer.albers)))
plot(AlbersSP)
#Load vegetation raster layer textfile clipped in ArcMap 
veg <-raster("cropnlcd.tif")

#Check to see all our layers are now in Albers projection
proj4string(veg)
proj4string(deer.albers)
proj4string(AlbersSP)

plot(veg)
points(deer.albers, col="red")

#Clip using the raster imported with "raster" package 
bbclip <- crop(veg, AlbersSP) 
cliproads <- gIntersection(roads, AlbersSP, byid=TRUE)
cliprivers <- gIntersection(rivers, AlbersSP, byid=TRUE)

plot(bbclip)
points(deer.albers, col="red")
plot(cliproads, add=T)
plot(cliprivers, col="blue",add=T)

###########################################################
#Code below will be for use with the spatstat package to
#convert segments of line layers (e.g., roads, rivers)
#to lines to enable distance to feature from deer locations
############################################################

#Most calculations with spatstat require 3 new classes so most
#code is created to achieve this goal:
#"owin" Observation windows
#"ppp" Planar point patterns
#"owin" Planar segment patterns

#Let's start with the road layer by converting a single line to a set of segments
#packaged as a function:
foo <- function(cliproads){
x <- cliproads@Lines[[1]]@coords
cbind(
head(x,-1),
tail(x,-1))}
#The function can be applied successively to each line in the list we
#extracted from roads. Results are output as a list, then converted to a
#matrix.
segs.lst <- lapply(cliproads@lines,foo)
segs <- do.call(rbind,segs.lst)

segs.x <- c(segs[,c(1,3)])
segs.y <- c(segs[,c(2,4)])
segs.owin <- as.owin(c(range(segs.x),range(segs.y)))

#The segments as a planar segment pattern:
segs.psp <- as.psp(segs, window=segs.owin)
plot(segs.psp)
points(deer.albers)
segs.psp[1:5]
lengths.psp(segs.psp[1:10])

#We can cut road segments into distances we control
dist <- pointsOnLines(segs.psp, eps=1000)

#Or using the pgirmess package (but some unknown errors)
#sppoints <- transLines2pix(clip, mindist=0.2)
#points <- thintrack(sppoints,mindist=100)
#plot(thintrack(mySPDF),pch=19,cex=0.7,col="red",add=TRUE)

#Convert deer.albers from SPDF back to a dataframe because we need xy coordinates to be in Albers.
#If all data is in UTM 12N then no need for this step.
deer2 <-as.data.frame(deer.albers)
deer2[1:5,]
newdeer <-cbind(deer2$x,deer2$y)
newdeer[1:5,]

#Create Spatial Polygons Data Frame of bounding polygon to set window for package spatstat
AlbersSPDF <- as(AlbersSP, "SpatialPolygonsDataFrame")
#bdy.owin <- gpc2owin(bdy.gpc)#new code using polyCub package
bdy.owin <- as.owin(AlbersSPDF)#NOTE:If 2 lines of new code does not work use this line with Spatial Polygons Data Frame
bdy <- as.polygonal(bdy.owin)
xy.ppp <- as.ppp(newdeer,W=bdy)
plot(xy.ppp)

#Let's buffer around the bounding box to be sure it encompasses all locations
buffSP <- gBuffer(AlbersSP,width=1000)
plot(buffSP)
points(deer.albers,col="red")
AlbersSPDF <- buffSP #Now re-run from bdy.gpc code from above to plot(xy.ppp)

is.owin(bdy.owin)
#[1] TRUE
is.ppp(xy.ppp)
#[1] TRUE
is.psp(segs.psp)
#[1] TRUE

#All is TRUE so now we can move forward with the analysis
roaddist <- nncross(xy.ppp, segs.psp)$dist
roaddist[1:5]
#Or identify segment number closest to each point
v <- nearestsegment(xy.ppp,segs.psp)#Identifies segment number not a distance
plot(segs.psp)
plot(xy.ppp[101], add=TRUE, col="red")
plot(segs.psp[v[101]], add=TRUE, lwd=5, col="red")

#Now we do the same to a river layer by converting a single line to a set of segments
#packaged as a function:
foo <- function(cliprivers){
x <- cliprivers@Lines[[1]]@coords
cbind(
head(x,-1),
tail(x,-1))}
#The function can be applied successively to each line in the list we
#extracted from roads. Results are output as a list, then converted to a
#matrix.
rivs.lst <- lapply(cliprivers@lines,foo)
rivs <- do.call(rbind,rivs.lst)

rivs.x <- c(rivs[,c(1,3)])
rivs.y <- c(rivs[,c(2,4)])
rivs.owin <- as.owin(c(range(rivs.x),range(rivs.y)))

#The segments as a planar segment pattern:
rivs.psp <- as.psp(rivs, window=rivs.owin)
plot(rivs.psp)
points(deer.albers)

#All is TRUE so now we can move forward with the analysis
rivdist <- nncross(xy.ppp, rivs.psp)$dist
rivdist[1:5]

#Or identify segment number closest to each point
riv <- nearestsegment(xy.ppp,rivs.psp)
plot(rivs.psp, lwd=1)
plot(xy.ppp[1], add=TRUE, col="red")
plot(rivs.psp[riv[1]], add=TRUE, lwd=5, col="red")
plot(xy.ppp[290], add=TRUE, col="blue")
plot(rivs.psp[riv[290]], add=TRUE, lwd=5, col="blue")
points(deer.albers)

#We can now summarize the distances in some meaning way
br <- seq(0,4000,500)
lbl <- paste(head(br,-1),tail(br,-1),sep="-")
road.tbl <- table(cut(roaddist,breaks=br,labels=lbl))
Rdresults <- road.tbl/sum(road.tbl)
Rdresults

br1 <- seq(0,4000,500)
lbl1 <- paste(head(br1,-1),tail(br1,-1),sep="-")
river.tbl <- table(cut(rivdist,breaks=br1,labels=lbl1))
Rivresults <- river.tbl/sum(river.tbl)
Rivresults

#Or we can place each distance into a category or Bin for each deer.
BinRoad <- bin.var(roaddist, bins=5, method='intervals', labels=c('1','2','3','4','5'))
BinRoad <- cut(roaddist, 5, method='intervals', include.lowest=TRUE, labels=c('1','2','3','4','5'))
BinRoad

BinRivers <- bin.var(rivdist, bins=5, method='intervals', labels=c('1','2','3','4','5'))
BinRivers <- cut(rivdist, 5, method='intervals', include.lowest=TRUE, labels=c('1','2','3','4','5'))
BinRivers 
deer.albers[1:5,]
Dist <- cbind(BinRoad,BinRivers)#end re-run here##########

muleys <- cbind(muleys, Dist)#This won't work if any locations are out of the box
str(muleys)

