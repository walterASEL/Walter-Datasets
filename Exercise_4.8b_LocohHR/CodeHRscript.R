library(adehabitatHR)
library(shapefiles)
library(rgeos)
library(rgdal)
library(maptools)

#Get input file
panther <- read.csv("pantherjitter2.csv")
str(panther)
panther$CatID <- as.factor(panther$CatID)

#Or explore with one panther with 381 relocations
cat159 <- subset(panther, CatID=="159")
str(cat159)
cat159$CatID <- factor(cat159$CatID)

#Get the relocation data from the source file
data.xy = cat159[c("x","y")]

xysp <- SpatialPoints(data.xy)

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)
#Creates a spatial data frame of ID
idsp<-data.frame(cat159[1])
#Adds ID and Date data frame with locations data frame
coordinates(idsp)<-sppt
proj4string(idsp) <- CRS("+proj=utm +zone=17 +ellps=WGS84")
locsdf <-as.data.frame(idsp)
head(locsdf)
## Shows the relocations
plot(data.xy, col="red")
locsdf[1:5,]
table(panther$CatID)

## Examinates the changes in home-range size for various values of k
## Be patient! the algorithm can be very long
ar <- LoCoH.k.area(idsp, k=c(15:25))
## 24 points seems to be a good choice (rough asymptote for all animals)
## the k-LoCoH method:
nn <- LoCoH.k(idsp[,1], k=19)
## Graphical display of the results
plot(nn, border=NA)
## the object nn is a list of objects of class
## SpatialPolygonsDataFrame
length(nn)
names(nn)
class(nn[[1]])
## shows the content of the object for the first animal
as.data.frame(nn[[1]])

## The 95% home range is the smallest area for which the
## proportion of relocations included is larger or equal
## to 95% In this case, it is the 339th row of the
## SpatialPolygonsDataFrame.
plot(nn[[1]][339,],lwd=2)

#The 50% home range code is on line 146
plot(nn[[1]][146,],add=TRUE)

#The 99% home range code is on line 359
plot(nn[[1]][359,],lwd=3, add=TRUE)

#Save shapefiles of resulting home range
ver <- getverticeshr(nn)
ver
plot(ver)
writeOGR(ver,dsn="FixedK",layer="FixedK24", driver = "ESRI Shapefile",
     overwrite=TRUE)
##Overwrite will not work so must edit path so "FixedK" folder is created with code below.
ver50 <-getverticeshr(nn, percent=50)
writeOGR(ver50,dsn="FixedK",layer="50FixedK24", driver = "ESRI Shapefile",overwrite=TRUE)
ver95 <-getverticeshr(nn, percent=95)
writeOGR(ver95,dsn="FixedK",layer="95FixedK24", driver = "ESRI Shapefile",overwrite=TRUE)
ver99 <-getverticeshr(nn, percent=99)
writeOGR(ver99,dsn="FixedK",layer="99FixedK24", driver = "ESRI Shapefile",overwrite=TRUE)
