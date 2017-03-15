library(adehabitatLT)
library(chron)
library(raster)
library(sp)
library(rgdal)

muleys <-read.csv("DCmuleysedited.csv", header=T)
str(muleys)

#CODE FOR AN INDIVIDUAL ANIMAL
muley15 <- subset(muleys, id=="D15")
muley15[1:10,]
muley15$id <- factor(muley15$id)
str(muley15)
summary <- table(muley15$UTM_Zone,muley15$id)
summary

#Sort data to address error in code and then look at first 10 records of data to confirm
muley15 <- muley15[order(muley15$GPSFixTime),]
muley15[1:10,]

######################################################
## Example of a trajectory of type II (time recorded)
### Conversion of the date to the format POSIX
#Needs to be done to get proper digits of date into R then POSIXct
#uses library(chron)
da <- as.character(muley15$GPSFixTime)
da <- as.POSIXct(strptime(muley15$GPSFixTime,format="%Y.%m.%d %H:%M:%S"))
#Attach da to muley15
muley15$da <- da

timediff <- diff(muley15$da)
muley15 <-muley15[-1,]
muley15$timediff <-as.numeric(abs(timediff)) 
str(muley15)

#Clean up muley15 for outliers
newmuleys <-subset(muley15, muley15$X > 599000 & muley15$X < 705000 & muley15$Y > 4167000 & muley15$timediff < 14401)
muley15 <- newmuleys
str(muley15)


data.xy = muley15[c("X","Y")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
proj4string(xysp) <- CRS("+proj=utm +zone=12 +ellps=WGS84")

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)
#Creates a spatial data frame of ID
idsp<-data.frame(muley15[2])
#Creates a spatial data frame of dt
dtsp<-data.frame(muley15[24])
#Creates a spatial data frame of Burst
busp<-data.frame(muley15[23])
#Merges ID and Date into the same spatial data frame
merge<-data.frame(idsp,dtsp,busp)
#Adds ID and Date data frame with locations data frame
coordinates(merge)<-sppt
plot(merge)
str(merge)

#Now let's have a little fun with these mule deer locations and
#explore
#Load vegetation raster layer textfile clipped in ArcMap 
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
proj4string(merge) <- CRS("+proj=utm +zone=12 +ellps=WGS84")
deer.albers <-spTransform(merge, CRS=Albers.crs)

ltr.albers <- as.ltraj(coordinates(deer.albers),merge$da,id=merge$id)

veg <-raster("cropnlcd.tif")
plot(veg)
points(deer.albers)
#Code below is used to just zoom in on grid using raster layer
e <- drawExtent()
#click on top left of crop box and bottom right of crop box create zoom
newclip <- crop(veg,e)
plot(newclip)
points(deer.albers, col="red")
vegspdf <- as(newclip,"SpatialPixelsDataFrame")
plot(ltr.albers, spixdf=vegspdf)

#Let's zoom in even closer
e2 <- drawExtent()
newclip2 <- crop(newclip,e2)
plot(newclip2)
points(deer.albers)

zoomspdf <- as(newclip2,"SpatialPixelsDataFrame")
zoom.ltr <- crop(deer.albers,zoomspdf)
ltr.zoom <- as.ltraj(coordinates(zoom.ltr),zoom.ltr$da,id=zoom.ltr$id)
plot(ltr.zoom, spixdf=zoomspdf)

windows() #NOTE: a new window is needed in new version of R studio
#Line of code below plots trajectory one location at a time
trajdyn(ltr.zoom, burst = attr(ltr.zoom[[1]], "burst"),spixdf=zoomspdf) 



