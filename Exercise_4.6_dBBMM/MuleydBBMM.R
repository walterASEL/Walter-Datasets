library(adehabitatLT)
library(move)
library(circular)
library(sp)
library(maptools)
library(adehabitatHR)
library(PBSmapping)

muleys <-read.csv("muleysexample.csv")
str(muleys)

#TIME DIFF ONLY NECESSARY AS A MEANS EXCLUDE POOR DATA LATER
muleys$Date <- as.numeric(muleys$GPSFixTime)
timediff <- diff(muleys$Date)*24*60
muleys <-muleys[-1,]
muleys$timediff <-as.numeric(abs(timediff))

#FOR ALL DEER
muleys$DT <-as.POSIXct(strptime(muleys$GPSFixTime, format='%Y.%m.%d %H:%M:%OS'))
muleys$DT

#Sort data to address error in code
muleys <- muleys[order(muleys$id,muleys$DT),]
summary(muleys$id)

#EXCLUDE OUTLIERS AND POOR DATA FIXES
newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.80)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys

#Create a move object for all deer using the Move package
loc <- move(x=muleys$X, y=muleys$Y, time=as.POSIXct(muleys$GPSFixTime, 
     format="%Y.%m.%d %H:%M:%S"), proj=CRS("+proj=utm +zone=12 +datum=NAD83"),data=muleys, 
     animal=muleys$id)

#Now create a dBBMM object
two_dbbmm <- brownian.bridge.dyn(object=loc, location.error=22, window.size=19, margin=7,  dimSize=100,time.step=180)
plot(two_dbbmm)

#or select a single deer
dataD8 <- subset(muleys, muleys$id == "D8")
dataD8$id <- droplevels.factor(dataD8$id)
#dataD8 <- muleys
d8 <- move(x=dataD8$X, y=dataD8$Y, time=as.POSIXct(dataD8$GPSFixTime, format="%Y.%m.%d %H:%M:%S"), proj=CRS("+proj=utm +zone=12 +datum=NAD83"),data=dataD8, animal=dataD8$id)
d8_dbbmm <- brownian.bridge.dyn(object=d8, location.error=22, window.size=19, margin=7,  dimSize=10,time.step=180)
plot(d8_dbbmm)
contour(d8_dbbmm, levels=c(.5,.9,.95,.99))
show(d8_dbbmm)

#Plot the movement of the animal
plot(d8, type="o", col=3, lwd=2, pch=20, xlab="location_east",ylab="location_north")

#Code below will get area of each isopleth in R
d8_cont <- getVolumeUD(d8_dbbmm)
d8_cont50 <- d8_cont<=.50
d8_cont95 <- d8_cont<=.95
area50 <- sum(values(d8_cont50))
area50
area95 <- sum(values(d8_cont95))
area95

##Cast the data over to an adehabitatHR estUD
dbbmm.px <- as(d8_dbbmm, "SpatialPixelsDataFrame")
image(dbbmm.px)
dbbmm.ud <- new("estUD",dbbmm.px)
dbbmm.ud@vol = FALSE
dbbmm.ud@h$meth = "dBBMM"
##Convert the raw UD values to volume
#udvol <- getvolumeUD(dbbmm.ud, standardize=FALSE)

shp50 <- getverticeshr(dbbmm.ud, percent=50, standardize=TRUE)
class(shp50)#Now is a SpatialPolygonsDataFrame
map.ps50 <- SpatialPolygons2PolySet(shp50)
#diss.map.50 <- joinPolys(map.ps50, operation = 'UNION')
diss.map.50 <- as.PolySet(map.ps50, projection = 'UTM', zone = '12')
diss.map.p50 <- PolySet2SpatialPolygons(diss.map.50, close_polys = TRUE)
data50 <- data.frame(PID = 1)
diss.map.p50 <- SpatialPolygonsDataFrame(diss.map.p50, data = data50)
writeOGR(diss.map.p50, dsn = ".", layer="d8contour50", driver = "ESRI Shapefile")
d8map.50 <- readOGR(dsn=".", layer="d8contour50")
plot(d8map.50)

shp95 <- getverticeshr(dbbmm.ud, percent=95, standardize=TRUE)
class(shp95)#Now is a SpatialPolygonsDataFrame
plot(shp95, add=TRUE)
map.ps95 <- SpatialPolygons2PolySet(shp95)
#diss.map.50 <- joinPolys(map.ps50, operation = 'UNION')
diss.map.95 <- as.PolySet(map.ps95, projection = 'UTM', zone = '12')
diss.map.p95 <- PolySet2SpatialPolygons(diss.map.95, close_polys = TRUE)
data95 <- data.frame(PID = 1)
diss.map.p95 <- SpatialPolygonsDataFrame(diss.map.p95, data = data95)
writeOGR(diss.map.p95, dsn = ".", layer="d8contour95", driver = "ESRI Shapefile")
d8map.95 <- readOGR(dsn=".", layer="d8contour95")
plot(d8map.95, add=T)

shp99 <- getverticeshr(dbbmm.ud, percent=99, standardize=TRUE)
class(shp99)#Now is a SpatialPolygonsDataFrame
plot(shp99, add=TRUE)
map.ps99 <- SpatialPolygons2PolySet(shp99)
#diss.map.50 <- joinPolys(map.ps50, operation = 'UNION')
diss.map.99 <- as.PolySet(map.ps99, projection = 'UTM', zone = '12')
diss.map.p99 <- PolySet2SpatialPolygons(diss.map.99, close_polys = TRUE)
data99 <- data.frame(PID = 1)
diss.map.p99 <- SpatialPolygonsDataFrame(diss.map.p99, data = data99)
writeOGR(diss.map.p99, dsn = ".", layer="d8contour99", driver = "ESRI Shapefile")
d8map.99 <- readOGR(dsn=".", layer="d8contour99")
plot(d8map.99)
