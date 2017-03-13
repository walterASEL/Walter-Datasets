library(rgdal)

study.states<-readOGR(dsn=".",layer="MDcounties")

plot(study.states, col="grey")

study.zoom<-readOGR(dsn=".",layer="MDzoom")

plot(study.zoom, col="grey")

muleys <-read.csv("muleysexample.csv", header=T)
str(muleys)

#Raw locations
coords<-data.frame(x = muleys$Long, y = muleys$Lat)
crs<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
coords
plot(coords)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.80)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys

#now remove NAs from excluding locations
muleys <- subset(muleys, !is.na(muleys$Lat))

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$Long, y = muleys$Lat)
coords
plot(coords)

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(crs))
deer.spdf[1:5,]
class(deer.spdf)
proj4string(deer.spdf)
points(deer.spdf)
points(deer.spdf, col="yellow")

# Now let's project both the mule deer locations and study site shapefile to NAD83 UTM Zone 12
new.crs <-CRS("+proj=utm +zone=12 +datum=WGS84")
MDzoomUTM12 <-spTransform(study.zoom, CRS=new.crs)
par(new=TRUE)
plot(MDzoomUTM12, col="bisque")
class(MDzoomUTM12)
proj4string(MDzoomUTM12)
summary(MDzoomUTM12)

deer.crs <-CRS("+proj=utm +zone=12 +datum=WGS84")
deerUTM12 <-spTransform(deer.spdf, CRS=deer.crs)
points(deerUTM12, col="red")
class(deerUTM12)
proj4string(deerUTM12)
deerUTM12[1:5,]

coordinates(deerUTM12)[1:5,]

       x       y
[1,] 677825.2 4192832
[2,] 677853.8 4192787
[3,] 677736.3 4192728
[4,] 677595.9 4192398
[5,] 677666.2 4192362

#plot coordinates in Lat Long over coordinates in UTM 12N
plot(coords)
par(new=TRUE)
plot(deerUTM12, col="red")

windows()
plot(study.zoom)
par(new=TRUE)
plot(deer.spdf, col="red")

windows()
plot(MDzoomUTM12,col="bisque")
par(new=TRUE)
plot(deerUTM12, col="red")

#another method to project data
require(rgdal)
EPSG <- make_EPSG()

nad83 <- EPSG[grep("NAD83",EPSG$note),]
nad83[grep("UTM zone 12N", nad83$note),]


