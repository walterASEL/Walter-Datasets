library(adehabitatLT)
library(chron)
library(sp)
library(rgdal)

muleys <-read.csv("DCmuleysedited.csv", header=T)
str(muleys)

#Code to look at number of relocations per animal
table(muleys$id)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.50 & muleys$Lat > 37.3 & muleys$Long < -107)
muleys <- newmuleys

######################################################
## Example of a trajectory of type II (time recorded)
### Conversion of the date to the format POSIX
#Needs to be done to get proper digits of date into R then POSIXct
da <- as.character(muleys$GPSFixTime)
da <- as.POSIXct(strptime(muleys$GPSFixTime,format="%Y.%m.%d %H:%M:%S"))
muleys$da <- da

timediff <- diff(muleys$da)*60*60
muleys <-muleys[-1,]
muleys$timediff <-as.numeric(abs(timediff)) 
str(muleys)

data.xy = muleys[c("X","Y")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
#proj4string(xysp) <- CRS("+proj=utm +zone=17 +ellps=WGS84")

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)
#Creates a spatial data frame of ID
idsp<-data.frame(muleys[2])
#Creates a spatial data frame of dt
dtsp<-data.frame(muleys[24])
#Creates a spatial data frame of Burst
busp<-data.frame(muleys[23])
#Merges ID and Date into the same spatial data frame
merge<-data.frame(idsp,dtsp,busp)
#Adds ID and Date data frame with locations data frame
coordinates(merge)<-sppt
plot(merge)
str(merge)

### Creation of an object of class "ltraj", with for ### example the first animal
ltraj <- as.ltraj(coordinates(merge),merge$da,id=merge$id)
plot(ltraj)
ltraj

#First Passage Time

plot(ltraj[1])
i1 <- fpt(ltraj[1], seq(300,1000, length=30))
plot(i1, scale = 200, warn = FALSE)

plot(ltraj[2])
i2 <- fpt(ltraj[2], seq(300,1000, length=30))
plot(i2, scale = 500, warn = FALSE)

toto2 <- meanfpt(i2)
toto2
attr(toto2, "radii")

toto2 <- varlogfpt(i2)
toto2
attr(toto2, "radii")

plot(ltraj[3])
i3 <- fpt(ltraj[3], seq(300,1000, length=30))
plot(i3, scale = 500, warn = FALSE)

toto3 <- meanfpt(i3)
toto3
attr(toto3, "radii")

toto3 <- varlogfpt(i3)
toto3
attr(toto3, "radii")

plot(ltraj[4])
i4 <- fpt(ltraj[4], seq(300,1000, length=30))
plot(i4, scale = 500, warn = FALSE)

toto4 <- meanfpt(i4)
toto4
attr(toto4, "radii")

toto4 <- varlogfpt(i4)
toto4
attr(toto4, "radii")

plot(ltraj[5])
i5 <- fpt(ltraj[5], seq(300,1000, length=30))
plot(i5, scale = 500, warn = FALSE)

toto5 <- meanfpt(i5)
toto5
attr(toto5, "radii")

toto5 <- varlogfpt(i5)
toto5
attr(toto5, "radii")

plot(ltraj[6])
i6 <- fpt(ltraj[6], seq(300,1000, length=30))
plot(i6, scale = 500, warn = FALSE)

plot(ltraj[7])
i7 <- fpt(ltraj[7], seq(300,1000, length=30))
plot(i7, scale = 500, warn = FALSE)

toto7 <- meanfpt(i7)
toto7
attr(toto7, "radii")

toto7 <- varlogfpt(i7)
toto7
attr(toto7, "radii")


is.regular(ltraj[1])
plotltr(ltraj[4], "dt")
windows()
plotltr(ltraj[1], "dist")

is.regular(ltraj[2])
plotltr(ltraj[7], "dt")
windows()
plotltr(ltraj[2], "dist")
ltraj[2]

#Code to export each trajectory as a shapefile if needed
##Define the projection of the coordinates
toto1 <-ltraj2sldf(ltraj[1])
plot(toto1)
writeOGR(toto1,dsn=".",layer="D12", driver = "ESRI Shapefile",overwrite=TRUE)
summary(toto1)

#Write lines and points as a shapefile
toto2lines <-ltraj2sldf(ltraj[2],byid=TRUE)
toto2pts <- ltraj2spdf(ltraj[2])

#If we want to define projection before making a shapefile
proj4string <- CRS("+proj=utm +zone=13N +ellps=WGS84")
toto2lines@proj4string <- proj4string
toto2pts@proj4string <- proj4string

plot(toto2pts)
plot(toto2lines, add=T)

writeOGR(toto2pts,dsn=".",layer="D15pts", driver = "ESRI Shapefile",overwrite_layer=TRUE)
writeOGR(toto2lines, dsn=".", paste("traj_line_",sep=""),driver = "ESRI Shapefile",overwrite=TRUE)


toto3 <-ltraj2sldf(ltraj[3])
plot(toto3)
writeOGR(toto3,dsn=".",layer="D16", driver = "ESRI Shapefile",overwrite=TRUE)

toto4 <-ltraj2sldf(ltraj[4])
plot(toto4)
writeOGR(toto4,dsn=".",layer="D19", driver = "ESRI Shapefile",overwrite=TRUE)

toto5 <-ltraj2sldf(ltraj[5])
plot(toto5)
writeOGR(toto5,dsn=".",layer="D4", driver = "ESRI Shapefile",overwrite=TRUE)

toto6 <-ltraj2sldf(ltraj[6])
plot(toto6)
writeOGR(toto6,dsn=".",layer="D6", driver = "ESRI Shapefile",overwrite=TRUE)

toto7 <-ltraj2sldf(ltraj[7])
plot(toto7)
writeOGR(toto7,dsn=".",layer="D8", driver = "ESRI Shapefile",overwrite=TRUE)

#Another look at the time and distance between relocations of each animal.
#The code can be used initially to inspect the data if there is concern 
#about consistancy in location fixes or distance between each location.
is.regular(ltraj[3])
plotltr(ltraj[3], "dt")
plotltr(ltraj[3], "dist")

