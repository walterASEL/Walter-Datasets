library(adehabitatLT)
library(chron)
library(spatstat)#for "duplicate" function

muleys <-read.csv("DCmuleysedited.csv", header=T)
str(muleys)

#Check for duplicate locations in dataset
summary(duplicated(muleys))

#Sort data to address error in code if needed
#muleys <- muleys[order(muleys$id),]

### Conversion of the date to the format POSIX
#Date <- as.character(muleys$GPSFixTime)
#Date <- as.POSIXct(strptime(as.character(muleys$GPSFixTime),"%Y.%m.%d %H:%M:%S"))
#muleys$Date <- Date


######################################################
##
## Example of a trajectory of type II (time recorded)
### Conversion of the date to the format POSIX
#Needs to be done to get proper digits of date into R then POSIXct
da <- as.POSIXct(strptime(muleys$GPSFixTime,format="%Y.%m.%d %H:%M:%S"))
da
muleys$da <- da
str(muleys)

timediff <- diff(muleys$da)*60
muleys <-muleys[-1,]
muleys$timediff <-as.numeric(abs(timediff)) 
str(muleys)#check to see timediff column was added to muleys

summary(muleys$id)

newmuleys <-subset(muleys, muleys$X > 599000 & muleys$X < 705000 & muleys$Y > 4167000 & muleys$timediff < 14401)
muleys <- newmuleys
str(muleys)

data.xy = muleys[c("X","Y")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
#proj4string(xysp) <- CRS("+proj=utm +zone=13 +ellps=WGS84")

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
head(ltraj[1])#Describes the trajectory
plot(ltraj[1])
plot(ltraj[2])
plot(ltraj[3])
plot(ltraj[4])
plot(ltraj[5])
plot(ltraj[6])
plot(ltraj[7])

hist(ltraj[1], "dt", freq = TRUE)
hist(ltraj[1], "dist", freq = TRUE)
hist(ltraj[2], "dt", freq = TRUE)
hist(ltraj[2], "dist", freq = TRUE)
hist(ltraj[3], "dt", freq = TRUE)
hist(ltraj[3], "dist", freq = TRUE)
hist(ltraj[4], "dt", freq = TRUE)
hist(ltraj[4], "dist", freq = TRUE)
hist(ltraj[5], "dt", freq = TRUE)
hist(ltraj[5], "dist", freq = TRUE)
hist(ltraj[6], "dt", freq = TRUE)
hist(ltraj[6], "dist", freq = TRUE)
hist(ltraj[7], "dt", freq = TRUE)
hist(ltraj[7], "dist", freq = TRUE)
