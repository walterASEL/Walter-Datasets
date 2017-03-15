library(adehabitatLT)
library(chron)
library(raster)
library(sp)

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
muley15[1:10,]#code displays the first 20 records to look at what sorting did to data

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

### Creation of an object of class "ltraj", with for ### example the first animal
ltraj <- as.ltraj(coordinates(merge),merge$da,id=merge$id)
plot(ltraj)
ltraj

## We want to study the trajectory of the day at the scale
## of the day. We define one trajectory per day. The trajectory should begin
## at 22H00
## The following function returns TRUE if the date is comprised between
## 06H00 and 23H00 (i.e. results in 3 locations/day bursts)
foo <- function(date) {
da <- as.POSIXlt(date)
ho <- da$hour + da$min
return(ho>18.0&ho<23.9)
}
deer <- cutltraj(ltraj, "foo(date)", nextr = TRUE)
deer

## Remove the first and last burst if needed?
#deer2 <- deer[-c(1,length(deer))]

#bind the trajectories
deer3 <- bindltraj(deer)
deer3
plot(deer3)
is.regular(deer)
plotltr(deer3, "dist")

## The relocations have been collected every 3 hours, and there are some
## missing data
## The reference date: the hour should be exact (i.e. minutes=0):
refda <- strptime("00:00", "%H:%M")
refda
## Set the missing values
deerset <- setNA(deer3, refda, 3, units = "hour")
## now, look at dt for the bursts:
plotltr(deerset, "dt")
## dt is nearly regular: round the date:
deerset1 <- sett0(deerset, refda, 3, units = "hour")
plotltr(deerset1, "dt/3600")
is.regular(deerset1)
## deerset1 is now regular

## Is the resulting object "sd" ?
is.sd(deerset1)
str(deerset1)

#Show the changes in the distance between successive relocations with the time
plotltr(deerset1, "dist")
## Segmentation of the trajectory based on these distances
lav <- lavielle(deerset1, Lmin=2, Kmax=20)
## Choose the number of segments
chooseseg(lav)
## 20 segments seem a good choice
## Show the partition
kk <- findpath(lav, 20)
kk

#Notice that the results show for each burst:
#(1) number of relocations
#(2) number of relocations removed (i.e., NA)
#(3) begin and end dates

#Now if we reduce the number of segments we get the following bursts:
## Segmentation of the trajectory based on these distances
lav <- lavielle(deerset1, Lmin=2, Kmax=10)
## Choose the number of segments
chooseseg(lav)
## 20 segments seem a good choice
##Show the partition
k <- findpath(lav, 10)
kk

#*********** List of class ltraj ***********
#Type of the traject: Type II (time recorded)
#Regular traject. Time lag between two locs: 10800 seconds

#Characteristics of the bursts:
plot(kk[1])
plot(kk[2])
plot(kk[3])
plot(kk[6])
plot(kk[7])
plot(kk[9])