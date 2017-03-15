library(SDMTools)
library(raster)
library(rgeos)
library(plyr)

rm(list=ls())

#Load vegetation raster layer textfile clipped in ArcMap 
crops <-raster("crop2012utm12.tif")
plot(crops)
class(crops)
as.matrix(table(values(crops)))
proj4string(crops)
crops

# reclassify the values into 9 groups
# all values between 0 and 20 equal 1, etc.
m <- c(-Inf,0,NA,2, 7, 2, 20, 60, 3, 60, 70, 4, 110, 132, 5, 133, 150, 6, 151, 191, 7, 192,Inf,NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(crops, rclmat)
plot(rc)
rc
as.matrix(table(values(rc)))

########################################################
#What if we wanted to compare difference in patch
#statistics among all deer with all locations combined?
########################################################

#Let's create buffers around individual locations
muleys <-read.csv("muleysexample.csv", header=T)
summary(muleys$id)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.80)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys
table(muleys$id)

muleys$GPSFixTime<-as.POSIXct(muleys$GPSFixTime, format="%Y.%m.%d%H:%M:%S")

#Only use the 8 lines of code below to subsample for demonstration purposes!
onemuley <-muleys[1:5,]
onemuley
twomuley <-muleys[102:106,]
twomuley
shortmd <- rbind(onemuley, twomuley)
str(shortmd)
shortmd$id <- factor(shortmd$id)#removed deer D12 because to few locations
muleys <- shortmd
muleys
##########################################################################
###########################################################################
buff3rd <- function(muleys) {
  coords<-data.frame(x = muleys$X, y = muleys$Y)
  deer.spdf <- SpatialPointsDataFrame(coords=coords, data = muleys, proj4string = CRS("+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs +datum=GRS80 +towgs84=0,0,0"))
  settbuff <- gBuffer(deer.spdf, width=1000, byid=FALSE)
  buffclip <- mask(rc, settbuff)
  buff.data <- PatchStat(buffclip)
  newline <- muleys$id
  bind <-cbind(newline[1], buff.data)
}

results <- ddply(muleys, .(id), buff3rd)
results
class(results)

###########################################################################
###########################################################################

################################################################
#Code above looks at patch and class metrics for each deer
#by combining all buffers into one polygon for each deer (i.e., 
#to define available habitat in 3rd order selection.   
#However, what if we wanted to compare difference in patch
#statistics among all deer by averaging metrics across buffers?
#################################################################

coords<-data.frame(x = muleys$X, y = muleys$Y)
deer.spdf <- SpatialPointsDataFrame(coords=coords, data = muleys, proj4string = CRS("+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs +datum=GRS80 +towgs84=0,0,0"))
setbuff <- gBuffer(deer.spdf, width=1000, byid=TRUE)
setbuff
muleys$newID <- paste(muleys$id, setbuff@plotOrder, sep="_")

buff3rdA <- function(muleys) {
  bufclip <- mask(rc, setbuff)
  buf.data <- PatchStat(bufclip)
}

results2 <- ddply(muleys, .(newID), buff3rdA)
results2
