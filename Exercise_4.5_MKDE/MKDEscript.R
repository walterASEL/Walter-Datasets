library(adehabitatHR)
library(adehabitatLT)
library(sp)
library(rgdal)
library(raster)
library(chron) 
library(rgeos)#for function "crop"

#We then need to create a Spatial Pixels Data Frame of the ascii grid we imported
habitat = as(readGDAL("beauzoom100.asc"), "SpatialPixelsDataFrame")
#Results: perimclipasc.asc has GDAL driver AAIGrid 
#and has 2801 rows and 2962 columns
str(habitat)
image(habitat)
#proj4string(habitat) <- CRS("+proj=utm +zone=17N +ellps=WGS84")
str(habitat)

#Creates a Spatial Points Data Frame for 2 animals by ID
twobirds <-read.csv("Twobirds.csv", header=T)
twobirds$id <-as.factor(twobirds$id)

#Needs to be done to get proper digits of date into R then POSIXct
xtime <- paste(twobirds$OldDate,twobirds$Hour)
twobirds$PosTime <- xtime

#Calculates time difference to use as dt
twobirds$date_time <- chron(as.character(twobirds$OrigDate),twobirds$Hour, format=c(dates="m/d/y", times="h:m:s"))
timediff <- diff(twobirds$date_time)*24*60
twobirds <-twobirds[-1,]
twobirds$timediff <-as.numeric(abs(timediff))

data.xy = twobirds[c("x","y")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)

#Creates a spatial data frame of ID
idsp<-data.frame(twobirds[2])
#Creates a spatial data frame of dt
dtsp<-data.frame(twobirds[17])
#Creates a spatial data frame of Burst
busp<-data.frame(twobirds[19])
#Merges ID and Date into the same spatial data frame
merge<-data.frame(idsp,dtsp,busp)
#Adds ID and Date data frame with locations data frame
coordinates(merge)<-sppt
proj4string(merge) <- CRS("+proj=utm +zone=17N +ellps=WGS84")

#Cast the Dates as POSIXct dates 
merge$DT <-as.POSIXct(strptime(merge$PosTime, format="%Y%m%d %H:%M:%S"))

#Create an ltraj trajectory object. 
ltraj <- as.ltraj(coordinates(merge), merge$DT, id = merge$id, burst = merge$id, typeII = TRUE)
plot(ltraj)

plot(ltraj, spixdf=habitat)

#Be sure to do this step after plotting ltraj onto spixdf or won't work!
#This step just builds a "fake" habitat map with habitat=1
fullgrid(habitat) <- TRUE
hab <- habitat
hab[[1]] <- as.numeric(!is.na(hab[[1]])) 
image(hab)
 
#This step is needed to convert SpatialGrid to SpatialPixels for use in "ud" estimation if needed
#"habitat" in "grid=habitat" must be of class SpatialPixels
fullgrid(hab) <- FALSE
class(hab)
image(hab)

#CODE TO CONDUCT BRB
#Assigne parameter values for BRB 
# Parameters for the Biased Random Bridge Kernel approach 
tmax <- 1*(24*60*60) + 1 ## set the maximum time between locations to be just more than 1 day
lmin <- 50 ## locations less than 50 meters apart are considered inactive.   
hmin <- 100 ## arbitrarily set to be same as hab grid cell resolution 

#Diffusion component for each habitat type using plug-in method
vv<- BRB.D(ltraj, Tmax = tmax, Lmin = lmin,  habitat = hab)
vv

ud <- BRB(ltraj, D = vv, Tmax = tmax, Lmin = lmin, hmin=hmin, grid = hab, b=TRUE, extent=0.1, tau = 300)
ud
image(ud)

#Address names in ud by assigning them to be the same as the ids in ltraj
#Must be done before using "getverticeshr" function
names(ud) <- id(ltraj) 

par(mfrow=c(2,3))

ver1_95 <- getverticeshr(ud, percent=95, standardize = TRUE, whi = id(ltraj[1]))
plot(ver1_95)

ver2_95 <- getverticeshr(ud, percent=95, standardize = TRUE, whi = id(ltraj[2]))
plot(ver2_95)

#Now let's create a new UD using an actual habitat layer that has more than
#"used/unused"
#Start by importing the habitat layer again and run the following
habitat2 = as(readGDAL("beauzoom100.asc"), "SpatialPixelsDataFrame")
plot(ltraj, spixdf=habitat2)

#CODE TO CONDUCT BRB
#Assigne parameter values for BRB 
# Parameters for the Biased Random Bridge Kernel approach 
tmax <- 1*(24*60*60) + 1 ## set the maximum time between locations to be just more than 1 day
lmin <- 50 ## locations less than 50 meters apart are considered inactive.   
hmin <- 100 ## arbitrarily set to be same as hab grid cell resolution 

#Diffusion component for each habitat type using plug-in method
vv2<- BRB.D(ltraj, Tmax = tmax, Lmin = lmin,  habitat = habitat2)
vv2

ud2 <- BRB(ltraj, D = vv2, Tmax = tmax, Lmin = lmin, hmin=hmin, habitat = habitat2, b=TRUE, extent=0.1, tau = 300, same4all=FALSE)
image(ud2)

names(ud2) <- id(ltraj) 

ver1a <- getverticeshr(ud2, percent=95, standardize = TRUE, whi = id(ltraj[1]))
plot(ver1a)

ver2a <- getverticeshr(ud2, percent=95, standardize = TRUE, whi = id(ltraj[2]))
plot(ver2a)


#Diffusion component for each habitat type using plug-in method
vv3<- BRB.D(ltraj, Tmax = tmax, Lmin = lmin,  habitat = NULL)
vv3

ud3 <- BRB(ltraj, D = vv3, Tmax = tmax, Lmin = lmin, hmin=hmin, habitat = NULL, b=TRUE, extent=0.1, tau = 300, same4all=FALSE)
image(ud3)

names(ud3) <- id(ltraj) 

ver1b <- getverticeshr(ud3, percent=95, standardize = TRUE, whi = id(ltraj[1]))
plot(ver1b)

ver2b <- getverticeshr(ud3, percent=95, standardize = TRUE, whi = id(ltraj[2]))
plot(ver2b)


#Some extraneous code to plot more detailed figures if needed.

#Plot bird 49
#udvol1<-getvolumeUD(ud[[1]], standardize = FALSE)
#myPal <- colorRampPalette( c("red","orange","yellow") )
#udvoltmp1<-udvol1
#udvoltmp1@data$n<- ifelse(udvoltmp1@data$n>=94.9,NA,udvoltmp1@data$n)
#image(udvoltmp1,col=myPal(64),frame.plot=FALSE)
#title(main=paste("MDKE","49",sep=" "),line=0,cex.main=1)

#windows()
#Plot bird 50
#udvol2<-getvolumeUD(ud[[2]], standardize = FALSE)
#myPal <- colorRampPalette( c("red","orange","yellow") )
#udvoltmp2<-udvol2
#udvoltmp2@data$n<- ifelse(udvoltmp2@data$n>=99.9,NA,udvoltmp2@data$n)
#image(udvoltmp2,col=myPal(64),frame.plot=FALSE)
#title(main=paste("MDKE","50",sep=" "),line=0,cex.main=1)