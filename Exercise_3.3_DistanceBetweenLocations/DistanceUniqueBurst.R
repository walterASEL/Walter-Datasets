library(adehabitatLT)
library(chron)
library(class)
library(Rcmdr)

muleys <-read.csv("DCmuleysedited.csv", header=T)
str(muleys)

#CODE FOR AN INDIVIDUAL ANIMAL
muley15 <- subset(muleys, id=="D15")
muley15[1:20,]
str(muley15)
summary <- table(muley15$UTM_Zone,muley15$id)
summary
muley15$id <- factor(muley15$id)

#Sort data to address error in code and then look at first 20 records of data to confirm
muley15 <- muley15[order(muley15$GPSFixTime),]
muley15[1:10,]#code displays the first 20 records to look at what sorting did to data

######################################################
## Example of a trajectory of type II (time recorded)
### Conversion of the date to the format POSIX
#Needs to be done to get proper digits of date into R then POSIXct
#uses library(chron)
da <- as.character(muley15$GPSFixTime)
da <- as.POSIXct(strptime(muley15$GPSFixTime,format="%Y.%m.%d %H:%M:%S"))
head(da)
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

#Now let's look at time differences between locations before moving forward
summary(muley15$timediff)

## We want to study the trajectory of the day at the scale
## of the day. We define one trajectory per day. The trajectory should begin
## at 22H00
## The following function returns TRUE if the date is comprised between
## 06H00 and 23H00 (i.e. results in 3 locations/day bursts)
foo <- function(date) {
da <- as.POSIXlt(date)
ho <- da$hour + da$min
return(ho>15.9&ho<23.9)
}
deer <- cutltraj(ltraj, "foo(date)", nextr = TRUE)

#Notice that the above code will remove 345 relocations that fall
#outside of your time criteria
#Warning message:
#In cutltraj(ltraj, "foo(date)", nextr = TRUE) :
#  At least 3 relocations are needed for a burst
# 345 relocations have been deleted


head(deer)

#Shows results of cutting the traj into individual bursts
#*********** List of class ltraj ***********

#Type of the traject: Type II (time recorded)
#Irregular traject. Variable time lag between two locs

#Characteristics of the bursts:
#     id   burst nb.reloc NAs          date.begin            date.end
#1   D15 D15.001        6   0 2011-10-12 03:00:52 2011-10-12 18:00:52
#2   D15 D15.003        7   0 2011-10-13 00:00:35 2011-10-13 18:00:35
#3   D15 D15.005        7   0 2011-10-14 00:00:42 2011-10-14 18:00:42
#4   D15 D15.007        7   0 2011-10-15 00:00:35 2011-10-15 18:00:45
#5   D15 D15.009        7   0 2011-10-16 00:00:39 2011-10-16 18:00:49
#6   D15 D15.011        6   0 2011-10-17 00:01:07 2011-10-17 15:01:03
#7   D15 D15.014        7   0 2011-10-18 00:00:34 2011-10-18 18:00:48
#8   D15 D15.016        7   0 2011-10-19 00:00:36 2011-10-19 18:00:40
#9   D15 D15.018        7   0 2011-10-20 00:00:53 2011-10-20 18:00:40
#10  D15 D15.020        7   0 2011-10-21 00:00:39 2011-10-21 18:00:37

#Code to change ltraj to a data.frame to summarize distance between 
#locations for each daily burst
dfdeer <- ld(deer)
head(dfdeer)
str(dfdeer)

#'data.frame':   2243 obs. of  13 variables:
# $ x        : num  677932 679037 679429 679750 679453 ...
# $ y        : num  4189551 4189493 4189406 4189053 4188461 ...
# $ date     : POSIXct, format: "2011-10-12 03:00:52" "2011-10-12 06:00:52" ...
# $ dx       : num  1105 392 321 -297 163 ...
# $ dy       : num  -58 -87 -353 -592 -89 NA -189 756 395 95 ...
# $ dist     : num  1107 402 477 662 186 ...
# $ dt       : num  10800 10786 10796 10808 10810 ...
# $ R2n      : num  0 1224389 2262034 3553128 3501541 ...
# $ abs.angle: num  -0.0524 -0.2184 -0.8328 -2.0358 -0.4998 ...
# $ rel.angle: num  NA -0.166 -0.614 -1.203 1.536 ...
# $ id       : Factor w/ 1 level "D15": 1 1 1 1 1 1 1 1 1 1 ...
# $ burst    : Factor w/ 325 levels "D15.001","D15.003",..: 1 1 1 1 1 1 2 2 2 2 ...
# $ pkey     : Factor w/ 2588 levels "D15.2011-10-12 03:00:52",..: 1

#Code to get mean distance moved for each burst
summary <- numSummary(dfdeer[,"dist"],groups=dfdeer$burst, statistics=c("mean","sd"))

#Convert matrix from data.frame to export into csv file
mean <- as.matrix(summary$table)
head(mean)
#Write.table gives csv output of Summary.  Be sure to specify the directory and the output files will be stored there 
write.table(mean, file = "Distance.csv", sep =",", row.names = TRUE, col.names = TRUE, qmethod ="double")


#Convert "mean" matrix from above to dataframe to summarize in R
meantable <- as.data.frame(mean)
BufferRadius <- numSummary(meantable, statistics=c("mean","sd"))
BufferRadius