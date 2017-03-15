library(adehabitatHR)
library(rgdal)
library(gstat)
library(plyr)

# Weather-Station-Code_Snowfall.R
#
# This code is designed to process data downloaded from Climate Date Online.
# <http://www.ncdc.noaa.gov/cdo-web/> This version looks at weather stations 
# that provide snowfall data in and around Pennyslvania. The code pulls out
# the desired data from the downloaded aggregate weather-station file and 
# calculates the mean annual snowfall per weather station for 12/1/94 - 3/31/05.
# The data is then exported to a text file for interpolation in ArcGIS.
# - Bill Kanapaux, PA Cooperative Fish & Wildlife Research Unit
#
# Last modified: Jan. 13, 2014
#

### remove other objects from the working session ###
rm(list=ls())

### read weather station data - note the use of stringsAsFactors=FALSE and na.strings='-9999' ###
WS <-read.table('Weather_Station_Data-Sep_01Dec1994-31March2005.txt',stringsAsFactors=FALSE, na.strings='-9999',header=T)

### check data ###
dim(WS)
head(WS)
summary (WS)

### Reformat DATE and create Year Month Day columns from NewDate column ###
WS$NewDate <- as.Date(as.character(WS$DATE), format("%Y%m%d"))
WS$Year = as.numeric(format(WS$NewDate, format = "%Y"))
WS$Month = as.numeric(format(WS$NewDate, format = "%m"))
WS$Day = as.numeric(format(WS$NewDate, format = "%d"))

head(WS)

### make a subset of WS that includes only the months of Dec-March ###
Winter <- WS[WS$Month %in% c(1,2,3,12), ]

### For December, add 1 to Year so that Year matches Jan-March in that season ###
Winter <- within(Winter, Year[Month==12] <- Year[Month==12] +1)

### check subset, including random row to make sure only selected months included ###
dim(Winter)
head(Winter)
Winter[699,]

### Create a matrix of unique STATION values (GHCND ) with Lat/Long values for later reference. ###
### Data contains some multiple versions of individual GHCND coordinates. Only want 1 set per.  ###
PulledCoords <- Winter[!duplicated(Winter[,1]),]
head(PulledCoords)
dim(PulledCoords)

CoordChart <- ddply(PulledCoords, c('STATION'), function(x) c(Lat=x$LATITUDE, Long=x$LONGITUDE))
head(CoordChart)

### Get the number of snowfall records for each STATION for each year and name it ###
### RecordTotal. Note that NA is omited from the length count ###
WinterRecords <- ddply(Winter, .(STATION,Year), summarize, RecordTotal = length(na.omit(SNOW)))
head(WinterRecords)
tail(WinterRecords)
dim(WinterRecords)

### Get the total amount of snowfall per STATION per year and name it YearlySnow ###
YearlySnow <- ddply(Winter, .(STATION,Year), summarize, Snow = sum(SNOW, na.rm=TRUE))
head(YearlySnow)
tail(YearlySnow)
dim(YearlySnow)

### Combine WinterRecords and YearlySnow into one matrix ###
AllWinters <- cbind(WinterRecords,YearlySnow)
AllWinters <- AllWinters[,-4:-5]
head(AllWinters)
tail(AllWinters)
dim(AllWinters)

### Only include years that have more than 75% of days recorded ###
WinterDays <- 121
FullWinters <- AllWinters[AllWinters$RecordTotal/WinterDays > 0.75, ]
head(FullWinters)
tail(FullWinters)
dim(FullWinters)

### Get the number of years with more than 75% of days recorded for each STATION ###
WinterYears <- ddply(FullWinters, c('STATION'), function(x) c(TotalYears=length(x$Year)))
head(WinterYears)
tail(WinterYears)
dim(WinterYears)

### Get the total amount of snow for each station for all years ###
TotalWinterSnow <- ddply(FullWinters, c('STATION'), function(x) c(TotalWinterSnow=sum(x$Snow)))
head(TotalWinterSnow)
dim(TotalWinterSnow)

### Combine WinterYears and TotalWinterSnow into one matrix ###
SnowCalc <- cbind(WinterYears,TotalWinterSnow)
SnowCalc <- SnowCalc[,-3]
head(SnowCalc)

### Get rid of the stations that don't have at least 10 years recorded at >75% of days ###
Complete.Records <- SnowCalc[SnowCalc$TotalYears > 9, ]
head(Complete.Records)
dim(Complete.Records)

### Calculate average annual snowfall and round to nearest mm ###
Complete.Records$MeanAnnualSnowfall <- Complete.Records$TotalWinterSnow/Complete.Records$TotalYears
Complete.Records$MeanAnnualSnowfall <- round (Complete.Records$MeanAnnualSnowfall, digits = 0)
head(Complete.Records)

### Convert SnowDepth from mm to cm
Complete.Records$MeanAnnualSnowfall <- Complete.Records$MeanAnnualSnowfall/10
head(Complete.Records)

### Add a column to CoordChart showing whether each row matches  a STATION in Complete.Records
### Use "NA" for value if no match, then delete rows with "NA" value. 
### Number of rows in CoordChart should now equal number of rows in Complete.Records
CoordChart$match <- match(CoordChart$STATION, Complete.Records$STATION, nomatch=NA)
CoordChart <- na.omit(CoordChart)
head(CoordChart)
dim(CoordChart)
dim(Complete.Records)

### Combine Complete.Records and CoordChart. Make sure each STATION matches in row
### Delete any rows that don't match. Shouldn't be any. If # of rows in Final.Values
### is less than # of rows in CoordChart, there is a problem (but note that # of cols does change).
Final.Values <- cbind(Complete.Records,CoordChart)
Final.Values$match2 <- match(Final.Values[  ,1], Final.Values[ ,5], nomatch=NA)
Final.Values <- na.omit(Final.Values)
dim(Final.Values)
dim(CoordChart)

head(Final.Values)

### Take out unnecssary rows (2nd STATION, match, and match2) and round MeanSnow to 2 decimal places
Final.Values[,5] <- Final.Values[,8] <- Final.Values[,9] <- NULL
head(Final.Values)

### Make data frame to get rid of lists (in R) so can export to text file
### Use text file to load weather station points into ArcGIS
#Final.Values <- as.data.frame(lapply(Final.Values,unlist))
#write.table(Final.Values, "MeanSnowData_95-05.txt", sep="\t", row.names=F)

#Alternatively we can conduct interpolation directly in R using the steps below.

#Need to convert factors to numeric
Final.Values$Longitude <- as.numeric(as.character(Final.Values$Long))
Final.Values$Latitude <- as.numeric(as.character(Final.Values$Lat))

#Here we need to create Spatial points, attach ID and Date Time sorted
#Then transform UTM to Lat Long because dBBMM will only work in AEQD projection
data.xy = Final.Values[c("Longitude","Latitude")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
#proj4string(xysp) <- CRS("+proj=longlat +ellps=WGS84")

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)

#Creates a spatial data frame of STATION
ID<-data.frame(Final.Values[1])
#Creates a spatial data frame of Mean Annual Snow Fall
MAS<-data.frame(Final.Values[4])
#Merges ID and Date into the same spatial data frame
merge<-data.frame(ID,MAS)
#Adds ID and Date data frame with locations data frame
coordinates(merge)<-sppt
proj4string(merge) <- CRS("+proj=longlat +ellps=WGS84")

#Import a county layer for study site and check projections
counties<-readOGR(dsn=".",layer="PACountiesAlbers")
proj4string(counties)
plot(counties)

#Project Weather Stations to match Counties
Albers.crs <-CRS("+proj=aea +lat_1=29.3 +lat_2=45.3 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
stations <- spTransform(merge, CRS=Albers.crs)

#Plot Weather Stations over counties
points(stations)

stations@data$MAS <- stations@data$MeanAnnualSnowfall

## create a grid onto which we will interpolate:
xx = spsample(counties, type="regular", cellsize=10000)
class(xx)

points(xx, bg="red", cex=.5,col="red")

#Convert to a SpatialPixels class
gridded(xx) <- TRUE
class(xx)

plot(xx)
points(stations, bg="red", cex=.5,col="red")
str(stations)

#Plot out the MAS across the study region
bubble(stations, zcol='MAS', fill=FALSE, do.sqrt=FALSE, maxsize=2, add=T)

## create a grid onto which we will interpolate:
## first get the range in data
 x.range <- as.integer(range(stations@coords[,1]))
 y.range <- as.integer(range(stations@coords[,2]))

## now expand to a grid with 500 meter spacing:
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=5000), y=seq(from=y.range[1], to=y.range[2], by=5000) )

## convert to SpatialPixel class
 coordinates(grd) <- ~ x+y
 gridded(grd) <- TRUE

## test it out:
plot(grd, cex=0.5)
points(stations, pch=1, col='red', cex=0.7)
title("Interpolation Grid and Sample Points")


x <- krige(stations@data$MAS~1, stations, xx)
class(x)
image(x)
points(stations, pch=1, col='blue', cex=0.7)