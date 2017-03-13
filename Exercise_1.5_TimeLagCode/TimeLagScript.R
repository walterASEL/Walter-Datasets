library(chron)
library(rgdal)
library(RAtmosphere)

# import file 
temp <- read.csv("Y2005_UTM_date.csv", header=T)
str(temp)

# modify time to include seconds 
temp$time <- paste(as.character(temp$LMT_TIME),"00",sep=":") 

# convert to chron date 
temp$date_time <- chron(as.character(temp$LMT_DATE),temp$time,format=c(dates="m/d/y",times="h:m:s")) 

# calculate difference in time in minutes 
timediff <- diff(temp$date_time)*24*60 
summary(timediff)

# remove first entry without any difference 
temp <- temp[-1,] 

# assign timediff column to original "temp" dataset for use later
temp$timediff <- as.numeric(timediff) 
str(temp)

#we can also export the file for use in excel or other programs
#write.table(temp, "TimeDiffdata.txt", row.names=TRUE,sep=" ", col.names=TRUE, quote=TRUE, na="NA")

#NOTE: Code below is to include night and day into datasets and to also to account for daylight savings.
#Package RAtmosphere will eliminate the need for the chunk of code below if working on an earlier version
#of R (i.e., code not available for R 3.2.1 plus).
#Thanks to Duane Diefenbach, PA Coop Unit Leader, for compiling all this from online sources.

##############################################################################################
#We first need to create a SPDF and transform to Lat Long then return to a Data Frame
#You only need this section of code if you need to acquire Lat Long coordinates for dataset
#############################################################################################
#utm.crs <-CRS("+proj=utm +zone=12 +datum=WGS84")
#dataSPDF<-data.frame(x = temp$UTMe, y = temp$UTMn)
#utm.spdf <- SpatialPointsDataFrame(coords = dataSPDF, data = temp, proj4string = utm.crs)

#ll.crs <- CRS("+proj=longlat +ellps=WGS84")
#datall <-spTransform(utm.spdf, CRS=ll.crs)
#str(datall)
#temp <- as.data.frame(datall)
#str(temp)
##############################################################################################
#Separate times into categories "Day" and "Night" based on sunrise-sunset table by running function
#below or simply using the RAtmosphere package
#############################################################################################
suncalc <- function(d,Lat=39.14133,Long=-106.7722){
  # code obtained from http://www.r-bloggers.com/approximate-sunrise-and-sunset-times/
  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)
  
  ##This method is copied from:
  ##Teets, D.A. 2003. Predicting sunrise and sunset times.
  ##  The College Mathematics Journal 34(4):317-321.
 
  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.

  ## Function to convert degrees to radians
  rad <- function(x)pi*x/180
  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)

  ##Convert observer's latitude to radians
  L=rad(Lat)

  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  ##NOTE: If working with UTC times use timezone = -4*(abs(Long)%%15)*sign(Long)
  timezone = -6*(abs(Long)%%15)*sign(Long)

  ## The earth's mean distance from the sun (km)
  r = 149598000

  theta = 2*pi/365.25*(d-80)

  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)

  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 

  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)

  ## now sunrise and after sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  suntime <- cbind(sunrise,sunset)  

  return(suntime)
}

############## Read in location data and retain lat, lon, and date
temp$Date <- paste((temp$Year),substr(temp$date_time, 2,3),substr(temp$date_time, 5,6),sep="-")
str(temp)

#calculate calendar day and center of locations
calday <- as.numeric(as.Date(temp$Date)-as.Date("2005-01-01"), units="days")
dat1 <- cbind(temp,calday)
moda <- format(as.Date(temp$Date),"%d-%b")
str(dat1)

dat1 <- cbind(dat1, suncalc(dat1$calday, Lat=dat1$LATITUDE, Long=dat1$LONGITUDE),moda)
hrchar <- as.character(substr(dat1$time,1,2))
hr <- as.numeric(as.character(substr(dat1$time,1,2)))
minchar <- as.character(substr(dat1$time,4,5))
min <- as.numeric(minchar)
localhr <- hr+min/60
dat1 <- cbind(dat1,hr,hrchar,minchar,localhr)
Diel <- ifelse(localhr<dat1$sunrise | localhr>dat1$sunset, 'Night', 'Day')
dat1 <- cbind(dat1,Diel)
str(dat1)
dat1[1:50,]

###############################################################################
############################################################################################



