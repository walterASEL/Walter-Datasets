library(trip)
library(stringr)
library(adehabitatHR)
library(lattice)
library(gmodels)
library(spatstat)
library(maptools)

muleys<-read.csv("DCmuleysedited.csv", header=T, sep=",")
str(muleys)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.50 & muleys$Lat > 37.3 & muleys$Long < -107)
muleys <- newmuleys

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$Long, y = muleys$Lat)
crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
head(coords)
plot(coords)

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(crs))
head(deer.spdf)
class(deer.spdf)
proj4string(deer.spdf)
plot(deer.spdf,col=deer.spdf$id)

muleys$NewDate<-as.POSIXct(muleys$GPSFixTime, format="%Y.%m.%d %H:%M:%S", origin="1970-01-01")

#TIME DIFF NECESSARY IN BBMM CODE
timediff <- diff(muleys$NewDate)*60*60
# remove first entry without any difference 
muleys <- muleys[-1,] 
muleys$timelag <-as.numeric(abs(timediff))
summary(muleys$timelag)
#Remove locations greater than 24 hours apart in time
muleys <- subset(muleys, muleys$timelag < 18000)

#However, this sample size represents multiple years of data so causes errors in running
#some home range estimators. Therefore, let's separate each deer into the years data
#are available with the name 048_2006 for example
muleys$Year <- format(muleys$NewDate, "%Y")
muleys <- subset(muleys, muleys$Year != "NA")
muleys$YearBurst <- c(paste(muleys$id,muleys$Year,sep="_"))
muleys$YearBurst <- as.factor(muleys$YearBurst)
str(muleys)
summary(muleys$YearBurst)

muleys <- subset(muleys, table(muleys$YearBurst)[muleys$YearBurst] > 100)
muleys$YearBurst <- factor(muleys$YearBurst)

range(muleys$NewDate)
#muleys$Year <- NULL
muleys$Year[muleys$NewDate > "2011-09-30 00:30:00" & muleys$NewDate < "2012-04-01 23:00:00"] <- 2011
muleys$Year[muleys$NewDate > "2012-03-31 00:30:00" & muleys$NewDate < "2012-10-01 23:00:00"] <- 2012
muleys$Year <- as.factor(muleys$Year)
muleys$YearBurst <- c(paste(muleys$id,muleys$Year,sep="_"))
muleys$YearBurst <- as.factor(muleys$YearBurst)
table(muleys$YearBurst)

muleys <- subset(muleys, table(muleys$YearBurst)[muleys$YearBurst] > 100)
muleys$YearBurst <- factor(muleys$YearBurst)

d1 <- muleys
str(d1)
#Code separate each animal into a shapefile or text file to use as "List" in Cumming and Cornelis 
# get input file
indata <- d1
innames <- unique(d1$YearBurst)
innames <- innames[1:10]#needs to be number of unique IDs
outnames <- innames
# begin loop to calculate home ranges
for (i in 1:length(innames)){
  data <- indata[which(indata$YearBurst==innames[i]),]
  if(dim(data)[1] != 0){
    #data <-data[c(-21)]
    # export the point data into a shp file
    data.xy = data[c("X", "Y")]
    coordinates(data.xy) <- ~X+Y
    sppt <- SpatialPointsDataFrame(coordinates(data.xy),data)
    proj4string(sppt) <- CRS("+proj=utm +zone=12 +datum=WGS84")
    #writePointsShape(sppt,fn=paste(outnames[i],sep="/"),factor2char=TRUE)
    #sppt <-data[c(-22,-23)] 
    write.table(sppt, paste(outnames[i],"txt",sep="."), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(paste(outnames[i],"txt",sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, "In_list.txt", append=TRUE)
    #The write.table line above should only be run once to create the In_list.txt file otherwise it rights all birds each time
}
}
##############################################################################
##############################################################################
#
#
#Code below to get NSD and best movement model for each bird
#
#
##############################################################################
##############################################################################

date()

# Reads the List file of GPS datasets

List<-read.table("In_list.txt",sep="\t",header=F)
head(List) #“List” contains the filenames of the 6 deer sets

# Generation of results vectors
ID <- rep(0,nrow(List))
LOCS <- rep(0,nrow(List))
MIGR <- rep(0,nrow(List))
MIXM <- rep(0,nrow(List))
DISP <- rep(0,nrow(List))
HORA <- rep(0,nrow(List))
NOMA <- rep(0,nrow(List))
ID <- rep(0,nrow(List))
LOCS <- rep(0,nrow(List))
MIGR <- rep(0,nrow(List))
MIXM <- rep(0,nrow(List))
DISP <- rep(0,nrow(List))
HORA <- rep(0,nrow(List))
NOMA <- rep(0,nrow(List))
AICC_1 <- rep(0,nrow(List))
AICC_2 <- rep(0,nrow(List))
AICC_3 <- rep(0,nrow(List))
AICC_4 <- rep(0,nrow(List))
AICC_5 <- rep(0,nrow(List))

minAIC <- rep(0,nrow(List))
d_AICC_1 <- rep(0,nrow(List))
d_AICC_2 <- rep(0,nrow(List))
d_AICC_3 <- rep(0,nrow(List))
d_AICC_4 <- rep(0,nrow(List))
d_AICC_5 <- rep(0,nrow(List))
LL_AICC_1 <- rep(0,nrow(List))
LL_AICC_2 <- rep(0,nrow(List))
LL_AICC_3 <- rep(0,nrow(List))
LL_AICC_4 <- rep(0,nrow(List))
LL_AICC_5 <- rep(0,nrow(List))
sumLL_AICC <- rep(0,nrow(List))
wi_AICC_1 <- rep(0,nrow(List))
wi_AICC_2 <- rep(0,nrow(List))
wi_AICC_3 <- rep(0,nrow(List))
wi_AICC_4 <- rep(0,nrow(List))
wi_AICC_5 <- rep(0,nrow(List))

#i=1 (use to test code before doing full run)
for(i in 1:nrow(List)) { 

coords<-read.table(as.character(List[i,]),sep="\t",header=T)
coords$DT<-as.POSIXct(coords$NewDate, format="%Y-%m-%d %H:%M:%S")

##make a data.frame of coordinates. Here the raw values are divided 
#by 1000 so that trajectories are calculated using km as the unit of measurement not meters
coord<-data.frame((coords$Y/1000),(coords$X/1000))    
# make ltraj: a trajectory of all the relocations
d2<-as.ltraj(coord,coords$DT,
coords$YearBurst,        #separate your data by individual.  
burst=coords$YearBurst, #burst is used to creat subdivisions within an individual.
typeII=TRUE)       #typeII can be TRUE: radio-track data, or FALSE: not time 
                   #recorded, such as tracks in the snow

#[1] "2007-06-05 16:00:00 EDT" "2015-03-12 14:00:00 EDT"
#you can now make your trajectory regular 
#firstly create a reference start time
#refda <- strptime("00:00:00", "%H:%M:%S")   #all relocations should be altered 
#to occur at 30 seconds past each minute

#you can now make your trajectory regular, as radio tracks tend to lose 
#a few seconds / minutes with each relocation
#firstly add "NA" for each missing location in your trajectory
#d3<-setNA(d2,refda,
#as.POSIXct("2007-06-01 06:00:00 EDT"), #any time before earliest timedate
#7200,            #stating there should be a location every 2 hours
#tol=7200,        #how many time units to search each side of expected location 
#units="sec")   #specifying the time units

#you can now make your trajectory regular 
#firstly create a reference start time
refda <- strptime("00:00:30", "%H:%M:%S")

##NOTE: The refda and d3 code above was not run as in Papworth because it results in too many
#relocations as "NA" that get removed below. Not quite sure the reason behind it being included?
#you can now make your trajectory regular 
d4<-sett0(d2, refda, 
10800,                         #stating the interval at which relocations should be
correction.xy =c("none"),   #if "cs" performs location correction based on the 
#assumption the individual moves at a constant speed 
tol=10800,       #how many time units to search either side of an expected location
units = "sec")  #specifying the time units
                              
#to view your regular trajectory of points with NA's
summary(d4)
#now calculating NSD for each point
datansd<-NULL
for(n in 1:length(summary(d4)[,1])) #stating that NSD should be 
#calculated separately for each burst
{
nsdall<-d4[[n]][,8]             #extracting the NSD for each location
nsdtimeall<-d4[[n]][,3]         #extracting the time for each location
nsdtimestartzero<-d4[[n]][,3]-d4[[n]][1,3]  
#extracting the time since trip start for each location
nsdid<-rep(as.vector(summary(d4)[n,1]),
length.out=summary(d4)[n,3])     
#extracting the individual associated with each location
nsdtrip<-rep(as.vector(summary(d4)[n,2]),length.out=summary(d4)[n,3])
#extracting the trip associated with each location
datansd1<-data.frame(nsdall,nsdtimeall,nsdtimestartzero,nsdid,nsdtrip)                  
#joining all these variables together in a data frame
datansd<-rbind(datansd,datansd1)                                                        
#joining all the data frames together
}
datansd$zero1<-as.numeric(unclass(datansd$nsdtimestartzero))                            
# making seconds since trip start numeric
datansd$zerostart<-datansd$zero1/60                                                     
#changing the time since trip start from seconds to minutes
datansd$minslitr2<-as.numeric(strftime(as.POSIXlt(datansd$nsdtimeall),
format="%M"))     
#making a vector of the hour of the day a location occured
datansd$hdaylitr2<-as.numeric(strftime(as.POSIXlt(datansd$nsdtimeall),
format="%H"))     
#making a vector of the minute in an hour a location occured
datansd$minsday<-((datansd$hdaylitr2*60)+datansd$minslitr2)                             
#calculating the minute in the day a location occured

summary(datansd)
datansd1<-na.omit(datansd)            #remove NA's


datansd1$coordinates<-coord           #add the coordinates for each point
#you now have the dataframe you need (datansd) to start analysis

#NSD 
#table(datansd1$nsdid)

#Now you can start modelling NSD using nlme. 
#Equations are from Bunnefeld at al (2011) A model-driven approach to quantify migration patterns: 
#individual, regional and yearly differences. 
#Journal of Animal Ecology 80: 466 - 476

#First we are going to model the data using nls, a least squares method,
#the simplest method and first method in Bunnefeld et al. 2011 (i.e., MIGRATION)
#that uses a double sigmoid or s-shaped function. 

###########################
##
##  MIGRATION
##
###########################

m1<-tryCatch(nls(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) + 
(-asym / (1 + exp((xmidB-zerostart)/scale2))), #Equation 1 in Bunnefeld et al. 2011
start = c(asym=15000000,xmidA=200000,xmidB=450000,scale1=1000,scale2=1000)                  
#these are the starting values for each parameter of the equation 
,data=na.omit(datansd1)),error=function(e)99999)   #this is the data
#summary(m1)        #this will print a summary of the converged model
#NOTE: The error function is simply to prevent the loop from crashing if model does not converge

###########################
##
##  MIXED MIGRATORY
##
###########################

m2 <-tryCatch(nls(nsdall ~  asymA /(1+exp((xmidA-zerostart)/scale1)) + 
(-asymB / (1 + exp((xmidB-zerostart)/scale2))), #Equation 2 in Bunnefeld et al. 2011
start = c(asymA=15000000,asymB=10000000, xmidA=200000,xmidB=450000,scale1=1000,scale2=1000)                  
#these are the starting values for each parameter of the equation 
,data=na.omit(datansd1)),error=function(e)99999)   #this is the data 
#summary(m2)        #this will print a summary of the converged model

###########################
##
##  DISPERSAL
##
###########################

m3 <-tryCatch(nls(nsdall ~  asym /(1+exp((xmid-zerostart)/scale)), #Equation 3 in Bunnefeld et al. 2011
start = c(asym=15000000,xmid=200000,scale=1000)                  
#these are the starting values for each parameter of the equation 
,data=na.omit(datansd1)),error=function(e)99999)   #this is the data
#summary(m3)        #this will print a summary of the converged model

###########################
##
## HOME RANGE
##
###########################

m4 <- tryCatch(nls(nsdall ~ intercept, data=na.omit(datansd1),start = list(intercept = 0)),error=function(e)99999) #Equation 4 in Bunnefeld et al. 2011
#where c is a constant
#summary(m4)        #this will print a summary of the converged model

###########################
##
## NOMADIC
##
###########################

m5 <- tryCatch(nls(nsdall ~ beta*zerostart,start=c(beta=1), data=na.omit(datansd1)),error=function(e)99999) #Equation 5 in Bunnefeld et al. 2011
#where beta is a constant and t the number of days since initial start date (i.e., 1 June of each year)
#summary(m5)        #this will print a summary of the converged model

#Below we are going to set up the AIC table 
ID[i] <- paste(unique(as.factor(datansd$nsdid)))
LOCS[i] <- nrow(coords)
MIGR[i] <- print(tryCatch(AIC(m1),error=function(e)0))
MIXM[i] <- print(tryCatch(AIC(m2),error=function(e)0))
DISP[i] <- print(tryCatch(AIC(m3),error=function(e)0))
HORA[i] <- print(tryCatch(AIC(m4),error=function(e)0))
NOMA[i] <- print(tryCatch(AIC(m5),error=function(e)0))


AICC_1[i] <- print(tryCatch(AIC(m1),error=function(e)99999))
AICC_2[i] <- print(tryCatch(AIC(m2),error=function(e)99999))
AICC_3[i] <- print(tryCatch(AIC(m3),error=function(e)99999))
AICC_4[i] <- print(tryCatch(AIC(m4),error=function(e)99999))
AICC_5[i] <- print(tryCatch(AIC(m5),error=function(e)99999))

minAIC[i] <- min(AICC_1[i],AICC_2[i],AICC_3[i],AICC_4[i],AICC_5[i])

d_AICC_1[i] <- (AICC_1[i] - minAIC[i])
d_AICC_2[i] <- (AICC_2[i] - minAIC[i])
d_AICC_3[i] <- (AICC_3[i] - minAIC[i])
d_AICC_4[i] <- (AICC_4[i] - minAIC[i])
d_AICC_5[i] <- (AICC_5[i] - minAIC[i])

LL_AICC_1[i] <- exp(-0.5*d_AICC_1[i])
LL_AICC_2[i] <- exp(-0.5*d_AICC_2[i])
LL_AICC_3[i] <- exp(-0.5*d_AICC_3[i])
LL_AICC_4[i] <- exp(-0.5*d_AICC_4[i])
LL_AICC_5[i] <- exp(-0.5*d_AICC_5[i])

sumLL_AICC[i] <- sum(LL_AICC_1[i],LL_AICC_2[i],LL_AICC_3[i],LL_AICC_4[i],LL_AICC_5[i])

wi_AICC_1[i] <- LL_AICC_1[i]/sumLL_AICC[i]
wi_AICC_2[i] <- LL_AICC_2[i]/sumLL_AICC[i]
wi_AICC_3[i] <- LL_AICC_3[i]/sumLL_AICC[i]
wi_AICC_4[i] <- LL_AICC_4[i]/sumLL_AICC[i]
wi_AICC_5[i] <- LL_AICC_5[i]/sumLL_AICC[i]

filename<-paste(substr(List[i,],1,8),"png",sep=".")
#NOTE:Numbers after "List[i,] need to encompass possible lengths of output name (i.e., D19.txt is 6 characters)
png(filename,height=20,width=30,units="cm",res=600)
#graphical exploration of the data will help you find sensible starting values 
#for each of the parameters asym, xmidA, xmidB, scale1 and scale2. 
#to graph nsd against time, use:
#xyplot(nsdall~zerostart|nsdtrip,data=datansd)
#str(nsdtest)
#now plot the data with the predicted curve  
nsdplot <- xyplot(nsdall ~ zerostart/3600, data=datansd1,
col="grey",    #color for the observed locations
type='b',      # 'b' shows the locations as dots, with a line connecting 
#successive locations. Can also be 'p' for just the locations, or 'l' for just 
#the line between locations
ylab=expression(paste('Net squared displacement ',' ', (km^2))), #y axis label
xlab="Hours after trip start")

plot(nsdplot)

dev.off()

}

#Create table of AIC values with lower AIC identifying best model
RESULT<-cbind(ID,LOCS,MIGR,MIXM,DISP,HORA,NOMA)
colnames(RESULT)<- c("ID","LOCS","MIGR","MIXM","DISP","HORA","NOMA")#
RESULT
#write.table(RESULT,"OUT_AUC_dBBMM.txt",sep="\t")

#Create table of raw values to calculate AICweights
Migratory <- rbind(ID,AICC_1,d_AICC_1,LL_AICC_1,wi_AICC_1)
MixedMig <- rbind(ID,AICC_2,d_AICC_2,LL_AICC_2,wi_AICC_2)
Disperser <- rbind(ID,AICC_3,d_AICC_3,LL_AICC_3,wi_AICC_3)
HomeRange <- rbind(ID,AICC_4,d_AICC_4,LL_AICC_4,wi_AICC_5)
Nomadic <- rbind(ID,AICC_5,d_AICC_5,LL_AICC_5,wi_AICC_5)

RESULT2 <- rbind(Migratory,MixedMig,Disperser,HomeRange,Nomadic)
RESULT2

date()