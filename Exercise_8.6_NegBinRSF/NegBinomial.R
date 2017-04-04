library(rgdal)
library(plyr)
library(raster)
library(rgeos)
library(MASS)#For nb models
#library(pscl)#For zero-inflated nb models

#rm(list=ls())

## 1.1.1. Reads and prepares the data 
muleys<-read.csv("FinalMD.csv", header=T, sep=",")
str(muleys)
muleys$ID <- paste(substr(muleys$fname,3,9))
muleys$ID <- as.factor(muleys$ID)

##Remove extra columns of data
newdata <- muleys[c(-1:-3,-13:-21)]
muleys<-newdata

summary(muleys$GPS.UTM.Easting)

muleys <- subset(muleys, muleys$GPS.UTM.Easting > 670000 & muleys$GPS.UTM.Easting != "NA")
str(muleys)

muleys$NewDate<-as.POSIXct(muleys$GPS.Fix.Time, format="%Y.%m.%d %H:%M:%S", origin="1970-01-01")

#Dates deer GPS collars were on the air
rangefunct <- function(muleys){
	dates=range(muleys$NewDate)
	print(dates)
	days=diff(range(muleys$NewDate))
	print(days)		
	}
collars <- dlply(muleys, .(ID), rangefunct)

##Sort Data
muleys <- muleys[order(muleys$ID, muleys$NewDate),]

##TIME DIFF NECESSARY IN BBMM CODE
timediff <- diff(muleys$NewDate)*60
## remove first entry without any difference 
muleys <- muleys[-1,] 
muleys$timelag <-as.numeric(abs(timediff))
##Remove locations greater than 5.5 hours apart in time
muleys <- subset(muleys, muleys$timelag < 19800)
summary(muleys$timelag)

##Make a spatial data frame of locations after removing outliers
muleysSPDF<-data.frame(x = muleys$GPS.UTM.Easting, y = muleys$GPS.UTM.Northing)
crs<-"+proj=utm +zone=12 +datum=WGS84"
head(muleysSPDF)
plot(muleysSPDF, axes=T)#To visualize all locations

str(muleysSPDF)
utm.spdf <- SpatialPointsDataFrame(coords = muleysSPDF, data = muleys, proj4string = CRS(crs))

##change muleysSPDF from UTM to Albers
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
muleys.spdf <-spTransform(utm.spdf, CRS=Albers.crs)

#Make data of albers x and y and replace with original muleys dataframe
muleys <- as.data.frame(muleys.spdf)
 
##We can create and clip with a rectangular grid
bbox(muleys.spdf)
bb1 <- cbind(x=c(-1138756,-1138756,-1109117,-1109117, -1138756), y=c(1692216, 1729463,1729463,1692216,1692216))
muleysAlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(muleys.spdf)))

##Let's buffer around the bounding box to be sure it encompasses all locations
muleysbuffSP <- gBuffer(muleysAlbersSP,width=2000)
plot(muleysbuffSP)
points(muleys.spdf,col="red", par(new=TRUE))
str(muleys.spdf@data)
summary(muleys.spdf@data$ID)

#######################################################################
##Subset locations by year for season-specific RSFs

range(muleys$NewDate)
#[1] "2011-10-13 03:00:34 EDT" "2013-09-14 21:00:52 EDT"

muleys$Date <- as.Date(muleys$GPS.Fix.Time, "%Y.%m.%d")
str(muleys)

winter2012 <- subset(muleys, Date > "2011-10-01" & Date < "2012-06-01")
str(winter2012)

win12.xy<-data.frame(x = winter2012$x, y = winter2012$y)
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
win12spdf <- SpatialPointsDataFrame(coords = win12.xy, data = winter2012, proj4string = CRS(albers))
plot(win12spdf, axes=T)
str(win12spdf)

################################################################################################## 
#Begin covariate preparate 

##Load crops raster layers
crop11 <- raster("crop11clip")
proj4string(crop11)
class(crop11)
plot(crop11)
points(win12spdf, col="red", cex=0.05)

##Clip out buff in crop11 layer
mulcrop11clip<-crop(crop11, muleysbuffSP)
as.matrix(table(values(mulcrop11clip)))

############################################################################################
#Reclassify clipped layers from above into 8 groups for each 2011, 2012,2013
#Crop layers. We then will create distance to rasters for each of 8 habitats:

#1 = Sunflower
#2 = Winter crops
#3 = Alfalfa
#4 = Summer crops
#5 = Random crops
#6 = Forest
#7 = Grassland
#8 = Shrubland
# 
mulcrop11clip
m11 <- c(-Inf,0,NA, 5.5, 6.5, 1, 22.5, 24.5, 2, 26.5, 27.5, 2, 29.5, 30.5, 2, 35.5, 36.5, 3, 3.5, 4.5, 4, .5, 1.5, 4, 11.5, 12.5, 4, 27.5, 28.5, 4, 31.5, 33.5, 4,  42.5, 43.5, 5, 47.5, 49.5, 5, 58.5, 59.5, 5, 60.5, 61.5, 5, 65.5, 69.5, 5, 76.5, 77.5, 5, 110.5, 111.5, 5, 120.5, 124.5, 5, 130.5, 131.5, 5, 189.5, 190.5, 5, 194.5, 195.5, 5, 228.5, 229.5, 5, 140.5, 143.5, 6, 170.5, 171.5, 7, 180.5, 181.5, 7, 36.5, 37.5, 7, 151.5, 152.5, 8, 41.5, 42.5, 4, 204.5, 205.5, 4, 230,Inf,NA)
rclmat11 <- matrix(m11, ncol=3, byrow=TRUE)
crop11rc <- reclassify(mulcrop11clip, rclmat11)
crop11rc
as.matrix(table(values(crop11rc)))#check on 8 categories
plot(crop11rc)

##########################################################################
##########################################################################
##NOTE: Code in this box was simply for demonstration purposes to reduce overall
## time for processing during class. Skip this section of code if using your
## own data and your computer has the appropriate processing capabilities.
#First we will clip out the raster layers by zooming into only a few locations
plot(crop11rc)
points(win12spdf, col="red", cex=0.05)
#Code below is used to just zoom in on grid using raster layer
e <- drawExtent()
#click on top left of crop box and bottom right of crop box create zoom
classclip <- crop(crop11rc,e)
plot(classclip)
points(win12spdf, col="red", cex=0.05)

#Now clip out locations within clipped raster
e2 <- drawExtent()
#click on top left of crop box and bottom right of crop box create zoom
win12clip <- crop(win12spdf,e2)
plot(classclip)
points(win12clip, col="red", cex=0.05)
table(win12clip$ID)
win12clip$ID <- droplevels(win12clip$ID)

crop11rc <- classclip #Code here not needed if this box of code is not run
muleys.spdf  <- win12clip#Code here not needed if this box of code is not run
#############################################################################################
#############################################################################################
#Create layers for CROP 2011
#Create Distance to Cover with cover being only forest
crop11df <- as.data.frame(as(crop11rc, "SpatialGridDataFrame"))
head(crop11df)
names(crop11df) <- c("crop", "x", "y")
head(crop11df)

#Create individual layers for each habitat category
sunflower11 <- subset(crop11df,crop11df$crop=="1")
head(sunflower11)
sun_only11 <- sunflower11[c(-1)]#need to remove first column 
d_sun11 <- distanceFromPoints(crop11rc,sun_only11)
plot(d_sun11)

wincrops11 <- subset(crop11df,crop11df$crop=="2")
head(wincrops11)
wintercrop_only11 <- wincrops11[c(-1)]#need to remove first column 
d_wincrop11 <- distanceFromPoints(crop11rc,wintercrop_only11)
plot(d_wincrop11)

alfalfa11 <- subset(crop11df,crop11df$crop=="3")
head(alfalfa11)
alf_only11 <- alfalfa11[c(-1)]#need to remove first column 
d_alf11 <- distanceFromPoints(crop11rc,alf_only11)
plot(d_alf11)

summercrop11 <- subset(crop11df,crop11df$crop=="4")
head(summercrop11)
summercrop_only11 <- summercrop11[c(-1)]#need to remove first column 
d_sumcrop11 <- distanceFromPoints(crop11rc,summercrop_only11)
plot(d_sumcrop11)

randomcrop11 <- subset(crop11df,crop11df$crop=="5")
head(randomcrop11)
randomcrop_only11 <- randomcrop11[c(-1)]#need to remove first column 
d_rancrop11 <- distanceFromPoints(crop11rc,randomcrop_only11)
plot(d_rancrop11)

forest11 <- subset(crop11df,crop11df$crop=="6")
head(forest11)
forest_only11 <- forest11[c(-1)]#need to remove first column 
d_forest11 <- distanceFromPoints(crop11rc,forest_only11)
plot(d_forest11)

grass11 <- subset(crop11df,crop11df$crop=="7")
head(grass11)
grass_only11 <- grass11[c(-1)]#need to remove first column 
d_grass11 <- distanceFromPoints(crop11rc,grass_only11)
plot(d_grass11)

shrub11 <- subset(crop11df,crop11df$crop=="8")
head(shrub11)
shrub_only11 <- shrub11[c(-1)]#need to remove first column 
d_shrub11 <- distanceFromPoints(crop11rc,shrub_only11)
plot(d_shrub11)

##############################################################################
#############################################################################
#Distance to Roads Raster start

##Bring in roads layer
roads<-readOGR(dsn=".",layer="AlbersRoads")
proj4string(roads)
plot(roads, pch=16)
plot(muleysbuffSP, col="green", add=TRUE)
points(muleys.spdf, col="red")

#Clip using muleysbuffSP
roadclip <- crop(roads, crop11rc)
#cliproads <- gIntersection(roads, muleysbuffSP, byid=TRUE)

plot(roadclip)
points(muleys.spdf, col="red")

#Rasterize function to create raster of road shapefile with crop data as mask
roadrast <- rasterize(roadclip,crop11rc, mask=TRUE)
image(roadrast)

#Now make a dataframe of all raster cell locations for each road segment
roadrastdf <- as.data.frame(as(roadrast, "SpatialGridDataFrame"))
str(roadrastdf)
#Remove the non-xy column of the data frame
roadrastdf <- roadrastdf[c(-1)]
#Use raster package to create distance from each raster cell to each 
#road layer raster cell.
d_roadrast <- distanceFromPoints(crop11rc, roadrastdf)
plot(d_roadrast)

#################################################################################################
##END Distance to Roads Code
###################################################################################################
#Create Distance to Cover with cover being only forest and same for all seasons
d_cover <- d_forest11
plot(d_cover)

##MERGE ALL RASTER LAYERS INTO A RASTER STACK IF USING EXTRACT FUNCTION
r11 <- stack(list(d_sun11=d_sun11,d_wincrop11=d_wincrop11,d_alf11=d_alf11, d_sumcrop11=d_sumcrop11, d_rancrop11=d_rancrop11, d_forest11=d_forest11, d_grass11=d_grass11, d_shrub11=d_shrub11,d_roads=d_roadrast))
plot(r11)

##OR MERGE ALL RASTER LAYER DATAFRAMES TO COMBINE INTO ONE DATASET IF USING GRAB VALUES

#Distance to habitat 2011
d_sun11df <- as.data.frame(as(d_sun11, "SpatialGridDataFrame"))
d_wincrop11df <- as.data.frame(as(d_wincrop11, "SpatialGridDataFrame"))
d_alf11df <- as.data.frame(as(d_alf11, "SpatialGridDataFrame"))
d_sumcrop11df <- as.data.frame(as(d_sumcrop11, "SpatialGridDataFrame"))
d_rancrop11df <- as.data.frame(as(d_rancrop11, "SpatialGridDataFrame"))
d_forest11df <- as.data.frame(as(d_forest11, "SpatialGridDataFrame"))
d_grass11df <- as.data.frame(as(d_grass11, "SpatialGridDataFrame"))
d_shrub11df <- as.data.frame(as(d_shrub11, "SpatialGridDataFrame"))

#Distance to cover
d_covdf <- as.data.frame(as(d_cover, "SpatialGridDataFrame"))
#Distance to roads
final_roaddf <- as.data.frame(as(d_roadrast, "SpatialGridDataFrame"))

#Combine dataframes for all layers for each year
layers1 = cbind(crop11df,d_sun11df,d_wincrop11df,d_alf11df,d_sumcrop11df,d_rancrop11df,d_forest11df,d_grass11df,d_shrub11df,d_covdf,final_roaddf)
head(layers1)
layers1 = layers1[,-c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)]
names(layers1) = c("crop","d_sunflower","d_wintercrop","d_alfalfa","d_summercrop","d_randomcrop","d_forest","d_grass","d_shrub","d_cover","d_roads","x", "y")
head(layers1)
#-----------------------------------------------------------------------------------
##############################################################################
######BEGIN CODE TO RUN NEGATIVE BINOMIAL REGRESSION FOR RSFs
#####
bbox(crop11rc)
#min      max
#s1 -1127205 -1118085
#s2  1712055  1720725


#628 m +- 262 (SD) mean move distance

## create vectors of the x and y points 
x <- seq(from = xmin(crop11rc), to = xmax(crop11rc), by = 1256)#628 m daily move distance times 2 
y <- seq(from = ymin(crop11rc), to = ymax(crop11rc), by = 1256) 

## create a grid of all pairs of coordinates (as a data.frame) 
#Also the restart point for later!!
xy <- expand.grid(x = x, y = y)
class(xy)
str(xy)

#Identifiy projection before creating Spatial Points Data Frame
crs2<-"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
grid.pts<-SpatialPointsDataFrame(coords= xy, data=xy, proj4string = CRS(crs2))
proj4string(grid.pts)
plot(grid.pts)
gridded(grid.pts)
class(grid.pts)
str(grid.pts)
griddata <- as.data.frame(grid.pts)
str(griddata)

#Create buffers around grid points
mbuffer <- gBuffer(grid.pts,width=628,byid=TRUE)
mbuff.spdf <- SpatialPolygonsDataFrame(mbuffer, data=data.frame(id=row.names(mbuffer), row.names=row.names(mbuffer)))
class(mbuff.spdf)
plot(grid.pts)
plot(mbuff.spdf,add=T)

plot(crop11rc)
plot(mbuff.spdf,add=T)
points(muleys.spdf,col=muleys.spdf@data$ID)

#Then we can get a mean for each covariate in each sample circle
#Extracts all rasters within each sample circle
date()
buff_win12 <- extract(r11, mbuff.spdf, weights=TRUE, fun = mean)
date()
head(buff_win12) 

#Now we need to convert to a data frame for nb modeling
# Read in animal_locations.txt or convert to data frame from above
#locations = read.table("deer_locations.txt", sep='\t', header=T)
locations.df = as.data.frame(muleys.spdf)
str(locations.df)
locations <- locations.df[-c(1:9,11:12,15:17)]

habitat_units_buffwin12 <- as.data.frame(buff_win12)
str(habitat_units_buffwin12)

#Add xy columns of circle centroids
mbuff.xy <- as.data.frame(grid.pts)
str(mbuff.xy)
habitat_units_buffwin12$x <- mbuff.xy$x
habitat_units_buffwin12$y <- mbuff.xy$y

# Source code file containing functions used below.
source("count_locations.R")

# Summary stats
names(locations)
head(locations)
tail(locations)
str(locations)
table(locations$ID)
plot(locations$x, locations$y, type='p')

# Summary stats
names(habitat_units_buffwin12)
head(habitat_units_buffwin12)
tail(habitat_units_buffwin12)
summary(habitat_units_buffwin12)

plot(habitat_units_buffwin12$x, habitat_units_buffwin12$y, type='p', cex=1.5)
points(locations$x, locations$y, col="red", cex=0.5, pch=19)
plot(mbuff.spdf, add=T)
#Test number of locations falling within actual buffers
test_buff <- crop(muleys.spdf, mbuff.spdf)
str(test_buff)
#points(test_buff,col="red")
#NB645572 <- subset(NB,NB$animal.ID=="647573A")
#plot(NB645572,NB645572$x,NB645572$y,col="red")

habitat_units_buffwin12[,1:9]=scale(habitat_units_buffwin12[,1:9],scale=TRUE)
#str(habitat_units_buffwin12)

# Calculate number of animal locations in each sampled habitat unit.
# (see "count_locations.R")
pooled.locations = locations
pooled.locations$ID = 1
NB = F.count.relocations(locations.df = pooled.locations, 
                    habitat.units.df = habitat_units_buffwin12, 
                    habitat.unit.size = 1256)
str(NB)
table(NB$animal.ID,NB$n.locations)

#habitat.units should match nrows of NB

nb = glm.nb(n.locations ~ offset(log(total)) + d_alf11 + d_forest11 + d_shrub11 + d_roads, data=NB)
summary(nb)

#----------------------------------------------------------------------------
# # Proportion of 0 counts in data
 sum(NB$n.locations == 0)/nrow(NB)
# 
 nb.density = structure(
 function # Probability mass function for NB2
# 
# # Description: This function gives the probability that a 
# # discrete random variable, X, is exactly equal to some value
# # according to a NB2 distribution.
# # Returns: Pr(X=k)
# 
 (k, 
# ### value at which to estimate probability
# ### Pr(X=k)
 mu, 
# ### NB2 estimate of mu
 theta
# ### NB2 estimate of theta
 ){
 
 	(gamma(theta+k)/(gamma(theta)*factorial(k)))*
 		(mu^k)*(theta^theta)/
 		((mu+theta)^(theta+k))
 		
 })
# # Expected proportion under NB2 model
 nb.density(k=0, mu=mean(NB$n.locations), theta=0.1136) 
# (Note: 0.1136 comes from the estimated theta of the model output (i.e., summary of nb))
# The value above can be interpreted as: "A NB2 distribution with theta=0.1136 and
# mu=136.8 should have an average of 45-% zero values"

editlayer4 <- layers1
names(editlayer4) = c("crop", "std_cover", "std_roads","x", "y")
head(editlayer4)

#Need to standardize the raw distance rasters first
editlayer4$std_sun <- editlayer4$d_sunflower
editlayer4$std_winter <- editlayer4$d_wintercrop
editlayer4$std_alf <- editlayer4$d_alfalfa
editlayer4$std_summer <- editlayer4$d_summercrop
editlayer4$std_rancrop <- editlayer4$d_randomcrop
editlayer4$std_forest <- editlayer4$d_forest
editlayer4$std_grass <- editlayer4$d_grass
editlayer4$std_shrub <- editlayer4$d_shrub
editlayer4$std_cover <- editlayer4$d_cover
editlayer4$std_roads <- editlayer4$d_roads

editlayer4[,14:23]=scale(editlayer4[,14:23],scale=TRUE)

editlayer4$cropnew <- 0
editlayer4$cropnew[editlayer4$crop == "1"] <- 1
editlayer4$cropnew[editlayer4$crop == "2"] <- 1
editlayer4$cropnew[editlayer4$crop == "3"] <- 1
editlayer4$cropnew[editlayer4$crop == "4"] <- 2
editlayer4$cropnew[editlayer4$crop == "5"] <- 3
editlayer4$cropnew <- as.factor(editlayer4$cropnew)


str(editlayer4)
editlayer4$crop <- as.numeric(editlayer4$crop)

#-----------------------------------------------------------------------------------
# predictions based on full model 
predictions = predict(nb, newdata=editlayer4, re.form=NA, type="response")# based on the scale of the linear predictors
predictions = exp(predictions)

range(predictions)

#-----------------------------------------------------------------------------------
# create Ascii grid of raw predictions
editlayer4$predictions = predictions
preds = editlayer4

preds = SpatialPixelsDataFrame(points=preds[c("x", "y")], data=preds)
preds = as(preds, "SpatialGridDataFrame")
names(preds)

#writeAsciiGrid(preds, "predictions.asc", attr=9) # attr should be column number for 'predictions'

#-----------------------------------------------------------------------------------
# assign each cell or habitat unit to a 'prediction class'.
# classes have (nearly) equal area, if the cells or habitat units have equal areas.
# output is a vector of class assignments (higher is better).
F.prediction.classes <- function(raw.prediction, n.classes){
    # raw.prediction = vector of raw (or scaled) RSF predictions
    # n.classes = number of prediction classes.
    pred.quantiles = quantile(raw.prediction, probs=seq(1/n.classes, 1-1/n.classes, by=1/n.classes))
    ans = rep(n.classes, length(raw.prediction))
    for(i in (n.classes-1):1){
        ans[raw.prediction < pred.quantiles[i]] = i
    }
    return(ans)
}

editlayer4$prediction.class = F.prediction.classes(editlayer4$predictions, 5)
table(editlayer4$prediction.class)
str(editlayer4)
##############################################
# create map of RSF prediction classes in R
m = SpatialPixelsDataFrame(points = editlayer4[c("x", "y")], data=editlayer4)
names(m)
par(mar=c(0,0,0,0))
image(m, attr=10, col=c("grey90", "grey70", "grey50", "grey30", "grey10")) # attr should be column number for 'predictions'
par(lend=1)
legend("bottomright", col=rev(c("grey90", "grey70", "grey50", "grey30", "grey10")),
	legend=c("High", "Medium-high", "Medium", "Medium-low", "Low"),
	title="Prediction Class", pch=15, cex=1.5,bty != "n", bg="white")
points(summer2013, col="red", cex=0.5)

# create Ascii grid of prediction classes
#m = as(m, "SpatialGridDataFrame")
#names(m)
#writeAsciiGrid(m, "prediction_classes.asc", attr=10) # attr = the column number of your binned prediction classes (i.e., prediction.class)


#-----------------------------------------------------------------------------------
# save the layers data frame for use later
#save(layers, file="layers")