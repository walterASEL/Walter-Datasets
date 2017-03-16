library(sp)
library(lattice)
library(rgdal)#readOGR
library(rgeos)#gIntersection
library(raster)#to use "raster" function
library(adehabitatHR)
library(maptools)#readAsciiGrid

muleys <-read.csv("muleysexample.csv", header=T)
str(muleys)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.80)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys

#Only use the line below for example exercise
#muleys <- muleys[sample(nrow(muleys), 100),]

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$X, y = muleys$Y)
utm.crs<-"+proj=utm +zone=12 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(utm.crs))

Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
deer.spdf <-spTransform(utm.spdf, CRS=Albers.crs)

#We can create and clip with a rectangular grid
bbox(deer.spdf)
bb1 <- cbind(x=c(-1127964,-1127964,-1115562,-1115562,-1127964), y=c(1718097,1724868,1724868,1718097,1718097))
AlbersSP <- SpatialPolygons(list(Polygons(list(Polygon(bb1)),"1")), proj4string=CRS(proj4string(deer.spdf)))
plot(AlbersSP)
points(deer.spdf, col="red")

#Let's buffer around the bounding box to be sure it encompasses all locations
buffSP <- gBuffer(AlbersSP,width=1000)
plot(buffSP)
points(deer.spdf,col="red")

elevation <- raster("dem_albers.txt")
image(elevation, col=terrain.colors(10))
contour(elevation, add=TRUE)
slope = terrain(elevation,opt='slope', unit='degrees')
aspect = terrain(elevation,opt='aspect', unit='degrees')
elevation
slope
aspect
demclip <- crop(elevation, buffSP)
sloclip <- crop(slope, buffSP)
aspclip <- crop(aspect, buffSP)

#Load crops raster layer 
crops <-raster("crop12clip.tif")
plot(crops)
class(crops)
as.matrix(table(values(crops)))
proj4string(crops)
crops

#Clip using the raster imported with "raster" package 
bbclip <- crop(crops, demclip)

# Reclassify crops raster from above into 9 groups
# all values between 0 and 20 equal 1, etc.
m <- c(-Inf,0,NA,2, 7, 2, 20, 60, 3, 60, 70, 4, 110, 132, 5, 133, 150, 6, 151, 172, 7, 180, 183, 8, 189, 191, 9,192,205,10)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(bbclip, rclmat)
plot(rc)
rc
as.matrix(table(values(rc)))

plot(bbclip)
points(deer.spdf,col="red")
plot(buffSP, add=T)

proj4string(bbclip)
proj4string(deer.spdf)
proj4string(buffSP)

nlcd <- as.data.frame(as(rc, "SpatialGridDataFrame")) 
elev <- as.data.frame(as(demclip, "SpatialGridDataFrame")) 
slo <- as.data.frame(as(sloclip, "SpatialGridDataFrame"))
asp <- as.data.frame(as(aspclip, "SpatialGridDataFrame")) 
str(elev)
str(slo)
str(asp)
str(nlcd)

layers = cbind(nlcd, elev, asp, slo)
head(layers)
layers = layers[,-c(2,3,5,6,8,9)]
names(layers) = c("nlcd", "elevation", "aspect", "slope","x", "y")
head(layers)

# turn aspect into categorical
aspect_categorical = rep(NA, nrow(layers))
aspect_categorical[layers$aspect < 45 | layers$aspect >= 315] = "N"
aspect_categorical[layers$aspect >= 45 & layers$aspect < 135] = "E"
aspect_categorical[layers$aspect >= 135 & layers$aspect < 225] = "S"
aspect_categorical[layers$aspect >= 225 & layers$aspect < 315] = "W"
table(aspect_categorical)
table(is.na(aspect_categorical))
layers$aspect_categorical = aspect_categorical
head(layers)
write.table(layers,"layer1.txt",sep=",",col.names=TRUE, quote=FALSE)

layers2 <- layers #change to layers2 simply to avoid confusion with "layer" term in function below

#-----------------------------------------------------------------------------------
# grab values for hexagonal sample of points (taken above)
grab.values = function(layer, x, y){
	# layer is data.frame of spatial layer, with values 'x', 'y', and ____?
	# x is a vector 
	# y is a vector
	if(length(x) != length(y)) stop("x and y lengths differ")
	z = NULL
	for(i in 1:length(x)){
		dist = sqrt((layer$x - x[i])^2 + (layer$y-y[i])^2) 
		#Could adjust this line or add another line to calculate moving window or distance to nearest feature
		z = rbind(z, layer[dist == min(dist),][1,])
	}
	return(z)
}

#Grab all values from muleys for each layer in r
test = grab.values(layers2, muleys$X, muleys$Y)
head(test)
str(muleys)

##NOTE that all values are the same but this is not correct.
##What is the problem here and how do we fix it?

#Need to grab Albers XY not UTM as in muleys above
muleys <- as.data.frame(deer.spdf)
str(muleys)
# grab all values for used and available points based on combined layer data set
# can take 5+ minutes
used = grab.values(layers2, muleys$x, muleys$y)
used$x = muleys$x
used$y = muleys$y
used$animal_id = muleys$id
used$use = 1
head(used)

#Create MCP for all locations (3nd order selection)
cp = mcp(deer.spdf[,2],percent=100)
as.data.frame(cp)
## Plot the home ranges
plot(cp)
## ... And the relocations
plot(deer.spdf, col=as.data.frame(deer.spdf)[,2], add=TRUE)

#Determine the habitat available using all code below
#First create random sample of points in each polygon
random <- sapply(slot(cp, 'polygons'), function(i) spsample(i, n=50, type='random', offset=c(0,0)))
plot(cp) ; points(random[[1]], col='red', pch=3, cex=.5)#The number in double brackets changes polygons
# stack into a single SpatialPoints object
random.merged <- do.call('rbind', random)
#Extract the original IDs
ids <- sapply(slot(cp, 'polygons'), function(i) slot(i, 'ID'))
#Determine the number of ACTUAL sample points generated for each polygon
newpts <- sapply(random, function(i) nrow(i@coords))
newpts #Nice check of how many points generated per polygon 
# generate a reconstituted vector of point IDs
pt_id <- rep(ids, newpts)
 
# promote to SpatialPointsDataFrame
random.final <- SpatialPointsDataFrame(random.merged, data=data.frame(poly_id=pt_id))
 
#Plot random final on MCPs
plot(cp) ; points(random.final, col=random.final$poly_id, pch=3, cex=0.5)

random.final

# make 'random.final' a data.frame
random.df = as.data.frame(random.final)
names(random.df) = c("ID","x","y")
# can take 5+ minutes
available = grab.values(layers2, random.df$x, random.df$y)
available$x = random.df$x
available$y = random.df$y
available$animal_id = pt_id
available$use = 0
head(available)

data = rbind(available, used)
str(data)
#''data.frame':	200 obs. of  9 variables:
#100 locations used +
#100 locations available (2 animals X 50 random locations)
#= 100 #Confirmed in code below
str(data)
# nlcd              : num  7 6 7 7 7 7 7 7 7 6 ...
# elevation         : int  2058 2058 2068 2068 2070 2072 2076 2062 2071 2071 ...
# aspect            : num  105 278 105 80 135 ...
# slope             : num  2.72 3.37 4.68 4.11 6.05 ...
# x                 : num  -1127639 -1127610 -1127864 -1127805 -1127862 ...
# y                 : num  1724257 1724271 1724091 1724218 1724174 ...
# aspect_categorical: chr  "E" "W" "E" "E" ...
# animal_id         : chr  "D12" "D12" "D12" "D12" ...
# use               : num  0 0 0 0 0 0 0 0 0 0 ...

#The above code is for 3rd order selection within home range of each deer
#What if we wanted to look at 3rd order selection around buffered circels for each deer location or set up
#data for Discrete Choice?
settbuff=gBuffer(deer.spdf, width=500, byid=TRUE)

#Determine the habitat available using all code below
#First create random sample of points in each polygon
ranbuff <- sapply(slot(settbuff, 'polygons'), function(i) spsample(i, n=3, type='random', offset=c(0,0)))
plot(settbuff) ; points(ranbuff[[100]], col='red', pch=3, cex=.5)#The number in double brackets changes polygons
# stack into a single SpatialPoints object
ranbuff.merged <- do.call('rbind', ranbuff)

#Extract the original IDs
buff_ids <- sapply(slot(settbuff, 'polygons'), function(i) slot(i, 'ID'))
buff_ids <- paste(settbuff$id, buff_ids, sep="_")
#Determine the number of ACTUAL sample points generated for each polygon
buffpts <- sapply(ranbuff, function(i) nrow(i@coords))
buffpts[1:20] #Nice check of how many points generated per polygon 
# generate a reconstituted vector of point IDs
buffpt_id <- rep(buff_ids, buffpts)
 
# promote to SpatialPointsDataFrame
buff.final <- SpatialPointsDataFrame(ranbuff.merged, data=data.frame(poly_id=buffpt_id))
 
#Plot buff.final on buffered circles
plot(settbuff) ; points(buff.final, col=random.final$poly_id, pch=3, cex=0.5)
buff.final

# make 'buff.final' a data.frame
buffer.df = as.data.frame(buff.final)
names(buffer.df) = c("ID","x", "y")
head(buffer.df)
str(random.df)
str(buffer.df)

# can take 5+ minutes
buff_avail = grab.values(layers2, buffer.df$x, buffer.df$y)
buff_avail$x = buffer.df$x
buff_avail$y = buffer.df$y
buff_avail$animal_id = buffpt_id
buff_avail$use = 0
buff_avail[1:10,]

data2 = rbind(buff_avail, used)
str(data2)

#Save workspace so all analysis are available
save.image("RSF_dataprep.RData")

#Before closing, let's save the "used" and available data set to use in the next exercise
write.table(used, "MD_used.txt")
write.table(available, "MD_avail.txt")
