library(sp)
library(rgdal)
library(raster)
library(adehabitatMA)

muleys <-read.csv("muleysexample.csv", header=T)
summary(muleys$id)
str(muleys)

#Remove outlier locations
newmuleys <-subset(muleys, muleys$Long > -110.50 & muleys$Lat > 37.8 & muleys$Long < -107)
muleys <- newmuleys

#Make a spatial data frame of locations after removing outliers
coords<-data.frame(x = muleys$Long, y = muleys$Lat)
crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
head(coords)
plot(coords)

deer.spdf <- SpatialPointsDataFrame(coords= coords, data = muleys, proj4string = CRS(crs))
deer.spdf[1:5,]
proj4string(deer.spdf)

#Project the deer.spdf to Albers as in previous exercise
Albers.crs <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
deer.albers <-spTransform(deer.spdf, CRS=Albers.crs)
proj4string(deer.albers)
bbox(deer.albers)
#       min      max
#x -1127964 -1115562
#y  1718097  1724867

plot(deer.albers)
## create vectors of the x and y points 
x <- seq(from = -1127964, to = -1115562, by = 1500) 
y <- seq(from = 1718097, to = 1724867, by = 1500) 

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

gridded(grid.pts) <- TRUE
gridded(grid.pts)
str(grid.pts)

grid <- as(grid.pts, "SpatialPolygons") 
plot(grid)
str(grid)
class(grid)
summary(grid)
gridspdf <- SpatialPolygonsDataFrame(grid, data=data.frame(id=row.names(grid), row.names=row.names(grid))) 
names.grd<-sapply(gridspdf@polygons, function(x) slot(x,"ID"))
text(coordinates(gridspdf), labels=sapply(slot(gridspdf, "polygons"), function(i) slot(i, "ID")), cex=0.5)
points(deer.albers, col="red") 
str(gridspdf@polygons)

o = over(deer.albers,gridspdf)
head(o)
new = cbind(deer.albers@data, o)
head(new)

#Run the x and y over to EXPAND the grid extent (by adding 300 units to max y) to encompass all locations then re-run the code over from 
#xy through new2 again
x <- seq(from = -1127964, to = -1115562, by = 1500) 
y <- seq(from = 1718097, to = 1725867, by = 1500) 

##BE SURE TO RUN CODE FROM XY CREATION THROUGH NEW2 AGAIN THEN LOOK AT DATA!!

o2 = over(deer.albers,gridspdf)
head(o2)
new2 = cbind(deer.albers@data, o2)#No more NAs causing errors!
new2[1:15,]

#Create sampling grid using extracted raster within locations of mule deer 8

#Load vegetation raster layer textfile clipped in ArcMap 
veg <-raster("cropnlcd.tif")
plot(veg)
class(veg)

#Clip the raster within the extent of the newly created grid
bbclip <- crop(veg, gridspdf)
plot(bbclip)
points(deer.albers, col="red")
plot(gridspdf, add=T)

#Cell size of raster layer
xres(bbclip)

#Create histogram of vegetation categories in bbclip
hist(bbclip)

#Calculate cell size in square meters
ii <- calcperimeter(gridspdf)#requires adehabitatMA package
as.data.frame(ii[1:5,])

#Extract
table = extract(bbclip,gridspdf)
str(table[1])

area = extract(bbclip,gridspdf)
combine=lapply(area,table)
combine[[1]]#Shows vegetation categories and numbers of cells in grid #1
combine[[27]]

