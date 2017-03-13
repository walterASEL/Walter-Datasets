#install.packages(c("adehabitat","adehabitatHR","maptools","raster", "rgdal"))

library(adehabitat)
library(adehabitatHR)
library(raster)
library(rgdal)
library(maptools)

#Import raster as text using the "raster" package 
r <-raster("polyascii2.txt")
proj4string(r)#Note: No projection information was imported with the raster
plot(r)
#Assign projection information for imported text file
proj4string(r) <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
r
proj4string(r)

##Import raster as an Ascii files (factor) using "adehabitat" package
##for file1, file2, and file3 in the 3 sections of code below:
## Path of the file to be imported
file1 <-  paste("polyextract2.asc", sep = "\\")
levelfile <- paste("TableExtract.txt", sep = "\\")
asp <- import.asc(file1, lev = levelfile, type = "factor")
image(asp)
asp
str(asp)

#Now let's look at the vegetation categories of the file
ta <- table(as.vector(asp))
names(ta) <- levels(asp)[as.numeric(names(ta))]
ta

file2 <-  paste("polyclip.asc", sep = "\\")
levelfile2 <- paste("TableClip.txt", sep = "\\")
asp2 <- import.asc(file2, lev = levelfile2, type = "factor")
image(asp2)
asp2
str(asp2)
#Shows 7 vegetation categories

#Now let's look at the vegetation categories of the file
ta2 <- table(as.vector(asp2))
names(ta2) <- levels(asp2)[as.numeric(names(ta2))]
ta2
 
#NOTE: R won't recognize double digit veg categories
## Import raster as an Ascii files (factor) using:
## Path of the file to be imported
file3 <-  paste("polyascii2.asc")
levelfile3 <- paste("TableCode.txt")
asp3 <- import.asc(file3, lev = levelfile3, type = "factor")
image(asp3)
asp3
str(asp3)

#Now let's look at the vegetation categories of the file
ta3 <- table(as.vector(asp3))
names(ta3) <- levels(asp3)[as.numeric(names(ta3))]
ta3

#Or we can also use the "adehabitat" package to import DEM
fileElev <-  paste("demascii.asc")
elev <- import.asc(fileElev)
image(elev)
elev

windows()
#We can also use the "rgdal" package to import an ascii grid as a Spatial Grid Data Frame
habitat <- readGDAL("polyascii2.asc")
proj4string(habitat)
proj4string(habitat) <-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
image(habitat)
str(habitat)

#CRS of shapefile layers
crs <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#CRS of raster layers
crs2 <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#Load county shapefile
county<-readOGR(dsn=".",layer="BeaufortCoAlbers")
proj4string(county)
#Now let's make county a SpatialPolygon class to simply data contained within it
polys <- as(county, "SpatialPolygons")
plot(polys,add=T,lwd=2)
polys
text(coordinates(polys), labels="Beaufort")
proj4string(polys)

#Load airport runway shapefile
run<-readOGR(dsn=".",layer="RunwayAlbers")
proj4string(run)
polys2 <- as(run, "SpatialPolygons")
plot(polys2,add=T,lwd=2)
polys2
text(coordinates(polys2), labels="Runway")
proj4string(polys2)

#Load aircraft flight pattern shapefile
path<-readOGR(dsn=".",layer="FlightImage")
proj4string(path)
polys3 <- as(path, "SpatialLines")
plot(polys3,add=T,lty="32", col="blue")
polys3
proj4string(polys3)

#Load aircraft flight pattern shapefile
road<-readOGR(dsn=".",layer="CountyRoadAlbers")
proj4string(road)
polys4 <- as(road, "SpatialLines")
plot(polys4,add=T,lty="22", col="green")
polys4
proj4string(polys4)

plot(county)
plot(road, add=T)
plot(run, col="red",add=T)
plot(path, col="blue",add=T)

#Clip using the raster imported with "raster" package 
clip <- crop(r, polys)
proj4string(clip)
plot(clip)
plot(polys,add=T,lwd=2)
plot(polys2,add=T,lwd=2, col="red")
plot(polys3,add=T,lty="62", col="blue")
plot(polys4,add=T,lty="22", col="yellow")

#Reclassify veg layer by importing it as a raster not a SpatialGridDataFrame 
veg <-raster("polydouble.txt")
plot(veg)
veg

# reclassify the values into 7 groups
# all values between 0 and 20 equal 1, etc.
m <- c(0, 19, 1, 20, 39, 2, 40, 50, 3, 51,68, 4, 69,79, 5, 80, 88, 6, 89, 99, 7)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(veg, rclmat)
plot(rc)
rc

#Now, let's remove water that is coded 11 and No Data that is coded as 127
m1 <- c(0, 19, NA, 20, 39, 1, 40, 50, 2, 51,68, 3, 69,79, 4, 80, 88, 5, 89, 99, 6, 100, 150, NA)
rclmat1 <- matrix(m1, ncol=3, byrow=TRUE)
rc1 <- reclassify(veg, rclmat1)
plot(rc1)
rc1

#Import bird 49 locations to R
bv49 <-read.csv("Bird49.csv", header=T)
str(bv49)#How many bird locations?
#Make a spatial points data frame of locations and convert to Albers
coords<-data.frame(x = bv49$x, y = bv49$y)
crs<-"+proj=utm +zone=17N +ellps=WGS84"
head(coords)

bvspdf <- SpatialPointsDataFrame(coords= coords, data = bv49, proj4string = CRS(crs))
str(bvspdf)
bvspdf[1:5,]
points(bvspdf, col="red")
#NOTE: Locations must be assigned to the UTM coordinate system prior to projection 
#to Albers so won't overlay on veg layer at this point because veg is in Albers
bv49Albers <-spTransform(bvspdf, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
class(bv49Albers)
proj4string(bv49Albers)
bv49Albers[1:5,]
points(bv49Albers, col="red")
#Determine which of those points lie within a cell that contains data by using the extract function. The extract function will
#extract covariate information from the raster at a particular point.
veg.survey<-extract(veg, bv49Albers)
veg.survey
veg.survey<-subset(bv49Albers,!is.na(veg.survey))
plot(veg.survey, col="black", add=T)

#We can also create a regular grid and do the same thing
Sample.points<-expand.grid(seq(veg@extent@xmin, veg@extent@xmax, by=1000), weight = seq(veg@extent@ymin, veg@extent@ymax, by=1000))
points(Sample.points, bg="red", cex=.5,col="red")

#Let's create some random points using the minimum and maximum coordinates of the raster to
#determine the range of points from which to select x and y
x.pts<-sample(seq(veg@extent@xmin, veg@extent@xmax, by=10),1000) ##generate
#x coordinates for random points
y.pts<-sample(seq(veg@extent@ymin, veg@extent@ymax, by=10),1000)
#Now create a spatial points file from the randomly generated points
coords2<-data.frame(x = x.pts, y = y.pts)
crs2<-"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
head(coords2)
points(coords2, bg="red", cex=.5,col="blue")

#Determine which of those points lie within a cell that contains data by using the extract function. The extract function will
#extract covariate information from the raster at a particular point.
str(coords2)
veg.sample<-extract(veg, coords2)
head(veg.sample)
veg.sample<-subset(coords2,!is.na(veg.sample))
str(veg.sample)
head(veg.sample)
points(veg.sample,col="red")

##CLIP VEGETATION AND DO SAME PROCESS
plot(clip)
#Determine which of those points lie within a cell that contains data by using the extract function. The extract function will
#extract covariate information from the raster at a particular point.
clip.survey<-extract(clip, bv49Albers)
clip.survey
clip.survey<-subset(bv49Albers,!is.na(clip.survey))
plot(clip.survey, col="black", add=T)

#We can also create a regular grid and do the same thing
Sample.points2<-expand.grid(seq(clip@extent@xmin, clip@extent@xmax, by=1500), weight = seq(clip@extent@ymin, clip@extent@ymax, by=1500))
points(Sample.points2, bg="red", cex=.5,col="red")

#Let's create some random points using the minimum and maximum coordinates of the raster to
#determine the range of points from which to select x and y
x.pts2<-sample(seq(clip@extent@xmin, clip@extent@xmax, by=10),500) ##generate
#x coordinates for random points
y.pts2<-sample(seq(clip@extent@ymin, clip@extent@ymax, by=10),500)
#Now create a spatial points file from the randomly generated points
coords3<-data.frame(x = x.pts2, y = y.pts2)
crs2<-"+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
points(coords3, bg="red", cex=.5,col="blue")

#Determine which of those points lie within a cell that contains data by using the extract function. The extract function will
#extract covariate information from the raster at a particular point.
str(coords3)
clip.sample<-extract(clip, coords3)
clip.sample
clip.sample<-subset(coords3,!is.na(clip.sample))
str(clip.sample)
points(clip.sample, cex=.5, col="red")
points(bv49Albers)





#NOTE: Code below this line must be done in updated version of R
# and on clipped data or a small dataset
#Reclassify
veg <-raster("polydouble.txt")
plot(veg)
proj4string(veg)
veg

county<-readOGR(dsn=".",layer="BeaufortCoAlbers")
proj4string(county)
spplot(county)
polys <- as(county, "SpatialPolygons")
plot(polys,add=T,lwd=2)


clip <- crop(veg, polys)

val <- getValues(clip)
xy <- as.data.frame(xyFromCell(clip,1:ncell(clip)))
xy <- cbind(xy,val)

library(ggplot2)
ggplot(na.omit(xy), aes(x=x, y=y, fill=val)) + geom_raster() + coord_equal()
p <- ggplot(na.omit(xy), aes(x=x, y=y, fill=factor(val))) + 
  geom_raster() + 
  coord_equal()

#To save layer if needed
try(ggsave(plot=p,coast,height=8,width=8))



