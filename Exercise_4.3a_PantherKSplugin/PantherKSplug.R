library(adehabitatHR)
library(ks)
library(rgdal)
library(raster)
library(PBSmapping)
library(sp)

panther<-read.csv("pantherjitter.csv",header=T)
str(panther)
panther$CatID <- as.factor(panther$CatID)

cat143 <- subset(panther, panther$CatID == "143")
##Get only the coordinates
loc <- data.frame("x"=cat143$X, "y"=cat143$Y)

##Define the projection of the coordinates
proj4string <- CRS("+proj=utm +zone=17N +ellps=WGS84")

##Make SpatialPointsDataFrame using the XY, attributes, and projection
spdf <- SpatialPointsDataFrame(loc, cat143, proj4string = proj4string)

#Calculate the bandwidth matrix to use later in creating the KDE
Hpi1 <- Hpi(x = loc)
Hpi1

##write out the bandwidth matrix to a file as you might want to refer to it later
#write.table(Hpi1, paste("hpivalue_", "143", ".txt", sep=""), row.names=FALSE,sep="\t")

##Create spatial points from just the xy’s
loc.pts <- SpatialPoints(loc, proj4string=proj4string)

##For home range calculations, ##some packages require evaluation points (ks) while others
##require grid as spatial pixels (adehabitatHR).

##Set the expansion value for the grid and get the bbox from the SpatialPointsDataFrame
expandValue <- 5000 #This value is the amount to add on each side of the bbox.
#Change to 5000 if error occurs at 99% ud
boundingVals <- spdf@bbox

##Get the change in x and y and adjust using expansion value
deltaX <- as.integer(((boundingVals[1,2]) - (boundingVals[1,1])) + (2*expandValue))
deltaY <- as.integer(((boundingVals[2,2]) - (boundingVals[2,1])) + (2*expandValue))

##100 meter grid for data in this exercise
gridRes <- 100
gridSizeX <- deltaX / gridRes
gridSizeY <- deltaY / gridRes
##Offset the bounding coordinates to account for the additional area
boundingVals[2,1] <- boundingVals[2,1] - expandValue
boundingVals[2,2] <- boundingVals[2,2] + expandValue
boundingVals[1,1] <- boundingVals[1,1] - expandValue
boundingVals[1,2] <- boundingVals[1,2] + expandValue

##Grid Topology object is basis for sampling grid (offset, cellsize, dim)
gridTopo <- GridTopology((boundingVals[,1]), c(gridRes,gridRes),c(gridSizeX,gridSizeY))

##Using the Grid Topology and projection create a SpatialGrid class
sampGrid <- SpatialGrid(gridTopo, proj4string = proj4string)

##Cast over to Spatial Pixels
sampSP <- as(sampGrid, "SpatialPixels")

##convert the SpatialGrid class to a raster
sampRaster <- raster(sampGrid)

##set all the raster values to 1 such as to make a data mask
sampRaster[] <- 1

##Get the center points of the mask raster with values set to 1
evalPoints <- xyFromCell(sampRaster, 1:ncell(sampRaster))

##Here we can see how grid has a buffer around the locations and trajectory. This will ensure that we
#project our home range estimates into a slightly larger extent than the original points extent (bbox) alone.
plot(sampRaster)
points(loc, pch=1, cex=0.5)

##Create the KDE using the evaluation points
hpikde <- kde(x=loc, H=Hpi1, eval.points=evalPoints)

##Create a template raster based upon the mask and then assign the values from the kde to the template
hpikde.raster <- raster(sampRaster)

hpikde.raster <- setValues(hpikde.raster,hpikde$estimate)

##Lets take this raster and put it back into an adehabitatHR object 
##This is convenient to use other adehabitatHR capabilities such as overlap indices or percent volume contours

##Cast over to SPxDF
hpikde.px <- as(hpikde.raster,"SpatialPixelsDataFrame")

##create new estUD using the SPxDF
hpikde.ud <- new("estUD", hpikde.px)

##Assign values to a couple slots of the estUD
hpikde.ud@vol = FALSE
hpikde.ud@h$meth = "Plug-in Bandwidth"

##Convert the UD values to volume using getvolumeUD from adehabitatHR and cast over to a raster
hpikde.ud.vol <- getvolumeUD(hpikde.ud, standardize=TRUE)
hpikde.ud.vol.raster <- raster(hpikde.ud.vol)

##Here we generate volume contours using the UD
hpikde.50vol <- getverticeshr(hpikde.ud, percent = 50,ida = NULL, unin = "m", unout = "ha", standardize=TRUE)
hpikde.80vol <- getverticeshr(hpikde.ud, percent = 80,ida = NULL, unin = "m", unout = "ha", standardize=TRUE)
hpikde.90vol <- getverticeshr(hpikde.ud, percent = 90,ida = NULL, unin = "m", unout = "ha", standardize=TRUE)
hpikde.95vol <- getverticeshr(hpikde.ud, percent = 95,ida = NULL, unin = "m", unout = "ha", standardize=TRUE)
hpikde.99vol <- getverticeshr(hpikde.ud, percent = 99,ida = NULL, unin = "m", unout = "ha", standardize=TRUE)


##Let’s put the HR, volume, volume contours, trajectory, and points on a plot
plot.new()
breaks <- c(0, 50, 80, 90, 95, 99)
plot(hpikde.ud.vol.raster, col=heat.colors(3), breaks=breaks,interpolate=TRUE, main="Kernel Density Estimation, Plug-in Bandwidth",xlab="Coords X", ylab="Coords Y", legend.shrink=0.80,legend.args=list(text="UD by Volume (%)",side=4, font=2, line=2.5, cex=0.8))
plot(hpikde.99vol, add=T)#col="grey",axes=T)
plot(hpikde.90vol, add=TRUE)
plot(hpikde.80vol, add=TRUE)
plot(hpikde.50vol, add=TRUE)
points(loc, pch=1, cex=0.5)

writeOGR(hpikde.50vol,dsn="output", layer="cat143plug50",driver="ESRI Shapefile",overwrite=TRUE)
writeOGR(hpikde.80vol,dsn="output", layer="cat143plug80",driver="ESRI Shapefile",overwrite=TRUE)
writeOGR(hpikde.90vol,dsn="output", layer="cat143plug90",driver="ESRI Shapefile",overwrite=TRUE)
writeOGR(hpikde.95vol,dsn="output", layer="cat143plug95",driver="ESRI Shapefile",overwrite=TRUE)
writeOGR(hpikde.99vol,dsn="output", layer="cat143plug99",driver="ESRI Shapefile",overwrite=TRUE)
