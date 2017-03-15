library(rasterVis)
library(raster)
library(rgl)

#We can simply look at a Digital Elevatiom=n Model (DEM)
dem <- raster("elev_dem.tif")
plot(dem)

#Plot in three dimensions
plot3D(dem)

#We can also create a variety of covariates using DEMs (i.e., slope, aspect)
highGround = dem > 2000
plot3D(highGround)

slope = terrain(dem,opt='slope')
aspect = terrain(dem,opt='aspect')
hill = hillShade(slope,aspect,40,270)
plot(hill,col=grey(0:100/100),legend=FALSE)
plot3D(hill)

plot(dem,col=rainbow(25,alpha=0.35))
plot(slope,col=rainbow(25,alpha=0.35))
plot(aspect,col=rainbow(25,alpha=0.35))

x <- terrain(dem, opt=c('slope', 'aspect'), unit='degrees')
plot(x)

# TPI for different neighborhood size:
tpiw <- function(x, w=5) {
m <- matrix(1/(w^2-1), nc=w, nr=w)
m[ceiling(0.5 * length(m))] <- 0
f <- focal(x, m)
x - f
}

tpi5 <- tpiw(dem)
plot3D(tpi5)

#Topographic Position Index - difference between mean 
#and surrounding 8 cells
tpi = terrain(dem,opt='tpi')
summary(tpi)
plot3D(tpi)

#Terrain Ruggedness Index - mean of absolute difference 
#between the cell and surrounding 8 cells
tri = terrain(dem,opt='tri')
summary(tri)
plot(tri)
plot3D(tri)

#Roughness Index - difference between mean 
#and surrounding 8 cells
ruf = terrain(dem,opt="roughness")
summary(ruf)
plot3D(ruf)

#Load vegetation raster layer textfile clipped in ArcMap 
crop <-raster("crop2012utm12.tif")
plot(crop)
plot3D(dem,drape=crop)

drape <- cut(crop^5, 5)
plot3D(dem, drape=drape)
plot3D(dem, drape=crop, col.regions=terrain.colors)
plot3D(dem, at=4)
crop

# reclassify the values into 9 groups all values between 0 and 20 equal 1, etc.
m <- c(-Inf,0,NA,2, 7, 2, 20, 60, 3, 60, 70, 4, 110, 132, 5, 133, 150, 6, 151, 172, 7, 180, 183, 8, 189, 191, 9,192,Inf,NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(crop, rclmat)
plot(rc)
rc
plot3D(dem,drape=rc)
#Not sure why this does not work as expected
