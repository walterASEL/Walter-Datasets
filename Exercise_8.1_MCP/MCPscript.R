library(adehabitatHR)
library(maptools)

muleys <-read.csv("muleysexample.csv", header=T)
muleys
str(muleys)

newmuleys <-subset(muleys, muleys$Long > -110.90 & muleys$Lat > 37.8)
muleys <- newmuleys
newmuleys <-subset(muleys, muleys$Long < -107)
muleys <- newmuleys
data.xy = muleys[c("X","Y")]
#Creates class Spatial Points for all locations
xysp <- SpatialPoints(data.xy)
proj4string(xysp) <- CRS("+proj=utm +zone=17 +ellps=WGS84")

#Creates a Spatial Data Frame from 
sppt<-data.frame(xysp)
#Creates a spatial data frame of ID
idsp<-data.frame(muleys[2])
#Merges ID and Date into the same spatial data frame
merge<-data.frame(idsp)
#Adds ID and Date data frame with locations data frame
coordinates(merge)<-sppt
plot(merge)
str(merge)

## estimates the MCP
cp <- mcp(merge[,1], percent=95)#(95% is the default)
## The home-range size
as.data.frame(cp)
## Plot the home ranges
plot(cp)
## ... And the relocations
plot(merge, col=as.data.frame(merge@data)[,1], add=TRUE)
plot(cp[2,])#only plot MCP for deer D8

#Using the maptools package write to a shapefile
writePolyShape(cp, "MCPhomerange")

#In this example, we have chosen to exclude 5% of the most extreme relocations,
#but we could have made another choice. We may compute the home-range size
#for various choices of the number of extreme relocations to be excluded, using
#the function mcp.area:
hrs <- mcp.area(merge[,1], percent=seq(50, 100, by = 5))
hrs
