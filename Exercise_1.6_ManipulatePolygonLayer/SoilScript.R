library(rgdal)
library(maptools)
library(foreign)
#Example using rgdal, rgdal automatically imports the projection file
###change dsn to the directory where your example files are stored
soils<-readOGR(dsn=".",layer="Soil_Properties")
soils@proj4string ###get projection
plot(soils)
names(soils) ###get attribute data

soils$Clay <- soils$SdvOutpu_1
soils$pH <- soils$SdvOutpu_2
soils$CEC <- soils$SdvOutpu_3

#shapefiles descriptions
#Shapefiles contain several slots which can be called with the "@" symbol or slot(object, "data")
soils@data[1:10,] #= a data frame with n observations associated with X covariates,
soils@polygons #=the number of polygons that the shapefile consists of
soils@plotOrder #= the order of the polygons
soils@bbox #= boundary box
soils@proj4string #= projection
#Within the slot
soils@polygons[[1]] ###will bring up the first polygon
soils@polygons[[1]]@area ###will bring up the area for the first polygon
soils@polygons[[1]]@ID ##will retrieve the ID of the first polygon
soils@polygons[[1]]@plotOrder ##will retrieve the order of the first polygon

#Select the areas that Percent Clay polygons are over 30%
plot(soils)
high.clay<- soils[soils$Clay>30,]
plot(high.clay, border="red", add=TRUE)

high.CEC<- soils[soils$CEC>14,]
plot(high.CEC, border="green", add=TRUE)

high.pH <- soils[soils$pH>8,]
plot(high.pH, border="yellow", add=TRUE)

#Import mule deer locations from harvested animals tested for CWD
#Note: Locations have been offset or altered so do not reflect actual 
#locations of samples
mule <- read.csv("MDjitterclip.csv", header=T)
str(mule)

coords<-data.frame(x = mule$x, y = mule$y)
crs<-"+proj=utm +zone=13 +datum=WGS84 +no_defs +towgs84=0,0,0"
coords

plot(coords, col="blue")
par(new=TRUE)

#Sampling points in a Spatial Object###type=?regular? will give a regular grid
samples<-spsample(soils, n=1000, type="random")
samples@proj4string

plot(soils, col="wheat")
points(coords, col="blue")
points(samples, col="red")

samples@bbox <- soils@bbox
samples@proj4string <- soils@proj4string

#Matches points with polygons:
soils.idx<- over(samples,soils)
locs <- SpatialPoints(coords)
locs@proj4string <- soils@proj4string
soils.locs<- over(locs, soils)

#Extracts and tallies Clay soil types:
obs.tbl <- table(soils.idx$Clay[soils.idx$Clay])
obs.tbl

#Converts the counts to proportions:
obs <- obs.tbl/sum(obs.tbl)
obs

obs.tbl2 <- table(soils.locs$Clay[soils.locs$Clay])
obs.tbl2

#Converts the counts to proportions:
obs2 <- obs.tbl2/sum(obs.tbl2)
obs2

