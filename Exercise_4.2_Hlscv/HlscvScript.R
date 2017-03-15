library(adehabitatHR)
library(sp)

panther <-read.csv("pantherjitter.csv", header=T)
str(panther)

loc <- data.frame("x"=panther$X,"y"=panther$Y)
proj4string <- CRS("+proj=utm +zone=17N +ellps=WGS84")
pantherspdf <- SpatialPointsDataFrame(loc,panther, proj4string = proj4string)
plot(pantherspdf, col=pantherspdf$CatID)

## Example of estimation using LSCV
udbis2 <- kernelUD(pantherspdf[,2], h = "href", hlim = c(10,50),extent=1)
image(udbis2)

#Now change h to href for comparison to LSCV later

## Compare the estimation with ad hoc and LSCV method
## for the smoothing parameter
cuicui2 <- kernel.area(udbis2)
cuicui2

#Note that regardless of change hlim or extent, LSCV will not converge for these animals so let's try 
#a trick here. I believe LSCV is a poor estimator with GPS locations being too numerous 
#and very close together compared to traditional VHF datasets which LSCV were originally evaluated.
#So lets jitter locations 50 meters from their original location and try again.

## Example of estimation using LSCV
panther$jitterX <- jitter(panther$X, factor=500)
panther$jitterY <- jitter(panther$Y, factor=500)
locjitter <- data.frame("x"=panther$jitterX,"y"=panther$jitterY)
proj4string <- CRS("+proj=utm +zone=17N +ellps=WGS84")
jitterspdf <- SpatialPointsDataFrame(locjitter,panther, proj4string = proj4string)
plot(jitterspdf, col=pantherspdf$id)
points(pantherspdf, col="blue")
udbis3 <- kernelUD(jitterspdf[,2], h = "LSCV")#, hlim = c(1, 5),extent=1)
image(udbis3)

#Now rerun with jitter factor = 100 then 500 instead of 50 and see what happens?

cuicui3 <- kernel.area(udbis3) ## LSCV
cuicui3

#Now rerun with jitter factor = 500 instead of 100 and see what happens?

iso <- cbind(cuicui2,cuicui3)
colnames(iso) <- c("FP121_lscv","FP143_lscv","FP121_Jitter","FP143_Jitter")
iso
